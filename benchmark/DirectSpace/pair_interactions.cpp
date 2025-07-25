#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/hilbert_sfc.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#ifdef STORMM_USE_HPC
#  include "../../src/Potential/hpc_pme_potential.h"
#endif
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Potential/pme_potential.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/tile_manager.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/hpc_phasespace_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// A class to manage particle groups when subdividing a set of coordinates.
//-------------------------------------------------------------------------------------------------
class ParticleGroup {
public:

  // The constructor takes no arguments and initializes a blank object.
  ParticleGroup(int system_index_in, int target_group_size_in,
                double cutoff_in = default_pme_cutoff, double margin_in = 1.5);

  // Get the number of particles currently held by the group.
  int getParticleCount() const ;
  
  // Get the number of particles in the group's halo.
  int getHaloParticleCount() const;

  // Compute the centroid of all particles in the group.  Because particles are already imaged when
  // they enter the group,  the centroid can be computed without imaging transforms.  If they are
  // provided, the best centroid will be determined.
  double3 computeCentroid(const double* umat = nullptr, const double* invu = nullptr) const;
  
  // Estimate the size of the minimal sphere enclosing all particles in the group.  Return the
  // value as a four-tuple containing the Cartesian X, Y, and Z coordinates of the sphere center
  // in its "x", "y", and "z" members, and the radius of the sphere in its "w" member.

  // Add a particle to the group.
  //
  // Arguments:
  //   topl_idx:  Index of the particle in the coordinate synthesis or synthesis of topologies
  //   reim_loc:  Re-imaged location of the particle
  //   cell_loc:  Neighbor list cell location of the particle
  void addParticle(int topl_idx, const double3 reim_loc, const int3 cell_loc);

  // Find the identities of halo atoms based on the known atoms of the group, their neighbor list
  // cells, and the geometry of the neighbor list cell grid.  This function should be called after
  // the group's list of particles is complete.
  //
  // Arguments:
  //   cg:  The neighbor list describing the associated synthesis of systems
  template <typename T, typename Tacc, typename Tcalc, typename T4>
  void findHaloAtoms(const CellGrid<T, Tacc, Tcalc, T4> &cg);

  // Determine the number of tiles needed to evaluate all interactions for a particle group, the
  // number of interaction tests (tiles with low numbers of sending atoms will get evaluated by
  // reduced iterations), and the overall number of valid interactions for the group.  Return these
  // results in the "x", "y", and "z" members of a double-precision tuple (the format is conducive
  // to adding very large numbers, and doing divisions or averaging later).  This function should
  // only be called after calling findHaloAtoms() on the group.
  //
  // Arguments:
  //   poly_ps:  The associated coordiante synthesis (this must match the CellGrid object given
  //             to the findHaloAtoms() member function)
  double3 evaluateInteractions(const PhaseSpaceSynthesis &poly_ps);

private:
  int system_index;                       ///< Index number of the system within the associated
                                          ///<   synthesis
  int target_group_size;                  ///< The number of atoms that the group is intended to
                                          ///<   hold.  For any given system, all but perhaps one
                                          ///<   group will hold the target number of atoms, unless
                                          ///<   other criteria are brought into play that limit
                                          ///<   group size based on locality.
  int atom_count;                         ///< The number of atoms in the current group
  int halo_count;                         ///< The number of halo atoms connected to the group
  double cutoff;                          ///< The particle-particle interaction cutoff used to
                                          ///<   determine halo atoms
  double margin;                          ///< Margin used to inflate the particle-particle
                                          ///<   interaction cutoff when determining halo atoms
  std::vector<int> atom_topl_idx;         ///< Topological indices of atoms in the group
  std::vector<int> halo_topl_idx;         ///< Topological indices of atoms in the halo region
  std::vector<int3> atom_cell_locations;  ///< Locations of cells containing group atoms, in the
                                          ///<   indexing system of the neighbor list cell grid (a
                                          ///<   system-specific offset may apply)
  std::vector<double3> imaged_locations;  ///< Locations of all atoms, re-imaged to the primary
                                          ///<   unit cell
};

//-------------------------------------------------------------------------------------------------
ParticleGroup::ParticleGroup(const int system_index_in, const int target_group_size_in,
                             const double cutoff_in, const double margin_in) :
    system_index{system_index_in}, target_group_size{target_group_size_in}, atom_count{0},
    halo_count{0}, cutoff{cutoff_in}, margin{margin_in}, atom_topl_idx{}, halo_topl_idx{},
    atom_cell_locations{}, imaged_locations{}
{}

//-------------------------------------------------------------------------------------------------
int ParticleGroup::getParticleCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int ParticleGroup::getHaloParticleCount() const {
  return halo_count;
}

//-------------------------------------------------------------------------------------------------
double3 ParticleGroup::computeCentroid(const double* umat, const double* invu) const {
  if (umat != nullptr && invu != nullptr) {

    // Determine the centroid in periodic boundary conditions
    std::vector<double> fx(atom_count), fy(atom_count), fz(atom_count);
    for (int i = 0; i < atom_count; i++) {

    }
  }

  // Average the locations of all particles
  double mean_x = 0.0;
  double mean_y = 0.0;
  double mean_z = 0.0;
  for (int i= 0; i < atom_count; i++) {
    mean_x += imaged_locations[i].x;
    mean_y += imaged_locations[i].y;
    mean_z += imaged_locations[i].z;
  }
  const double inv_atom_count = 1.0 / static_cast<double>(atom_count);
  mean_x *= inv_atom_count;
  mean_y *= inv_atom_count;
  mean_z *= inv_atom_count;
  return { mean_x, mean_y, mean_z }; 
}

//-------------------------------------------------------------------------------------------------
void ParticleGroup::addParticle(const int topl_idx, const double3 reim_loc, const int3 cell_loc) {
  atom_topl_idx.push_back(topl_idx);
  imaged_locations.push_back(reim_loc);
  atom_cell_locations.push_back(cell_loc);
  atom_count++;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void ParticleGroup::findHaloAtoms(const CellGrid<T, Tacc, Tcalc, T4> &cg) {
  const CellGridReader cgr = cg.data();
  const ullint sys_dims = cgr.system_cell_grids[system_index];
  const int cell_na = ((sys_dims >> 28) & 0xfff);
  const int cell_nb = ((sys_dims >> 40) & 0xfff);
  const int cell_nc = (sys_dims >> 52);
  const int cell_start = (sys_dims & 0xfffffff);
  const int xfrm_stride = roundUp(9, 32);
  const T* cell_invu = &cgr.system_cell_invu[xfrm_stride * system_index];
  const bool tcoord_is_integral = isSignedIntegralScalarType<T>();
  const double sq_bound = (cutoff + margin) * (cutoff + margin);
  std::vector<int> halo_atom_list;
  halo_atom_list.reserve(1024);
  for (int i = 0; i < atom_count; i++) {

    // Find the atom within its home cell and take its location from the cell grid.
    const uint iatom_img_idx = cgr.img_atom_idx[atom_topl_idx[i]];
    double iatom_x, iatom_y, iatom_z;
    const T4 i_crdq = cgr.image[iatom_img_idx];
    if (tcoord_is_integral) {
      iatom_x = static_cast<double>(i_crdq.x) * cgr.inv_lpos_scale;
      iatom_y = static_cast<double>(i_crdq.y) * cgr.inv_lpos_scale;
      iatom_z = static_cast<double>(i_crdq.z) * cgr.inv_lpos_scale;
    }
    else {
      iatom_x = i_crdq.x;
      iatom_y = i_crdq.y;
      iatom_z = i_crdq.z;
    }
    
    // Find the neighbor list cell, then search within the cell (for atoms of lower neighbor list
    // index) and in the half shell of the surrounding cells, +/- 3 cell to account for whatever
    // pair list margin.
    const int icell_a = atom_cell_locations[i].x;
    const int icell_b = atom_cell_locations[i].y;
    const int icell_c = atom_cell_locations[i].z;
    const int icell = cell_start + (((icell_c * cell_nb) + icell_b) * cell_na) + icell_a;
    for (int acon = -3; acon <= 3; acon++) {
      int jcell_a = icell_a + acon;
      jcell_a += ((jcell_a < 0) - (jcell_a >= cell_na)) * cell_na;
      for (int bcon = -3; bcon <= 3; bcon++) {
        int jcell_b = icell_b + bcon;
        jcell_b += ((jcell_b < 0) - (jcell_b >= cell_nb)) * cell_nb;
        for (int ccon = -3; ccon <= 0; ccon++) {
          int jcell_c = icell_c + ccon;
          jcell_c += ((jcell_c < 0) - (jcell_c >= cell_nc)) * cell_nc;
          const int jcell = cell_start + (((jcell_c * cell_nb) + jcell_b) * cell_na) + jcell_a;
          
          // Avoid double-counting in the grand list of cell-to-cell interactions, as with the
          // tower-plate scheme.
          if (ccon == 0 && (bcon > 0 || (bcon == 0 && acon > 0))) {
            continue;
          }
          const double da = acon;
          const double db = bcon;
          const double dc = ccon;
          const double da_x = (cell_invu[0] * da) + (cell_invu[3] * db) + (cell_invu[6] * dc);
          const double da_y =                       (cell_invu[4] * db) + (cell_invu[7] * dc);
          const double da_z =                                             (cell_invu[8] * dc);

          // Loop over all atoms in the cell, place them relative to the atom of interest, and
          // detemine whether they are in range of the ith atom in the group.  The list of all
          // applicable atoms will continue to grow, eventually to be sorted and then culled of
          // duplicates.  This process is obviously inefficient, but the purpose is to get a
          // measure of the tile enrichment.  Some notion with the work required to make the pair
          // list in this way may come along with the experience.
          const uint2 jcell_lims = cgr.cell_limits[jcell];
          const uint jlim = (jcell == icell) ? iatom_img_idx : jcell_lims.x + (jcell_lims.y >> 16);
          for (uint j = jcell_lims.x; j < jlim; j++) {
            const T4 j_crdq = cgr.image[j];
            double jatom_x, jatom_y, jatom_z;
            if (tcoord_is_integral) {
              jatom_x = static_cast<double>(j_crdq.x) * cgr.inv_lpos_scale;
              jatom_y = static_cast<double>(j_crdq.y) * cgr.inv_lpos_scale;
              jatom_z = static_cast<double>(j_crdq.z) * cgr.inv_lpos_scale;
            }
            else {
              jatom_x = j_crdq.x;
              jatom_y = j_crdq.y;
              jatom_z = j_crdq.z;
            }
            jatom_x += da_x;
            jatom_y += da_y;
            jatom_z += da_z;
            const double dx = jatom_x - iatom_x;
            const double dy = jatom_y - iatom_y;
            const double dz = jatom_z - iatom_z;
            const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            if (r2 < sq_bound) {
              halo_atom_list.push_back(cgr.nonimg_atom_idx[j]);
            }
          }
        }
      }
    }
  }
  halo_topl_idx = reduceUniqueValues(halo_atom_list);
  halo_count = halo_topl_idx.size();
  int nboot = 0;
  std::vector<bool> keep_atoms(halo_topl_idx.size(), true);
  for (int i = 0; i < atom_count; i++) {
    for (int j = 0; j < halo_count; j++) {
      if (halo_topl_idx[j] == atom_topl_idx[i]) {
        keep_atoms[j] = false;
        nboot++;
      }
    }
  }
  int keepcon = 0;
  for (int i = 0; i < halo_count; i++) {
    if (keep_atoms[i]) {
      halo_topl_idx[keepcon] = halo_topl_idx[i];
      keepcon++;
    }
  }
  halo_count -= nboot;
  halo_topl_idx.resize(halo_count);
}

//-------------------------------------------------------------------------------------------------
double3 ParticleGroup::evaluateInteractions(const PhaseSpaceSynthesis &poly_ps) {
  const PsSynthesisReader poly_psr = poly_ps.data();
  const int xfrm_stride = roundUp(9, warp_size_int);
  const double* umat = &poly_psr.umat[xfrm_stride * system_index];
  const double* invu = &poly_psr.invu[xfrm_stride * system_index];
  double3 result = { 0.0, 0.0, 0.0 };
  const int tile_count = (halo_count + warp_size_int - 1) / warp_size_int;
  result.x = tile_count;
  const int last_tile_halo_atoms = (halo_count % warp_size_int);
  int total_tests = (tile_count - 1) * target_group_size * warp_size_int;
  if (last_tile_halo_atoms == 1) {
    total_tests += warp_size_int;
  }
  else {
    total_tests += target_group_size * warp_size_int;
  }
  total_tests += target_group_size * target_group_size;
  result.y = total_tests;
  int halo_base = 0;
  const double sq_bound = cutoff * cutoff;
  int success = 0;
  while (halo_base < halo_count) {
    const int jlim = std::min(halo_base + warp_size_int, halo_count);
    for (int i = 0; i < atom_count; i++) {
      const int iatom = atom_topl_idx[i];

      // Pull out the group's atom.  Save coordinate scaling for the displacement calculation.
      double group_x, group_y, group_z;
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        group_x = llround(poly_psr.xcrd[iatom]);
        group_y = llround(poly_psr.ycrd[iatom]);
        group_z = llround(poly_psr.zcrd[iatom]);
      }
      else {
        group_x = hostInt95ToDouble(poly_psr.xcrd[iatom], poly_psr.xcrd_ovrf[iatom]);
        group_y = hostInt95ToDouble(poly_psr.ycrd[iatom], poly_psr.ycrd_ovrf[iatom]);
        group_z = hostInt95ToDouble(poly_psr.zcrd[iatom], poly_psr.zcrd_ovrf[iatom]);
      }
      for (int j = halo_base; j < jlim; j++) {
        const int jatom = halo_topl_idx[j];
        double halo_x, halo_y, halo_z;
        if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
          halo_x = llround(poly_psr.xcrd[jatom]);
          halo_y = llround(poly_psr.ycrd[jatom]);
          halo_z = llround(poly_psr.zcrd[jatom]);
        }
        else {
          halo_x = hostInt95ToDouble(poly_psr.xcrd[jatom], poly_psr.xcrd_ovrf[jatom]);
          halo_y = hostInt95ToDouble(poly_psr.ycrd[jatom], poly_psr.ycrd_ovrf[jatom]);
          halo_z = hostInt95ToDouble(poly_psr.zcrd[jatom], poly_psr.zcrd_ovrf[jatom]);
        }
        double dx = (halo_x - group_x) * poly_psr.inv_gpos_scale;
        double dy = (halo_y - group_y) * poly_psr.inv_gpos_scale;
        double dz = (halo_z - group_z) * poly_psr.inv_gpos_scale;
        imageCoordinates<double, double>(&dx, &dy, &dz, umat, invu, poly_psr.unit_cell,
                                         ImagingMethod::MINIMUM_IMAGE);
        if ((dx * dx) + (dy * dy) + (dz * dz) < sq_bound) {
          success++;
        }
      }
    }
    halo_base += warp_size_int;
  }
  success += atom_count * (atom_count - 1) / 2;
  result.z = success;
  return result;
}

//-------------------------------------------------------------------------------------------------
// Determine the number and sizes of atom batches needed to span a given total number of atoms.
//
// Arguments:
//   total_atoms:       The total number of atoms to cover
//   max_batch:         The maximum batch size to use in covering all atoms
//   min_batch:         The minimum batch size to use in covering all atoms
//   trim_final_batch:  Indicate that the final batch should be trimmed to the smallest possible
//                      size, conserving store operations while trimming the total tile area
//-------------------------------------------------------------------------------------------------
std::vector<int> batchList(const int total_atoms, const int max_batch, const int min_batch,
                           const bool trim_final_batch) {
  int ncovered = 0;
  int bsize = max_batch;
  int bmask = min_batch - 1;
  std::vector<int> result;
  const int twice_min_batch = min_batch * 2;
  while (ncovered < total_atoms) {
    const int atoms_left = total_atoms - ncovered;
    if (atoms_left >= bsize || (bsize == min_batch * 4 && atoms_left < bsize &&
                                atoms_left > 3 * min_batch) || bsize == min_batch ||
        (bsize == twice_min_batch && atoms_left > min_batch)) {
      if (atoms_left >= bsize || trim_final_batch == false) {
        result.push_back(bsize);
      }
      else {
        int best_final_batch = bsize;
        while (best_final_batch / 2 >= atoms_left) {
          best_final_batch >>= 1;
        }
        result.push_back(best_final_batch);
      }
      ncovered += bsize;
    }
    else if (bsize >= min_batch) {
      bsize >>= 1;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Analyze the theoretical occupancy of the tower-plate kernel.
//
// Arguments:
//   cg:                The cell grid to analyze
//   sys_idx:           Index of the system of interest within the cell grid
//   cutoff:            The particle-particle cutoff determining valid interactions
//   min_batch:         The minimum batch size to consider when making tiles
//   max_batch:         The minimum batch size to consider when making tiles
//   trim_final_batch:  Indicate whether to trim the final batch along any portion of the
//                      tower-plate decomposition
//   cull_irrelevant:   Indicate whether to cull atoms in the tower and plate which could not
//                      possibly be within range of those they would need to interact with
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void reportTowerPlateOccupancy(const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg,
                               const double cutoff, const int min_batch, const int max_batch,
                               const bool trim_final_batch, const bool cull_irrelevant) {

  // Construct the tower-plate arrangement
  const std::vector<int> tw_rel_a = {  0,  0,  0,  0,  0 };
  const std::vector<int> pl_rel_a = { -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1 };
  const std::vector<int> tw_rel_b = {  0,  0,  0,  0,  0 };
  const std::vector<int> pl_rel_b = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0 };
  const std::vector<int> tw_rel_c = { -2, -1,  0,  1,  2 };
  const std::vector<int> pl_rel_c = {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };

  // Loop over all systems
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cgr = cg.data();
  int overall_max_pop = 0;
  std::vector<double> running_sum, squared_sum;
  std::vector<double> tp_ideal_work(cgr.system_count, 0.0), tp_total_work(cgr.system_count, 0.0);
  std::vector<double> tp_successful(cgr.system_count, 0.0), tp_stores(cgr.system_count, 0.0);
  for (int sys_idx = 0; sys_idx < cgr.system_count; sys_idx++) {
    const ullint sys_dims = cgr.system_cell_grids[sys_idx];
    const int na_cell = ((sys_dims >> 28) & 0xfff);
    const int nb_cell = ((sys_dims >> 40) & 0xfff);
    const int nc_cell = ((sys_dims >> 52) & 0xfff);
    const int cell_start = (sys_dims & 0xfffffff);
    const int nabc_cells = na_cell * nb_cell * nc_cell;
    std::vector<int> tower_cells(5), plate_cells(12), tower_offsets(5), plate_offsets(12);
    std::vector<int> tower_prefix(6), plate_prefix(13), tower_;
    const int xfrm_offset = sys_idx * roundUp(9, warp_size_int);
    std::vector<double> cg_invu(9);
    for (int i = 0; i < 9; i++) {
      cg_invu[i] = cgr.system_cell_invu[xfrm_offset + i];
    }
    const double cut_sq = cutoff * cutoff;

    // Gauge the maximum cell population
    int max_pop = 0;
    for (int i = 0; i < nabc_cells; i++) {
      max_pop = std::max(max_pop, static_cast<int>(cgr.cell_limits[cell_start + i].y >> 16));
    }
    if (max_pop > overall_max_pop) {
      overall_max_pop = max_pop;
      running_sum.resize(overall_max_pop + 1, 0.0);
      squared_sum.resize(overall_max_pop + 1, 0.0);
    }
    std::vector<int> populations(max_pop + 1, 0);
    for (int i = 0; i < na_cell; i++) {
      for (int j = 0; j < nb_cell; j++) {
        for (int k = 0; k < nc_cell; k++) {
          
          // Determine the cell identities
          for (int m = 0; m < 5; m++) {
            int c_a = i + tw_rel_a[m];
            int c_b = j + tw_rel_b[m];
            int c_c = k + tw_rel_c[m];
            c_a += ((c_a < 0) - (c_a >= na_cell)) * na_cell;
            c_b += ((c_b < 0) - (c_b >= nb_cell)) * nb_cell;
            c_c += ((c_c < 0) - (c_c >= nc_cell)) * nc_cell;
            tower_cells[m] = (((c_c * nb_cell) + c_b) * na_cell) + c_a + cell_start;
            const uint2 cldat = cgr.cell_limits[tower_cells[m]];
            tower_offsets[m] = cldat.x;
            if (cull_irrelevant) {
              const uint nlim = cldat.x + (cldat.y >> 16);
              if (m == 0) {
                int nrel = 0;
                for (uint n = cldat.x; n < nlim; n++) {
                  nrel += ((cgr.relevance[n] >> 12) & 0x1);
                }
                tower_prefix[m] = nrel;
              }
              else if (m == 4) {
                int nrel = 0;
                for (uint n = cldat.x; n < nlim; n++) {
                  nrel += ((cgr.relevance[n] >> 15) & 0x1);
                }
                tower_prefix[m] = nrel;
              }
              else {
                tower_prefix[m] = (cldat.y >> 16);
              }
            }
            else {
              tower_prefix[m] = (cldat.y >> 16);
            }
          }
          tower_prefix[5] = 0;
          prefixSumInPlace<int>(&tower_prefix, PrefixSumType::EXCLUSIVE);
          populations[tower_prefix[3] - tower_prefix[2]] += 1;
          for (int m = 0; m < 12; m++) {
            int c_a = i + pl_rel_a[m];
            int c_b = j + pl_rel_b[m];
            int c_c = k + pl_rel_c[m];
            c_a += ((c_a < 0) - (c_a >= na_cell)) * na_cell;
            c_b += ((c_b < 0) - (c_b >= nb_cell)) * nb_cell;
            c_c += ((c_c < 0) - (c_c >= nc_cell)) * nc_cell;
            plate_cells[m] = (((c_c * nb_cell) + c_b) * na_cell) + c_a + cell_start;
            const uint2 cldat = cgr.cell_limits[plate_cells[m]];
            plate_offsets[m] = cldat.x;
            if (cull_irrelevant) {
              int nrel = 0;
              const uint nlim = cldat.x + (cldat.y >> 16);
              for (uint n = cldat.x; n < nlim; n++) {
                nrel += ((cgr.relevance[n] >> m) & 0x1);
              }
              plate_prefix[m] = nrel;
            }
            else {
              plate_prefix[m] = (cldat.y >> 16);
            }
          }
          plate_prefix[12] = 0;
          prefixSumInPlace<int>(&plate_prefix, PrefixSumType::EXCLUSIVE);
          tp_ideal_work[sys_idx] += static_cast<double>(tower_prefix[5] * plate_prefix[12]);

          // Determine the number and sizes of tower sets for tower-plate interactions
          const std::vector<int> tower_sets = batchList(tower_prefix[5], max_batch, min_batch,
                                                        trim_final_batch);
          const std::vector<int> plate_sets = batchList(plate_prefix[12], max_batch, min_batch,
                                                        trim_final_batch);
          const std::vector<int> centr_sets = batchList(tower_prefix[3] - tower_prefix[2],
                                                        max_batch, min_batch, trim_final_batch);
          const std::vector<int> lower_sets = batchList(tower_prefix[2], max_batch, min_batch,
                                                        trim_final_batch);
          const double tower_atoms_tiled = sum<int>(tower_sets);
          const double plate_atoms_tiled = sum<int>(plate_sets);
          const double centr_atoms_tiled = sum<int>(centr_sets);
          const double lower_atoms_tiled = sum<int>(lower_sets);
          tp_total_work[sys_idx] += (tower_atoms_tiled * plate_atoms_tiled) +
                                    (centr_atoms_tiled * ((centr_atoms_tiled / 2) +
                                                          lower_atoms_tiled));
          tp_stores[sys_idx] += (tower_sets.size() * (plate_sets.size() + 1)) +
                                (centr_sets.size() * (lower_sets.size() + 1)) +
                                (centr_sets.size() * (centr_sets.size() + 1) / 2);

          // Stage the tower and plate atoms
          std::vector<double> tower_x(tower_prefix[5]);
          std::vector<double> tower_y(tower_prefix[5]);
          std::vector<double> tower_z(tower_prefix[5]);
          std::vector<double> plate_x(plate_prefix[12]);
          std::vector<double> plate_y(plate_prefix[12]);
          std::vector<double> plate_z(plate_prefix[12]);
          size_t npt = 0;
          for (int m = 0; m < 5; m++) {
            const uint2 cldat = cgr.cell_limits[tower_cells[m]];
            const uint poslim = cldat.x + (cldat.y >> 16);
            const double stack_mult = m - 2;
            const double x_del = stack_mult * cg_invu[6];
            const double y_del = stack_mult * cg_invu[7];
            const double z_del = stack_mult * cg_invu[8];
            for (uint pos = cldat.x; pos < poslim; pos++) {
              if (cull_irrelevant == false || (m != 0 && m != 4) ||
                  (m == 0 && ((cgr.relevance[pos] >> 12) & 0x1)) ||
                  (m == 4 && ((cgr.relevance[pos] >> 15) & 0x1))) {
                const Tcoord4 atom_img = cgr.image[pos];
                tower_x[npt] = atom_img.x + x_del;
                tower_y[npt] = atom_img.y + y_del;
                tower_z[npt] = atom_img.z + z_del;
                npt++;
              }
            }
          }
          npt = 0;
          for (int m = 0; m < 12; m++) {
            const uint2 cldat = cgr.cell_limits[plate_cells[m]];
            const uint poslim = cldat.x + (cldat.y >> 16);
            const double stack_mult_x = static_cast<double>(m - ((m / 5) * 5) - 2);
            const double stack_mult_y = static_cast<double>((m / 5) - 2);
            const double x_del = (stack_mult_x * cg_invu[0]) + (stack_mult_y * cg_invu[3]);
            const double y_del =                               (stack_mult_y * cg_invu[4]);
            for (uint pos = cldat.x; pos < poslim; pos++) {
              if (cull_irrelevant == false || ((cgr.relevance[pos] >> m) & 0x1)) {
                const Tcoord4 atom_img = cgr.image[pos];
                plate_x[npt] = atom_img.x + x_del;
                plate_y[npt] = atom_img.y + y_del;
                plate_z[npt] = atom_img.z;
                npt++;
              }
            }
          }

          // Calculate the total number of interactions that are in range
          int n_success = 0;
          for (int m = 0; m < tower_prefix[5]; m++) {
            for (int  n = 0; n < plate_prefix[12]; n++) {
              const double dx = plate_x[n] - tower_x[m];
              const double dy = plate_y[n] - tower_y[m];
              const double dz = plate_z[n] - tower_z[m];
              const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
              if (r2 < cut_sq) {
                n_success++;
              }
            }
          }
          for (int m = tower_prefix[2]; m < tower_prefix[3]; m++) {
            for (int n = tower_prefix[2]; n < m; n++) {
              const double dx = tower_x[n] - tower_x[m];
              const double dy = tower_y[n] - tower_y[m];
              const double dz = tower_z[n] - tower_z[m];
              const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
              if (r2 < cut_sq) {
                n_success++;
              }
            }
            for (int n = tower_prefix[0]; n < tower_prefix[2]; n++) {
              const double dx = tower_x[n] - tower_x[m];
              const double dy = tower_y[n] - tower_y[m];
              const double dz = tower_z[n] - tower_z[m];
              const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
              if (r2 < cut_sq) {
                n_success++;
              }
            }
          }
          tp_successful[sys_idx] += static_cast<double>(n_success);

          // Calculate the number of trivial exclusions based on distance to line and distance to
          // plane.  The first plane of interest is defined by the unit cell A and C vectors, the
          // second by the unit cell B and C vectors.  The AB plane (the Cartesian XY plane) is
          // also of interest and most simple to compute in the general case (just consider the
          // particle's Cartesian Z coordinate), but as a culling criterion for the tower itself
          // it could be useful.  The tower atoms would need to be re-ordered so that the small
          // percentage of tower atoms that fall into this category of things that could never
          // interact with anything in the plate are more likely to be concentrated into a single
          // batch.  Otherwise, the distances to six critical lines will define whether a particle
          // is within range of anything in the central column.
        }
      }
    }

    // Contribute the populations to running sums and standard deviations
    for (int i = 0; i < max_pop; i++) {
      const double dpi = populations[i];
      running_sum[i] += dpi;
      squared_sum[i] += dpi * dpi;
    }
  }
  std::vector<double> ideal_oo_total(cgr.system_count), success_oo_ideal(cgr.system_count);
  std::vector<double> success_oo_total(cgr.system_count);
  for (int i = 0; i < cgr.system_count; i++) {
    ideal_oo_total[i] = tp_ideal_work[i] / tp_total_work[i];
    success_oo_ideal[i] = tp_successful[i] / tp_ideal_work[i];
    success_oo_total[i] = tp_successful[i] / tp_total_work[i];
  }
  const VarianceMethod stdev = VarianceMethod::STANDARD_DEVIATION;
  const double v_idt = (cgr.system_count > 2) ? variance(ideal_oo_total, stdev) : 0.0;
  const double v_sid = (cgr.system_count > 2) ? variance(success_oo_ideal, stdev) : 0.0;
  const double v_it  = (cgr.system_count > 2) ? variance(success_oo_total, stdev) : 0.0;
  const double v_str = (cgr.system_count > 2) ? variance(tp_stores, stdev) : 0.0;
  printf("  Cutoff %9.4lf : Real Atom Content       %7.2lf %% +/- %7.2lf %%\n"
         "                     Success in Real Pairs   %7.2lf %% +/- %7.2lf %%\n"
         "                     Success in All Pairs    %7.2lf %% +/- %7.2lf %%\n"
         "                     Total store events    %9.1lf %% +/- %7.0lf\n\n", cutoff,
         mean(ideal_oo_total) * 100.0, v_idt * 100.0, mean(success_oo_ideal) * 100.0,
         v_sid * 100.0, mean(success_oo_total) * 100.0, v_it * 100.0, mean(tp_stores),
         v_str);
  printf("  Cell Population  Count\n  ---------------  -----\n");
  const double dnsys = cgr.system_count;
  for (int i = 0; i <= overall_max_pop; i++) {
    printf("        %3d    %8.2lf +/- %7.2lf\n", i, running_sum[i] / dnsys,
           ((dnsys * squared_sum[i]) - (running_sum[i] * running_sum[i])) / (dnsys * dnsys));
  }
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// Analyze the occupancy that might be obtained with a method focused on optimizing the
// co-localization of atom groups of a particular size using a Hilbert space-filling curve.
//
// Arguments:
//   cg:          The cell grid to analyze
//   cutoff:      The particle-particle cutoff determining valid interactions
//   margin:      A hypothetical pair list margin to impose.  This value inflates the cutoff such
//                that all particles which are within the cutoff plus the margin of any particle
//                in the group are considered interacting.
//   group_size:  The size of the groups to develop
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void reportHilbertSpaceOccupancy(const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg,
                                 const double cutoff, const double margin, const int group_size,
                                 const int hs_grid_x = 256, const int hs_grid_y = 256,
                                 const int hs_grid_z = 256) {
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cgr = cg.data();
  const PhaseSpaceSynthesis *poly_ps = cg.getCoordinateSynthesisPointer();
  const PsSynthesisReader poly_psr = poly_ps->data();
  HilbertSFC curve(hs_grid_x, hs_grid_y, hs_grid_z, HilbertCurveMode::STRETCH);
  const int xfrm_stride = roundUp(9, warp_size_int);
  const double dhs_grid_x = hs_grid_x;
  const double dhs_grid_y = hs_grid_y;
  const double dhs_grid_z = hs_grid_z;
  std::vector<double3> all_contrib(cgr.system_count, { 0.0, 0.0, 0.0 });
  for (int sys_idx = 0; sys_idx < cgr.system_count; sys_idx++) {
    const ullint sys_dims = cgr.system_cell_grids[sys_idx];
    const int na_cell = ((sys_dims >> 28) & 0xfff);
    const int nb_cell = ((sys_dims >> 40) & 0xfff);
    const int nc_cell = ((sys_dims >> 52) & 0xfff);
    const double dna_cell = na_cell;
    const double dnb_cell = nb_cell;
    const double dnc_cell = nc_cell;
    const int cell_start = (sys_dims & 0xfffffff);
    const int nabc_cells = na_cell * nb_cell * nc_cell;
    
    // Assign all atoms with a key into the Hilbert space-filling curve.
    const int atom_llim = poly_psr.atom_starts[sys_idx];
    const int natom = poly_psr.atom_counts[sys_idx];
    const int atom_hlim = atom_llim + natom;
    std::vector<int2> hs_keys(natom);
    const double* umat = &poly_psr.umat[xfrm_stride * sys_idx];
    const double* invu = &poly_psr.invu[xfrm_stride * sys_idx];
    std::vector<double3> reim_loc(natom);
    std::vector<int3> cell_loc(natom);
    for (int i = atom_llim; i < atom_hlim; i++) {
      double dx, dy, dz;
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        dx = llround(poly_psr.xcrd[i]) * poly_psr.inv_gpos_scale;
        dy = llround(poly_psr.ycrd[i]) * poly_psr.inv_gpos_scale;
        dz = llround(poly_psr.zcrd[i]) * poly_psr.inv_gpos_scale;
      }
      else {
        dx = hostInt95ToDouble(poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]) * poly_psr.inv_gpos_scale;
        dy = hostInt95ToDouble(poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]) * poly_psr.inv_gpos_scale;
        dz = hostInt95ToDouble(poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]) * poly_psr.inv_gpos_scale;
      }
      double ndx = (umat[0] * dx) + (umat[3] * dy) + (umat[6] * dz);
      double ndy =                  (umat[4] * dy) + (umat[7] * dz);
      double ndz =                                   (umat[8] * dz);
      ndx -= floor(ndx);
      ndy -= floor(ndy);
      ndz -= floor(ndz);
      reim_loc[i - atom_llim].x = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
      reim_loc[i - atom_llim].y =                   (invu[4] * ndy) + (invu[7] * ndz);
      reim_loc[i - atom_llim].z =                                     (invu[8] * ndz);
      cell_loc[i - atom_llim].x = ndx * dna_cell;
      cell_loc[i - atom_llim].y = ndy * dnb_cell;
      cell_loc[i - atom_llim].z = ndz * dnc_cell;
      ndx *= dhs_grid_x;
      ndy *= dhs_grid_y;
      ndz *= dhs_grid_z;
      int indx = ndx;
      int indy = ndy;
      int indz = ndz;
      indx -= (indx == hs_grid_x);
      indy -= (indy == hs_grid_y);
      indz -= (indz == hs_grid_z);
      hs_keys[i - atom_llim] = { curve.getCurveIndex(indx, indy, indz), i };
    }

    // Sort the keys to obtain the list of atoms, indexed by their appearance in the topology, into
    // the Hilbert space-filling curve.
    std::sort(hs_keys.begin(), hs_keys.end(), [](int2 a, int2 b) { return (a.x < b.x); });

    // With the colocalized ordering, select atoms one batch at a time.
    int grp_base = 0;
    int grp_idx = 0;
    std::vector<ParticleGroup> pgroups;
    pgroups.reserve((natom + group_size - 1) / group_size);
    while (grp_base < natom) {
      pgroups.emplace_back(sys_idx, group_size, cutoff, margin);

      // Add particles to the group
      for (int j = 0; j < group_size; j++) {
        if (grp_base + j < natom) {
          const int topl_idx = hs_keys[grp_base + j].y;
          pgroups[grp_idx].addParticle(topl_idx, reim_loc[topl_idx - atom_llim],
                                       cell_loc[topl_idx - atom_llim]);
        }
      }
      pgroups[grp_idx].findHaloAtoms(cg);
      const double3 contrib = pgroups[grp_idx].evaluateInteractions(*poly_ps);
      
      all_contrib[sys_idx].x += contrib.x;
      all_contrib[sys_idx].y += contrib.y;
      all_contrib[sys_idx].z += contrib.z;
      
      // Increment the group and atom counters
      grp_idx++;
      grp_base += group_size;
    }
  }
  std::vector<double> all_proportion(cgr.system_count);
  for (int i = 0; i < cgr.system_count; i++) {
    all_proportion[i] = 100.0 * all_contrib[i].z / all_contrib[i].y;
  }
  const double hsfc_mean = mean(all_proportion);
  const double hsfc_var = (cgr.system_count > 2) ? variance(all_proportion,
                                                            VarianceMethod::STANDARD_DEVIATION) :
                                                   0.0;
  printf("Hilbert space-filling curve:  %7.2lf +/- %7.2lf %% success\n\n", hsfc_mean, hsfc_var);
}

//-------------------------------------------------------------------------------------------------
// Check the forces computed by the GPU kernel against CPU results.
//
// Arguments:
//   poly_ps:    The synthesis of coordinates, including forces on all particles
//   eval_frc:   Indicate whether forces were evaluated as part of the run
//   timer:      Time tracking object, to absorb the time running the prior CPU calculation as well
//               as the test to compare the results
//   time_test:  The integer code for the test time tracking
//-------------------------------------------------------------------------------------------------
void checkPairBenchmarkForces(const PhaseSpaceSynthesis &poly_ps, const EvaluateForce eval_frc,
                              StopWatch *timer, const int time_test) {
#ifdef STORMM_USE_HPC
  const UnitCellType system_uc = poly_ps.getUnitCellType();
  CoordinateFrame cpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  CoordinateFrame gpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  const CoordinateFrameReader cpu_cfr = cpu_frc.data();
  const CoordinateFrameReader gpu_cfr = gpu_frc.data();
  std::vector<double> cpu_chk(cpu_cfr.natom, 0.0);
  std::vector<double> gpu_chk(cpu_cfr.natom, 0.0);
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      coordCopy(&cpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::HOST);
      coordCopy(&gpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::DEVICE);
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.xcrd[j];
        gpu_chk[j] = gpu_cfr.xcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian X axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.ycrd[j];
        gpu_chk[j] = gpu_cfr.ycrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian Y axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.zcrd[j];
        gpu_chk[j] = gpu_cfr.zcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian Z axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      break;
    case EvaluateForce::NO:
      break;
    }
  }
#endif
  timer->assignTime(time_test);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Baseline variables
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif
  
  // Lay out time categories for program profiling
  const int time_load  = timer.addCategory("Load input files");
  const int time_build = timer.addCategory("Class object setup");
  const int time_test  = timer.addCategory("CPU testing");
  const int time_init  = timer.addCategory("Initialization");
  const int time_pairs = timer.addCategory("Pairs kernel with initialization");
  
  // Take in command line inputs
  CommandLineParser clip("pair_interactions", "A benchmarking program for measuring the rate at "
                         "which various GPU kernels can compute all pairwise particle-particle "
                         "interactions within one or more systems.", { "-timings" });
  clip.addStandardAmberInputs("-c", "-p", "-ig_seed");
  clip.addStandardBenchmarkingInputs({ "-iter", "-trials", "-replicas", "-skip_cpu_check",
                                       "-cutoff", "-elec_cutoff", "-vdw_cutoff", "-pad",
                                       "-eval_nrg", "-omit_frc", "-fp_bits", "-blur" });

  // Custom inputs for this benchmarking program
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-dual_cg", NamelistType::BOOLEAN);
  t_nml->addHelp("-dual_cg", "Force the use of dual neighbor lists, if the electrostatic and "
                 "van-der Waals cutoffs are not already distinct.");
  t_nml->addKeyword("-warp_mult", NamelistType::INTEGER, std::to_string(1));
  t_nml->addHelp("-warp_mult", "The number of warps to devote to processing each neighbor list "
                 "cell's assigned pair interactions.");
  t_nml->addKeyword("-occupancy", NamelistType::BOOLEAN);
  t_nml->addHelp("-occupancy", "Compute the occupancy of GPU warps during the pairs calculation.");
  t_nml->addKeyword("-min_batch", NamelistType::INTEGER, std::to_string(8));
  t_nml->addHelp("-min_batch", "The minimum batch size to consider when creating mock tiles for "
                 "the STORMM pairs calculation based on interactions among a tower-plate neutral "
                 "territory decomposition");
  t_nml->addKeyword("-max_batch", NamelistType::INTEGER, std::to_string(32));
  t_nml->addHelp("-max_batch", "The minimum batch size to consider when creating mock tiles for "
                 "the STORMM pairs calculation based on interactions among a tower-plate neutral "
                 "territory decomposition");
  t_nml->addKeyword("-trim_last", NamelistType::BOOLEAN);
  t_nml->addHelp("-trim_last", "Trim the final batch of atoms in any sector of the tower-plate "
                 "decomposition, conserving store operations while reducing total tile area.");
  t_nml->addKeyword("-cull_distal", NamelistType::BOOLEAN);
  t_nml->addHelp("-cull_distal", "Trim the distal atoms in the tower and plate which will never "
                 "be within range of any they would need to interact with.");
  t_nml->addKeyword("-hs_margin", NamelistType::REAL, std::to_string(1.5));
  t_nml->addHelp("-hs_margin", "The width of the margin to include in neighbor list building "
                 "based on a Hilbert space-fillin curve, in units of Angstroms.  This margin will "
                 "be added to whatever applicable cutoff.  The Hilbert space-filling curve "
                 "implementation is for quantification of the tile enrichment only and does not "
                 "guide pair interaction computations in practice.");
  t_nml->addKeyword("-hs_group", NamelistType::INTEGER, std::to_string(16));
  t_nml->addHelp("-hs_group", "The number of 'receiving' atoms to include in each batch when "
                 "building a hypothetical neighbor list based on a Hilbert space-filling curve.");
  t_nml->addKeyword("-hs_grid", NamelistType::INTEGER, std::to_string(256));
  t_nml->addHelp("-hs_grid", "The level of detail in the Hilbert space-filling curve which will "
                 "be stretched across each simulation unit cell in order to build a hypothetical "
                 "neighbor list.");
  
  // Initialize the testing environment such that it cooperates with this program's own
  // CommandLineParser to read user input.
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Read command line instructions.
  clip.parseUserInput(argc, argv);
  const std::string inpcrd_file = t_nml->getStringValue("-c");
  const std::string topology_file = t_nml->getStringValue("-p");
  bool use_dual_cg = t_nml->getBoolValue("-dual_cg");
  const bool skip_cpu_check = t_nml->getBoolValue("-skip_cpu_check");
  const bool est_occupancy = t_nml->getBoolValue("-occupancy");
  const int n_trials = t_nml->getIntValue("-trials");
  const int n_repeats = t_nml->getIntValue("-iter");
  const int n_replicas = t_nml->getIntValue("-replicas");
  const int fp_bits = t_nml->getIntValue("-fp_bits");
  const int warp_mult = t_nml->getIntValue("-warp_mult");
  const int ig_seed = t_nml->getIntValue("-ig_seed");
  bool have_user_cutoff = false;
  double elec_cutoff, vdw_cutoff;
  if (t_nml->getKeywordStatus("-cutoff") == InputStatus::USER_SPECIFIED) {
    vdw_cutoff = t_nml->getRealValue("-cutoff");
    elec_cutoff = vdw_cutoff;
    have_user_cutoff = true;
  }
  if (t_nml->getKeywordStatus("-elec_cutoff") == InputStatus::USER_SPECIFIED) {
    elec_cutoff = t_nml->getRealValue("-elec_cutoff");
    have_user_cutoff = true;
  }
  if (t_nml->getKeywordStatus("-vdw_cutoff") == InputStatus::USER_SPECIFIED) {
    vdw_cutoff = t_nml->getRealValue("-vdw_cutoff");
    have_user_cutoff = true;
  }
  if (have_user_cutoff == false) {
    vdw_cutoff = t_nml->getRealValue("-cutoff");
    elec_cutoff = vdw_cutoff;
  }
  double cutoff_pad = t_nml->getRealValue("-pad");
  EvaluateEnergy eval_nrg = (t_nml->getBoolValue("-eval_nrg")) ? EvaluateEnergy::YES :
                                                                 EvaluateEnergy::NO;
  EvaluateForce eval_frc;
  if (t_nml->getBoolValue("-omit_frc")) {

    // Force the evaluation of energy if forces are to be omitted.
    eval_nrg = EvaluateEnergy::YES;
    eval_frc = EvaluateForce::NO;
  }
  else {
    eval_frc = EvaluateForce::YES;
  }
  const double nl_margin = t_nml->getRealValue("-hs_margin");
  const int hs_group_size = t_nml->getIntValue("-hs_group");
  const int hs_grid_size = t_nml->getIntValue("-hs_grid");
  const bool trim_final_batch = t_nml->getBoolValue("-trim_last");
  const bool cull_irrelevant = t_nml->getBoolValue("-cull_distal");
  const int min_batch = t_nml->getIntValue("-min_batch");
  const int max_batch = t_nml->getIntValue("-max_batch");

  // Input checks
  if (n_replicas <= 0) {
    rtErr("A replica count of " + std::to_string(n_replicas) + " is invalid.\n",
          "pair_interactions");
  }
  if (nl_margin > std::min(elec_cutoff, vdw_cutoff) * 0.5) {
    rtErr("A neighbor list margin of " + minimalRealFormat(nl_margin, 1.0e-4, true) +
          " is invalid for cutoffs of " + minimalRealFormat(elec_cutoff, 1.0e-4, true) +
          " (electrostatic) and " + minimalRealFormat(vdw_cutoff, 1.0e-4, true) +
          " (van-der Waals).", "pair_interactions");
  }
  
  // A Hybrid object was created to engage the GPU.  Absorb any bootup time into "miscellaneous."
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
  }
  
  // Check that a system has been provided.  Load the basic topology and input coordinates.
#ifdef STORMM_USE_HPC
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), ExceptionResponse::DIE);
  timer.assignTime(time_load);
  
  // Stage the input parameters: cutoffs for each interaction.
  MolecularMechanicsControls mmctrl(0.01, 1, 1, warp_mult, elec_cutoff, vdw_cutoff);

  // Create the basic objects
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(n_replicas, 0),
                                                              0.0, ig_seed, 40, 44, fp_bits);
  PsSynthesisWriter host_psw = poly_ps.data();
  Xoshiro256ppGenerator xrs(ig_seed);
  for (int i = 1; i < host_psw.system_count; i++) {
    const size_t aoff = host_psw.atom_starts[i];
    addRandomNoise(&xrs, &host_psw.xcrd[aoff], &host_psw.xcrd_ovrf[aoff], &host_psw.ycrd[aoff],
                   &host_psw.ycrd_ovrf[aoff], &host_psw.zcrd[aoff], &host_psw.zcrd_ovrf[aoff],
                   host_psw.atom_counts[i], 0.01, host_psw.gpos_scale);
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(n_replicas, 0));
  ScoreCard sc(poly_ps.getSystemCount(), 1, 32);
  use_dual_cg = (use_dual_cg || (fabs(elec_cutoff - vdw_cutoff) > 1.0e-6));
  if (use_dual_cg == false) {
    elec_cutoff = vdw_cutoff;
  }
  const NeighborListKind layout = (use_dual_cg) ? NeighborListKind::DUAL : NeighborListKind::MONO;
  const CoreKlManager launcher(gpu, poly_ag);
  PPITable direct_space_table(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                              TableIndexing::SQUARED_ARG, elec_cutoff);
  const double ew_coeff = direct_space_table.getEwaldCoefficient();
  timer.assignTime(time_build);
  
  // Upload the basic objects
  poly_ps.upload();
  poly_ag.upload();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  SyNonbondedKit<float, float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  MMControlKit<float> ctrl = mmctrl.spData(devc);
  PPIKit<float, float4> nrg_tab = direct_space_table.spData(devc);
  ScoreCardWriter scw = sc.data(devc);
  timer.assignTime(time_build);

  // Create the neighbor list cell grids and run calculations
  if (use_dual_cg) {
    CellGrid<float, int, float, float4> cg_qq(&poly_ps, &poly_ag, elec_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::ELECTROSTATIC);
    CellGrid<float, int, float, float4> cg_lj(&poly_ps, &poly_ag, vdw_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::VAN_DER_WAALS);
    const TinyBoxPresence has_tiny_box = (cg_qq.getTinyBoxPresence() == TinyBoxPresence::YES ||
                                          cg_lj.getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                         TinyBoxPresence::YES : TinyBoxPresence::NO;
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL, has_tiny_box,
                                                      eval_frc, eval_nrg, ClashResponse::NONE);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout, has_tiny_box,
                                 poly_ag);
    mmctrl.upload();
    cg_qq.upload();
    cg_lj.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cg_qqw = cg_qq.data(devc);
    CellGridWriter<float, int, float, float4> cg_ljw = cg_lj.data(devc);
    const PsSynthesisBorders sysbrd = cg_qq.getUnitCellTransforms(devc);
    timer.assignTime(time_build);
    
    // Run the initialization of forces multiple times to gain a sense of the overhead paid in the
    // test.
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
        switch (has_tiny_box) {
        case TinyBoxPresence::NO:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, &cg_qqw, &cg_ljw, &tlpn, &scw, &ctrl, eval_frc,
                         eval_nrg, tp_bt, 0.0, 0.0);
          break;
        case TinyBoxPresence::YES:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cg_qqw, &cg_ljw, &tlpn, &scw, &ctrl,
                         eval_frc, eval_nrg, tp_bt, 0.0, 0.0);
          break;
        }          
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg_qq.contributeForces(devc, gpu);
    cg_lj.contributeForces(devc, gpu);

    // Estimate the thread occupancy of the GPU during the calculation
    if (est_occupancy) {
      reportTowerPlateOccupancy(cg_qq, elec_cutoff, min_batch, max_batch, trim_final_batch,
                                cull_irrelevant);
      reportTowerPlateOccupancy(cg_lj, vdw_cutoff, min_batch, max_batch, trim_final_batch,
                                cull_irrelevant);
      reportHilbertSpaceOccupancy(cg_qq, elec_cutoff, nl_margin, hs_group_size, hs_grid_size,
                                  hs_grid_size, hs_grid_size);
      reportHilbertSpaceOccupancy(cg_lj, vdw_cutoff, nl_margin, hs_group_size, hs_grid_size,
                                  hs_grid_size, hs_grid_size);
    }
  }
  else {
    CellGrid<float, int, float, float4> cg(&poly_ps, &poly_ag, vdw_cutoff, cutoff_pad, 4,
                                           NonbondedTheme::ALL);
    const TinyBoxPresence has_tiny_box = cg.getTinyBoxPresence();
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL,
                                                      cg.getTinyBoxPresence(), eval_frc, eval_nrg,
                                                      ClashResponse::NONE);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout,
                                 cg.getTinyBoxPresence(), poly_ag);
    mmctrl.upload();
    cg.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cgw = cg.data(devc);
    const PsSynthesisBorders sysbrd = cg.getUnitCellTransforms(devc);
    timer.assignTime(time_build);

    // Test the basic force initialization
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
        switch (cg.getTinyBoxPresence()) {
        case TinyBoxPresence::NO:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                         tp_bt, 0.0, 0.0);
          break;
        case TinyBoxPresence::YES:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                         eval_nrg, tp_bt, 0.0, 0.0);
          break;
        }
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg.contributeForces(devc, gpu);

    // Estimate the thread occupancy of the GPU during the calculation
    if (est_occupancy) {
      reportTowerPlateOccupancy(cg, vdw_cutoff, min_batch, max_batch, trim_final_batch,
                                cull_irrelevant);
      reportHilbertSpaceOccupancy(cg, vdw_cutoff, nl_margin, hs_group_size, hs_grid_size,
                                  hs_grid_size, hs_grid_size);
    }
  }

  // Run a separate check on the validity of the forces
  if (skip_cpu_check == false) {
    for (int i = 0; i < poly_ps.getSystemCount(); i++) {
      PhaseSpace ps_i(poly_ps.getAtomCount(i), poly_ps.getUnitCellType(),
                      HybridFormat::HOST_ONLY);
      coordCopy(&ps_i, poly_ps, i);
      ps_i.initializeForces();
      const AtomGraph *ag_i = poly_ps.getSystemTopologyPointer(i);
      const LocalExclusionMask lem_i(ag_i);
      evaluateParticleParticleEnergy(&ps_i, ag_i, lem_i, PrecisionModel::SINGLE, elec_cutoff,
                                     vdw_cutoff, ew_coeff);
      coordCopy(&poly_ps, i, ps_i);
    }
    checkPairBenchmarkForces(poly_ps, eval_frc, &timer, time_test);
  }
#else // STORMM_USE_HPC
  rtWarn("This benchmarking program requires GPU support to run.", "pair_interactions");
#endif // STORMM_USE_HPC

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return 0;
}
