#include <algorithm>
#include "../../../src/copyright.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Chemistry/atommask.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/Math/summation.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "mmgbsa_carveout.h"

namespace mmgbsa {

using stormm::constants::CaseSensitivity;
using stormm::chemistry::AtomMask;
using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::IsomerPlan;
using stormm::parse::strcmpCased;
using stormm::stmath::indexingArray;
using stormm::stmath::sum;
using stormm::trajectory::PhaseSpaceReader;
using stormm::topology::ChemicalDetailsKit;
  
//-------------------------------------------------------------------------------------------------
std::vector<int> receptorCarveOut(const MMGBSAControls &gbsacon, const SystemCache &sc) {
  std::vector<int> result;

  // Unpack the systems cache based on ligands and the receptor.
  const std::vector<int> rec_cache_idx = sc.getMatchingSystemIndices("receptor");
  const std::vector<int> lig_cache_idx = sc.getMatchingSystemIndices("ligand");
  const int nrcpt = rec_cache_idx.size();
  const int nlgnd = lig_cache_idx.size();
  std::vector<PhaseSpaceReader> receptor_psr;
  std::vector<PhaseSpaceReader> ligand_psr;
  receptor_psr.reserve(nrcpt);
  ligand_psr.reserve(nlgnd);
  for (int i = 0; i < nrcpt; i++) {
    receptor_psr.push_back(sc.getCoordinatePointer(rec_cache_idx[i])->data());
  }
  std::vector<const AtomGraph*> ligand_agv;
  ligand_agv.reserve(nlgnd);
  for (int i = 0; i < nlgnd; i++) {
    ligand_psr.push_back(sc.getCoordinatePointer(lig_cache_idx[i])->data());
    ligand_agv.push_back(sc.getSystemTopologyPointer(lig_cache_idx[i]));
  }
  const AtomGraph* receptor_ag = sc.getSystemTopologyPointer(rec_cache_idx[0]);
  std::vector<bool> active_atom(receptor_ag->getAtomCount(), false);
  if (gbsacon.proximityCarveout()) {

    // Define a grid for all receptor coordinate sets to make neighbor lists for each receptor
    // state.
    const double mask_cutoff = gbsacon.getCarveoutCutoff();
    const double sq_cutoff = mask_cutoff * mask_cutoff;
    double g_xmin, g_ymin, g_zmin, g_xmax, g_ymax, g_zmax;
    for (int i = 0; i < nrcpt; i++) {
      for (int j = 0; j < receptor_psr[i].natom; j++) {
        if (i == 0 && j == 0) {
          g_xmin = receptor_psr[i].xcrd[j];
          g_ymin = receptor_psr[i].ycrd[j];
          g_zmin = receptor_psr[i].zcrd[j];
          g_xmax = receptor_psr[i].xcrd[j];
          g_ymax = receptor_psr[i].ycrd[j];
          g_zmax = receptor_psr[i].zcrd[j];
        }
        else {
          g_xmin = std::min(g_xmin, receptor_psr[i].xcrd[j]);
          g_ymin = std::min(g_ymin, receptor_psr[i].ycrd[j]);
          g_zmin = std::min(g_zmin, receptor_psr[i].zcrd[j]);
          g_xmax = std::max(g_xmax, receptor_psr[i].xcrd[j]);
          g_ymax = std::max(g_ymax, receptor_psr[i].ycrd[j]);
          g_zmax = std::max(g_zmax, receptor_psr[i].zcrd[j]);
        }
      }
    }
    const int g_nx = ceil((g_xmax - g_xmin) / mask_cutoff);
    const int g_ny = ceil((g_ymax - g_ymin) / mask_cutoff);
    const int g_nz = ceil((g_zmax - g_zmin) / mask_cutoff);
    const double gx_cell_len = (g_xmax - g_xmin) / static_cast<double>(g_nx);
    const double gy_cell_len = (g_ymax - g_ymin) / static_cast<double>(g_ny);
    const double gz_cell_len = (g_zmax - g_zmin) / static_cast<double>(g_nz);
    const int g_nxyz = g_nx * g_ny * g_nz;
    std::vector<std::vector<int>> ngbr_asgn(nrcpt), ngbr_bounds(nrcpt), ngbr_list(nrcpt);
    for (int i = 0; i < nrcpt; i++) {
      ngbr_asgn[i].resize(receptor_psr[i].natom);
      ngbr_list[i].resize(receptor_psr[i].natom);
      ngbr_bounds[i].resize(g_nxyz + 1);
      for (int j = 0; j < receptor_psr[i].natom; j++) {
        int blk_x = (receptor_psr[i].xcrd[j] - g_xmin) / gx_cell_len;
        int blk_y = (receptor_psr[i].ycrd[j] - g_ymin) / gy_cell_len;
        int blk_z = (receptor_psr[i].zcrd[j] - g_zmin) / gz_cell_len;
        blk_x = std::max(blk_x, 0);
        blk_y = std::max(blk_y, 0);
        blk_z = std::max(blk_z, 0);
        blk_x = std::min(blk_x, g_nx - 1);
        blk_y = std::min(blk_y, g_ny - 1);
        blk_z = std::min(blk_z, g_nz - 1);
        ngbr_asgn[i][j] = (((blk_z * g_ny) + blk_y) * g_nx) + blk_x;
      }
    }
    for (int i = 0; i < nrcpt; i++) {
      indexingArray(ngbr_asgn[i], &ngbr_list[i], &ngbr_bounds[i]);
    }
    
    // Create a list of ligand coordinate abstracts and begin checking off receptor atoms that
    // may lie within range.
    const bool use_all_ligands = strcmpCased(gbsacon.getCarveoutLigandReference(), "all",
                                             CaseSensitivity::NO);
    for (int lg_idx = 0; lg_idx < nlgnd; lg_idx++) {
      const std::string& fname = ligand_agv[lg_idx]->getFileName();
      if (use_all_ligands || fname.find(gbsacon.getCarveoutLigandReference()) < fname.size()) {
        for (int i = 0; i < ligand_psr[lg_idx].natom; i++) {

          // Determine the hypothetical neighbor list cell, absent periodic boundary conditions,
          // and then loop over all atoms within that cell or neighboring cells.  If the ligand
          // atom exists in a cell that is not an actual part of the neighbor list, then there
          // will be no relevant receptor atoms in that cell but there may be receptor atoms in
          // adjacent cells.
          const double atmx = ligand_psr[lg_idx].xcrd[i];
          const double atmy = ligand_psr[lg_idx].ycrd[i];
          const double atmz = ligand_psr[lg_idx].zcrd[i];
          const int blk_min_x = floor((atmx - g_xmin - mask_cutoff) / gx_cell_len);
          const int blk_min_y = floor((atmy - g_ymin - mask_cutoff) / gy_cell_len);
          const int blk_min_z = floor((atmz - g_zmin - mask_cutoff) / gz_cell_len);
          const int blk_max_x = ceil((atmx - g_xmin + mask_cutoff) / gx_cell_len);
          const int blk_max_y = ceil((atmy - g_ymin + mask_cutoff) / gy_cell_len);
          const int blk_max_z = ceil((atmz - g_zmin + mask_cutoff) / gz_cell_len);
          const int search_min_x = std::max(0, blk_min_x - 1);
          const int search_min_y = std::max(0, blk_min_y - 1);
          const int search_min_z = std::max(0, blk_min_z - 1);
          const int search_max_x = std::min(g_nx, blk_max_x + 1);
          const int search_max_y = std::min(g_ny, blk_max_y + 1);
          const int search_max_z = std::min(g_nz, blk_max_z + 1);
          for (int gi = search_min_x; gi < search_max_x; gi++) {
            for (int gj = search_min_y; gj < search_max_y; gj++) {
              for (int gk = search_min_z; gk < search_max_z; gk++) {
                const int gijk = (((gk * g_ny) + gj) * g_nx) + gi;
                for (int rc_idx = 0; rc_idx < nrcpt; rc_idx++) {
                  for (int j = ngbr_bounds[rc_idx][gijk]; j < ngbr_bounds[rc_idx][gijk + 1]; j++) {
                    const int jatom = ngbr_list[rc_idx][j];
                    if (! active_atom[jatom]) {
                      const double dx = receptor_psr[rc_idx].xcrd[jatom] - atmx;
                      const double dy = receptor_psr[rc_idx].ycrd[jatom] - atmy;
                      const double dz = receptor_psr[rc_idx].zcrd[jatom] - atmz;
                      const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
                      if (r2 < sq_cutoff) {
                        active_atom[jatom] = true;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  // Expand the list of relevant atoms by looping over rotatable groups in the receptor, then
  // selecting the minimal rotatable groups that subsume all atoms in the carveout.  The carveout
  // expands to comprise all atoms in these rotatable groups.
  const ChemicalFeatures& rcpt_cf = sc.getFeatures(rec_cache_idx[0]);
  const std::vector<IsomerPlan> rcpt_rg = rcpt_cf.getRotatableBondGroups();
  const int n_rg = rcpt_rg.size();
  for (int i = 0; i < n_rg; i++) {
    const int nrot = rcpt_rg[i].getMovingAtomCount();
    if (rcpt_rg[i].getRealMovingAtomCount() <= gbsacon.getRotatingGroupLimit()) {
      bool touched = false;
      for (int j = 0; j < nrot; j++) {
        touched = (touched || active_atom[rcpt_rg[i].getMovingAtom(j)]);
      }
      if (touched) {
        for (int j = 0; j < nrot; j++) {
          active_atom[rcpt_rg[i].getMovingAtom(j)] = true;
        }
      }
    }
  }

  // Expand the list of relevant atoms by looping over affected residues.
  if (gbsacon.completeResidueCarveout()) {
    const ChemicalDetailsKit cdk = receptor_ag->getChemicalDetailsKit();
    for (int i = 0; i < cdk.natom; i++) {
      if (active_atom[i]) {
        const int res_idx = receptor_ag->getResidueIndex(i);
        for (int j = cdk.res_limits[res_idx]; j < cdk.res_limits[res_idx + 1]; j++) {
          active_atom[j] = true;
        }

        // Advance to the end of this residue, which is expected to consist of a contiguous set
        // of atoms that are now all active.
        i = cdk.res_limits[res_idx + 1] - 1;
      }
    }
  }

  // Expand the list of relevant atoms using a user-specified atom mask.
  const int natom = receptor_ag->getAtomCount();
  if (gbsacon.getCarveoutExtraMask() != std::string(default_carveout_mask)) {
    for (int i = 0; i < nrcpt; i++) {
      const AtomMask xmask(gbsacon.getCarveoutExtraMask(), receptor_ag, rcpt_cf,
                           sc.getCoordinates(rec_cache_idx[i]));
      for (int j = 0; j < natom; j++) {
        if (xmask.isAtomInMask(j)) {
          active_atom[j] = true;
        }
      }
    }
  }

  // Transcribe the mask of active atoms into a list of indices
  const int n_active = sum<int>(active_atom);
  result.reserve(n_active);
  for (int i = 0; i < natom; i++) {
    if (active_atom[i]) {
      result.push_back(i);
    }
  }
  return result;
}
 
} // namespace mmgbsa
