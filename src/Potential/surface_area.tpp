// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc3>
void surfaceArea(PsSynthesisWriter *poly_psw, ScoreCard *sc,
                 const SyNonbondedKit<Tcalc, Tcalc2> &synbk, const uint* sasa_mask,
                 const std::vector<Tcalc3> &sphere_pts, const double probe_radius,
                 const double weight, const EvaluateForce eval_frc,
                 const SasaReference radii_source) {
  std::vector<int> sh_near_neighbors(sasa_neighbor_list_buffer_size);
  for (int i = 0; i < poly_psw->system_count; i++) {

    // For the GPU kernel, each block will take a segment of 32 atoms, then find all atoms which
    // are close enough to one of them that they might affect the exposed surface area.  A buffer
    // of atom indices kept by the kernel (up to 8192) will be held in __shared__ memory in order
    // to avoid allocating arrays in a separate object stored in __global__ memory.  In the
    // extreme circumstance that 32 atoms have more than 8192 near neighbors, the entirety of
    // other atoms in the molecular system will be scanned.
    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const std::vector<bool> &sasa_mask,
                   const std::vector<Tcalc3> &sphere_pts, const Tcalc* atomic_radii,
                   const double probe_radius, const double weight, const EvaluateForce eval_frc,
                   const int* radii_idx) {
  const double3 grid_origin = { minValue<double>(psw->xcrd, psw->natom),
                                minValue<double>(psw->ycrd, psw->natom),
                                minValue<double>(psw->zcrd, psw->natom) };
  const double3 grid_lengths = { maxValue<double>(psw->xcrd, psw->natom) - grid_origin.x,
                                 maxValue<double>(psw->ycrd, psw->natom) - grid_origin.y,
                                 maxValue<double>(psw->zcrd, psw->natom) - grid_origin.z };
  double max_rad = 0.0;
  if (radii_idx == nullptr) {
    for (int i = 0; i < psw->natom; i++) {
      max_rad = std::max(max_rad, static_cast<double>(atomic_radii[i]));
    }
  }
  else {
    for (int i = 0; i < psw->natom; i++) {
      max_rad = std::max(max_rad, static_cast<double>(atomic_radii[radii_idx[i]]));
    }
  }
  const double cell_min_width = 2.0 * (max_rad + probe_radius);
  const int3 cell_count = { static_cast<int>(ceil(grid_lengths.x / cell_min_width)),
                            static_cast<int>(ceil(grid_lengths.y / cell_min_width)),
                            static_cast<int>(ceil(grid_lengths.z / cell_min_width)) };
  const double3 cell_dims = { grid_lengths.x / static_cast<double>(cell_count.x),
                              grid_lengths.y / static_cast<double>(cell_count.y),
                              grid_lengths.z / static_cast<double>(cell_count.z) };
  int nmask_atom = 0;
  for (int i = 0; i < psw->natom; i++) {
    nmask_atom += sasa_mask[i];
  }
  std::vector<int> cell_homes(nmask_atom), cell_contents(nmask_atom);
  std::vector<int> cell_content_bounds((cell_count.x * cell_count.y * cell_count.z) + 1);
  for (int i = 0; i < psw->natom; i++) {
    if (sasa_mask[i]) {
      int cx = (psw->xcrd[i] - grid_origin.x) / cell_dims.x;
      if (cx < 0) {
        cx = 0;
      }
      else if (cx >= cell_count.x) {
        cx = cell_count.x - 1;
      }
      int cy = (psw->ycrd[i] - grid_origin.y) / cell_dims.y;
      if (cy < 0) {
        cy = 0;
      }
      else if (cy >= cell_count.y) {
        cy = cell_count.y - 1;
      }
      int cz = (psw->zcrd[i] - grid_origin.z) / cell_dims.z;
      if (cz < 0) {
        cz = 0;
      }
      else if (cz >= cell_count.z) {
        cz = cell_count.z - 1;
      }
      cell_homes[i] = (((cz * cell_count.y) + cy) * cell_count.x) + cx;
    }
  }
  indexingArray(cell_homes, &cell_contents, &cell_content_bounds);
  const size_t nsph_pts = sphere_pts.size();
  const Tcalc value_patch = 4.0 * pi / static_cast<Tcalc>(nsph_pts);
  double result = 0.0;
  for (int i = 0; i < cell_count.x; i++) {
    for (int j = 0; j < cell_count.y; j++) {
      for (int k = 0; k < cell_count.z; k++) {
        const int cell_ijk = (((k * cell_count.y) + j) * cell_count.x) + i;
        for (int m = cell_content_bounds[cell_ijk]; m < cell_content_bounds[cell_ijk + 1]; m++) {
          const int m_atom_idx = cell_contents[m];
          const Tcalc atom_x = psw->xcrd[m_atom_idx];
          const Tcalc atom_y = psw->ycrd[m_atom_idx];
          const Tcalc atom_z = psw->zcrd[m_atom_idx];
          const Tcalc r_eff = (radii_idx == nullptr) ?
                              atomic_radii[m_atom_idx] + probe_radius :
                              atomic_radii[radii_idx[m_atom_idx]] + probe_radius;
          const Tcalc pt_area = value_patch * r_eff * r_eff;
          for (int ipt = 0; ipt < nsph_pts; ipt++) {
            const Tcalc pt_x = atom_x + (sphere_pts[ipt].x * r_eff);
            const Tcalc pt_y = atom_y + (sphere_pts[ipt].y * r_eff);
            const Tcalc pt_z = atom_z + (sphere_pts[ipt].z * r_eff);

            // Another triple-nested loop over all neighboring cells' contents
            const int k_llim = std::max(k - 1, 0);
            const int k_hlim = std::min(k + 1, cell_count.z);
            bool survives = true;
            for (int kv = k_llim; kv < k_hlim; kv++) {
              const int j_llim = std::max(j - 1, 0);
              const int j_hlim = std::min(j + 1, cell_count.y); 
              for (int jv = j_llim; jv < j_hlim; jv++) {
                const int i_llim = std::max(i - 1, 0);
                const int i_hlim = std::min(i + 1, cell_count.x);
                for (int iv = i_llim; iv < i_hlim; iv++) {
                  const int cell_ijkv = (((kv * cell_count.y) + jv) * cell_count.x) + iv;
                  for (int mv = cell_content_bounds[cell_ijkv];
                       mv < cell_content_bounds[cell_ijkv + 1]; mv++) {
                    const int mv_atom_idx = cell_contents[mv];
                    if (mv_atom_idx == m_atom_idx) {
                      continue;
                    }
                    if (mv_atom_idx == m_atom_idx) {
                      continue;
                    }
                    const Tcalc dx = pt_x - psw->xcrd[mv_atom_idx];
                    const Tcalc dy = pt_y - psw->ycrd[mv_atom_idx];
                    const Tcalc dz = pt_z - psw->zcrd[mv_atom_idx];
                    const Tcalc rv_eff = (radii_idx == nullptr) ?
                                         atomic_radii[mv_atom_idx] + probe_radius :
                                         atomic_radii[radii_idx[mv_atom_idx]] + probe_radius;
                    if ((dx * dx) + (dy * dy) + (dz * dz) < rv_eff * rv_eff) {
                      survives = false;

                      // Bail out of all nested loops
                      mv = cell_content_bounds[cell_ijkv + 1];
                      iv = i_hlim;
                      jv = j_hlim;
                      kv = k_hlim;
                    }
                  }
                }
              }
            }
            if (survives) {
              result += pt_area;
            }
          }
        }
      }
    }
  }
  return result * weight;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const std::vector<bool> &sasa_mask,
                   const NonbondedKit<Tcalc> &nbk, const std::vector<Tcalc3> &sphere_pts,
                   const double probe_radius, const double weight, const EvaluateForce eval_frc) {
  std::vector<Tcalc> atomic_radii(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    atomic_radii[i] = 0.5 * nbk.lj_sigma[i];
  }
  return surfaceArea(psw, sasa_mask, sphere_pts, atomic_radii.data(), probe_radius, weight,
                     eval_frc, nbk.lj_idx);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const std::vector<bool> &sasa_mask,
                   const ImplicitSolventKit<Tcalc> &isk, const std::vector<Tcalc3> &sphere_pts,
                   const double probe_radius, const double weight, const EvaluateForce eval_frc) {
  return surfaceArea(psw, sasa_mask, sphere_pts, isk.pb_radii, probe_radius, weight, eval_frc);
}

} // namespace energy
} // namespace stormm
