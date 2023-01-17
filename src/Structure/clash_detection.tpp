// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc maxClashingDistance(const NonbondedKit<Tcalc> &nbk, const Tcalc elec_limit,
                          const Tcalc vdw_ratio) {
  Tcalc result = elec_limit;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    for (int j = 0; j <= i; j++) {
      result = std::max(result, nbk.lj_sigma[(i * nbk.n_lj_types) + j]);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool trivialClashCheck(const std::vector<Tcalc> &cachi_xcrd, const std::vector<Tcalc> &cachj_xcrd,
                       const std::vector<Tcalc> &cachi_ycrd, const std::vector<Tcalc> &cachj_ycrd,
                       const std::vector<Tcalc> &cachi_zcrd, const std::vector<Tcalc> &cachj_zcrd,
                       const int ni_atoms, const int nj_atoms, const Tcalc max_clash) {
  Tcalc min_ix = cachi_xcrd[0];
  Tcalc min_iy = cachi_ycrd[0];
  Tcalc min_iz = cachi_zcrd[0];
  Tcalc max_ix = cachi_xcrd[0];
  Tcalc max_iy = cachi_ycrd[0];
  Tcalc max_iz = cachi_zcrd[0];
  for (int i = 1; i < ni_atoms; i++) {
    min_ix = std::min(cachi_xcrd[i], min_ix);
    min_iy = std::min(cachi_ycrd[i], min_iy);
    min_iz = std::min(cachi_zcrd[i], min_iz);
    max_ix = std::max(cachi_xcrd[i], max_ix);
    max_iy = std::max(cachi_ycrd[i], max_iy);
    max_iz = std::max(cachi_zcrd[i], max_iz);
  }
  Tcalc min_jx = cachj_xcrd[0];
  Tcalc min_jy = cachj_ycrd[0];
  Tcalc min_jz = cachj_zcrd[0];
  Tcalc max_jx = cachj_xcrd[0];
  Tcalc max_jy = cachj_ycrd[0];
  Tcalc max_jz = cachj_zcrd[0];
  for (int j = 1; j < nj_atoms; j++) {
    min_jx = std::min(cachj_xcrd[j], min_jx);
    min_jy = std::min(cachj_ycrd[j], min_jy);
    min_jz = std::min(cachj_zcrd[j], min_jz);
    max_jx = std::max(cachj_xcrd[j], max_jx);
    max_jy = std::max(cachj_ycrd[j], max_jy);
    max_jz = std::max(cachj_zcrd[j], max_jz);
  }
  min_ix -= max_clash;
  min_iy -= max_clash;
  min_iz -= max_clash;
  max_ix += max_clash;
  max_iy += max_clash;
  max_iz += max_clash;
  return (max_ix > min_jx && min_ix < max_jx && max_iy > min_jy && min_iy < max_jy &&
          max_iz > min_jz && min_iz < max_jz);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool directClashTesting(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &maskr, const Tcalc elec_limit,
                        const Tcalc vdw_ratio, const Tcalc inv_scale) {

  // Determine the calculation type and a maximum distance parameter
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc max_clash = maxClashingDistance<Tcalc>(nbk, elec_limit, vdw_ratio);

  // Lay out a workspace and begin searching
  bool clash_found = false;
  std::vector<Tcalc> cachi_xcrd(tile_length), cachi_ycrd(tile_length), cachi_zcrd(tile_length);
  std::vector<Tcalc> cachj_xcrd(tile_length), cachj_ycrd(tile_length), cachj_zcrd(tile_length);
  std::vector<int> cachi_ljidx(tile_length), cachj_ljidx(tile_length);
  for (int sti = 0; sti < maskr.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {
      const int stni_atoms = std::min(maskr.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(maskr.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;
      const int diag_supertile = (sti == stj);
      const int stij_map_index = maskr.supertile_map_idx[(stj * maskr.supertile_stride_count) +
                                                         sti];
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);
          const int tij_map_index = maskr.tile_map_idx[stij_map_index +
                                                       (tj * tile_lengths_per_supertile) + ti];
          if (isFloatingPointScalarType<Tcoord>()) {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              cachi_xcrd[i]  = xcrd[atom_i];
              cachi_ycrd[i]  = ycrd[atom_i];
              cachi_zcrd[i]  = zcrd[atom_i];
              cachi_ljidx[i] = nbk.lj_idx[atom_i];
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              cachj_xcrd[j]  = xcrd[atom_j];
              cachj_ycrd[j]  = ycrd[atom_j];
              cachj_zcrd[j]  = zcrd[atom_j];
              cachj_ljidx[j] = nbk.lj_idx[atom_j];
            }
          }
          else {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              cachi_xcrd[i]  = xcrd[atom_i] * inv_scale;
              cachi_ycrd[i]  = ycrd[atom_i] * inv_scale;
              cachi_zcrd[i]  = zcrd[atom_i] * inv_scale;
              cachi_ljidx[i] = nbk.lj_idx[atom_i];
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              cachj_xcrd[j]  = xcrd[atom_j] * inv_scale;
              cachj_ycrd[j]  = ycrd[atom_j] * inv_scale;
              cachj_zcrd[j]  = zcrd[atom_j] * inv_scale;
              cachj_ljidx[j] = nbk.lj_idx[atom_j];
            }
          }

          // Perform a trivial clash check
          if (diag_tile || (diag_supertile && ti - tj < 3) ||
              trivialClashCheck<Tcalc>(cachi_xcrd, cachj_xcrd, cachi_ycrd, cachj_ycrd, cachi_zcrd,
                                       cachj_zcrd, ni_atoms, nj_atoms, max_clash)) {
            for (int i = 0; i < ni_atoms; i++) {
              const int ljt_i = cachi_ljidx[i];
              const uint mask_i = maskr.mask_data[tij_map_index + i];
              const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
              if (tcalc_is_double) {
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const int ljt_ij = (nbk.n_lj_types * ljt_i) + cachj_ljidx[j];
                  const Tcalc dx = cachj_xcrd[j] - cachi_xcrd[i];
                  const Tcalc dy = cachj_ycrd[j] - cachi_ycrd[i];
                  const Tcalc dz = cachj_zcrd[j] - cachi_zcrd[i];
                  const Tcalc r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                  if (r < nbk.lj_sigma[ljt_ij] * vdw_ratio || r < elec_limit) {
                    clash_found = true;
                  }
                }
              }
              else {
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const int ljt_ij = (nbk.n_lj_types * ljt_i) + cachj_ljidx[j];
                  const Tcalc dx = cachj_xcrd[j] - cachi_xcrd[i];
                  const Tcalc dy = cachj_ycrd[j] - cachi_ycrd[i];
                  const Tcalc dz = cachj_zcrd[j] - cachi_zcrd[i];
                  const Tcalc r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                  if (r < nbk.lj_sigma[ljt_ij] * vdw_ratio || r < elec_limit) {
                    clash_found = true;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return clash_found;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMaskReader &maskr, const Tcalc elec_limit,
                 const Tcalc vdw_ratio, const Tcalc inv_scale) {

  // Determine the calculation type
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  bool clash_found = false;
  
  // Do the direct calculation if the size is very small
  if (nbk.natom < clash_direct_calculation_size_limit) {
    clash_found = directClashTesting(xcrd, ycrd, zcrd, vk, nbk, maskr, vdw_ratio, elec_limit,
                                     inv_scale);
  }
  else {

    // Find the largest clashing distance.
    const Tcalc max_clash = maxClashingDistance<Tcalc>(nbk, elec_limit, vdw_ratio);
    if (max_clash >= constants::tiny) {

      // Design a grid based on the maximum clash distance.  In order to be worthwhile, the grid
      // cells must be at least max_clash thick between all parallel faces, and more than five
      // grid cells in all directions
    }
  }

  // Loop over all 1:4 interactions to clean up possible clashes among these pairs.
  for (int pos = 0; pos < vk.ndihe; pos++) {
    if (vk.dihe14_param_idx[pos] == 0) {
      continue;
    }
    const int atom_i = vk.dihe_i_atoms[pos];
    const int atom_l = vk.dihe_l_atoms[pos];
    const int ljt_il = (nbk.n_lj_types * nbk.lj_idx[atom_i]) + nbk.lj_idx[atom_l];
    const Tcalc dx = xcrd[atom_l] - xcrd[atom_i];
    const Tcalc dy = ycrd[atom_l] - ycrd[atom_i];
    const Tcalc dz = zcrd[atom_l] - zcrd[atom_i];
    const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                        sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    if (r < nbk.lj_14_sigma[ljt_il] * vdw_ratio) {
      return true;
    }
  }
  for (int pos = 0; pos < vk.ninfr14; pos++) {
    if (vk.infr14_param_idx[pos] == 0) {
      continue;
    }
    const int atom_i = vk.infr14_i_atoms[pos];
    const int atom_l = vk.infr14_l_atoms[pos];
    const int ljt_il = (nbk.n_lj_types * nbk.lj_idx[atom_i]) + nbk.lj_idx[atom_l];
    const Tcalc dx = xcrd[atom_l] - xcrd[atom_i];
    const Tcalc dy = ycrd[atom_l] - ycrd[atom_i];
    const Tcalc dz = zcrd[atom_l] - zcrd[atom_i];
    const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                        sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    if (r < nbk.lj_14_sigma[ljt_il] * vdw_ratio) {
      return true;
    }
  }
  
  // If this point has been reached, no clash was detected.
  return false;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame,
                            const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                            const StaticExclusionMaskReader &maskr, const double elec_limit,
                            const Tcalc vdw_ratio) {
  const size_t padded_atoms = roundUp(csr.natom, warp_size_int);
  const size_t atom_start = frame * padded_atoms;
  const size_t xfrm_start = frame * roundUp<size_t>(9, warp_size_zu);
  return detectVanDerWaalsClash(&csr.xcrd[atom_start], &csr.ycrd[atom_start],
                                &csr.zcrd[atom_start], vk, nbk, maskr, elec_limit, vdw_ratio,
                                csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> *cs, const int frame,
                            const AtomGraph *ag, const StaticExclusionMask &mask,
                            const Tcalc elec_limit, const Tcalc vdw_ratio) {
  const size_t ct = std::type_index(typeid(Tcalc)).hash_code();
  if (ct == double_type_index) {
    return detectClash<Tcoord, double>(cs->data(), frame, ag->getDoublePrecisionValenceKit(),
                                       ag->getDoublePrecisionNonbondedKit(), mask.data(),
                                       elec_limit, vdw_ratio);
  }
  else if (ct == float_type_index) {
    return detectClash<Tcoord, float>(cs->data(), frame, ag->getSinglePrecisionValenceKit(),
                                      ag->getSinglePrecisionNonbondedKit(), mask.data(),
                                      elec_limit, vdw_ratio);
  }
  else {
    rtErr("The clash detection must be performed in double or float.  " +
          getStormmScalarTypeName<Tcalc> + " is invalid.", "detectVanDerWaalsClash");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
bool detectVanDerWaalsClash(const CoordinateSeries<Tcoord> &cs, const int frame,
                            const AtomGraph &ag, const StaticExclusionMask &mask,
                            const Tcalc elec_limit, const Tcalc vdw_ratio) {
  return detectVanDerWaalsClash<Tcoord, Tcalc>(cs.getSelfPointer(), frame, ag.getSelfPointer(),
                                               mask, elec_limit, vdw_ratio);
}

} // namespace structure
} // namespace stormm
