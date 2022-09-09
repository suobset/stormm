// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const Tcoord* xcrd_a, const Tcoord* ycrd_a, const Tcoord* zcrd_a, const Tcoord* xcrd_b,
           const Tcoord* ycrd_b, const Tcoord* zcrd_b, const Tcalc* masses,
           const RmsdMethod method, const int lower_limit, const int upper_limit,
           const Tcalc inv_gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  Tcalc result = 0.0;
  Tcalc inv_mass_divisor;
  switch (method) {
  case RmsdMethod::ALIGN_MASS:
  case RmsdMethod::NO_ALIGN_MASS:
    inv_mass_divisor = value_one / sum<Tcalc>(&masses[lower_limit], upper_limit - lower_limit);
    break;
  case RmsdMethod::ALIGN_GEOM:
  case RmsdMethod::NO_ALIGN_GEOM:
    inv_mass_divisor = value_one / static_cast<Tcalc>(upper_limit - lower_limit);
    break;
  }
  switch (method) {
  case RmsdMethod::ALIGN_MASS:
  case RmsdMethod::ALIGN_GEOM:
    {
      // In order to compute the aligned RMSD without moving the coordinates, the movement must
      // be computed in temporary variables only.
      Tcalc sab_xx = 0.0;
      Tcalc sab_xy = 0.0;
      Tcalc sab_xz = 0.0;
      Tcalc sab_yx = 0.0;
      Tcalc sab_yy = 0.0;
      Tcalc sab_yz = 0.0;
      Tcalc sab_zx = 0.0;
      Tcalc sab_zy = 0.0;
      Tcalc sab_zz = 0.0;
      Tcalc sa_x = 0.0;
      Tcalc sa_y = 0.0;
      Tcalc sa_z = 0.0;
      Tcalc sb_x = 0.0;
      Tcalc sb_y = 0.0;
      Tcalc sb_z = 0.0;
      const bool use_mass = (method == RmsdMethod::ALIGN_MASS);
      if (tcoord_is_sgnint) {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc xca = static_cast<Tcalc>(xcrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc yca = static_cast<Tcalc>(ycrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc zca = static_cast<Tcalc>(zcrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc xcb = static_cast<Tcalc>(xcrd_b[i]) * inv_gpos_scale_factor;
          const Tcalc ycb = static_cast<Tcalc>(ycrd_b[i]) * inv_gpos_scale_factor;
          const Tcalc zcb = static_cast<Tcalc>(zcrd_b[i]) * inv_gpos_scale_factor;
          sab_xx += tmass * xca * xcb;
          sab_xy += tmass * xca * ycb;
          sab_xz += tmass * xca * zcb;
          sab_yx += tmass * yca * xcb;
          sab_yy += tmass * yca * ycb;
          sab_yz += tmass * yca * zcb;
          sab_zx += tmass * zca * xcb;
          sab_zy += tmass * zca * ycb;
          sab_zz += tmass * zca * zcb;
          sa_x += tmass * xca;
          sa_y += tmass * yca;
          sa_z += tmass * zca;
          sb_x += tmass * xcb;
          sb_y += tmass * ycb;
          sb_z += tmass * zcb;
        }
      }
      else {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          sab_xx += tmass * xcrd_a[i] * xcrd_b[i];
          sab_xy += tmass * xcrd_a[i] * ycrd_b[i];
          sab_xz += tmass * xcrd_a[i] * zcrd_b[i];
          sab_yx += tmass * ycrd_a[i] * xcrd_b[i];
          sab_yy += tmass * ycrd_a[i] * ycrd_b[i];
          sab_yz += tmass * ycrd_a[i] * zcrd_b[i];
          sab_zx += tmass * zcrd_a[i] * xcrd_b[i];
          sab_zy += tmass * zcrd_a[i] * ycrd_b[i];
          sab_zz += tmass * zcrd_a[i] * zcrd_b[i];
          sa_x += tmass * xcrd_a[i];
          sa_y += tmass * ycrd_a[i];
          sa_z += tmass * zcrd_a[i];
          sb_x += tmass * xcrd_b[i];
          sb_y += tmass * ycrd_b[i];
          sb_z += tmass * zcrd_b[i];
        }
      }
      const Tcalc coma_x = sa_x * inv_mass_divisor;
      const Tcalc coma_y = sa_y * inv_mass_divisor;
      const Tcalc coma_z = sa_z * inv_mass_divisor;
      const Tcalc comb_x = sb_x * inv_mass_divisor;
      const Tcalc comb_y = sb_y * inv_mass_divisor;
      const Tcalc comb_z = sb_z * inv_mass_divisor;

      // Assemble the Kabsch matrix and diagonalize it.  The weighted product of the centers of
      // mass is equal to the weighted sum sa_(x,y,z) or sb_(x,y,z) times the center of mass
      // comb_(x,y,z) or coma_(x,y,z), so rather than evaluate the full F.O.I.L. to complete
      // mass * ((x,y,z)crd_a - coma_(x,y,z)) * ((x,y,z)crd_b - comb_(x,y,z)), the final two
      // terms cancel.
      const Tcalc aa = sab_xx - (sb_x * coma_x);
      const Tcalc ab = sab_xy - (sb_y * coma_x);
      const Tcalc ac = sab_xz - (sb_z * coma_x);
      const Tcalc ba = sab_yx - (sb_x * coma_y);
      const Tcalc bb = sab_yy - (sb_y * coma_y);
      const Tcalc bc = sab_yz - (sb_z * coma_y);
      const Tcalc ca = sab_zx - (sb_x * coma_z);
      const Tcalc cb = sab_zy - (sb_y * coma_z);
      const Tcalc cc = sab_zz - (sb_z * coma_z);
      std::vector<Tcalc> rmat(16);
      rmat[ 0] = aa + bb + cc;
      rmat[ 1] = cb - bc;
      rmat[ 2] = ac - ca;
      rmat[ 3] = ba - ab;
      rmat[ 5] = aa - bb - cc;
      rmat[ 6] = ab + ba;
      rmat[ 7] = ca + ac;
      rmat[10] = bb - aa - cc;
      rmat[11] = bc + cb;
      rmat[15] = cc - aa - bb;
      rmat[ 4] = rmat[ 1];
      rmat[ 8] = rmat[ 2];
      rmat[12] = rmat[ 3];
      rmat[ 9] = rmat[ 6];
      rmat[13] = rmat[ 7];
      rmat[14] = rmat[11];
      for (int i = 0; i < 16; i++) {
        rmat[i] *= inv_mass_divisor;
      }
      std::vector<Tcalc> vmat(16, 0.0), eigval(4, 0.0);
      jacobiEigensolver(&rmat, &vmat, &eigval, 4);
      int max_eig_loc = 0;
      for (int i = 1; i < 4; i++) {
        if (eigval[i] > eigval[max_eig_loc]) {
          max_eig_loc = i;
        }
      }

      // Form the rotation matrix
      const Tcalc a = vmat[(4 * max_eig_loc)    ];
      const Tcalc x = vmat[(4 * max_eig_loc) + 1];
      const Tcalc y = vmat[(4 * max_eig_loc) + 2];
      const Tcalc z = vmat[(4 * max_eig_loc) + 3];
      std::vector<Tcalc> umat(9);
      umat[0] = (a * a) + (x * x) - (y * y) - (z * z);
      umat[3] = 2.0 * ((x * y) + (a * z));
      umat[6] = 2.0 * ((x * z) - (a * y));
      umat[1] = 2.0 * ((x * y) - (a * z));
      umat[4] = (a * a) - (x * x) + (y * y) - (z * z);
      umat[7] = 2.0 * ((y * z) + (a * x));
      umat[2] = 2.0 * ((x * z) + (a * y));
      umat[5] = 2.0 * ((y * z) - (a * x));
      umat[8] = (a * a) - (x * x) - (y * y) + (z * z);

      // Shift and rotate the coordinates of the first frame (in temporary variables only) and
      // compare them to the shifted, unrotated coordinates of the second frame.
      if (tcoord_is_sgnint) {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc shfta_x = (static_cast<Tcalc>(xcrd_a[i]) * inv_gpos_scale_factor) - coma_x;
          const Tcalc shfta_y = (static_cast<Tcalc>(ycrd_a[i]) * inv_gpos_scale_factor) - coma_y;
          const Tcalc shfta_z = (static_cast<Tcalc>(zcrd_a[i]) * inv_gpos_scale_factor) - coma_z;
          const Tcalc shftb_x = (static_cast<Tcalc>(xcrd_b[i]) * inv_gpos_scale_factor) - comb_x;
          const Tcalc shftb_y = (static_cast<Tcalc>(ycrd_b[i]) * inv_gpos_scale_factor) - comb_y;
          const Tcalc shftb_z = (static_cast<Tcalc>(zcrd_b[i]) * inv_gpos_scale_factor) - comb_z;
          const Tcalc rota_x = (umat[0] * shfta_x) + (umat[3] * shfta_y) + (umat[6] * shfta_z);
          const Tcalc rota_y = (umat[1] * shfta_x) + (umat[4] * shfta_y) + (umat[7] * shfta_z);
          const Tcalc rota_z = (umat[2] * shfta_x) + (umat[5] * shfta_y) + (umat[8] * shfta_z);
          const Tcalc dx = shftb_x - rota_x;
          const Tcalc dy = shftb_y - rota_y;
          const Tcalc dz = shftb_z - rota_z;
          result += tmass * ((dx * dx) + (dy * dy) + (dz * dz));
        }
      }
      else {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc shfta_x = xcrd_a[i] - coma_x;
          const Tcalc shfta_y = ycrd_a[i] - coma_y;
          const Tcalc shfta_z = zcrd_a[i] - coma_z;
          const Tcalc shftb_x = xcrd_b[i] - comb_x;
          const Tcalc shftb_y = ycrd_b[i] - comb_y;
          const Tcalc shftb_z = zcrd_b[i] - comb_z;
          const Tcalc rota_x = (umat[0] * shfta_x) + (umat[3] * shfta_y) + (umat[6] * shfta_z);
          const Tcalc rota_y = (umat[1] * shfta_x) + (umat[4] * shfta_y) + (umat[7] * shfta_z);
          const Tcalc rota_z = (umat[2] * shfta_x) + (umat[5] * shfta_y) + (umat[8] * shfta_z);
          const Tcalc dx = shftb_x - rota_x;
          const Tcalc dy = shftb_y - rota_y;
          const Tcalc dz = shftb_z - rota_z;
          result += tmass * ((dx * dx) + (dy * dy) + (dz * dz));
        }
      }
    }
    break;
  case RmsdMethod::NO_ALIGN_MASS:
    if (tcoord_is_sgnint) {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = static_cast<Tcalc>(xcrd_b[i] - xcrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dy = static_cast<Tcalc>(ycrd_b[i] - ycrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dz = static_cast<Tcalc>(zcrd_b[i] - zcrd_a[i]) * inv_gpos_scale_factor;
        result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
      }      
    }
    else {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = xcrd_b[i] - xcrd_a[i];
        const Tcalc dy = ycrd_b[i] - ycrd_a[i];
        const Tcalc dz = zcrd_b[i] - zcrd_a[i];
        result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
      }
    }
    break;
  case RmsdMethod::NO_ALIGN_GEOM:
    if (tcoord_is_sgnint) {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = static_cast<Tcalc>(xcrd_b[i] - xcrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dy = static_cast<Tcalc>(ycrd_b[i] - ycrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dz = static_cast<Tcalc>(zcrd_b[i] - zcrd_a[i]) * inv_gpos_scale_factor;
        result += (dx * dx) + (dy * dy) + (dz * dz);
      }      
    }
    else {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = xcrd_b[i] - xcrd_a[i];
        const Tcalc dy = ycrd_b[i] - ycrd_a[i];
        const Tcalc dz = zcrd_b[i] - zcrd_a[i];
        result += (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
    break;
  }
  return (tcalc_is_double) ? sqrt(result * inv_mass_divisor) : sqrtf(result * inv_mass_divisor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
           const int upper_limit) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const size_t padded_natom_zu = roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t fa_offset = frame_a * padded_natom_zu;
  const size_t fb_offset = frame_b * padded_natom_zu;
  const Tcalc* mass_ptr = (tcalc_is_double) ? (Tcalc*)(cdk.masses) : (Tcalc*)(cdk.sp_masses);
  return rmsd<Tcoord, Tcalc>(&csr.xcrd[fa_offset], &csr.ycrd[fa_offset], &csr.zcrd[fa_offset],
                             &csr.xcrd[fb_offset], &csr.ycrd[fb_offset], &csr.zcrd[fb_offset],
                             mass_ptr, method, lower_limit, upper_limit, csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesWriter<Tcoord> &csw, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
           const int upper_limit) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const size_t padded_natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  const size_t fa_offset = frame_a * padded_natom_zu;
  const size_t fb_offset = frame_b * padded_natom_zu;
  const Tcalc* mass_ptr = (tcalc_is_double) ? (Tcalc*)(cdk.masses) : (Tcalc*)(cdk.sp_masses);
  return rmsd<Tcoord, Tcalc>(&csw.xcrd[fa_offset], &csw.ycrd[fa_offset], &csw.zcrd[fa_offset],
                             &csw.xcrd[fb_offset], &csw.ycrd[fb_offset], &csw.zcrd[fb_offset],
                             mass_ptr, method, lower_limit, upper_limit, csw.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeries<Tcoord> &cs, const int frame_a, const int frame_b,
           const AtomGraph &ag, const RmsdMethod method, const int lower_limit,
           const int upper_limit) {
  const int nframes = cs.getFrameCount();
  if (frame_a < 0 || frame_a >= nframes || frame_b < 0 || frame_b >= nframes) {
    rtErr("Frames " + std::to_string(frame_a) + " and " + std::to_string(frame_b) + " are not "
          "accessible in a series of " + std::to_string(nframes) + " frames.", "rmsd");
  }
  return rmsd(cs.data(), frame_a, frame_b, ag.getChemicalDetailsKit(), method, lower_limit,
              upper_limit);
}

} // namespace structure
} // namespace stormm
