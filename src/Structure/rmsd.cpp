#include "Math/vector_ops.h"
#include "rmsd.h"

namespace omni {
namespace structure {

using math::sum;

//-------------------------------------------------------------------------------------------------
double rmsd(const double* xcrd_a, const double* ycrd_a, const double* zcrd_a, const double* xcrd_b,
            const double* ycrd_b, const double* zcrd_b, const double* masses,
            const RmsdMethod method,  const int lower_limit, const int upper_limit) {
  double result = 0.0;
  double inv_mass_divisor;
  switch (method) {
  case RmsdMethod::ALIGN_MASS:
  case RmsdMethod::NO_ALIGN_MASS:
    inv_mass_divisor = 1.0 / sum<double>(&masses[lower_limit], upper_limit - lower_limit);
    break;
  case RmsdMethod::ALIGN_GEOM:
  case RmsdMethod::NO_ALIGN_GEOM:
    inv_mass_divisor = 1.0 / static_cast<double>(upper_limit - lower_limit);
    break;
  }
  switch (method) {
  case RmsdMethod::ALIGN_MASS:
  case RmsdMethod::ALIGN_GEOM:
    {
      // In order to compute the aligned RMSD without moving the coordinates, the movement must
      // be computed in temporary variables only. 
      double sab_xx = 0.0;
      double sab_xy = 0.0;
      double sab_xz = 0.0;
      double sab_yx = 0.0;
      double sab_yy = 0.0;
      double sab_yz = 0.0;
      double sab_zx = 0.0;
      double sab_zy = 0.0;
      double sab_zz = 0.0;
      double sa_x = 0.0;
      double sa_y = 0.0;
      double sa_z = 0.0;
      double sb_x = 0.0;
      double sb_y = 0.0;
      double sb_z = 0.0;
      const bool use_mass = (method == RmsdMethod::ALIGN_MASS);
      for (int i = lower_limit; i < upper_limit; i++) {
        const double tmass = (use_mass) ? masses[i] : 1.0;
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
      const double coma_x = sa_x * inv_mass_divisor;
      const double coma_y = sa_y * inv_mass_divisor;
      const double coma_z = sa_z * inv_mass_divisor;
      const double comb_x = sb_x * inv_mass_divisor;
      const double comb_y = sb_y * inv_mass_divisor;
      const double comb_z = sb_z * inv_mass_divisor;

      // Assemble the Kabsch matrix and diagonalize it
      const double aa = sab_xx - (sb_x * coma_x) - (sa_x * comb_x) + (coma_x * comb_x);
      const double ab = sab_xy - (sb_y * coma_x) - (sa_x * comb_y) + (coma_x * comb_y);
      const double ac = sab_xz - (sb_z * coma_x) - (sa_x * comb_z) + (coma_x * comb_z);
      const double ba = sab_yx - (sb_x * coma_y) - (sa_y * comb_x) + (coma_y * comb_x);
      const double bb = sab_yy - (sb_y * coma_y) - (sa_y * comb_y) + (coma_y * comb_y);
      const double bc = sab_yz - (sb_z * coma_y) - (sa_y * comb_z) + (coma_y * comb_z);
      const double ca = sab_zx - (sb_x * coma_z) - (sa_z * comb_x) + (coma_z * comb_x);
      const double cb = sab_zy - (sb_y * coma_z) - (sa_z * comb_y) + (coma_z * comb_y);
      const double cc = sab_zz - (sb_z * coma_z) - (sa_z * comb_z) + (coma_z * comb_z);
      std::vector<double> rmat(16);
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
      std::vector<double> vmat(16, 0.0), eigval(4, 0.0);
      jacobiEigensolver(&rmat, &vmat, &eigval, 4);
      eigenSort(&eigval, &vmat, 4);

      // Form the rotation matrix
      const double a = vmat[0];
      const double x = vmat[1];
      const double y = vmat[2];
      const double z = vmat[3];
      std::vector<double> umat(9);
      u[0] = (a * a) + (x * x) - (y * y) - (z * z);
      u[3] = 2.0 * ((x * y) + (a * z));
      u[6] = 2.0 * ((x * z) - (a * y));
      u[1] = 2.0 * ((x * y) - (a * z));
      u[4] = (a * a) - (x * x) + (y * y) - (z * z);
      u[7] = 2.0 * ((y * z) + (a * x));
      u[2] = 2.0 * ((x * z) + (a * y));
      u[5] = 2.0 * ((y * z) - (a * x));
      u[8] = (a * a) - (x * x) - (y * y) + (z * z);

      // Shift and rotate the coordinates of the first frame (in temporary variables only) and
      // compare them to the shifted, unrotated coordinates of the second frame.
      for (int i = lower_limit; i < upper_limit; i++) {
        const double shfta_x = xcrd_a[i] - coma_x;
        const double shfta_y = ycrd_a[i] - coma_y;
        const double shfta_z = zcrd_a[i] - coma_z;
        const double shftb_x = xcrd_b[i] - comb_x;
        const double shftb_y = ycrd_b[i] - comb_y;
        const double shftb_z = zcrd_b[i] - comb_z;
        const double rota_x = (u[0] * shfta_x) + (u[3] * shfta_y) + (u[6] * shfta_z);
        const double rota_y = (u[1] * shfta_x) + (u[4] * shfta_y) + (u[7] * shfta_z);
        const double rota_z = (u[2] * shfta_x) + (u[5] * shfta_y) + (u[8] * shfta_z);
        const double dx = shftb_x - rota_x;
        const double dy = shftb_y - rota_y;
        const double dz = shftb_z - rota_z;
        result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
      } 
    }
    break;
  case RmsdMethod::NO_ALIGN_MASS:
    for (int i = lower_limit; i < upper_limit; i++) {
      const double dx = xcrd_b[i] - xcrd_a[i];
      const double dy = ycrd_b[i] - ycrd_a[i];
      const double dz = zcrd_b[i] - zcrd_a[i];
      result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
    }
    break;
  case RmsdMethod::NO_ALIGN_GEOM:
    for (int i = lower_limit; i < upper_limit; i++) {
      const double dx = xcrd_b[i] - xcrd_a[i];
      const double dy = ycrd_b[i] - ycrd_a[i];
      const double dz = zcrd_b[i] - zcrd_a[i];
      result += (dx * dx) + (dy * dy) + (dz * dz);
    }
    break;
  }
  return sqrt(result * inv_mass_divisor);
}
  
//-------------------------------------------------------------------------------------------------
double rmsd(const PhaseSpaceReader &psr_a, const PhaseSpaceReader &psr_b,
            const ChemicalDetailsKit &cdk, const RmsdMethod method, const int lower_limit,
            const int upper_limit) {
}

//-------------------------------------------------------------------------------------------------
double rmsd(const PhaseSpace &ps_a, const PhaseSpace &ps_b, const AtomGraph &ag,
            const RmsdMethod method, const int lower_limit, const int upper_limit) {

}

//-------------------------------------------------------------------------------------------------
double massWeightedRmsd(const PhaseSpace &ps_a, const PhaseSpace &ps_b const ) {

}

} // namespace structure
} // namespace omni
