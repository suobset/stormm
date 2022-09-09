// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> beardRotationMatrix(const T om_x, const T om_y, const T om_z) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);

  std::vector<T> result(9);

  // Convenient quantities
  const T om2_x = om_x * om_x;
  const T om2_y = om_y * om_y;
  const T om2_z = om_z * om_z;
  const T om2 = om2_x + om2_y + om2_z;
  const T om = (t_is_double) ? sqrt(om2) : sqrtf(om2);
  const T cos_om = (t_is_double) ? cos(om) : cosf(om);
  const T sin_om = (t_is_double) ? sin(om) : sinf(om);
  const T inv_om = 1.0 / om;
  const T inv_om2 = inv_om * inv_om;

  // Compute rotation matrix
  result[0] = (((om2_y + om2_z) * cos_om) + om2_x) * inv_om2;
  result[3] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) + ((om_z * inv_om) * sin_om);
  result[6] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) - ((om_y * inv_om) * sin_om);
  result[1] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) - ((om_z * inv_om) * sin_om);
  result[4] = (((om2_x + om2_z) * cos_om) + om2_y) * inv_om2;
  result[7] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) + ((om_x * inv_om) * sin_om);
  result[2] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) + ((om_y * inv_om) * sin_om);
  result[5] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) - ((om_x * inv_om) * sin_om);
  result[8] = (((om2_x + om2_y) * cos_om) + om2_z) * inv_om2;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tcalc alpha,
                       const Tcalc beta, const Tcalc gamma, const int lower_limit,
                       const int upper_limit, const VirtualSiteKit<Tcalc> *vsk,
                       const Tcalc globalpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const std::vector<Tcalc> rmat = beardRotationMatrix(alpha, beta, gamma);
  if (isSignedIntegralScalarType<Tcoord>()) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;
    for (int i = lower_limit; i < upper_limit; i++) {
      const Tcalc xc = static_cast<Tcalc>(xcrd[i]) * inv_globalpos_scale_factor;
      const Tcalc yc = static_cast<Tcalc>(ycrd[i]) * inv_globalpos_scale_factor;
      const Tcalc zc = static_cast<Tcalc>(zcrd[i]) * inv_globalpos_scale_factor;
      xcrd[i] = llround(((rmat[0] * xc) + (rmat[3] * yc) + (rmat[6] * zc)) *
                        globalpos_scale_factor);
      ycrd[i] = llround(((rmat[1] * xc) + (rmat[4] * yc) + (rmat[7] * zc)) *
                        globalpos_scale_factor);
      zcrd[i] = llround(((rmat[2] * xc) + (rmat[5] * yc) + (rmat[8] * zc)) *
                        globalpos_scale_factor);
    }
  }
  else {
    for (int i = lower_limit; i < upper_limit; i++) {
      const Tcoord xc = xcrd[i];
      const Tcoord yc = ycrd[i];
      const Tcoord zc = zcrd[i];
      xcrd[i] = (rmat[0] * xc) + (rmat[3] * yc) + (rmat[6] * zc);
      ycrd[i] = (rmat[1] * xc) + (rmat[4] * yc) + (rmat[7] * zc);
      zcrd[i] = (rmat[2] * xc) + (rmat[5] * yc) + (rmat[8] * zc);
    }
  }

  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites<double, double>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                      *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void translateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit, const VirtualSiteKit<Tcalc> *vsk,
                          const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<Tcoord>()) {
    for (int i = lower_limit; i < upper_limit; i++) {
      xcrd[i] += llround(xmove * globalpos_scale_factor);
      ycrd[i] += llround(ymove * globalpos_scale_factor);
      zcrd[i] += llround(zmove * globalpos_scale_factor);
    }
  }
  else {
    for (int i = lower_limit; i < upper_limit; i++) {
      xcrd[i] += xmove;
      ycrd[i] += ymove;
      zcrd[i] += zmove;
    }
  }
  
  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites<double, double>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                      *vsk);
  }
}

} // namespace structure
} // namespace stormm
