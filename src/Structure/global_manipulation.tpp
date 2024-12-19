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
template <typename T, typename Tcalc>
void rotateCoordinates(T *x, T *y, T *z, const Tcalc* axis_a, const Tcalc* axis_b,
                       const Tcalc* axis_c, const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<T>()) {
    const Tcalc xval = static_cast<Tcalc>(*x) / globalpos_scale_factor;
    const Tcalc yval = static_cast<Tcalc>(*y) / globalpos_scale_factor;
    const Tcalc zval = static_cast<Tcalc>(*z) / globalpos_scale_factor;
    const Tcalc tmp_x = (axis_a[0] * xval) + (axis_a[1] * yval) + (axis_a[2] * zval);
    const Tcalc tmp_y = (axis_b[0] * xval) + (axis_b[1] * yval) + (axis_b[2] * zval);
    const Tcalc tmp_z = (axis_c[0] * xval) + (axis_c[1] * yval) + (axis_c[2] * zval);
    *x = llround(tmp_x * globalpos_scale_factor);
    *y = llround(tmp_y * globalpos_scale_factor);
    *z = llround(tmp_z * globalpos_scale_factor);
  }
  else {
    const Tcalc xval = *x;
    const Tcalc yval = *y;
    const Tcalc zval = *z;
    const Tcalc tmp_x = (axis_a[0] * xval) + (axis_a[1] * yval) + (axis_a[2] * zval);
    const Tcalc tmp_y = (axis_b[0] * xval) + (axis_b[1] * yval) + (axis_b[2] * zval);
    const Tcalc tmp_z = (axis_c[0] * xval) + (axis_c[1] * yval) + (axis_c[2] * zval);
    *x = tmp_x;
    *y = tmp_y;
    *z = tmp_z;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(T *x, T *y, T *z, const Tcalc* umat, const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<T>()) {
    const Tcalc xval = static_cast<Tcalc>(*x) / globalpos_scale_factor;
    const Tcalc yval = static_cast<Tcalc>(*y) / globalpos_scale_factor;
    const Tcalc zval = static_cast<Tcalc>(*z) / globalpos_scale_factor;
    const Tcalc tmp_x = (umat[0] * xval) + (umat[3] * yval) + (umat[6] * zval);
    const Tcalc tmp_y = (umat[1] * xval) + (umat[4] * yval) + (umat[7] * zval);
    const Tcalc tmp_z = (umat[2] * xval) + (umat[5] * yval) + (umat[8] * zval);
    *x = llround(tmp_x * globalpos_scale_factor);
    *y = llround(tmp_y * globalpos_scale_factor);
    *z = llround(tmp_z * globalpos_scale_factor);
  }
  else {
    const Tcalc xval = *x;
    const Tcalc yval = *y;
    const Tcalc zval = *z;
    const Tcalc tmp_x = (umat[0] * xval) + (umat[3] * yval) + (umat[6] * zval);
    const Tcalc tmp_y = (umat[1] * xval) + (umat[4] * yval) + (umat[7] * zval);
    const Tcalc tmp_z = (umat[2] * xval) + (umat[5] * yval) + (umat[8] * zval);
    *x = tmp_x;
    *y = tmp_y;
    *z = tmp_z;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(T* x, T* y, T* z, const size_t n, const Tcalc* axis_a, const Tcalc* axis_b,
                       const Tcalc* axis_c, const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < n; i++) {
      const Tcalc xval = static_cast<Tcalc>(x[i]) / globalpos_scale_factor;
      const Tcalc yval = static_cast<Tcalc>(y[i]) / globalpos_scale_factor;
      const Tcalc zval = static_cast<Tcalc>(z[i]) / globalpos_scale_factor;
      const Tcalc tmp_x = (axis_a[0] * xval) + (axis_a[1] * yval) + (axis_a[2] * zval);
      const Tcalc tmp_y = (axis_b[0] * xval) + (axis_b[1] * yval) + (axis_b[2] * zval);
      const Tcalc tmp_z = (axis_c[0] * xval) + (axis_c[1] * yval) + (axis_c[2] * zval);
      x[i] = llround(tmp_x * globalpos_scale_factor);
      y[i] = llround(tmp_y * globalpos_scale_factor);
      z[i] = llround(tmp_z * globalpos_scale_factor);
    }
  }
  else {
    for (size_t i = 0; i < n; i++) {
      const Tcalc xval = x[i];
      const Tcalc yval = y[i];
      const Tcalc zval = z[i];
      const Tcalc tmp_x = (axis_a[0] * xval) + (axis_a[1] * yval) + (axis_a[2] * zval);
      const Tcalc tmp_y = (axis_b[0] * xval) + (axis_b[1] * yval) + (axis_b[2] * zval);
      const Tcalc tmp_z = (axis_c[0] * xval) + (axis_c[1] * yval) + (axis_c[2] * zval);
      x[i] = tmp_x;
      y[i] = tmp_y;
      z[i] = tmp_z;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(T* x, T* y, T* z, const size_t n, const Tcalc* umat,
                       const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < n; i++) {
      const Tcalc xval = static_cast<Tcalc>(x[i]) / globalpos_scale_factor;
      const Tcalc yval = static_cast<Tcalc>(y[i]) / globalpos_scale_factor;
      const Tcalc zval = static_cast<Tcalc>(z[i]) / globalpos_scale_factor;
      const Tcalc tmp_x = (umat[0] * xval) + (umat[3] * yval) + (umat[6] * zval);
      const Tcalc tmp_y = (umat[1] * xval) + (umat[4] * yval) + (umat[7] * zval);
      const Tcalc tmp_z = (umat[2] * xval) + (umat[5] * yval) + (umat[8] * zval);
      x[i] = llround(tmp_x * globalpos_scale_factor);
      y[i] = llround(tmp_y * globalpos_scale_factor);
      z[i] = llround(tmp_z * globalpos_scale_factor);
    }
  }
  else {
    for (size_t i = 0; i < n; i++) {
      const Tcalc xval = x[i];
      const Tcalc yval = y[i];
      const Tcalc zval = z[i];
      const Tcalc tmp_x = (umat[0] * xval) + (umat[3] * yval) + (umat[6] * zval);
      const Tcalc tmp_y = (umat[1] * xval) + (umat[4] * yval) + (umat[7] * zval);
      const Tcalc tmp_z = (umat[2] * xval) + (umat[5] * yval) + (umat[8] * zval);
      x[i] = tmp_x;
      y[i] = tmp_y;
      z[i] = tmp_z;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                       const std::vector<Tcalc> &axis_a, const std::vector<Tcalc> &axis_b,
                       const std::vector<Tcalc> &axis_c, const Tcalc globalpos_scale_factor) {
  const size_t target_size = x->size();
  if (target_size != y->size() || target_size != z->size()) {
    rtErr("Standard Template Library vectors of sizes " + std::to_string(target_size) + ", " +
          std::to_string(y->size()) + ", and " + std::to_string(z->size()) + " are incompatible.",
          "rotateCoordinates");
  }
  rotateCoordinates(x->data(), y->data(), z->data(), target_size, axis_a.data(), axis_b.data(),
                    axis_c.data(), globalpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                       const std::vector<Tcalc> &umat, const Tcalc globalpos_scale_factor) {
  const size_t target_size = x->size();
  if (target_size != y->size() || target_size != z->size()) {
    rtErr("Standard Template Library vectors of sizes " + std::to_string(target_size) + ", " +
          std::to_string(y->size()) + ", and " + std::to_string(z->size()) + " are incompatible.",
          "rotateCoordinates");
  }
  rotateCoordinates(x->data(), y->data(), z->data(), target_size, umat.data(),
                    globalpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z, const Hybrid<Tcalc> &axis_a,
                       const Hybrid<Tcalc> &axis_b, const Hybrid<Tcalc> &axis_c,
                       const Tcalc globalpos_scale_factor) {
  const size_t target_size = x->size();
  if (target_size != y->size() || target_size != z->size()) {
    rtErr("Hybrid objects of sizes " + std::to_string(target_size), + ", " +
          std::to_string(y->size()) + ", and " + std::to_string(z->size()) + " are incompatible.",
          "rotateCoordinates");
  }
  rotateCoordinates(x->data(), y->data(), z->data(), target_size, axis_a.data(), axis_b.data(),
                    axis_c.data(), globalpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc>
void rotateCoordinates(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z, const Hybrid<Tcalc> &umat,
                       const Tcalc globalpos_scale_factor) {
  const size_t target_size = x->size();
  if (target_size != y->size() || target_size != z->size()) {
    rtErr("Hybrid objects of sizes " + std::to_string(target_size), + ", " +
          std::to_string(y->size()) + ", and " + std::to_string(z->size()) + " are incompatible.",
          "rotateCoordinates");
  }
  rotateCoordinates(x->data(), y->data(), z->data(), target_size, umat.data(),
                    globalpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename T3, typename Tcalc>
void rotateCoordinates(T3 *c, const T3 axis_a, const T3 axis_b, const T3 axis_c,
                       const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralHpcVectorType<T3>()) {
    const Tcalc xval = static_cast<Tcalc>(c->x) / globalpos_scale_factor;
    const Tcalc yval = static_cast<Tcalc>(c->y) / globalpos_scale_factor;
    const Tcalc zval = static_cast<Tcalc>(c->z) / globalpos_scale_factor;
    const Tcalc tmp_x = (axis_a.x * xval) + (axis_a.y * yval) + (axis_a.z * zval);
    const Tcalc tmp_y = (axis_b.x * xval) + (axis_b.y * yval) + (axis_b.z * zval);
    const Tcalc tmp_z = (axis_c.x * xval) + (axis_c.y * yval) + (axis_c.z * zval);
    c->x = llround(tmp_x * globalpos_scale_factor);
    c->y = llround(tmp_y * globalpos_scale_factor);
    c->z = llround(tmp_z * globalpos_scale_factor);
  }
  else {
    const Tcalc xval = c->x;
    const Tcalc yval = c->y;
    const Tcalc zval = c->z;
    const Tcalc tmp_x = (axis_a.x * xval) + (axis_a.y * yval) + (axis_a.z * zval);
    const Tcalc tmp_y = (axis_b.x * xval) + (axis_b.y * yval) + (axis_b.z * zval);
    const Tcalc tmp_z = (axis_c.x * xval) + (axis_c.y * yval) + (axis_c.z * zval);
    c->x = tmp_x;
    c->y = tmp_y;
    c->z = tmp_z;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tcalc alpha,
                       const Tcalc beta, const Tcalc gamma, const int lower_limit,
                       const int upper_limit, const VirtualSiteKit<Tcalc> *vsk,
                       const Tcalc globalpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const std::vector<Tcalc> rmat = beardRotationMatrix(alpha, beta, gamma);
  rotateCoordinates<Tcoord, Tcalc>(&xcrd[lower_limit], &ycrd[lower_limit], &zcrd[lower_limit],
                                   upper_limit - lower_limit, rmat.data(), globalpos_scale_factor);

  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites<Tcoord, Tcalc>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                     *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeriesWriter<Tcoord> *csw, const size_t frame_index,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc alpha, const Tcalc beta,
                       const Tcalc gamma, const int lower_limit, const int upper_limit) {
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : csw->natom;
  const size_t atom_start = static_cast<size_t>(roundUp(csw->natom, warp_size_int)) * frame_index;
  rotateCoordinates<Tcoord, Tcalc>(&csw->xcrd[atom_start], &csw->ycrd[atom_start],
                                   &csw->zcrd[atom_start], alpha, beta, gamma, lower_limit,
                                   actual_upper_limit);
  placeVirtualSites<Tcoord, Tcalc>(&csw->xcrd[atom_start], &csw->ycrd[atom_start],
                                   &csw->zcrd[atom_start], nullptr, nullptr, UnitCellType::NONE,
                                   vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeries<Tcoord> *cs, const size_t frame_index, const AtomGraph &ag,
                       const Tcalc alpha, const Tcalc beta, const Tcalc gamma,
                       const int lower_limit, const int upper_limit) {

  // The calculation type is either single- or double-precision
  CoordinateSeriesWriter<Tcoord> csw = cs->data();
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
    rotateCoordinates<Tcoord, double>(&csw, frame_index, vsk, alpha, beta, gamma, lower_limit,
                                      upper_limit);
  }
  else {
    const VirtualSiteKit<float> vsk = ag.getSinglePrecisionVirtualSiteKit();
    rotateCoordinates<Tcoord, float>(&csw, frame_index, vsk, alpha, beta, gamma, lower_limit,
                                     upper_limit);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateCoordinates(PhaseSpaceSynthesis *poly_ps, const int system_index, const Tcalc alpha,
                       const Tcalc beta, const Tcalc gamma, const int lower_limit,
                       const int upper_limit) {
  const AtomGraph *ag = poly_ps->getSystemTopologyPointer(system_index);
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
    CoordinateFrame cf = poly_ps->exportCoordinates(system_index);
    CoordinateFrameWriter cfr = cf.data();
    const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfr.natom;
    rotateCoordinates<double, double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, alpha, beta, gamma,
                                      lower_limit, actual_upper_limit, &vsk);
    poly_ps->importSystem(cf, system_index, TrajectoryKind::POSITIONS);
  }
  else {
    const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
    CoordinateFrame cf = poly_ps->exportCoordinates(system_index);
    CoordinateFrameWriter cfr = cf.data();
    const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfr.natom;
    rotateCoordinates<double, float>(cfr.xcrd, cfr.ycrd, cfr.zcrd, alpha, beta, gamma,
                                      lower_limit, actual_upper_limit, &vsk);
    poly_ps->importSystem(cf, system_index, TrajectoryKind::POSITIONS);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateCoordinates(CondensateWriter *cdnsw, const int system_index,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc alpha, const Tcalc beta,
                       const Tcalc gamma, const int lower_limit, const int upper_limit) {
  const size_t atom_start = cdnsw->atom_starts[system_index];
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               cdnsw->atom_counts[system_index];
  switch (cdnsw->mode) {
  case PrecisionModel::DOUBLE:
    rotateCoordinates<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], alpha, beta, gamma, lower_limit,
                                     actual_upper_limit);
    placeVirtualSites<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], nullptr, nullptr,
                                     UnitCellType::NONE, vsk);
    break;
  case PrecisionModel::SINGLE:
    rotateCoordinates<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
    placeVirtualSites<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
    break;
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
    placeVirtualSites<Tcoord, Tcalc>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                     *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(PsSynthesisWriter *poly_psw, const int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit) {
  const int llim = poly_psw->atom_starts[system_index] + lower_limit;
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               poly_psw->atom_counts[system_index];
  const int hlim = poly_psw->atom_starts[system_index] + actual_upper_limit;

  // Translate the coordinates for the selected atoms
  const int95_t ixmove = hostDoubleToInt95(xmove * poly_psw->gpos_scale);
  const int95_t iymove = hostDoubleToInt95(ymove * poly_psw->gpos_scale);
  const int95_t izmove = hostDoubleToInt95(zmove * poly_psw->gpos_scale);
  for (int i = llim; i < hlim; i++) {
    const int95_t xnew = hostSplitFPSum(ixmove, poly_psw->xcrd[i], poly_psw->xcrd_ovrf[i]);
    const int95_t ynew = hostSplitFPSum(iymove, poly_psw->ycrd[i], poly_psw->ycrd_ovrf[i]);
    const int95_t znew = hostSplitFPSum(izmove, poly_psw->zcrd[i], poly_psw->zcrd_ovrf[i]);
    poly_psw->xcrd[i]      = xnew.x;
    poly_psw->xcrd_ovrf[i] = xnew.y;
    poly_psw->ycrd[i]      = ynew.x;
    poly_psw->ycrd_ovrf[i] = ynew.y;
    poly_psw->zcrd[i]      = znew.x;
    poly_psw->zcrd_ovrf[i] = znew.y;
  }

  // Reposition the virtual sites
  const int natom = poly_psw->atom_counts[system_index];
  CoordinateFrame cf(natom);
  CoordinateFrameWriter cfw = cf.data();
  for (int i = 0; i < natom; i++) {
    cfw.xcrd[i] = hostInt95ToDouble(poly_psw->xcrd[i], poly_psw->xcrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
    cfw.ycrd[i] = hostInt95ToDouble(poly_psw->ycrd[i], poly_psw->ycrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
    cfw.zcrd[i] = hostInt95ToDouble(poly_psw->zcrd[i], poly_psw->zcrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
  }
  placeVirtualSites<double, Tcalc>(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr,
                                   UnitCellType::NONE, vsk);
  for (int i = 0; i < vsk.nsite; i++) {
    bool frame_moved;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[i])) {
    case VirtualSiteKind::NONE:
      frame_moved = false;
      break;
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      frame_moved = ((vsk.frame1_idx[i] >= llim && vsk.frame1_idx[i] < hlim) ||
                     (vsk.frame2_idx[i] >= llim && vsk.frame2_idx[i] < hlim));
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      frame_moved = ((vsk.frame1_idx[i] >= llim && vsk.frame1_idx[i] < hlim) ||
                     (vsk.frame2_idx[i] >= llim && vsk.frame2_idx[i] < hlim) ||
                     (vsk.frame3_idx[i] >= llim && vsk.frame3_idx[i] < hlim));
      break;
    case VirtualSiteKind::FIXED_4:
      frame_moved = ((vsk.frame1_idx[i] >= llim && vsk.frame1_idx[i] < hlim) ||
                     (vsk.frame2_idx[i] >= llim && vsk.frame2_idx[i] < hlim) ||
                     (vsk.frame3_idx[i] >= llim && vsk.frame3_idx[i] < hlim) ||
                     (vsk.frame4_idx[i] >= llim && vsk.frame4_idx[i] < hlim));
      break;
    }
    if (frame_moved) {
      const int vs_idx = vsk.vs_atoms[i];
      const int95_t xvs_new = hostDoubleToInt95(cfw.xcrd[vs_idx] * poly_psw->gpos_scale);
      const int95_t yvs_new = hostDoubleToInt95(cfw.ycrd[vs_idx] * poly_psw->gpos_scale);
      const int95_t zvs_new = hostDoubleToInt95(cfw.zcrd[vs_idx] * poly_psw->gpos_scale);
      poly_psw->xcrd[i]      = xvs_new.x;
      poly_psw->xcrd_ovrf[i] = xvs_new.y;
      poly_psw->ycrd[i]      = yvs_new.x;
      poly_psw->ycrd_ovrf[i] = yvs_new.y;
      poly_psw->zcrd[i]      = zvs_new.x;
      poly_psw->zcrd_ovrf[i] = zvs_new.y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(PhaseSpaceSynthesis *poly_ps, const int system_index, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit) {
  const AtomGraph *ag = poly_ps->getSystemTopologyPointer(system_index);
  PsSynthesisWriter poly_psw = poly_ps->data();
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
    translateCoordinates(&poly_psw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                         upper_limit);
  }
  else {
    const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
    translateCoordinates(&poly_psw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                         upper_limit);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(CondensateWriter *cdnsw, const int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, const Tcalc xmove, const Tcalc ymove,
                          const Tcalc zmove, const int lower_limit, const int upper_limit) {
  const size_t atom_start = cdnsw->atom_starts[system_index];
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               cdnsw->atom_counts[system_index];
  switch (cdnsw->mode) {
  case PrecisionModel::DOUBLE:
    translateCoordinates<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                        &cdnsw->zcrd[atom_start], xmove, ymove, zmove, lower_limit,
                                        actual_upper_limit);
    placeVirtualSites<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], nullptr, nullptr,
                                     UnitCellType::NONE, vsk);
    break;
  case PrecisionModel::SINGLE:
    translateCoordinates<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                       &cdnsw->zcrd_sp[atom_start], xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
    placeVirtualSites<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
    break;
  }
}

} // namespace structure
} // namespace stormm
