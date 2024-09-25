#include "copyright.h"
#include "Numerics/host_bit_counting.h"
#include "formulas.h"
#include "hilbert_sfc.h"
#include "rounding.h"

namespace stormm {
namespace stmath {

using numerics::hostFfs;
  
//-------------------------------------------------------------------------------------------------
HilbertSFC::HilbertSFC(const int x_span_in, const int y_span_in, const int z_span_in,
                       const HilbertCurveMode method) :
  x_span{x_span_in}, y_span{y_span_in}, z_span{z_span_in},
  curve{},
  curve_indices{}
{
  switch (method) {
  case HilbertCurveMode::STRETCH:
    {
      const int max_pow = std::max(std::max(std::max(hostFfs(x_span),
                                                     hostFfs(y_span)), hostFfs(z_span)), 1);
      x_span = ipow(2, max_pow);
      y_span = x_span;
      z_span = x_span;
      curve = drawHilbertSpaceCurve(max_pow);
    }
    break;
  case HilbertCurveMode::OVERSPAN:
  case HilbertCurveMode::EXACT:
    break;
  }

  // Step along the space-filling curve, noting the indices that each element of the grid space
  // corresponds to.
  const int total_pts = x_span * y_span * z_span;
  curve_indices.resize(total_pts);
  for (int i = 0; i < total_pts; i++) {
    curve_indices[(((curve[i].z * y_span) + curve[i].y) * x_span) + curve[i].x] = i;
  }
}

//-------------------------------------------------------------------------------------------------
int3 HilbertSFC::getDimensions() const {
  return { x_span, y_span, z_span };
}

//-------------------------------------------------------------------------------------------------
int HilbertSFC::getDimension(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return x_span;
  case CartesianDimension::Y:
    return y_span;
  case CartesianDimension::Z:
    return z_span;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int HilbertSFC::getDimension(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return x_span;
  case UnitCellAxis::B:
    return y_span;
  case UnitCellAxis::C:
    return z_span;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::vector<int3>& HilbertSFC::getCurve() const {
  return curve;
}

//-------------------------------------------------------------------------------------------------
int HilbertSFC::getCurveIndex(const int grid_x, const int grid_y, const int grid_z) const {
  return curve_indices[(((grid_z * y_span) + grid_y) * x_span) + grid_x];
}
  
//-------------------------------------------------------------------------------------------------
void propagateHilbertCurve(const int s, std::vector<int3> *curve, const int x, const int y,
                           const int z, const int dx, const int dy, const int dz, const int dx2,
                           const int dy2, const int dz2, const int dx3, const int dy3,
                           const int dz3) {
  if (s == 1) { 
    curve->push_back({ x, y, z });
  }
  else {
    const int next_s = s / 2;
    const int next_x = x - (next_s  * (((dx < 0) * dx) + ((dx2 < 0) * dx2) + ((dx3 < 0) * dx3)));
    const int next_y = y - (next_s  * (((dy < 0) * dy) + ((dy2 < 0) * dy2) + ((dy3 < 0) * dy3)));
    const int next_z = z - (next_s  * (((dz < 0) * dz) + ((dz2 < 0) * dz2) + ((dz3 < 0) * dz3)));
    propagateHilbertCurve(next_s, curve, next_x, next_y, next_z, dx2, dy2, dz2, dx3, dy3, dz3, dx,
                          dy, dz);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx), next_y + (next_s * dy),
                          next_z + (next_s * dz), dx3, dy3, dz3, dx, dy, dz, dx2, dy2, dz2);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx) + (next_s * dx2),
                          next_y + (next_s * dy) + (next_s * dy2),
                          next_z + (next_s * dz) + (next_s * dz2), dx3, dy3, dz3, dx, dy, dz, dx2,
                          dy2, dz2);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx2), next_y + (next_s * dy2),
                          next_z + (next_s * dz2), -dx, -dy, -dz, -dx2, -dy2, -dz2, dx3, dy3, dz3);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx2) + (next_s * dx3),
                          next_y + (next_s * dy2) + (next_s * dy3),
                          next_z + (next_s * dz2) + (next_s * dz3), -dx, -dy, -dz, -dx2, -dy2,
                          -dz2, dx3, dy3, dz3);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx) + (next_s * dx2) + (next_s * dx3),
                          next_y + (next_s * dy) + (next_s * dy2) + (next_s * dy3),
                          next_z + (next_s * dz) + (next_s * dz2) + (next_s * dz3), -dx3, -dy3,
                          -dz3, dx, dy, dz, -dx2, -dy2, -dz2);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx) + (next_s * dx3),
                          next_y + (next_s * dy) + (next_s * dy3),
                          next_z + (next_s * dz) + (next_s * dz3), -dx3, -dy3, -dz3, dx, dy, dz,
                          -dx2, -dy2, -dz2);
    propagateHilbertCurve(next_s, curve, next_x + (next_s * dx3), next_y + (next_s * dy3),
                          next_z + (next_s * dz3), dx2, dy2, dz2, -dx3, -dy3, -dz3, -dx, -dy, -dz);
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int3> drawHilbertSpaceCurve(const int order) {
  if (order > 10) {
    rtErr("A Hilbert Space Curve of order " + std::to_string(order) + " would be excessively "
          "large.", "drawHilbertSpaceCurve");
  }
  if (order < 0) {
    rtErr("A Hilbert Space Curve of order " + std::to_string(order) + " is invalid.",
          "drawHilbertSpaceCurve");
  }
  std::vector<int3> result;
  const int span = ipow(2, order);
  result.reserve(span * span * span);
  propagateHilbertCurve(span, &result);
  return result;
}

} // namespace stmath
} // namespace stormm
