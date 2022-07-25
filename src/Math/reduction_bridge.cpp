#include "reduction_bridge.h"

namespace omni {
namespace math {

using card::HybridKind;

//-------------------------------------------------------------------------------------------------
ReductionBridge::ReductionBridge(const size_t n_values) :
  x_buffer{HybridKind::POINTER, "bridge_xbuff"},
  y_buffer{HybridKind::POINTER, "bridge_ybuff"},
  z_buffer{HybridKind::POINTER, "bridge_zbuff"},
  storage{3LLU * roundUp(n_values, warp_size_zu)}
{
  const size_t padded_nval = roundUp(n_values, warp_size_zu);
  x_buffer.setPointer(&storage,                  0, n_values);
  y_buffer.setPointer(&storage,        padded_nval, n_values);
  z_buffer.setPointer(&storage, 2LLU * padded_nval, n_values);
}

//-------------------------------------------------------------------------------------------------
const double* ReductionBridge::getPointer(const CartesianDimension cdim,
                                          const HybridTargetLevel tier) const {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* ReductionBridge::getPointer(const CartesianDimension cdim, const HybridTargetLevel tier) {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);
  }
  __builtin_unreachable();
}

} // namespace math
} // namespace omni
