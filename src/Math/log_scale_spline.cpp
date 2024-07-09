#include "copyright.h"
#include "log_scale_spline.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
ullint doublePrecisionSplineDetailMask(const int mantissa_bits) {
  ullint result = 0LLU;
  for (int i = 0; i < 52 - mantissa_bits; i++) {
    result |= (0x1LLU << i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
uint singlePrecisionSplineDetailMask(const int mantissa_bits) {
  uint result = 0U;
  for (int i = 0; i < 23 - mantissa_bits; i++) {
    result |= (0x1U << i);
  }
  return result;
}

} // namespace stmath
} // namespace stormm
