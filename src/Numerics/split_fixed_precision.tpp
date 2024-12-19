// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace numerics {

//-------------------------------------------------------------------------------------------------
template <typename T> void hostSplitFPSum(T *a_x, int* a_y, const int95_t b) {
  const int95_t result = hostInt95Sum(*a_x, *a_y, b.x, b.y);
  *a_x = result.x;
  *a_y = result.y;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void hostSplitFPSum(T *a_x, int* a_y, const int2 b) {
  const int2 result = hostInt63Sum(*a_x, *a_y, b.x, b.y);
  *a_x = result.x;
  *a_y = result.y;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void hostSplitFPSubtract(T *a_x, int* a_y, const int95_t b) {
  const int95_t result = hostInt95Subtract(*a_x, *a_y, b.x, b.y);
  *a_x = result.x;
  *a_y = result.y;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void hostSplitFPSubtract(T *a_x, int* a_y, const int2 b) {
  const int2 result = hostInt63Subtract(*a_x, *a_y, b.x, b.y);
  *a_x = result.x;
  *a_y = result.y;
}


} // namespace numerics
} // namespace stormm
