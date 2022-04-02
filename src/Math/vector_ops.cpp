#include "vector_ops.h"

namespace omni {
namespace math {

//-------------------------------------------------------------------------------------------------
void accumulateBitmask(uint* va, const size_t pos) {
  const size_t acc_elem = pos / uint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * uint_bit_count_zu);
  va[acc_elem] |= (0x1 << acc_bit);
}
  
//-------------------------------------------------------------------------------------------------
void accumulateBitmask(std::vector<uint> *va, const size_t pos) {
  if (pos >= va->size() * uint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned integer bit mask of "
          "length " + std::to_string(va->size()) + ".", "accumulateBitmask");
  }
  accumulateBitmask(va->data(), pos);
}

//-------------------------------------------------------------------------------------------------
void accumulateBitmask(ullint* va, const size_t pos) {
  const size_t acc_elem = pos / ullint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ullint_bit_count_zu);
  va[acc_elem] |= (0x1 << acc_bit);
}
  
//-------------------------------------------------------------------------------------------------
void accumulateBitmask(std::vector<ullint> *va, const size_t pos) {
  if (pos >= va->size() * ullint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned long long integer bit "
          "mask of length " + std::to_string(va->size()) + ".", "accumulateBitmask");
  }
  accumulateBitmask(va->data(), pos);
}

//-------------------------------------------------------------------------------------------------
void accumulateBitmask(ushort* va, const size_t pos) {
  const size_t acc_elem = pos / ushort_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ushort_bit_count_zu);
  va[acc_elem] |= (0x1 << acc_bit);
}
  
//-------------------------------------------------------------------------------------------------
void accumulateBitmask(std::vector<ushort> *va, const size_t pos) {
  if (pos >= va->size() * ushort_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned short integer bit "
          "mask of length " + std::to_string(va->size()) + ".", "accumulateBitmask");
  }
  accumulateBitmask(va->data(), pos);
}

} // namespace math
} // namespace omni
