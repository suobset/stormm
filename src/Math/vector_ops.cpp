#include "copyright.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
bool foundIn(const std::vector<bool> &va, const std::vector<bool> &vb, const DataOrder format) {
  const size_t n_a = va.size();
  const size_t n_b = vb.size();
  bool a_has_true  = false;
  bool a_has_false = false;
  for (size_t i = 0; i < n_a; i++) {
    a_has_true  = (a_has_true || va[i]);
    a_has_false = (a_has_false || (! va[i]));
  }
  bool b_has_true  = false;
  bool b_has_false = false;
  for (size_t i = 0; i < n_b; i++) {
    b_has_true  = (b_has_true || vb[i]);
    b_has_false = (b_has_false || (! vb[i]));
  }
  bool result = true;
  if (a_has_true && (! b_has_true)) {
    result = false;
  }
  if (a_has_false && (! b_has_false)) {
    result = false;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> reduceUniqueValues(const std::vector<bool> &va, const std::vector<bool> &vb) {
  const size_t n_b = vb.size();
  std::vector<bool> result = va;
  for (size_t i = 0; i < n_b; i++) {
    bool found = false;
    const size_t n_r = result.size();
    for (size_t j = 0; j < n_r; j++) {
      found = (found || result[j] == vb[i]);
    }
    if (! found) {
      result.push_back(vb[i]);
    }
  }
  return result;
}

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
  va[acc_elem] |= (0x1LLU << acc_bit);
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

//-------------------------------------------------------------------------------------------------
void unsetBitInMask(uint* va, const size_t pos) {
  const size_t acc_elem = pos / uint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * uint_bit_count_zu);
  va[acc_elem] &= ~(0x1 << acc_bit);
}

//-------------------------------------------------------------------------------------------------
void unsetBitInMask(std::vector<uint> *va, const size_t pos) {
  if (pos >= va->size() * uint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned integer bit mask of "
          "length " + std::to_string(va->size()) + ".", "unsetBitInMask");
  }
  unsetBitInMask(va->data(), pos);
}

//-------------------------------------------------------------------------------------------------
void unsetBitInMask(ullint* va, const size_t pos) {
  const size_t acc_elem = pos / ullint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ullint_bit_count_zu);
  va[acc_elem] &= ~(0x1LLU << acc_bit);
}
  
//-------------------------------------------------------------------------------------------------
void unsetBitInMask(std::vector<ullint> *va, const size_t pos) {
  if (pos >= va->size() * ullint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned long long integer bit "
          "mask of length " + std::to_string(va->size()) + ".", "unsetBitInMask");
  }
  unsetBitInMask(va->data(), pos);
}

//-------------------------------------------------------------------------------------------------
void unsetBitInMask(ushort* va, const size_t pos) {
  const size_t acc_elem = pos / ushort_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ushort_bit_count_zu);
  va[acc_elem] &= ~(0x1 << acc_bit);
}
  
//-------------------------------------------------------------------------------------------------
void unsetBitInMask(std::vector<ushort> *va, const size_t pos) {
  if (pos >= va->size() * ushort_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned short integer bit "
          "mask of length " + std::to_string(va->size()) + ".", "unsetBitInMask");
  }
  unsetBitInMask(va->data(), pos);
}

//-------------------------------------------------------------------------------------------------
int readBitFromMask(const uint* va, const size_t pos) {
  const size_t acc_elem = pos / uint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * uint_bit_count_zu);
  return ((va[acc_elem] >> acc_bit) & 0x1);
}
  
//-------------------------------------------------------------------------------------------------
int readBitFromMask(const std::vector<uint> &va, const size_t pos) {
  if (pos >= va.size() * uint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned integer bit mask of "
          "length " + std::to_string(va.size()) + ".", "readBitFromMask");
  }
  return readBitFromMask(va.data(), pos);
}

//-------------------------------------------------------------------------------------------------
int readBitFromMask(const ullint* va, const size_t pos) {
  const size_t acc_elem = pos / ullint_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ullint_bit_count_zu);
  return ((va[acc_elem] >> acc_bit) & 0x1LLU);
}
  
//-------------------------------------------------------------------------------------------------
int readBitFromMask(const std::vector<ullint> &va, const size_t pos) {
  if (pos >= va.size() * ullint_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned long long integer bit "
          "mask of length " + std::to_string(va.size()) + ".", "readBitFromMask");
  }
  return readBitFromMask(va.data(), pos);
}

//-------------------------------------------------------------------------------------------------
int readBitFromMask(const ushort* va, const size_t pos) {
  const size_t acc_elem = pos / ushort_bit_count_zu;
  const int acc_bit  = pos - (acc_elem * ushort_bit_count_zu);
  return ((va[acc_elem] >> acc_bit) & 0x1);
}
  
//-------------------------------------------------------------------------------------------------
int readBitFromMask(const std::vector<ushort> &va, const size_t pos) {
  if (pos >= va.size() * ushort_bit_count_zu) {
    rtErr("Position " + std::to_string(pos) + " is invalid for an unsigned short integer bit mask "
          "of length " + std::to_string(va.size()) + ".", "readBitFromMask");
  }
  return readBitFromMask(va.data(), pos);
}

} // namespace stmath
} // namespace stormm
