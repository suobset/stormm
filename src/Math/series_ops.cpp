#include <algorithm>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "series_ops.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const int* number_series, const size_t length,
                                        const int output_size) {
  const int actual_size = (output_size == 0LLU) ? std::max(maxValue(number_series, length), 1) :
                                                  output_size;
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (actual_size + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);

  // Loop over all polar hydrogens and check the appropriate bits
  for (size_t i = 0; i < length; i++) {
    const int mask_idx = number_series[i] / n_bits;
    const int bit_idx = number_series[i] - (mask_idx * n_bits);
    result[mask_idx] |= (0x1 << bit_idx);
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const std::vector<int> &number_series,
                                        const int output_size) {
  return numberSeriesToBitMask(number_series.data(), number_series.size(), output_size);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const Hybrid<int> &number_series, const int output_size) {
  return numberSeriesToBitMask(number_series.data(), number_series.size(), output_size);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> getSubsetIndexPattern(std::vector<int> *x_subset, const int x_size,
                                       int *n_subset_indices) {

  // Protect against empty data
  const size_t nidx = x_subset->size();
  if (nidx == 0) {
    *n_subset_indices = 0;
    return std::vector<int>();
  }
  
  // Obtain the maximum index
  int* xsbs_ptr = x_subset->data();
  const size_t subset_size = x_subset->size();
  std::vector<bool> index_found(x_size, false);
  for (size_t i = 0; i < nidx; i++) {
    index_found[xsbs_ptr[i]] = true;
  }
  int new_index_counter = 0;
  std::vector<int> result(x_size, -1);
  for (size_t i = 0; i < x_size; i++) {
    if (index_found[i]) {
      result[i] = new_index_counter;     
      new_index_counter++;
    }
  }
  if (n_subset_indices != nullptr) {
    *n_subset_indices = new_index_counter;
  }
  for (size_t i = 0; i < nidx; i++) {
    xsbs_ptr[i] = result[xsbs_ptr[i]];
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<int> enumerateMask(const std::vector<bool> &bitmask) {
  std::vector<int> result;
  int nchk = 0;
  const size_t nval = bitmask.size();
  for (size_t i = 0; i < nval; i++) {
    nchk += static_cast<int>(bitmask[i]);
  }
  result.resize(nchk);
  nchk = 0;
  for (size_t i = 0; i < nval; i++) {
    if (bitmask[i]) {
      result[nchk] = i;
      nchk++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> enumerateMask(const uint* bitmask, const int length) {
  std::vector<int> result;
  int nchk = 0;
  for (int i = 0; i < length; i++) {
    const int seg = (i >> 5);
    const int pos = i - (seg << 5);
    nchk += ((bitmask[seg] >> pos) & 0x1);
  }
  result.resize(nchk);
  nchk = 0;
  for (int i = 0; i < length; i++) {
    const int seg = (i >> 5);
    const int pos = i - (seg << 5);
    if ((bitmask[seg] >> pos) & 0x1) {
      result[nchk] = i;
      nchk++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> enumerateMask(const std::vector<uint> &bitmask) {
  return enumerateMask(bitmask.data(), bitmask.size() * 32);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> enumerateMask(const Hybrid<uint> &bitmask) {
  return enumerateMask(bitmask.data(), bitmask.size() * 32);
}

} // namespace stmath
} // namespace stormm
