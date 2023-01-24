// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> incrementingSeries(const T start_value, const T end_value,
                                                        const T increment) {
  const T zero = 0;
  const T actual_increment = (end_value - start_value > zero) ?  std::abs(increment) :
                                                                -std::abs(increment);
  std::vector<T> result((end_value - start_value) / increment);
  T v = start_value;
  const size_t nval = result.size();
  for (size_t i = 0; i < nval; i++) {
    result[i] = v;
    v += actual_increment;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const T* original_values,
                                                          const size_t values_length,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  const int actual_length = (reduced_length == 0) ? maxValue(reduced_indices) : reduced_length;
  const size_t ridx_size = reduced_indices.size();
  if (values_length < ridx_size) {
    rtErr("The array of original values is smaller than the number of unique indices.",
          "extractIndexedValues");    
  }
  std::vector<T> result;
  for (size_t i = 0; i < ridx_size; i++) {
    if (reduced_indices[i] >= 0) {
      result[reduced_indices[i]] = original_values[i];
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const std::vector<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  return extractIndexedValues(original_values.data(), original_values.size(), reduced_indices,
                              reduced_length);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const Hybrid<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  return extractIndexedValues(original_values.data(), original_values.size(), reduced_indices,
                              reduced_length);
}

} // namespace math
} // namespace stormm
