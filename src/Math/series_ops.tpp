// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace math {

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

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const T* values, const size_t length) {
  double esum = 0.0;
  for (size_t i = 0LLU; i < length; i++) {
    if (values[i] <= (T)(0)) {
      rtErr("The logarithm of a negative number, or zero, is undefined.  " +
            realToString(values[i], 11, 4, NumberFormat::STANDARD_REAL) + " was encountered in "
            "position " + std::to_string(i) + ".", "logProduct");
    }
    esum += log(values[i]);
  }
  return esum;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const std::vector<T> &values) {
  return logProduct(values.data(), values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> double logProduct(const Hybrid<T> &values) {
  return logProduct(values.data(), values.size());
}

} // namespace math
} // namespace stormm
