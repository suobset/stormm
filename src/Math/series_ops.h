// -*-c++-*-
#ifndef OMNI_SET_OPS_H
#define OMNI_SET_OPS_H

#include <vector>
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "Math/vector_ops.h"

namespace omni {
namespace math {

using card::Hybrid;

/// \brief Convert a number series into a bitmask, such that for a number series consisting of
///        integers { a, b, c, d, e } bits a, b, c, d, and e are checked off in a mask broken into
///        segments of 32 bits (or whatever the size of an unsigned integer happens to be).  The
///        number series is trusted to have no values less than zero.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a std::vector of original values
///   - Operate on a Hybrid object of original values
///
/// \param number_series  The series of numbers to map
/// \param length         Length of the number_series array (if C-style)
/// \param output_size    Output size to prepare (if set to zero or unspecified, the number series
///                       will first be scanned to determine the necessary size from its largest
///                       value)
/// \{
std::vector<uint> numberSeriesToBitMask(const int* number_series, const size_t length,
                                        int output_size = 0);

std::vector<uint> numberSeriesToBitMask(const std::vector<int> &number_series,
                                        int output_size = 0);

std::vector<uint> numberSeriesToBitMask(const Hybrid<int> &number_series, int output_size = 0);
/// \}

/// \brief Let a vector x_subset contain indices into a table of values, but (probably) not making
///        use of all indices, i.e. 0, 2, 3, 1, 1, 0, 5, 6, 9 with ten indices in all.  Create a
///        key for extracting from the original table the values which x_subset does index, which
///        can then guide the construction of a smaller table for which x_subset does use all of
///        the values.
///
/// \param x_subset          A collection of indices, probably with incomplete usage of some
///                          larger list of possible indices
/// \param n_subset_indices  Pointer for returning the number of unique indices found in the
///                          subset
std::vector<int> getSubsetIndexPattern(const std::vector<int> &x_subset,
                                       int *n_subset_indices = nullptr);

/// \brief Extract values from a table into a (smaller) table that fits with a reduced indexing
///        scheme.  For example:
///
///        Original       = { 0.1, 0.5, 0.9, 0.2, 0.6, 0.8, 0.7, 0.4 }
///        indexed_values = {  -1,   0,  -1,  -1,   1,   2,   3,  -1 }
///        result         = { 0.5, 0.6, 0.8, 0.7 }
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a std::vector of original values
///   - Operate on a Hybrid object of original values
///
/// \param original_values  The original table of indexed values
/// \param values_length    Length of the original_values array (if C-style)
/// \param reduced_indices  Table of size no larger than original_values indicating the indices
///                         that each value should have in a reduced table, which is then returned
/// \param reduced_length   Length of the reduced parameter table (if not supplied, it will be
///                         calculated based on the maximum value of the reduced_indices array)
/// \{
template <typename T> std::vector<T> extractIndexedValues(const T* original_values,
                                                          size_t values_length,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length = 0);

template <typename T> std::vector<T> extractIndexedValues(const std::vector<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length = 0);

template <typename T> std::vector<T> extractIndexedValues(const Hybrid<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length = 0);
/// \}
  
} // namespace math
} // namespace omni

#include "series_ops.tpp"

#endif
