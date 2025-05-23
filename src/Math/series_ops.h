// -*-c++-*-
#ifndef STORMM_SET_OPS_H
#define STORMM_SET_OPS_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "summation.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using parse::NumberFormat;
using parse::realToString;

/// \brief Create a number series based on a starting number, a final value, and an increment.
///        The result is inclusive of the starting value but not the ending value.
///
/// \param start_value  The first value in the list
/// \param end_value    The final value in the list
/// \param increment    The increment to the list (this will be reversed in sign as necessary to
///                     ensure that the series can progress from the starting value to the ending
///                     value.
template <typename T> std::vector<T> incrementingSeries(const T start_value, const T end_value,
                                                        const T increment = 1);
  
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
///        the values.  Upon completion, the original indexing will have been replaced by the new
///        indexing and a key indicating how the modified index values map to the original index
///        values is returned.
///
/// \param x_subset          A collection of indices, probably with incomplete usage of some
///                          larger list of possible indices.  Modified and returned.
/// \param x_size            The true size of the original array.  The subset may not include the
///                          maximum index, hence this size must be provided.
/// \param n_subset_indices  Pointer for returning the number of unique indices found in the
///                          subset
std::vector<int> getSubsetIndexPattern(std::vector<int> *x_subset, int x_size,
                                       int *n_subset_indices = nullptr);

/// \brief Extract values from a table into a (smaller) table that fits with a reduced indexing
///        scheme.  Indices of -1 are considered invalid and removed.  For example:
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
  
/// \brief Given a series of integers in the range [ 0, ... N ), count the quantities of each
///        integer.  Produce a new series with the indices of all elements in the original series
///        whose value is 0, followed by the indices of all elements whose value is 1, ..., N-1.
///        Fill a series of bounds [ 0, N ] (N + 1 elements in length) with the bounds of the
///        indexing array for locating 0, 1, ..., N-1 in the original array.
///
/// Overloaded:
///   - Provide C-style arrays with a trused length and maximum value 
///   - Provide Standard Template Library vectors
///   - Provide Hybrid objects
///
/// \param raw_values       Array of values, each element in the range [ 0, ..., N )
/// \param value_locations  Array of indices in raw_values wherein instances of each value 0,
///                         1, ... N-1 can be found, filled and returned
/// \param value_bounds     Bounds array for value_locations, filled and returned
/// \param value_count      Trusted length of the raw_values array (if C-style arrays are provided)
/// \param value_limit      One less than the trusted length of the value_bounds array (this is for
///                         checking purposes, or in special cases limiting the organization to
///                         values in a specific range).  This is the upper limit of values that
///                         will be catalogged: the range becomes [ 0 ... value_limit ).
/// \{
template <typename Tdata, typename Tloc>
void indexingArray(const Tdata* raw_values, Tloc* value_locations, Tloc* value_bounds,
                   size_t value_count, size_t value_limit);

template <typename Tdata, typename Tloc>
void indexingArray(const std::vector<Tdata> &raw_values, std::vector<Tloc> *value_locations,
                   std::vector<Tloc> *value_bounds, size_t value_limit = 0);

template <typename Tdata, typename Tloc>
void indexingArray(const Hybrid<Tdata> &raw_values, Hybrid<Tloc> *value_locations,
                   Hybrid<Tloc> *value_bounds, size_t value_limit = 0);
/// \}

/// \brief Translate a bitmask into a series of integers indicating which indices are checked 1 or
///        'true'.  Theis reverses the operation of numberSeriesToBitMask, above, in a more general
///        way.
///
/// Overloaded:
///   - Provide a Standard Template Library bool vector
///   - Provide a vector of unsigned integer bit masks, either as C-style array, a Standard
///     Template Library bool vector, or as a Hybrid object
///
/// \param bitmask  The original bit mask
/// \param length   The trusted length of the bit mask, if the bit mask is provided as a C-style
///                 array (this option is only available for unsigned integer bit masks, not
///                 bool* arrays).  The length should be the number of valid bits, each considered
///                 an element of the array.  The length should not be submitted as the number of
///                 unsigned integers in the array.
/// \{
std::vector<int> enumerateMask(const std::vector<bool> &bitmask);

std::vector<int> enumerateMask(const uint* bitmask, const int length);

std::vector<int> enumerateMask(const std::vector<uint> &bitmask);

std::vector<int> enumerateMask(const Hybrid<uint> &bitmask);
/// \}

  
} // namespace stmath
} // namespace stormm

#include "series_ops.tpp"

#endif
