// -*-c++-*-
#ifndef OMNI_VECTOR_OPS_H
#define OMNI_VECTOR_OPS_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "Cuda/hybrid.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Reporting/error_format.h"
#include "statistics.h"

namespace omni {
namespace math {

using cuda::Hybrid;
using data_types::getHpcVectorTypeSize;
using data_types::isHpcVectorType;
using data_types::isFloatingPointScalarType;
using data_types::isScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;

/// \brief Check that two vectors of the same data type are compatible for various arithmetic
///        comparisons.  Throws and exception if the vectors cannot be compared.
///
/// Overloaded:
///   - Check two Standard Template Library vectors
///   - Check two Hybrid objects
///
/// \param va  The first vector
/// \param vb  The second vector
/// \{
template <typename T> void vectorComparisonCheck(const std::vector<T> &va,
                                                 const std::vector<T> &vb, const char* caller);
template <typename T> void vectorComparisonCheck(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                 const char* caller);
/// \}

/// \brief Maximum absolute deviation between two vectors
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> T maxAbsoluteDifference(const T* va, const T* vb, size_t length);
template <typename T> T maxAbsoluteDifference(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> T maxAbsoluteDifference(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Maximum relative deviation between two vectors.  The signed deviation with the largest
///        absolute value will be returned.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double maxRelativeDifference(const T* va, const T* vb, size_t length);
template <typename T> double maxRelativeDifference(const std::vector<T> &va,
                                                   const std::vector<T> &vb);
template <typename T> double maxRelativeDifference(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the mean unsigned error between two vectors.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double meanUnsignedError(const T* va, const T* vb, size_t length);
template <typename T> double meanUnsignedError(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double meanUnsignedError(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the root mean squared error between two vectors, normalized by the mean value
///        of the reference vector.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double relativeRmsError(const T* va, const T* vb, size_t length);
template <typename T> double relativeRmsError(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double relativeRmsError(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the mean value of a vector of familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double mean(const T* va, const size_t length);
template <typename T> double mean(const std::vector<T> &va);
template <typename T> double mean(const Hybrid<T> &va);
/// \}

/// \brief Compute the mean value of a vector of familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param method  Statistical method for computing a variance-related quantity
/// \{
template <typename T> double variance(const T* va, size_t length, VarianceMethod method);
template <typename T> double variance(const std::vector<T> &va, VarianceMethod method);
template <typename T> double variance(const Hybrid<T> &va, VarianceMethod method);
/// \}

/// \brief Compute the mean value of a vector of two-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double2 bivariateMean(const T* va, size_t length);
template <typename T> double2 bivariateMean(const std::vector<T> &va);
template <typename T> double2 bivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Compute the mean value of a vector of three-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double3 trivariateMean(const T* va, size_t length);
template <typename T> double3 trivariateMean(const std::vector<T> &va);
template <typename T> double3 trivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Compute the mean value of a vector of three-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double4 quadrivariateMean(const T* va, size_t length);
template <typename T> double4 quadrivariateMean(const std::vector<T> &va);
template <typename T> double4 quadrivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Obtain the maximum value of a vector of scalar numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T maxValue(const T* va, size_t length);
template <typename T> T maxValue(const std::vector<T> &va);
template <typename T> T maxValue(const Hybrid<T> &va);
/// \}

/// \brief Obtain the minimum value of a vector of scalar numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T minValue(const T* va, size_t length);
template <typename T> T minValue(const std::vector<T> &va);
template <typename T> T minValue(const Hybrid<T> &va);
/// \}

/// \brief Obtain the largest absolute value of a vector of scalar numbers
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T maxAbsValue(const T* va, size_t length);
template <typename T> T maxAbsValue(const std::vector<T> &va);
template <typename T> T maxAbsValue(const Hybrid<T> &va);
/// \}
  
/// \brief Compute the Pearson correlation coefficient for two vectors of numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The first vector
/// \param vb      The second vector
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double pearson(const T* va, const T* vb, size_t length);
template <typename T> double pearson(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double pearson(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Add a constant to all elements of a vector.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param inc     The amount by which to alter the vector (passed in by value with an implicit
///                conversion to the same type as the vector itself)
/// \{
template <typename T> void addScalarToVector(T* va, size_t length, T inc);
template <typename T> void addScalarToVector(std::vector<T> *va, T inc);
template <typename T> void addScalarToVector(Hybrid<T> *va, T inc);
/// \}

/// \brief Compute the cross product of two vectors.  Each should be of dimension 3, and if any of
///        the vectors have more than 3 elements, only the first three will be used.  If supplying
///        vectors, no check is made to ensure the float or double nature of the type (that's up to
///        the developer, and allows for corner cases in which two vectors of ints could produce a
///        vector of other ints, provided that the format is not overflowed).  If supplying HPC 
///        tuples, the compiler will catch cases of bad input where two-element tuples are
///        supplied, but four-element tuples will act like vectors longer than three elements.
///
/// Overloaded:
///   - Operate on C-style arrays of trusted length
///   - Operate on Standard Template Library vectors
///   - Operate on three-element HPC tuples
///
/// \param va  Compute va x vb = vc.
/// \param vb  Compute va x vb = vc.
/// \param vc  Pre-allocated third vector to hold the result.  Compute va x vb = vc.
/// \{
template <typename T> void crossProduct(const T* va, const T* vb, T* vc);
template <typename T> void crossProduct(const std::vector<T> &va, const std::vector<T> &vb,
                                        std::vector<T> *vc);
template <typename T> T crossProduct(const T va, const T vb);
/// \}

/// \brief Compute the magnitude of a vector.  The result is always returned as a double-precision
///        value, and accumulated as such, but individual elements are not converted to double
///        prior to taking their squares.  This can overflow integer data types if used improperly.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  Length of va
/// \{
template <typename T> double magnitude(const T* va, const size_t length);
template <typename T> double magnitude(const std::vector<T> &va);
template <typename T> double magnitude(const Hybrid<T> &va);
/// \}

/// \brief Compute the dot product of two vectors.  The result is always returned as a
///        double-precision value, and accumulated as such, but individual elements are not
///        converted to double prior to multiplication.  This can overflow integer data types if
///        used improperly.
///
/// Overloaded:
///   - Operate on a pair of C-style arrays of trusted length
///   - Operate on a pair of Standard Template Library vectors
///   - Operate on a pair of Hybrid objects
///
/// \param va      The first vector of interest
/// \param vb      The second vector of interest
/// \param length  The first vector of interest
/// \{
template <typename T> double dot(const T* va, const T* vb, const size_t length);
template <typename T> double dot(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double dot(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}
  
/// \brief Compute the projection of one vector onto another.
///
/// Overloaded:
///   - Operate on C-style arrays of trusted length
///   - Operate on Standard Template Library vectors
///   - Operate on a pair of Hybrid objects
///   - Return a new array or work with a provided, pre-existing one (only works with vector
///     objects, not raw pointers)
///
/// \param va  The projection of va onto vb will be computed 
/// \param vb  The projection of va onto vb will be computed 
/// \param vc  Store the result in a pre-allocated space vc
/// \{
template <typename T> void project(const T* va, const T* vb, T* vc, size_t length);
template <typename T> void project(const std::vector<T> &va, const std::vector<T> &vb,
                                   std::vector<T> *vc);
template <typename T> void project(const Hybrid<T> &va, const Hybrid<T> &vb, Hybrid<T> *vc);
template <typename T> std::vector<T> project(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> Hybrid<T> project(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

} // namespace math
} // namespace omni

#include "vector_ops.tpp"

#endif
