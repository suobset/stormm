// -*-c++-*-
#ifndef OMNI_SUMMATION_H
#define OMNI_SUMMATION_H

#include <vector>
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "UnitTesting/approx.h"

namespace omni {
namespace math {

using card::Hybrid;
using card::HybridTargetLevel;
using errors::rtErr;
using testing::Approx;
using testing::ComparisonType;

/// \brief The type of prefix sum to compute.
enum class PrefixSumType {
  EXCLUSIVE, INCLUSIVE
};

// \brief Compute an in-place prefix sum over a vector.  For "capped" prefix sums, which are
///       exclusive sums followed by one extra entry containing the total sum over the entire
///       array, provide an array that is one element larger than the actual number of indices of
///       interest with a zero in the final element.
///
/// \param v       The vector to accumulate
/// \param style   Method for computing the sum
/// \param caller  Name of the calling function (optional)
/// \{
template <typename TSum, typename TBase>
void prefixSumInPlace(std::vector<TBase> *v, const PrefixSumType style,
                      const char* caller = nullptr);
/// \}

/// \brief Sum the elements of a standard template library vector for elements of known types.
///        Developers must take care that the sum type have a sufficient range, particularly in
///        cases with integral types.
///
/// Overloaded:
///   - Signed integer sums, all returning long long int
///   - Unsigned integer sums, all returning unsigned long long int
///   - Real number vector sums, all returning a double precision floating point number
///   - Summations for HPC vector types following the above conventions for each component
///
/// \param v       The vector to sum
/// \param vlen    The length of the vector
/// \param buffer  [HPC only] Accumulator for individual thread blocks
/// \param tier    [HPC only] Indicator of whether to perform the summation over data on the host
///                or on the device
/// \param gpu     [HPC only] Details of the GPU where the data resides
/// \{
template <typename TSum, typename TBase>
TSum sum(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sum(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sum(const Hybrid<TBase> &hb);

template <typename TSum, typename TBase>
TSum sumTuple2(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple2(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple2(const Hybrid<TBase> &hb, const HybridTargetLevel tier = HybridTargetLevel::HOST);

template <typename TSum, typename TBase>
TSum sumTuple3(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple3(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple3(const Hybrid<TBase> &hb, const HybridTargetLevel tier = HybridTargetLevel::HOST);

template <typename TSum, typename TBase>
TSum sumTuple4(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple4(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple4(const Hybrid<TBase> &hb, const HybridTargetLevel tier = HybridTargetLevel::HOST);
/// \}

} // namespace math
} // namespace omni

#include "summation.tpp"

#endif

