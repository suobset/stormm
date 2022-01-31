// -*-c++-*-
#ifndef OMNI_STATISTICS_H
#define OMNI_STATISTICS_H

namespace omni {
namespace math {

/// \brief Enumerate various statistical methods for computing the variance of a sequence of
///        n numbers (i.e. as presented in a std::vector).
enum class VarianceMethod {
  VARIANCE,                    // Basic variance, sum_i((x_i - <x>)^2)
  STANDARD_DEVIATION,          // Normalized root variance, sqrt(sum_i((x_i - <x>)^2) / (n - 1))
  ROOT_MEAN_SQUARED_DEVIATION, // Unnormalized root variance, sqrt(sum_i((x_i - <x>)^2) / n)
  COEFFICIENT_OF_VARIATION,    // Standard deviation divided by mean absolute value <|x|>
  NORMALIZED_RMSD              // Root mean squared deviation divided by mean absolute value <|x|>
};

/// \brief Compute the variance of a set of data as a running quantity based on just the power
///        series of the one- and two-moments of the data.  Different methods fed to this function
///        can also produce standard deviation, root mean squared deviation, and either quantity
///        normalized by the mean of the data.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
/// \param method          Statistical analysis method
double running_variance(double sum_of_squares, double sum_of_values, int sample_count,
                        VarianceMethod method = VarianceMethod::VARIANCE);

/// \brief Compute the standard deviation of a set of data based on its first two moments.  Returns
///        zero if there are less than two data points.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_stdev(double sum_of_squares, double sum_of_values, int sample_count);

/// \brief Compute the root mean squared deviation of a set of data based on its first two moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_rmsd(double sum_of_squares, double sum_of_values, int sample_count);

/// \brief Compute the mean-normalized standard deviation of a set of data based on its first two
///        moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_coefficient_of_variation(double sum_of_squares, double sum_of_values,
                                        int sample_count);

/// \brief Compute the normalized root mean squared deviation of a set of data based on its first
///        two moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_normalized_rmsd(double sum_of_squares, double sum_of_values, int sample_count);

} // namespace math
} // namespace omni

#endif
