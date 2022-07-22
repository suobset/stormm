// -*-c++-*-
#ifndef OMNI_REDUCTION_H
#define OMNI_REDUCTION_H

#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "reduction_workunit.h"

namespace omni {
namespace math {

using math::addScalarToVector;
using math::sum;
using numerics::max_llint_accumulation;
using numerics::splitRealAccumulation;
using numerics::splitRealConversion;
  
/// \brief Define the abstract length for a ReductionWorkUnit's abstract
constexpr int rdwu_abstract_length = 8;

/// \brief Enumerate stages of a reduction operation, which can be taken piecemeal or all at once.
enum class ReductionStage {
  GATHER,     ///< Gather information from many particles or other sources.
  SCATTER,    ///< Scatter the accumulated information (this can also include a secondary gather
              ///<   on block-wide accumulated results).
  ALL_REDUCE  ///< Perform the entire reduction in one step, when no intermediate accumulation is
              ///<   needed.
};

/// \brief Indicate the action to take when scattering the fully gathered sum (which is always a
///        sum) to the writeable array positions.  If further manipulations are needed, i.e. the
///        gathered sum is to be multiplied by two and then used as a divisor, those must be done
///        explicitly on the buffered values at the end of a gathering operation, so that the
///        scatter can then be applied in a subsequent step.
enum class ReductionGoal {
  NORMALIZE,      ///< Compute the norm of the data, whether the norm of the vector of scalar
                  ///<   values present at all atoms or the norm of the collection of forces
                  ///<   acting on each atom.
  CENTER_ON_ZERO  ///< Compute the means of the data in up to three dimensions and subtract these
                  ///<   means from all values.
};
  
/// \brief Indicate whether parallel reductions can be performed with one work unit per system.
enum class RdwuPerSystem {
  ONE,      ///< Each system's gather and scatter operations are performed by a single work unit,
            ///<   making one-kernel all-reduce possible.
  MULTIPLE  ///< One or more systems' gather and scatter operations are shared across multiple
            ///<   work units, forcing separate gather and scatter kernels.
};

/// \brief Enumerate the features of a reduction work unit abstract, as encoded in the synthesis
///        of topologies.
enum class RdwuAbstractMap {
  ATOM_START,    ///< Absolute starting point of atoms in various (concatenated, in the case of a
                 ///<   synthesis of many systems) arrays
  ATOM_END,      ///< Absolute upper limit of atoms to read or write in various arrays
  RESULT_INDEX,  ///< Location in one or more holding arrays for intermediate results of gathering
  DEPN_START,    ///< Start of dependencies in one or more holding arrays for making a final
                 ///<   assembly of the results prior to scattering
  DEPN_END,      ///< Upper limit of dependencies--if this is only one greater than the start, it
                 ///<   means that the reduction can be done by a single work unit
  SYSTEM_ID      ///< Index of the system to which this reduction pertains
};

/// \brief Collect the simple components needed to guide reductions across all systems in a
///        topology synthesis (or the corresponding compilation of coordinates): the number of
///        work units, the strategy, and the list of work unit abstracts.
struct ReductionKit {

  /// \brief The constructor takes a straight list of values for all member variables.
  ReductionKit(int nrdwu_in, RdwuPerSystem rps_in, const int* rdwu_abstracts_in,
               const int* atom_counts_in);

  /// \brief Take the typical copy and move constructors for an abstract with constants.
  /// \{
  ReductionKit(const ReductionKit &original) = default;
  ReductionKit(ReductionKit &&original) = default;
  /// \}
  
  const int nrdwu;            ///< The number of reduction work units in the synthesis
  const RdwuPerSystem rps;    ///< Are there one or multiple work units serving each system?
  const int* rdwu_abstracts;  ///< Reduction work unit abstracts (strides of rdwu_abstract_length)
  const int* atom_counts;     ///< The number of atoms per system (for normalization purposes)
};

/// \brief Perform gathering operations for a normalization reduction.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
double gatherNormalization(const ReductionSubstrate<llint> rsbs, const int start_pos,
                           const int end_pos);

template <typename T>
double gatherNormalization(const ReductionSubstrate<T> rsbs, int start_pos, int end_pos);
/// \}
  
/// \brief Perform gathering operations for a centering-on-zero reduction.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
double3 gatherCenterOnZero(const ReductionSubstrate<llint> rsbs, int start_pos, int end_pos);

template <typename T>
double3 gatherCenterOnZero(const ReductionSubstrate<T> rsbs, int start_pos, int end_pos);
/// \}

/// \brief Perform scattering operations for a normalization reduction: divide the values found
///        in the writeable arrays by the magnitude of the vector represented by all values in
///        the readable arrays (computed gatherNormalization).
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param tsum       Sum of the vector components' squares (presumably the sum of the squares of
///                   values about to be divided by its square root)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
void scatterNormalization(ReductionSubstrate<llint> rsbs, double tsum, int start_pos, int end_pos);

template <typename T>
void scatterNormalization(ReductionSubstrate<T> rsbs, double tsum, int start_pos, int end_pos);
/// \}

/// \brief Perform scattering operations for a center-on-zero reduction of data in up to three
///        dimensions.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param tsum_x     Sum of the vector components in the 1st dimension
/// \param tsum_y     Sum of the vector components in the 2nd dimension
/// \param tsum_z     Sum of the vector components in the 3rd dimension
/// \param natom      Total number of atoms in the system (for taking the system-wide average)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
void scatterCenterOnZero(ReductionSubstrate<llint> rsbs, double tsum_x, double tsum_y,
                         double tsum_z, int natom, int start_pos, int end_pos);

template <typename T>
void scatterCenterOnZero(ReductionSubstrate<T> rsbs, double tsum_x, double tsum_y, double tsum_z,
                         int natom, int start_pos, int end_pos);
/// \}

/// \brief Evaluate a reduction.  In CPU code, this takes advantage of the fact that the
///        ReductionKit contains information about whether storage of the gathered results is
///        necessary and will recursively call itself to accomplish the all-reduce.  By default,
///        this function will perform normalizations across each system in all provided data
///        vectors.
///
/// \param rsbs     Data to reduce
/// \param redk     Instructions for the reduction operations to perform
/// \param process  The stage of reduction to perform (if an all-reduce is requested but must be
///                 carried out in stages, the function can call itself recursively)
template <typename T>
void evalReduction(ReductionSubstrate<T> *rsbs, const ReductionKit &redk,
                   const ReductionStage process = ReductionStage::ALL_REDUCE,
                   const ReductionGoal purpose = ReductionGoal::NORMALIZE);

} // namespace math
} // namespace omni

#include "reduction.tpp"

#endif
