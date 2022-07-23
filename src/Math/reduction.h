// -*-c++-*-
#ifndef OMNI_REDUCTION_H
#define OMNI_REDUCTION_H

#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "reduction_enumerators.h"
#include "reduction_workunit.h"

namespace omni {
namespace math {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using math::addScalarToVector;
using math::sum;
using numerics::max_llint_accumulation;
using numerics::splitRealAccumulation;
using numerics::splitRealConversion;
using synthesis::AtomGraphSynthesis;

/// \brief Allocate space for reduction operations to store temporary accumulations, bridging the
///        gap between gathering and scattering operations.
class ReductionBridge {
public:
  /// The constructor allocates space in all three buffers for a set amount of data.
  ///
  /// \param n_values  The number of intermediate values to store, most likely determined by the
  ///                  number of reduction work units
  ReductionBridge(size_t n_values);

  /// \brief Take the default copy and move constructors as well as assignment operators.
  /// \{
  ReductionBridge(const ReductionBridge &original) = default;
  ReductionBridge(ReductionBridge &&original) = default;
  ReductionBridge& operator=(const ReductionBridge &original) = default;
  ReductionBridge& operator=(ReductionBridge &&original) = default;
  /// \}
  
  /// Get pointers to one of the buffers.
  ///
  /// Overloaded:
  ///   - Get a const pointer to a const form of this object's data
  ///   - Get a non-const pointer to a non-const form of this object's data
  ///
  /// \param cdim  The "dimension" to obtain a pointer for
  /// \param tier  Obtain pointers to host or device data
  /// \{
  const double* getPointer(CartesianDimension cdim,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getPointer(CartesianDimension cdim, HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

private:
  Hybrid<double> x_buffer;  ///< Buffer for the first type of data (it could be data pertaining to
                            ///<   the Cartesian X dimension, or in another setting something like
                            ///<   the squared magnitude of all forces for Conjugate Gradient
                            ///<   energy minimization).
  Hybrid<double> y_buffer;  ///< Buffer for the second type of data
  Hybrid<double> z_buffer;  ///< Buffer for the third type of data
  Hybrid<double> storage;   ///< ARRAY-kind Hybrid object targeted by all of the preceding objects
};

/// \brief Collect the simple components needed to guide reductions across all systems in a
///        topology synthesis (or the corresponding compilation of coordinates): the number of
///        work units, the strategy, and the list of work unit abstracts.
struct ReductionKit {

  /// \brief The constructor takes a straight list of values for all member variables, or a
  ///        combination of familiar objects.
  ///
  /// \param poly_ag  Compilation of topologies with reduction work unit counts and abstracts
  /// \param tier     Level at which to obtain the data arrays (needed for passing to subsequent
  ///                 getter functions in the topology and coordinate syntheses)
  /// \{
  ReductionKit(int nrdwu_in, RdwuPerSystem rps_in, const int* rdwu_abstracts_in,
               const int* atom_counts_in);

  ReductionKit(const AtomGraphSynthesis &poly_ag,
               HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

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
