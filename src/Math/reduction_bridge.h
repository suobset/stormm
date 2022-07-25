// -*-c++-*-
#ifndef OMNI_REDUCTION_BRIDGE_H
#define OMNI_REDUCTION_BRIDGE_H

#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"

namespace omni {
namespace math {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;

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

} // namespace math
} // namespace omni

#endif
