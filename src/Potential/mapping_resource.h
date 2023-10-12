// -*-c++-*-
#ifndef STORMM_MAPPING_RESOURCE_H
#define STORMM_MAPPING_RESOURCE_H

#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

using card::CoreKlManager;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::ExceptionResponse;

/// \brief The abstract of the MappingResource is intended to be always writeable; otherwise, like
///        other resource objects, it would have no practical use.
struct MappingResourceKit {

  /// \brief The constructor takes a list of settings for all member variables as inputs.
  MappingResourceKit(int max_blocks_in, int max_grid_points_in, PrecisionModel mode_in,
                     QMapMethod approach_in, int* overflow_in, double* dstash_in,
                     float* fstash_in);

  /// \brief As with other abstracts, the presence of const members implicitly deletes the copy and
  ///        move assignment operators, but with copy and move initialization these constructors
  ///        are still valid and can take their default forms.
  /// \{
  MappingResourceKit(const MappingResourceKit &original) = default;
  MappingResourceKit(MappingResourceKit &&original) = default;
  /// \}

  const int max_blocks;       ///< Maximum number of blocks that can use the resource at one time
  const int max_grid_points;  ///< The maximum number of grid points per block
  const PrecisionModel mode;  ///< The precision mode for which real memory space is prepared.
                              ///<   This indicates whether fstash or dstash is valid.  If 
                              ///<   the 
  const QMapMethod approach;  ///< The manner of density mapping for which this object is prepared.
                              ///<   only ACC_REGISTER and ACC_SHARED will successfully construct
                              ///<   the underlying MappingResource.
  int* overflow;              ///< Overflow space for block-specific __shared__ accumulators.
                              ///<   This is analogous to, but not the same as, the overflow space
                              ///<   allocated in a PMIGrid object.  This can be useful to the
                              ///<   kernels associated with ACC_SHARED, if overflow is occurring.
  double* dstash;             ///< Space to temporarily stash float64_t numbers.  This is
                              ///<   useful to the double-precision kernels associated with
                              ///<   ACC_REGISTER.
  float* fstash;              ///< Space to temporarily stash float32_t numbers.  This is
                              ///<   useful to the single-precision kernels associated with
                              ///<   ACC_REGISTER.
};
  
/// \brief A cognate of the CacheResource kit for valence interactions, this provides space for
///        each of the blocks to utilize a private segment of __global__ memory for temporary
///        caching needs in a write > __syncthreads() > read manner.
class MappingResource {
public:

  /// \brief The MappingResource is tailored to serve a particular PMIGrid.  Its precision model
  ///        will be set based on that reference object and its allocations based off of the work
  ///        unit configuration.
  ///
  /// \param pm      The particle-mesh interaction grids which kernels will accumulate
  /// \param gpu     Specifications of the GPU that will perform the calculations
  /// \param policy  Indicate the path to take if the resource is created when unneeded
  /// \{
  MappingResource(const PMIGrid *pm, const CoreKlManager &launcher,
                  ExceptionResponse policy = ExceptionResponse::WARN);
  
  MappingResource(const PMIGrid &pm, const CoreKlManager &launcher,
                  ExceptionResponse policy = ExceptionResponse::WARN);
  /// \}

  /// \brief Get the number of blocks that the object is prepared to serve
  int getBlockCount() const;

  /// \brief Get the maximum number of grid points that each block can manage at any one time
  int getGridPointsPerBlock() const;

  /// \brief Get the object's precision mode
  PrecisionModel getMode() const;

  /// \brief Get the density mapping approach that this resource enables
  QMapMethod getApproach() const;

  /// \brief Get the object's writeable abstract.
  ///
  /// \param tier  Target pointers to data on the CPU host or GPU device
  MappingResourceKit data(HybridTargetLevel tier = HybridTargetLevel::HOST);
#ifdef STORMM_USE_HPC
  /// \brief Upload all data in the object to the GPU device (for experimental or testing purposes)
  void upload();

  /// \brief Download all data in the object from the GPU device, i.e. for debugging purposes
  void download();
#endif
private:
  int block_limit;               ///< Number of thread blocks for which the resource is prepared
  int grid_points_per_block;     ///< Upper limit of grid points per thread block.  If any block
                                 ///<  works on more than this number of grid points, the memory
                                 ///<  for different blocks will begin to overlap or possibly
                                 ///<  overrun the length of data arrays (see below).
  PrecisionModel mode;           ///< The precision model for which resources are prepared.  This
                                 ///<   determines whether double_buffer or float_buffer is
                                 ///<   allocated, but is relevant only when the approach is
                                 ///<   ACC_REGISTER.
  QMapMethod approach;           ///< Density mapping approach for which resources are prepared
  Hybrid<int> overflow;          ///< Space for overflow bits in fixed-precision accumulation
                                 ///<   required by kernels launched under the ACC_SHARED mapping
                                 ///<   protocol
  Hybrid<double> double_buffer;  ///< Temporary storage space for buffering double-precision values
                                 ///<   accumulated under the ACC_REGISTER mapping protocol
  Hybrid<float> float_buffer;    ///< Temporary storage space for buffering single-precision values
                                 ///<   accumulated under the ACC_REGISTER mapping protocol
};

} // namespace energy
} // namespace stormm

#endif
