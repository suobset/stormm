// -*-c++-*-
#include <vector>
#ifdef OMNI_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "Reporting/error_format.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_phasespace_synthesis.h"
#include "hpc_phasespace_synthesis.cuh"

namespace omni {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSystemUploader(PsSynthesisWriter devc_view, PsSynthesisWriter host_view, const int low_index,
                const int high_index) {

  // Assume that the bounds arrays for each system have been copied to the device, where they are
  // visible with the device view pointers.
  
}
  
//-------------------------------------------------------------------------------------------------
extern void systemUploader(PsSynthesisWriter *devc_view, PsSynthesisWriter *host_view,
                           const TrajectoryKind kind, const int low_index,
                           const int high_index, const GpuDetails &gpu) {
  kSystemUploader<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*devc_view, *host_view,
                                                                      low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSystemDownloader(PsSynthesisWriter devc_view, PsSynthesisWriter host_view, const int low_index,
                  const int high_index) {

  // Assume that the bounds arrays for each system have been copied to the device, where they are
  // visible with the device view pointers.
  
}
  
//-------------------------------------------------------------------------------------------------
extern void systemDownloader(PsSynthesisWriter *devc_view, PsSynthesisWriter *host_view,
                             const TrajectoryKind kind, const int low_index,
                             const int high_index, const GpuDetails &gpu) {
  kSystemDownloader<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*devc_view, *host_view,
                                                                        low_index, high_index);
}

} // namespace synthesis
} // namespace omni
