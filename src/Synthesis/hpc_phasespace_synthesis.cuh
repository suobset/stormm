// -*-c++-*-
#ifndef OMNI_HPC_PHASE_SPACE_SYNTHESIS_CUH
#define OMNI_HPC_PHASE_SPACE_SYNTHESIS_CUH

#include <vector>
#ifdef OMNI_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "Accelerator/gpu_details.h"
#include "Trajectory/trajectory_enumerators.h"
#include "Trajectory/phasespace_synthesis.h"

namespace omni {
namespace synthesis {

using card::GpuDetails;
using trajectory::PsSynthesisWriter;

/// \brief Upload scattered data from specific systems within a PhaseSpaceSynthesis object.  The
///        PhaseSpaceSynthesis is a central object within OMNI and can be huge (multiple gigabytes
///        worth of coordinate, velocity, and force data.  As such, it deserves dedicated kernels
///        for managing data transfer to and from the host.
///
/// \param devc_view   Collection of pointers to PhaseSpaceSynthesis data on the device
/// \param host_view   Collection of pointers to host mapped data, but visible on the device
/// \param low_index   Lower bound of systems to upload
/// \param high_index  Upper bound of systems to upload (the range is [low_index, high_index))
__global__ void __launch_bounds__(large_block_size, 1)
kSystemUploader(PsSynthesisWriter devc_view, PsSynthesisWriter host_view, int low_index,
                int high_index);

/// \brief Launch the above kernel to upload scattered data for specific systems within a
///        PhaseSpaceSynthesis object.
///
/// \param devc_view   Collection of pointers to PhaseSpaceSynthesis data on the device
/// \param host_view   Collection of pointers to host mapped data, but visible on the device
/// \param low_index   Lower bound of systems to upload
/// \param high_index  Upper bound of systems to upload (the range is [low_index, high_index))
void systemUploader(PsSynthesisWriter *devc_view, PsSynthesisWriter *host_view,
                    TrajectoryKind kind, int low_index, int high_index, const GpuDetails &gpu);

/// \brief Download scattered data from specific systems within a PhaseSpaceSynthesis object.  The
///        PhaseSpaceSynthesis is a central object within OMNI and can be huge (multiple gigabytes
///        worth of coordinate, velocity, and force data.  As such, it deserves dedicated kernels
///        for managing data transfer to and from the host.
///
/// \param devc_view   Collection of pointers to PhaseSpaceSynthesis data on the device
/// \param host_view   Collection of pointers to host mapped data, but visible on the device
/// \param low_index   Lower bound of systems to download
/// \param high_index  Upper bound of systems to download (the range is [low_index, high_index))
__global__ void __launch_bounds__(large_block_size, 1)
kSystemDownloader(PsSynthesisWriter devc_view, PsSynthesisWriter host_view, int low_index,
                  int high_index);

/// \brief Launch the above kernel to download scattered data for specific systems within a
///        PhaseSpaceSynthesis object.
///
/// \param devc_view   Collection of pointers to PhaseSpaceSynthesis data on the device
/// \param host_view   Collection of pointers to host mapped data, but visible on the device
/// \param low_index   Lower bound of systems to download
/// \param high_index  Upper bound of systems to download (the range is [low_index, high_index))
void systemDownloader(PsSynthesisWriter *devc_view, PsSynthesisWriter *host_view,
                      TrajectoryKind kind, int low_index, int high_index, const GpuDetails &gpu);

} // namespace synthesis
} // namespace omni

#endif
