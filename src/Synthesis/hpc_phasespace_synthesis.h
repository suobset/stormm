// -*-c++-*-
#ifndef OMNI_HPC_PHASE_SPACE_SYNTHESIS_H
#define OMNI_HPC_PHASE_SPACE_SYNTHESIS_H

#include "Accelerator/gpu_details.h"
#include "Trajectory/trajectory_enumerators.h"
#include "phasespace_synthesis.h"

namespace omni {
namespace synthesis {

using card::GpuDetails;
using trajectory::TrajectoryKind;

/// \brief Launch the above kernel to upload scattered data for specific systems within a
///        PhaseSpaceSynthesis object.
///
/// \param devc_view   Collection of pointers to PhaseSpaceSynthesis data on the device
/// \param host_view   Collection of pointers to host mapped data, but visible on the device
/// \param low_index   Lower bound of systems to upload
/// \param high_index  Upper bound of systems to upload (the range is [low_index, high_index))
void systemUploader(PsSynthesisWriter *devc_view, PsSynthesisWriter *host_view,
                    TrajectoryKind kind, int low_index, int high_index, const GpuDetails &gpu);

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
