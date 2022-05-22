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

/// \brief Launch the eponymous kernel to upload scattered data for specific systems within a
///        PhaseSpaceSynthesis object.
///
/// \param destination  Collection of pointers to PhaseSpaceSynthesis data on the host or device.
///                     If on the host, these poiniters must be to host-mapped data visible by the
///                     device.
/// \param source       Collection of pointers to PhaseSpaceSynthesis data on the host or device
///                     If on the host, these poiniters must be to host-mapped data visible by the
///                     device.
/// \param low_index    Lower bound of systems to upload
/// \param high_index   Upper bound of systems to upload (the range is [low_index, high_index))
/// \param gpu    Details of the GPU in use
void systemTransfer(PsSynthesisWriter *destination, PsSynthesisWriter *source,
                    TrajectoryKind kind, int low_index, int high_index, const GpuDetails &gpu);

/// \brief Launch the eponymous kernel to initialize forces on the GPU in a PhaseSpaceSynthesis
///        object.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis, containing system limits and
///               pointers to all coordinates and forces
/// \param index  Index of the system to initialize; if negative, all systems will be initialized.
/// \param gpu    Details of the GPU in use
void psyInitializeForces(PsSynthesisWriter *psyw, int index, const GpuDetails &gpu);
  
} // namespace synthesis
} // namespace omni

#endif
