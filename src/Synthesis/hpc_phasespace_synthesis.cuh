// -*-c++-*-
#ifndef STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH
#define STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH

#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "phasespace_synthesis.h"

namespace stormm {
namespace synthesis {

using constants::large_block_size;

/// \brief Transfer a subset of system coordinates (and perhaps box dimensions) data from specific
///        systems within a PhaseSpaceSynthesis object betwee the host and the accelerator device.
///        PhaseSpaceSynthesis is a central object within STORMM and can be huge (multiple
///        gigabytes worth of coordinate, velocity, and force data.  As such, it deserves a
///        dedicated kernel for managing data transfer to and from the host.
///
/// \param destination  Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param source       Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param low_index    Lower bound of systems to upload
/// \param high_index   Upper bound of systems to upload (the range is [low_index, high_index))
  __global__ void __launch_bounds__(large_block_size, 1)
kSystemTransfer(PsSynthesisWriter destination, PsSynthesisWriter source, int low_index,
                int high_index, const TrajectoryKind material);

/// \brief Initialize forces for one or all systems of a PhaseSpaceSynthesis on the GPU.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis, containing system limits and
///               pointers to all coordinates and forces
/// \param index  Index of the system to initialize; if negative, all systems will be initialized.
__global__ void __launch_bounds__(large_block_size, 1)
kPsyInitializeForces(PsSynthesisWriter psyw, const int index);

/// \brief Initialize critical buffers in the phase space (specifically, prior coordinates and
///        velocities) which would otherwise not be used in energy minimization calculations.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis
__global__ void __launch_bounds__(large_block_size, 1)
kPsyPrimeConjugateGradient(PsSynthesisWriter psyw);
  
} // namespace synthesis
} // namespace stormm

#endif
