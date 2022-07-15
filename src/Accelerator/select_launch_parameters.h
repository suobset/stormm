// -*-c++-*-
#ifndef OMNI_SELECT_LAUNCH_PARAMETERS_H
#define OMNI_SELECT_LAUNCH_PARAMETERS_H

#include "gpu_details.h"
#include "Potential/hpc_valence_potential.cuh"
#include "Potential/hpc_nonbonded_potential.cuh"
#include "Synthesis/atomgraph_synthesis.h"

namespace omni {
namespace card {

using synthesis::AtomGraphSynthesis;
  
/// \brief Given details of the GPU and the set of systems to simulate (the workload), query the
///        specifications of each version of the compiled kernels to determine the optimal launch
///        configurations.  Store the results in a KernelManager object as acquired wisdom.
///        This process amounts to constructing the object in a practical sense, but is done here,
///        with a free function in a separate implementation file, to allow the KernelManager
///        object to be reciprocally included in various HPC units to provide launch instructions
///        without generating circular dependencies.
///
/// \param gpu      Details of the GPU chosen for running calculations
KernelManager selectLaunchParameters(const GpuDetails &gpu);

} // namespace card
} // namespace omni

#endif
