#include "select_launch_parameters.h"

namespace omni {
namespace card {

using energy::queryValenceKernelRequirements;
using energy::queryNonbondedKernelRequirements;

//-------------------------------------------------------------------------------------------------
KernelManager selectLaunchParameters(const GpuDetails &gpu) {
  KernelManager result(gpu);
#ifdef OMNI_USE_HPC
  queryValenceKernelRequirements(&result);
  queryNonbondedKernelRequirements(&result);
#endif
  return result;
}

} // namespace card
} // namespace omni
