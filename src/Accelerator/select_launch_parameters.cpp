#include "select_launch_parameters.h"

namespace omni {
namespace card {

using energy::queryValenceKernelRequirements;
using energy::queryNonbondedKernelRequirements;

//-------------------------------------------------------------------------------------------------
KernelManager selectLaunchParameters(const GpuDetails &gpu, const AtomGraphSynthesis &poly_ag) {
  KernelManager result(gpu);
#ifdef OMNI_USE_HPC
  queryValenceKernelRequirements(&result, poly_ag);
  queryNonbondedKernelRequirements(&result);
#endif
  return result;
}

} // namespace card
} // namespace omni
