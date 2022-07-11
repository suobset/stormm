#include "select_launch_parameters.h"

namespace omni {
namespace card {

using energy::queryValenceKernelRequirements;
using energy::queryNonbondedKernelRequirements;

//-------------------------------------------------------------------------------------------------
KernelManager selectLaunchParameters() {
  KernelManager result;
#ifdef OMNI_USE_HPC
  queryValenceKernelRequirements(&result);
#endif
  return result;
}

} // namespace card
} // namespace omni
