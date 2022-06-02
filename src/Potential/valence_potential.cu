A// -*-c++-*-
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph_abstracts.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/valence_workunit.h"

namespace omni {
namespace energy {

using synthesis::maximum_valence_work_unit_atoms;

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ void splitForceContribution(const float fval, const float max_incr,
                                                       const int pos, int* sh_primary,
                                                       int* sh_overflow_active, int* overflow) {
  int ival;
  if (fabsf(fval) >= max_incr) {
    const int spillover = flocal / max_incr;
    ival = fval - (spillover * max_incr);
    atomicAdd(&overflow[pos], spillover);

    // No atomic needed as this just sets the value, which started as zero, to one
    sh_overflow_active[pos / warp_size_int] = 1;
  }
  else {
    ival = fval;
  }
  const int prim_old = atomicAdd(&primary[pos], ival);
  const int rollover = 2 * ((ival > 0 && prim_old + ival < prim_old) +
                            (ival < 0 && prim_old + ival > prim_old));
  if (rollover) {
    atomicAdd(&overflow[pos], rollover);
    sh_overflow_active[pos / warp_size_int] = 1;
  }
}
  
// Single-precision floating point definitions
#define TCALC float
#  define SQRT_FUNC sqrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define ABS_FUNC  fabsf

#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define KERNEL_NAME kfsValenceForceAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME  
#      define COMPUTE_ENERGY
#        define KERNEL_NAME kfsValenceForceEnergyAccumulation
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define KERNEL_NAME kfValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kfValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kfValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear single-precision floating point definitions
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef ABS_FUNC
#undef TCALC

// Double-precision floating point definitions
#define TCALC double
#  define SQRT_FUNC sqrt
#  define ACOS_FUNC acos
#  define COS_FUNC  cos
#  define ABS_FUNC  fabs

#  define COMPUTE_FORCE
#    define KERNEL_NAME kdValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kdValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear double-precision floating point definitions
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef ABS_FUNC
#undef TCALC

} // namespace energy
} // namespace omni

#endif
