// -*-c++-*-
#include "Constants/hpc_bounds.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph_abstracts.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/valence_workunit.h"

namespace omni {
namespace energy {

using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using synthesis::maximum_valence_work_unit_atoms;

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ void splitForceContribution(const float fval, const int pos,
                                                       int* sh_primary, int* sh_overflow_active,
                                                       int* overflow) {
  int ival;
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    ival = __float2int_rn(fval - ((float)(spillover) * max_int_accumulation_f));
    atomicAdd(&overflow[pos], spillover);

    // No atomic needed as this just sets the value, which started as zero, to one
    sh_overflow_active[pos / warp_size_int] = 1;
  }
  else {
    ival = __float2int_rn(fval);
  }
  const int prim_old = atomicAdd(&sh_primary[pos], ival);
  if ((ival > 0 && prim_old + ival < prim_old) || (ival < 0 && prim_old + ival > prim_old)) {
    atomicAdd(&overflow[pos],
              (1 - (2 * (ival < 0))) * 2 * ((ival > 0 && prim_old + ival < prim_old) +
                                            (ival < 0 && prim_old + ival > prim_old)));
    sh_overflow_active[pos / warp_size_int] = 1;
  }
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ int2 convertSplitFixedPrecision(const float fval) {
  int2 result = { 0, 0 };
  if (fabsf(fval) >= max_int_accumulation_f) {
    const int spillover = fval / max_int_accumulation_f;
    result.x = __float2int_rn(fval - ((float)(spillover) * max_int_accumulation_f));
    result.y = spillover;
  }
  else {
    result.x = __float2int_rn(fval);
    result.y = 0;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ void addSplitFixedPrecision(const int2 ival, const int pos,
                                                       int* sh_primary, int* sh_overflow_active,
                                                       int* overflow) {
  const int prim_old = atomicAdd(&sh_primary[pos], ival.x);
  if ((ival.x > 0 && prim_old + ival.x < prim_old) || (ival.x < 0 && prim_old + ival > prim_old)) {
    ival.y += (1 - (2 * (ival.x < 0))) * 2 * ((ival.x > 0 && prim_old + ival.x < prim_old) +
                                              (ival.x < 0 && prim_old + ival.x > prim_old));
  }
  if (ival.y) {
    atomicAdd(&overflow[pos], ival.y);
  }
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ int2 combineSplitFixedPrecision(const int2 a, const int2 b) {
  int2 result = { a.x + b.x, a.y + b.y };
  if ((b.x > 0 && result.x < a.x) || b.x < 0 && result.x > a.x) {
    result.y += (1 - (2 * (b.x < 0))) * 2 * ((b.x > 0 && result.x < a.x) +
                                             (b.x < 0 && result.x > a.x));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ int2 antiCombineSplitFixedPrecision(const int2 a, const int2 b) {
  int2 result = { a.x + b.x, a.y + b.y };
  if ((b.x > 0 && result.x < a.x) || b.x < 0 && result.x > a.x) {
    result.y += (1 - (2 * (b.x < 0))) * 2 * ((b.x > 0 && result.x < a.x) +
                                             (b.x < 0 && result.x > a.x));
  }
  result.x = -result.x;
  result.y = -result.y;
  return result;
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ double3 crossProduct(const double3 a, const double3 b) {
  return { (a.y * b.z) - (a.z - b.y), (a.z * b.x) - (a.x - b.z), (a.x * b.y) - (a.y * b.x) };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ float3 crossProduct(const float3 a, const float3 b) {
  return { (a.y * b.z) - (a.z - b.y), (a.z * b.x) - (a.x - b.z), (a.x * b.y) - (a.y * b.x) };
}

// Single-precision floating point definitions
#define TCALC float
#  define TCALC3 float3
#  define CONV_FUNC __float2int_rn
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define SIN_FUNC  sinf
#  define ABS_FUNC  fabsf
#  define CHECK_COSARG

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
#    define UPDATE_ATOMS
#      define KERNEL_NAME kfValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kfValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kfValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kfValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear single-precision floating point definitions
#  undef TCALC3
#  undef CONV_FUNC
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef CHECK_COSARG
#undef TCALC

// Double-precision floating point definitions
#define TCALC double
#  define TCALC3 double3
#  define CONV_FUNC __double2ll_rn
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define ACOS_FUNC acos
#  define COS_FUNC  cos
#  define SIN_FUNC  sin
#  define ABS_FUNC  fabs

#  define COMPUTE_FORCE
#    define KERNEL_NAME kdValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define UPDATE_ATOMS
#      define KERNEL_NAME kdValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kdValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear double-precision floating point definitions
#  undef TCALC3
#  undef CONV_FUNC
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#undef TCALC

} // namespace energy
} // namespace omni

#endif
