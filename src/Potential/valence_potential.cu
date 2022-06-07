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

#include accumulation.i

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
