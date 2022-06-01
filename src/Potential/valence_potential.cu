A// -*-c++-*-
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph_abstracts.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/valence_workunit.h"

namespace omni {
namespace energy {

using synthesis::maximum_valence_work_unit_atoms;

// Single-precision floating point definitions
#define TCALC float
#define SQRT_FUNC sqrtf
#define ACOS_FUNC acosf
#define COS_FUNC  cosf
#define ABS_FUNC  fabsf

//-------------------------------------------------------------------------------------------------
#define COMPUTE_FORCE
#define KERNEL_NAME kfValenceForceAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME  
#define COMPUTE_ENERGY
#define KERNEL_NAME kfValenceForceEnergyAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME  
#undef COMPUTE_FORCE
#define KERNEL_NAME kfValenceEnergyAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME
#undef  COMPUTE_FORCE

// Clear single-precision floating point definitions
#undef TCALC
#undef SQRT_FUNC
#undef ACOS_FUNC
#undef COS_FUNC
#undef ABS_FUNC

// Double-precision floating point definitions
#define TCALC double
#define SQRT_FUNC sqrt
#define ACOS_FUNC acos
#define COS_FUNC  cos
#define ABS_FUNC  fabs

//-------------------------------------------------------------------------------------------------
#define COMPUTE_FORCE
#define KERNEL_NAME kdValenceForceAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME  
#define COMPUTE_ENERGY
#define KERNEL_NAME kdValenceForceEnergyAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME
#undef COMPUTE_FORCE
#define KERNEL_NAME kdValenceEnergyAccumulation
#include "valence_potential.cui"
#undef KERNEL_NAME
#undef  COMPUTE_FORCE

// Clear double-precision floating point definitions
#undef TCALC
#undef SQRT_FUNC
#undef ACOS_FUNC
#undef COS_FUNC
#undef ABS_FUNC

} // namespace energy
} // namespace omni

#endif
