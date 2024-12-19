// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Accelerator/gpu_details.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Numerics/numeric_enumerators.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/cellgrid.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_valence_potential.h"

namespace stormm {
namespace energy {

using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using numerics::chooseAccumulationMethod;
using numerics::getEnumerationName;
using stmath::roundUp;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::boltzmann_constant_f;
using symbols::gafs_to_kcal_f;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::inverse_twopi_f;
using symbols::kcal_to_gafs_f;
using symbols::near_to_one_f;
using symbols::near_to_one_lf;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::VwuAbstractMap;
using synthesis::vwu_abstract_length;
using trajectory::ThermostatKind;
using trajectory::ThermostatPartition;
using topology::TorsionKind;
using topology::VirtualSiteKind;

#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"
#include "Trajectory/thermostat_utilities.cui"
#include "valence_util.cui"

// Single-precision floating point definitions
#define TCALC float
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define LLCONV_FUNC __float2ll_rn
#  define SPLITCONV_FUNC floatToInt63
#  define SPLIT_TYPE int2
#  define SQRT_FUNC sqrtf
#  define CBRT_FUNC cbrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define SIN_FUNC  sinf
#  define ABS_FUNC  fabsf
#  define MIX_FUNC  computeRestraintMixtureF
#  define TCALC_IS_SINGLE

// Compile the standard kernels with all combinations of energy and force accumulation methods.
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#          define KERNEL_NAME kfsValenceForceAccumulationSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 12
#            define KERNEL_NAME kfsValenceAtomUpdateSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 14
#            define KERNEL_NAME kfsValenceForceEnergyAccumulationSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 10
#              define KERNEL_NAME kfsValenceEnergyAtomUpdateSM
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 64
#      define VALENCE_BLOCK_MULTIPLICITY 16
#        define KERNEL_NAME kfValenceForceAccumulationSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define UPDATE_ATOMS
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 12
#          define KERNEL_NAME kfValenceAtomUpdateSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 14
#          define KERNEL_NAME kfValenceForceEnergyAccumulationSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 10
#            define KERNEL_NAME kfValenceEnergyAtomUpdateSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#    undef COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define VALENCE_KERNEL_THREAD_COUNT 64
#      define VALENCE_BLOCK_MULTIPLICITY 16
#        define KERNEL_NAME kfValenceEnergyAccumulationSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef COMPUTE_ENERGY

// Make new kernels with a clash forgiveness check.
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 16
#            define KERNEL_NAME kfsValenceForceAccumulationNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 12
#              define KERNEL_NAME kfsValenceAtomUpdateNonClashSM
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 12
#              define KERNEL_NAME kfsValenceForceEnergyAccumulationNonClashSM
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define UPDATE_ATOMS
#            define VALENCE_KERNEL_THREAD_COUNT 64
#              define VALENCE_BLOCK_MULTIPLICITY 10
#                define KERNEL_NAME kfsValenceEnergyAtomUpdateNonClashSM
#                  include "valence_potential.cui"
#                undef KERNEL_NAME
#              undef VALENCE_BLOCK_MULTIPLICITY
#            undef VALENCE_KERNEL_THREAD_COUNT
#          undef UPDATE_ATOMS
#        undef COMPUTE_ENERGY
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#          define KERNEL_NAME kfValenceForceAccumulationNonClashSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME  
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 12
#            define KERNEL_NAME kfValenceAtomUpdateNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 96
#          define VALENCE_BLOCK_MULTIPLICITY 16
#            define KERNEL_NAME kfValenceForceEnergyAccumulationNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 10
#              define KERNEL_NAME kfValenceEnergyAtomUpdateNonClashSM
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef COMPUTE_ENERGY
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#          define KERNEL_NAME kfValenceEnergyAccumulationNonClashSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
#  undef CLASH_FORGIVENES

// Clear single-precision floating point definitions
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SPLITCONV_FUNC
#  undef SPLIT_TYPE
#  undef SQRT_FUNC
#  undef CBRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#  undef TCALC_IS_SINGLE
#undef TCALC

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_CUDA
extern cudaFuncAttributes
queryValenceKernelRequirementsSM(const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                                 const AccumulationMethod acc_meth, const VwuGoal purpose,
                                 const ClashResponse collision_handling) {

  
  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes result;
  cudaError_t cfa = cudaErrorInvalidValue;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateSM);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateSM);
            break;
          }
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateSM);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateSM);
            break;
          }
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationSM);
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationNonClashSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateNonClashSM);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationNonClashSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateNonClashSM);
            break;
          }
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationNonClashSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateNonClashSM);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationNonClashSM);
            break;
          case VwuGoal::MOVE_PARTICLES:
            cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateNonClashSM);
            break;
          }
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationNonClashSM);
      break;
    }
    break;
  }
  
  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes for kernel k");
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      error_message += "fs";
      break;
    case AccumulationMethod::WHOLE:
      error_message += "f";
      break;
    case AccumulationMethod::AUTOMATIC:
      rtErr("Kernels do not accept " + getEnumerationName(acc_meth) + " accumulation.",
            "queryValenceKernelRequirements");
      break;
    }
    error_message += "Valence";
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          error_message += "ForceEnergyAccumulation";
          break;
        case VwuGoal::MOVE_PARTICLES:
          error_message += "EnergyAtomUpdate";
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          error_message += "ForceAccumulation";
          break;
        case VwuGoal::MOVE_PARTICLES:
          error_message += "AtomUpdate";
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      error_message += "EnergyAccumulation";
      break;
    }
    error_message += "SM.";

    // Report the error
    rtErr(error_message, "queryValenceKernelRequirementsSM");
  }
  
  return result;
}
#endif // STORMM_USE_CUDA

//-------------------------------------------------------------------------------------------------
extern void launchValenceSM(const SyValenceKit<float> &poly_vk,
                            const SyRestraintKit<float, float2, float4> &poly_rk,
                            MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                            const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                            ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                            CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const VwuGoal purpose,
                            const AccumulationMethod refined_force_sum, const int2 bt,
                            const float clash_distance, const float clash_ratio) {
  
  // Rather than a switch over cases of the ClashResponse enumerator, just use the nonzero values
  // of either parameter to indicate that clash damping has been requested.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
    
      // When the goal is to accumulate energies, forces, or both, the force accumulation method
      // becomes a critical detail when choosing the kernel.
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (refined_force_sum) {
        case AccumulationMethod::SPLIT:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfsValenceForceEnergyAccumulationNonClashSM<<<bt.x,
                                                          bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, clash_distance,
                                                                  clash_ratio, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfsValenceForceAccumulationNonClashSM<<<bt.x,
                                                    bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            clash_distance, clash_ratio,
                                                            *gmem_r);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfValenceForceEnergyAccumulationNonClashSM<<<bt.x,
                                                         bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, clash_distance,
                                                                 clash_ratio, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfValenceForceAccumulationNonClashSM<<<bt.x,
                                                   bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, *gmem_r);
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateForce::NO:
        kfValenceEnergyAccumulationNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio, *scw,
                                                              *gmem_r);
        break;
      }
      break;
    case VwuGoal::MOVE_PARTICLES:
    
      // When the goal is to move particles, evaluating the force is obligatory, but the manner in
      // which forces are accumulated is still important.  Whether to accumulate energies while
      // evaluating forces and moving the particles remains a consideration in choosing the proper
      // kernel.
      switch (refined_force_sum) {
      case AccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceEnergyAtomUpdateNonClashSM<<<bt.x,
                                                 bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         clash_distance, clash_ratio, poly_auk,
                                                         *tstw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceAtomUpdateNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         clash_distance, clash_ratio, poly_auk,
                                                         *tstw, *gmem_r);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceEnergyAtomUpdateNonClashSM<<<bt.x,
                                                bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        clash_distance, clash_ratio, poly_auk,
                                                        *tstw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceAtomUpdateNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        clash_distance, clash_ratio, poly_auk,
                                                        *tstw,*gmem_r);
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
  }
  else {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:

      // See above for the rationale on whether forces or energies are evaluated in each context.
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (refined_force_sum) {
        case AccumulationMethod::SPLIT:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfsValenceForceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                *poly_psw, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfsValenceForceAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          *gmem_r);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfValenceForceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                               *poly_psw, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfValenceForceAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         *gmem_r);
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateForce::NO:
        kfValenceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                      *gmem_r);
        break;
      }
      break;
    case VwuGoal::MOVE_PARTICLES:
    
      // See above for the rationale on the choice of each kernel.
      switch (refined_force_sum) {
      case AccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceEnergyAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                       poly_auk, *tstw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                 *tstw, *gmem_r);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceEnergyAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                      poly_auk, *tstw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                *tstw, *gmem_r);
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
  }
}

} // namespace energy
} // namespace stormm
