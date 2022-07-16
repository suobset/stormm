// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/synthesis_enumerators.h"
#include "hpc_nonbonded_potential.cuh"

namespace omni {
namespace energy {

using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using data_types::int95_t;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
using numerics::max_llint_accumulation_f;
using numerics::PrecisionModel;
using synthesis::NbwuKind;
using synthesis::small_block_max_imports;
using synthesis::small_block_max_atoms;
using synthesis::tile_groups_wu_abstract_length;

#include "accumulation.cui"
#include "Math/rounding.cui"

// Single-precision floating point definitions
#define NONBOND_KERNEL_THREAD_COUNT 256
#define TCALC float
#  define TCALC_IS_SINGLE
#  if (__CUDA_ARCH__ >= 750) && (__CUDA_ARCH__ < 800)
#    define NONBOND_KERNEL_BLOCKS_MULTIPLIER 4
#  else
#    define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#  endif
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define COMPUTE_ENERGY
#        define KERNEL_NAME ktgfsNonbondedForceEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#      undef COMPUTE_ENERGY
#      define KERNEL_NAME ktgfsNonbondedForce
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef SPLIT_FORCE_ACCUMULATION
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgfNonbondedForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgfNonbondedForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgfNonbondedEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  undef TCALC_IS_SINGLE
#undef TCALC

#define TCALC double
#  define SPLIT_FORCE_ACCUMULATION
#  define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgdsNonbondedForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgdsNonbondedForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgdNonbondedEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  undef SPLIT_FORCE_ACCUMULATION
#undef TCALC
#undef NONBOND_KERNEL_THREAD_COUNT

//-------------------------------------------------------------------------------------------------
extern void nonbondedKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedForce, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedForce __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedForceEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedForceEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdsNonbondedForce, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdsNonbondedForce __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdNonbondedEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdNonbondedEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdsNonbondedForceEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdsNonbondedForceEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern void queryNonbondedKernelRequirements(KernelManager *wisdom) {

  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes attr;
  if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForce) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  const GpuDetails wgpu = wisdom->getGpu();
  int arch_block_multiplier = (wgpu.getArchMajor() == 7 && wgpu.getArchMinor() >= 5) ? 4 : 5;
  wisdom->catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::NO, ForceAccumulationMethod::SPLIT,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForceEnergy) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgfNonbondedForce) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::NO, ForceAccumulationMethod::WHOLE,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgfNonbondedForceEnergy) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::YES, ForceAccumulationMethod::WHOLE,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgfNonbondedEnergy) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                                 EvaluateEnergy::YES, ForceAccumulationMethod::WHOLE,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForce) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  arch_block_multiplier = 3;
  wisdom->catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::NO, ForceAccumulationMethod::SPLIT,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForceEnergy) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                                 EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, ktgdNonbondedEnergy) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                                 EvaluateEnergy::YES, ForceAccumulationMethod::WHOLE,
                                 attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                                 attr.sharedSizeBytes);
}
  
//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsDp(const SyNonbondedKit<double> &poly_nbk,
                                        const SeMaskSynthesisReader &poly_ser,
                                        MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                                        ScoreCardWriter *scw, CacheResourceKit<double> *gmem_r,
                                        const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy,
                                        const KernelManager &launcher) {
  const int2 bt = launcher.getNonbondedKernelDims(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS,
                                                  eval_force, eval_energy,
                                                  ForceAccumulationMethod::SPLIT);
  switch (eval_force) {
  case EvaluateForce::YES:
    switch (eval_energy) {
    case EvaluateEnergy::YES:
      ktgdsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                *gmem_r);
      break;
    case EvaluateEnergy::NO:
      ktgdsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
      break;
    }
    break;
  case EvaluateForce::NO:
    ktgdNonbondedEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *gmem_r);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsDp(const AtomGraphSynthesis &poly_ag,
                                        const StaticExclusionMaskSynthesis &poly_se,
                                        MolecularMechanicsControls *mmctrl,
                                        PhaseSpaceSynthesis *poly_ps, ScoreCard *sc,
                                        CacheResource *tb_space, const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy,
                                        const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  MMControlKit<double> ctrl = mmctrl->dpData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
  launchNonbondedTileGroupsDp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r, eval_force,
                              eval_energy, launcher);
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsSp(const SyNonbondedKit<float> &poly_nbk,
                                        const SeMaskSynthesisReader &poly_ser,
                                        MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                                        ScoreCardWriter *scw, CacheResourceKit<float> *gmem_r,
                                        const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy,
                                        const ForceAccumulationMethod force_sum,
                                        const KernelManager &launcher) {
  const int2 bt = launcher.getNonbondedKernelDims(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                                  eval_force, eval_energy,
                                                  ForceAccumulationMethod::SPLIT);
  switch (eval_force) {
  case EvaluateForce::YES:
    switch (force_sum) {
    case ForceAccumulationMethod::SPLIT:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        ktgfsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                  *gmem_r);
        break;
      case EvaluateEnergy::NO:
        ktgfsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::WHOLE:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        ktgfNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                 *gmem_r);
        break;
      case EvaluateEnergy::NO:
        ktgfNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::AUTOMATIC:
      if (poly_psw->frc_bits <= 23) {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                    *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      else {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                   *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      break;
    }
    break;
  case EvaluateForce::NO:
    ktgfNonbondedEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *gmem_r);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsSp(const AtomGraphSynthesis &poly_ag,
                                        const StaticExclusionMaskSynthesis &poly_se,
                                        MolecularMechanicsControls *mmctrl,
                                        PhaseSpaceSynthesis *poly_ps, ScoreCard *sc,
                                        CacheResource *tb_space, const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy,
                                        const ForceAccumulationMethod force_sum,
                                        const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  MMControlKit<float> ctrl = mmctrl->spData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<float> gmem_r = tb_space->spData(tier);
  launchNonbondedTileGroupsSp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r, eval_force,
                              eval_energy, force_sum, launcher);
}

} // namespace energy
} // namespace omni
