// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "Synthesis/nonbonded_workunit.h"
#include "hpc_nonbonded_potential.cuh"

namespace omni {
namespace energy {

using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using synthesis::small_block_max_imports;
using synthesis::small_block_max_atoms;
using synthesis::tile_groups_wu_abstract_length;

#include "accumulation.cui"
#include "Math/rounding.cui"

// Single-precision floating point definitions
#define NONBOND_KERNEL_THREAD_COUNT 256
#define TCALC float
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
#undef TCALC

#define TCALC double
#  define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgdNonbondedForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgdNonbondedForce
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
#undef TCALC
#undef NONBOND_KERNEL_THREAD_COUNT

//-------------------------------------------------------------------------------------------------
extern void nonbondedKernelSetup() {
  cudaFuncSetSharedMemConfig(ktgfNonbondedForce,       cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(ktgfNonbondedEnergy,      cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(ktgfNonbondedForceEnergy, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(ktgdNonbondedForce,       cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(ktgdNonbondedEnergy,      cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(ktgdNonbondedForceEnergy, cudaSharedMemBankSizeEightByte);
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsDp(const SyNonbondedKit<double> &poly_nbk,
                                        const SeMaskSynthesisReader &poly_ser,
                                        MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                                        ScoreCardWriter *scw, CacheResourceKit<double> *gmem_r,
                                        const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy, const GpuDetails &gpu) {
  const int blocks_multiplier = (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  const int nblocks = gpu.getSMPCount() * blocks_multiplier;
  const int nthreads = 256;
  switch (eval_force) {
  case EvaluateForce::YES:
    switch (eval_energy) {
    case EvaluateEnergy::YES:
      ktgdNonbondedForceEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                      *gmem_r);
      break;
    case EvaluateEnergy::NO:
      ktgdNonbondedForce<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
      break;
    }
    break;
  case EvaluateForce::NO:
    ktgdNonbondedEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                               *gmem_r);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsDp(const AtomGraphSynthesis &poly_ag,
                                        const StaticExclusionMaskSynthesis &poly_se,
                                        MolecularMechanicsControls *mmctrl,
                                        PhaseSpaceSynthesis *poly_ps, ScoreCard *sc,
                                        CacheResource *tb_space, const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy, const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  MMControlKit<double> ctrl = mmctrl->dpData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
  launchNonbondedTileGroupsDp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r, eval_force,
                              eval_energy, gpu);
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbondedTileGroupsSp(const SyNonbondedKit<float> &poly_nbk,
                                        const SeMaskSynthesisReader &poly_ser,
                                        MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                                        ScoreCardWriter *scw, CacheResourceKit<float> *gmem_r,
                                        const EvaluateForce eval_force,
                                        const EvaluateEnergy eval_energy,
                                        const ForceAccumulationMethod force_sum,
                                        const GpuDetails &gpu) {
  const int blocks_multiplier = (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  const int nblocks = gpu.getSMPCount() * blocks_multiplier;
  const int nthreads = 256;
  switch (eval_force) {
  case EvaluateForce::YES:
    switch (force_sum) {
    case ForceAccumulationMethod::SPLIT:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        ktgfsNonbondedForceEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                         *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        ktgfsNonbondedForce<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::WHOLE:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        ktgfNonbondedForceEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                        *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        ktgfNonbondedForce<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::AUTOMATIC:
      if (poly_psw->frc_bits <= 23) {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfsNonbondedForceEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                           *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfsNonbondedForce<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                     *gmem_r);
          break;
        }
      }
      else {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfNonbondedForceEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                          *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfNonbondedForce<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      break;
    }
    break;
  case EvaluateForce::NO:
    ktgfNonbondedEnergy<<<nblocks, nthreads>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                               *gmem_r);
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
                                        const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  MMControlKit<float> ctrl = mmctrl->spData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<float> gmem_r = tb_space->spData(tier);
  launchNonbondedTileGroupsSp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r, eval_force,
                              eval_energy, force_sum, gpu);
}

} // namespace energy
} // namespace omni
