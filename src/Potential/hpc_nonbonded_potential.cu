// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/synthesis_enumerators.h"
#include "hpc_nonbonded_potential.h"

namespace stormm {
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

#include "Numerics/accumulation.cui"
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
extern cudaFuncAttributes
queryNonbondedKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                 const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                                 const ForceAccumulationMethod acc_meth) {

  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes attr;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForceEnergy) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdsNonbondedForceEnergy.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case EvaluateEnergy::NO:
        if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForce) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdsNonbondedForce.",
                "queryNonbondedKernelRequirements");
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      if (cudaFuncGetAttributes(&attr, ktgdNonbondedEnergy) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ktgdNonbondedEnergy.",
              "queryValenceKernelRequirements");
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (kind) {
      case NbwuKind::TILE_GROUPS:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (acc_meth) {
          case ForceAccumulationMethod::SPLIT:
            if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfsNonbondedForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case ForceAccumulationMethod::WHOLE:
            if (cudaFuncGetAttributes(&attr, ktgfNonbondedForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfNonbondedForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (acc_meth) {
          case ForceAccumulationMethod::SPLIT:
            if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfsNonbondedForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case ForceAccumulationMethod::WHOLE:
            if (cudaFuncGetAttributes(&attr, ktgfNonbondedForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfNonbondedForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:
        if (cudaFuncGetAttributes(&attr, ktgfNonbondedEnergy) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfNonbondedEnergy.",
                "queryValenceKernelRequirements");
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
      break;
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  }
  return attr;
}
  
//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<double> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<double> *ctrl,
                            PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                            CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
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
    break;
  case NbwuKind::SUPERTILES:
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<float> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<float> *ctrl,
                            PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                            CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy,
                            const ForceAccumulationMethod force_sum, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
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
        if (poly_psw->frc_bits <= 24) {
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
    break;
  case NbwuKind::SUPERTILES:
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            ScoreCard *sc, CacheResource *tb_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy,
                            const ForceAccumulationMethod force_sum,
                            const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const int2 bt = launcher.getNonbondedKernelDims(prec, nb_work_type, eval_force, eval_energy,
                                                  force_sum);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r,
                      eval_force, eval_energy, bt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r,
                      eval_force, eval_energy, force_sum, bt);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            ScoreCard *sc, CacheResource *tb_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const KernelManager &launcher) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, sc, tb_space, eval_force, eval_energy,
                    ForceAccumulationMethod::SPLIT, launcher);
  }
  else {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, sc, tb_space, eval_force, eval_energy,
                    ForceAccumulationMethod::WHOLE, launcher);
  }
}

} // namespace energy
} // namespace stormm
