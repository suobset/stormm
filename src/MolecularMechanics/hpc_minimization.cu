// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/omni_vector_types.h"
#include "Math/hpc_reduction.h"
#include "Math/reduction_abstracts.h"
#include "Math/reduction_enumerators.h"
#include "Potential/energy_enumerators.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Synthesis/synthesis_enumerators.h"
#include "hpc_minimization.h"

namespace omni {
namespace mm {

using constants::verytiny;
using data_types::int95_t;
using energy::CacheResourceKit;
using energy::EvaluateEnergy;
using energy::EvaluateForce;
using energy::StateVariable;
using math::ConjGradSubstrate;
using math::ReductionGoal;
using math::ReductionStage;
using math::RdwuAbstractMap;
using math::rdwu_abstract_length;
using numerics::max_int_accumulation_f;
using numerics::max_llint_accumulation;
using synthesis::NbwuKind;
using synthesis::SeMaskSynthesisReader;
using synthesis::SyNonbondedKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using synthesis::VwuGoal;

#include "../Potential/accumulation.cui"

// Conjugate gradient particle advancement
#define TCALC_IS_DOUBLE
#  define KERNEL_NAME kdLineAdvance
#    include "line_movement.cui"
#  undef KERNEL_NAME
#undef TCALC_IS_DOUBLE
#define KERNEL_NAME kfLineAdvance
#  include "line_movement.cui"
#undef KERNEL_NAME

//-------------------------------------------------------------------------------------------------
extern void minimizationKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kdLineAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdLineAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfLineAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfLineAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (cudaFuncGetAttributes(&result, kdLineAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kdLineAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  case PrecisionModel::SINGLE:
    if (cudaFuncGetAttributes(&result, kfLineAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kfLineAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchLineAdvance(const PrecisionModel prec, PsSynthesisWriter *poly_psw,
                              const ReductionKit &redk, const ScoreCardWriter &scw,
                              LinMinWriter *lmw, const int move_number, const int2 redu_lp) {
  switch (prec) {
  case PrecisionModel::DOUBLE:
    kdLineAdvance<<<redu_lp.x, redu_lp.y>>>(*poly_psw, redk, scw, *lmw, move_number);
    break;
  case PrecisionModel::SINGLE:
    kfLineAdvance<<<redu_lp.x, redu_lp.y>>>(*poly_psw, redk, scw, *lmw, move_number);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchLineAdvance(const PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                              const AtomGraphSynthesis &poly_ag, ScoreCard *sc,
                              LineMinimization *line_record, const int move_number,
                              const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  const ReductionKit redk(poly_ag, tier);
  ScoreCardWriter scw = sc->data();
  LinMinWriter lmw = line_record->data();
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, move_number, redu_lp);
}
  
//-------------------------------------------------------------------------------------------------
extern void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                               const StaticExclusionMaskSynthesis &poly_se,
                               PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                               ScoreCard *sc, CacheResource *vale_cache, CacheResource *nonb_cache,
                               ReductionBridge *rbg, LineMinimization *line_record,
                               const ForceAccumulationMethod acc_meth,
                               const GpuDetails &gpu, const KernelManager &launcher) {

  // Obtain abstracts of critical objects, re-using them throughout the inner loop.  Some abstracts
  // are needed by both branches.
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data(tier);
  MMControlKit<double> dctrl = mmctrl->dpData(tier);
  MMControlKit<float>  fctrl = mmctrl->spData(tier);
  ScoreCardWriter scw = sc->data(tier);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const int2 nonb_fe_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                          EvaluateEnergy::YES,
                                                          ForceAccumulationMethod::SPLIT);
  const int2 nonb_xe_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::NO,
                                                          EvaluateEnergy::YES,
                                                          ForceAccumulationMethod::SPLIT);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  LinMinWriter lmw = line_record->data(tier);
  const ReductionKit redk(poly_ag, tier);
  ConjGradSubstrate cgsbs(poly_psw, rbg, tier);
  
  // Progress through minimization cycles will not be measured with the step counter in the
  // molecular mechanics control object--that will increment at four times the rate in the case
  // of conjugate gradient line minimizations.
  const int total_steps = mmctrl->getTotalCycles();
  mmctrl->primeWorkUnitCounters(launcher, prec, poly_ag);
  poly_ps->primeConjugateGradientCalculation(gpu, tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit();
      const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit();
      const SyRestraintKit<double,
                           double2, double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit();
      CacheResourceKit<double> vale_tbr = vale_cache->dpData();
      CacheResourceKit<double> nonb_tbr = nonb_cache->dpData();
      for (int i = 0; i < total_steps; i++) {
        poly_ps->initializeForces(gpu, tier);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &dctrl, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, nonb_fe_lp);
        launchValence(poly_vk, poly_rk, &dctrl, &poly_psw, &scw, &vale_tbr, EvaluateForce::YES,
                      EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_fe_lp);
        launchConjugateGradient(redk, &cgsbs, &dctrl, redu_lp);
        kfLineAdvance<<<redu_lp.x, redu_lp.y>>>(poly_psw, redk, scw, lmw, 0);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit();
      const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit();
      const SyRestraintKit<float,
                           float2, float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit();
      CacheResourceKit<float> vale_tbr = vale_cache->spData();
      CacheResourceKit<float> nonb_tbr = nonb_cache->spData();
      for (int i = 0; i < total_steps; i++) {
        poly_ps->initializeForces(gpu, tier);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &fctrl, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, acc_meth, nonb_fe_lp);
        launchValence(poly_vk, poly_rk, &fctrl, &poly_psw, &scw, &vale_tbr, EvaluateForce::YES,
                      EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth, vale_fe_lp);
        launchConjugateGradient(redk, &cgsbs, &fctrl, redu_lp);
        kfLineAdvance<<<redu_lp.x, redu_lp.y>>>(poly_psw, redk, scw, lmw, 0);
      }
    }
    break;
  }
}

} // namespace mm
} // namespace omni
