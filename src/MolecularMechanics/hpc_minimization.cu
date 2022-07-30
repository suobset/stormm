// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/scaling.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"

namespace omni {
namespace mm {

using card::GpuDetails;
using card::KernelManager;
using constants::verytiny;
  
#include "../Potential/accumulation.cui"

// Conjugate gradient particle advancement
#define TCALC_IS_DOUBLE
#  define KERNEL_NAME kdConjGradAdvance
#    include "line_movement.cui"
#  undef KERNEL_NAME
#undef TCALC_IS_DOUBLE
#define KERNEL_NAME kfConjGradAdvance
#  include "line_movement.cui"
#undef KERNEL_NAME

//-------------------------------------------------------------------------------------------------
extern void minimizationKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kdConjGradAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdConjGradAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfConjGradAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfConjGradAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (cudaFuncGetAttributes(&result, kdConjGradAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kdConjGradAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  case PrecisionModel::SINGLE:
    if (cudaFuncGetAttributes(&result, kfConjGradAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kfConjGradAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                               const StaticExclusionMaskSynthesis &poly_se,
                               PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                               ScoreCard *sc, CacheResource *vale_cache, CacheResource *nonb_cache,
                               const ForceAccumulationMethod acc_meth,
                               const KernelManager &launcher) {

  // Obtain abstracts of critical objects, re-using them throughout the inner loop.  Some abstracts
  // are needed by both branches.
  PsSynthesisWriter poly_psw = poly_ps->data();
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  MMControlKit<double> dctrl = mmctrl->dpData();
  MMControlKit<float>  fctrl = mmctrl->spData();
  ScoreCardWriter scw = sc->data();
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

  // Progress through minimization cycles will not be measured with the step counter in the
  // molecular mechanics control object--that will increment at four times the rate in the case
  // of conjugate gradient line minimizations.
  const int total_steps = mmctrl->getTotalCycles();
  mmctrl->primeWorkUnitConters(launcher, prec, poly_ag);
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
        poly_ps->initializeForces()
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &dctrl, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                        nonb_lp);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit();
      const SyRestraintKit<float,
                           float2, float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit();
      CacheResourceKit<float> vale_tbr = vale_cache->spData();
      CacheResourceKit<float> nonb_tbr = nonb_cache->spData();
      for (int i = 0; i < total_steps; i++) {
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &dctrl, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, acc_meth, nonb_lp);
      }
    }
    break;
  }
}

} // namespace mm
} // namespace omni
