// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/hpc_reduction.h"
#include "Math/reduction_abstracts.h"
#include "Math/reduction_enumerators.h"
#include "Potential/energy_enumerators.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Synthesis/synthesis_enumerators.h"
#include "hpc_minimization.h"

namespace stormm {
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

// Conjugate gradient particle advancement: both single- and double-precision kernels do
// calculations in double-precision, but the double-precision form of the kernel respects the
// extended fixed-precision format.
#define TCALC double
#define FABS_FUNC fabs
#define SQRT_FUNC sqrt
#define LLCONV_FUNC __double2ll_rn
#  define TCALC_IS_DOUBLE
#    define KERNEL_NAME kdLineAdvance
#      include "line_movement.cui"
#    undef KERNEL_NAME
#  undef TCALC_IS_DOUBLE
#  define KERNEL_NAME kfLineAdvance
#    include "line_movement.cui"
#  undef KERNEL_NAME
#undef FABS_FUNC
#undef SQRT_FUNC
#undef LLCONV_FUNC
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void minimizationKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kdLineAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdLineAdvance __shared__ memory bank size to eight bytes.",
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
                               PhaseSpaceSynthesis *poly_ps, const MinimizeControl &mincon,
                               MolecularMechanicsControls *mmctrl_fe,
                               MolecularMechanicsControls *mmctrl_xe, ScoreCard *sc,
                               CacheResource *vale_fe_cache, CacheResource *vale_xe_cache,
                               CacheResource *nonb_cache, ReductionBridge *rbg,
                               LineMinimization *line_record,
                               const ForceAccumulationMethod acc_meth, const GpuDetails &gpu,
                               const KernelManager &launcher) {

  // Obtain abstracts of critical objects, re-using them throughout the inner loop.  Some abstracts
  // are needed by both branches.
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data(devc_tier);
  if (sc->getSystemCount() != poly_ag.getSystemCount()) {
    rtErr("The energy tracking object is not prepared for the number of systems contained in the "
          "topology synthesis (" + std::to_string(sc->getSystemCount()) + " vs. " +
          std::to_string(poly_ag.getSystemCount()) + ").", "launchMinimization");
  }

  // Make sure that the energy tracking object holds space for storing all energy components in
  // all snapshots that might be generated over the course of the energy minimization, plus the
  // final energies of all minimized snapshots.
  const int ntpr = mincon.getDiagnosticPrintFrequency();
  const int total_nrg_snapshots = roundUp(mincon.getTotalCycles(), ntpr) + 1;
  if (sc->getSampleCapacity() != total_nrg_snapshots) {
    sc->reserve(total_nrg_snapshots);
  }
  ScoreCardWriter scw = sc->data(devc_tier);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                       EvaluateEnergy::YES,
                                                       ForceAccumulationMethod::SPLIT);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  LinMinWriter lmw = line_record->data(devc_tier);
  const ReductionKit redk(poly_ag, devc_tier);
  ConjGradSubstrate cgsbs(poly_psw, rbg, devc_tier);
  
  // Progress through minimization cycles will not be measured with the step counter in the
  // molecular mechanics control object--that will increment at four times the rate in the case
  // of conjugate gradient line minimizations.
  const int total_steps = mincon.getTotalCycles();
  mmctrl_fe->primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES, prec,
                                   poly_ag);
  mmctrl_xe->primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES, prec,
                                   poly_ag);
  poly_ps->primeConjugateGradientCalculation(gpu, devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<double, double2, double4> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(devc_tier);
      CacheResourceKit<double> vale_fe_tbr = vale_fe_cache->dpData(devc_tier);
      CacheResourceKit<double> vale_xe_tbr = vale_xe_cache->dpData(devc_tier);
      CacheResourceKit<double> nonb_tbr = nonb_cache->dpData(devc_tier);
      MMControlKit<double> ctrl_fe = mmctrl_fe->dpData(devc_tier);
      MMControlKit<double> ctrl_xe = mmctrl_xe->dpData(devc_tier);
      for (int i = 0; i < total_steps; i++) {

        // First stage of the cycle: compute forces and obtain the conjugate gradient move.
        poly_ps->initializeForces(gpu, devc_tier);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_fe, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_fe, &poly_psw, &scw, &vale_fe_tbr,
                      EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_fe_lp);
        if (i % ntpr == 0) {
          sc->commit(devc_tier);
          sc->incrementStep();
          scw = sc->data(devc_tier);
        }
        ctrl_fe.step += 1;
        launchConjugateGradient(redk, &cgsbs, &ctrl_fe, redu_lp);

        // Second stage of the cycle: advance once along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 0, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        ctrl_xe.step += 1;

        // Third stage of the cycle: advance once more along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 1, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        ctrl_xe.step += 1;

        // Final stage of the cycle: advance a final time along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 2, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        ctrl_xe.step += 1;

        // Fit a cubic polynomial to guess the best overall advancement.  Place the system there.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 3, redu_lp);
      }

      // One additional energy calculation to get the final energy
      sc->initialize(devc_tier, gpu);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                      &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
      launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
      ctrl_xe.step += 1;
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(devc_tier);
      CacheResourceKit<float> vale_fe_tbr = vale_fe_cache->spData(devc_tier);
      CacheResourceKit<float> vale_xe_tbr = vale_xe_cache->spData(devc_tier);
      CacheResourceKit<float> nonb_tbr = nonb_cache->spData(devc_tier);
      MMControlKit<float> ctrl_fe = mmctrl_fe->spData(devc_tier);
      MMControlKit<float> ctrl_xe = mmctrl_xe->spData(devc_tier);
      for (int i = 0; i < total_steps; i++) {

        // First stage of the cycle: compute forces and obtain the conjugate gradient move.
        poly_ps->initializeForces(gpu, devc_tier);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_fe, &poly_psw, &scw, &nonb_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, acc_meth, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_fe, &poly_psw, &scw, &vale_fe_tbr,
                      EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                      vale_fe_lp);
        if (i % ntpr == 0) {
          sc->commit(devc_tier);
          sc->incrementStep();
          scw = sc->data(devc_tier);
        }
        ctrl_fe.step += 1;
        launchConjugateGradient(redk, &cgsbs, &ctrl_fe, redu_lp);

        // Second stage of the cycle: advance once along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 0, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                      vale_xe_lp);
        ctrl_xe.step += 1;

        // Third stage of the cycle: advance once more along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 1, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                      vale_xe_lp);
        ctrl_xe.step += 1;

        // Final stage of the cycle: advance a final time along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 2, redu_lp);
        sc->initialize(devc_tier, gpu);
        launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                        &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth, nonb_lp);
        launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                      EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                      vale_xe_lp);
        ctrl_xe.step += 1;

        // Fit a cubic polynomial to guess the best overall advancement.  Place the system there.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 3, redu_lp);
      }

      // One additional energy calculation to get the final energy
      sc->initialize(devc_tier, gpu);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &scw,
                      &nonb_tbr, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth, nonb_lp);
      launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, &scw, &vale_xe_tbr,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                    vale_xe_lp);
      ctrl_xe.step += 1;
    }
    break;
  }

  // Advance the energy tracking history counter to log the final energy results
  sc->commit(devc_tier);
  sc->incrementStep();
  scw = sc->data(devc_tier);
}

//-------------------------------------------------------------------------------------------------
extern ScoreCard launchMinimization(const AtomGraphSynthesis &poly_ag,
                                    const StaticExclusionMaskSynthesis poly_se,
                                    PhaseSpaceSynthesis *poly_ps, MinimizeControls mincon,
                                    const GpuDetails &gpu, const PrecisionModel prec) {

  // Prepare to track the energies of the structures as they undergo geometry optimization.
  ScoreCard result(poly_ps->getSystemCount(),
                   mincon.getTotalCycles() / mincon.getDiagnosticPrintFrequency(), 32);

  // Map out all kernels.  Only a few are needed but this is not a lot of work.
  KernelManager launcher(gpu, poly_ag);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, NbwuKind::TILE_GROUPS,
                                                       EvaluateForce::YES, EvaluateEnergy::YES,
                                                       ForceAccumulationMethod::SPLIT);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);

  // Prepare progress tracking objects.
  MolecularMechanicsControls mmctrl_fe(mincon);
  MolecularMechanicsControls mmctrl_xe(mincon);

  // Prepare cache space for each kernel
  CacheResource vale_fe_cache(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_xe_cache(vale_xe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonb_cache(nonb_lp.x, small_block_max_atoms);
  ReductionBridge poly_rbg(poly_ag.getReductionWorkUnitCount());
  LineMinimization line_record(poly_ag.getSystemCount());
  launchMinimization(prec, poly_ag, poly_se, poly_ps, &mmctrl_fe, &mmctrl_xe, result,
                     &vale_fe_cache, &vale_xe_cache, &nonb_cache, &rbg, &line_record,
                     chooseForceAccumulationMethod(poly_ps->getgetForceAccumulationBits()),
                     gpu, launcher);
  return result;
}

//-------------------------------------------------------------------------------------------------
extern ScoreCard launchMinimization(AtomGraphSynthesis *poly_ag,
                                    PhaseSpaceSynthesis *poly_ps, MinimizeControls mincon,
                                    const GpuDetails &gpu, const PrecisionModel prec) {
  switch (poly_ag.getNonbondedWorkType()) {
  case NbwuKind::TILE_GROUPS:
  case NbwuKind::SUPERTILES:
    {
      const StaticExclusionMaskSynthesis poly_se(poly_ag.getTopologyPointers(),
                                                 poly_ag.getTopologyIndices());
      poly_ag.loadNonbondedWorkUnits(poly_se);
      return launchMinimization(poly_ag, poly_se, poly_ps, mincon, gpu, prec);
    }
    break;
  case NbwuKind::HONEYCOMB:
    break;
  }
  __builtin_unreachable();
}

} // namespace mm
} // namespace stormm
