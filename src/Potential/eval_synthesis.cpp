// -*-c++-*-
#include "copyright.h"
#include "eval_synthesis.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
void evalSyNonbondedEnergy(const AtomGraphSynthesis &poly_ag,
                           const StaticExclusionMaskSynthesis &poly_se,
                           PhaseSpaceSynthesis *poly_ps, ScoreCard *ecard,
                           const NonbondedTask task, const PrecisionModel prec,
                           const EvaluateForce eval_elec_force, const EvaluateForce eval_vdw_force,
                           const double clash_minimum_distance, const double clash_ratio) {
  PsSynthesisWriter psyw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    evalSyNonbondedTileGroups<double>(poly_ag.getDoublePrecisionNonbondedKit(), poly_se.data(),
                                      &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                      clash_minimum_distance, clash_ratio);
    break;
  case PrecisionModel::SINGLE:
    evalSyNonbondedTileGroups<float>(poly_ag.getSinglePrecisionNonbondedKit(), poly_se.data(),
                                     &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                     clash_minimum_distance, clash_ratio);
    break;
  }
}

} // namespace energy
} // namespace stormm
