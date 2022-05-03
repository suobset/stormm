#include "nonbonded_potential.h"

namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser, PhaseSpaceWriter psw,
                                ScoreCard *ecard, const EvaluateForce eval_elec_force,
                                const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateNonbondedEnergy<double, double, double>(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd,
                                                         psw.umat, psw.invu, psw.unit_cell,
                                                         psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                         eval_elec_force, eval_vdw_force,
                                                         system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard, const EvaluateForce eval_elec_force,
                                const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateNonbondedEnergy(ag.getDoublePrecisionNonbondedKit(), se.data(), ps->data(), ecard,
                                 eval_elec_force, eval_vdw_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se, PhaseSpace *ps,
                                ScoreCard *ecard, const EvaluateForce eval_elec_force,
                                const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateNonbondedEnergy(ag->getDoublePrecisionNonbondedKit(), se.data(), ps->data(),
                                 ecard, eval_elec_force, eval_vdw_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser, CoordinateFrameReader cfr,
                                ScoreCard *ecard, const int system_index) {
  return evaluateNonbondedEnergy<double, double, double>(nbk, ser, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                         cfr.umat, cfr.invu, cfr.unit_cell,
                                                         nullptr, nullptr, nullptr, ecard,
                                                         EvaluateForce::NO, EvaluateForce::NO,
                                                         system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const NonbondedKit<double> nbk,
                                const StaticExclusionMaskReader ser,
                                const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                const int system_index) {
  return evaluateNonbondedEnergy(nbk, ser, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                const int system_index) {
  return evaluateNonbondedEnergy(ag.getDoublePrecisionNonbondedKit(), se.data(), cfr, ecard,
                                 system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                const int system_index) {
  return evaluateNonbondedEnergy(ag->getDoublePrecisionNonbondedKit(), se.data(), cfr, ecard,
                                 system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     PhaseSpaceWriter psw, ScoreCard *ecard,
                                     const EvaluateForce eval_force, const int system_index) {
  return evaluateGeneralizedBornEnergy(nbk, isk, ngb_tables, psw.xcrd, psw.ycrd, psw.zcrd,
                                       psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                       psw.zfrc, psw.xprv, psw.yprv, psw.zprv, ecard, eval_force,
                                       system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph &ag,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag.getDoublePrecisionNonbondedKit(),
                                       ag.getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables, ps->data(), ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph *ag,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag->getDoublePrecisionNonbondedKit(),
                                       ag->getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables, ps->data(), ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader cfr, ScoreCard *ecard,
                                     const int system_index) {
  std::vector<double>effective_gb_radii(nbk.natom);
  std::vector<double>psi(nbk.natom);
  std::vector<double>sumdeijda(nbk.natom);
  return evaluateGeneralizedBornEnergy<double, double, double>(nbk, isk, ngb_tables, cfr.xcrd,
                                                               cfr.ycrd, cfr.zcrd, cfr.umat,
                                                               cfr.invu, cfr.unit_cell, nullptr,
                                                               nullptr, nullptr,
                                                               effective_gb_radii.data(),
                                                               psi.data(), sumdeijda.data(), ecard,
                                                               EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(nbk, isk, ngb_tables, CoordinateFrameReader(cfw), ecard,
                                       system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph &ag,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag.getDoublePrecisionNonbondedKit(),
                                       ag.getDoublePrecisionImplicitSolventKit(), ngb_tables,
                                       cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph *ag,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrameReader &cfr, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag->getDoublePrecisionNonbondedKit(),
                                       ag->getDoublePrecisionImplicitSolventKit(), ngb_tables,
                                       cfr, ecard, system_index);
}

} // namespace energy
} // namespace omni
