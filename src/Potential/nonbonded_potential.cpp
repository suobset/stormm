#include "nonbonded_potential.h"

namespace stormm {
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
                                const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateNonbondedEnergy(ag.getDoublePrecisionNonbondedKit(), se.data(), cf.data(), ecard,
                                 system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateNonbondedEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateNonbondedEnergy(ag->getDoublePrecisionNonbondedKit(), se.data(), cf.data(), ecard,
                                 system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     PhaseSpaceWriter psw, ScoreCard *ecard,
                                     const EvaluateForce eval_force, const int system_index) {
  return evaluateGeneralizedBornEnergy<double, double, double>(nbk, ser, isk, ngb_kit, psw.xcrd,
                                                               psw.ycrd, psw.zcrd, psw.xfrc,
                                                               psw.yfrc, psw.zfrc, psw.xprv,
                                                               psw.yprv, psw.zprv, ecard,
                                                               eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag.getDoublePrecisionNonbondedKit(), se.data(),
                                       ag.getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables.getDoublePrecisionAbstract(), ps->data(), ecard,
                                       eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables, PhaseSpace *ps,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag->getDoublePrecisionNonbondedKit(), se.data(),
                                       ag->getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables.getDoublePrecisionAbstract(), ps->data(), ecard,
                                       eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     const CoordinateFrameReader cfr, ScoreCard *ecard,
                                     const int system_index) {
  std::vector<double>effective_gb_radii(nbk.natom);
  std::vector<double>psi(nbk.natom);
  std::vector<double>sumdeijda(nbk.natom);
  return evaluateGeneralizedBornEnergy<double, double, double>(nbk, ser, isk, ngb_kit, cfr.xcrd,
                                                               cfr.ycrd, cfr.zcrd, nullptr,
                                                               nullptr, nullptr,
                                                               effective_gb_radii.data(),
                                                               psi.data(), sumdeijda.data(), ecard,
                                                               EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornKit<double> ngb_kit,
                                     const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(nbk, ser, isk, ngb_kit, CoordinateFrameReader(cfw), ecard,
                                       system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph &ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrame &cf, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag.getDoublePrecisionNonbondedKit(), se.data(),
                                       ag.getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables.getDoublePrecisionAbstract(),
                                       cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateGeneralizedBornEnergy(const AtomGraph *ag, const StaticExclusionMask &se,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const CoordinateFrame &cf, ScoreCard *ecard,
                                     const int system_index) {
  return evaluateGeneralizedBornEnergy(ag->getDoublePrecisionNonbondedKit(), se.data(),
                                       ag->getDoublePrecisionImplicitSolventKit(),
                                       ngb_tables.getDoublePrecisionAbstract(),
                                       cf.data(), ecard, system_index);
}

} // namespace energy
} // namespace stormm
