#include "mm_evaluation.h"

namespace omni {
namespace mm {

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                const NonbondedKit<double> &nbk, const EvaluateForce eval_force,
                const int system_index) {
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                const EvaluateForce eval_force, const int system_index) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(),
                                     ag.getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                const EvaluateForce eval_force, const int system_index) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(),
                                     ag->getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                    const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                    const int step, const int system_index) {
  evalValeMM<double, double, double>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                     vk, nbk, eval_force, system_index);
  evaluateRestraints(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                     eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const RestraintApparatusDpReader &rar,
                    const EvaluateForce eval_force, const int step, const int system_index) {
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index);
  evaluateRestraints(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force,
                    const int system_index, const int step) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(),
                                     ag.getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force,
                    const int system_index, const int step) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc,
                                     psw.zfrc, sc, ag->getDoublePrecisionValenceKit(),
                                     ag->getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    const EvaluateForce eval_force, const int system_index) {
  evaluateNonbondedEnergy(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                          eval_force, eval_force, system_index);
  evalValeMM<double, double, double>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                     vk, nbk, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    const EvaluateForce eval_force, const int system_index) {
  evaluateNonbondedEnergy(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const StaticExclusionMask &se, const EvaluateForce eval_force,
                    const int system_index) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const StaticExclusionMask &se, const EvaluateForce eval_force,
                    const int system_index) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                        const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                        ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                        const int system_index, const int step) {
  evaluateNonbondedEnergy(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                          eval_force, eval_force, system_index);
  evalValeMM<double, double, double>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                     vk, nbk, eval_force, system_index);
  evaluateRestraints(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                     eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                        const int system_index, const int step) {
  evaluateNonbondedEnergy(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index);
  evaluateRestraints(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int system_index, const int step) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int system_index, const int step) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                          system_index);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, system_index, step);
}

} // namespace mm
} // namespace omni
