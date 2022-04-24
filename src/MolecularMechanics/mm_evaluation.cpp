#include "Potential/nonbonded_potential.h"
#include "Potential/valence_potential.h"
#include "mm_evaluation.h"

namespace omni {
namespace mm {

//-------------------------------------------------------------------------------------------------
void evalValeMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                ScoreCard *sc, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                const EvaluateForce eval_force) {
  evaluateBondTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, eval_force);
  evaluateAngleTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                     eval_force);
  evaluateDihedralTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                        eval_force);
  evaluateUreyBradleyTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                           eval_force);
  evaluateCharmmImproperTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                              eval_force);
  evaluateCmapTerms(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, eval_force);
  evaluateAttenuated14Terms(vk, nbk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                            sc, eval_force, eval_force);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                const NonbondedKit<double> &nbk, const EvaluateForce eval_force) {
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, vk, nbk, eval_force);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                const EvaluateForce eval_force) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionNonbondedKit(),
             eval_force);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                const EvaluateForce eval_force) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag->getDoublePrecisionValenceKit(),
             ag->getDoublePrecisionNonbondedKit(), eval_force);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                    const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                    ScoreCard *sc, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                    const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                    const int step) {
  evalValeMM(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, vk, nbk, eval_force);
  evaluateRestraints(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                     eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const RestraintApparatusDpReader &rar,
                    const EvaluateForce eval_force, const int step) {
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, vk, nbk, eval_force);
  evaluateRestraints(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force, const int step) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionNonbondedKit(),
             eval_force);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force, const int step) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag->getDoublePrecisionValenceKit(),
             ag->getDoublePrecisionNonbondedKit(), eval_force);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
                        const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                        ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                        const int step) {
  evaluateNonbondedEnergy(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                          eval_force, eval_force);
  evalValeMM(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, vk, nbk, eval_force);
  evaluateRestraints(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                     eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintApparatusDpReader &rar, const EvaluateForce eval_force,
                        const int step) {
  evaluateNonbondedEnergy(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force);
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, vk, nbk, eval_force);
  evaluateRestraints(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int step) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force);
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag.getDoublePrecisionValenceKit(), nbk, eval_force);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int step) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                          psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, eval_force);
  evalValeMM(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
             psw.zfrc, sc, ag->getDoublePrecisionValenceKit(), nbk, eval_force);
  evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                     psw.xfrc, psw.yfrc, psw.zfrc, sc, eval_force, 0, step);
}

} // namespace mm
} // namespace omni
