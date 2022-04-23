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

} // namespace mm
} // namespace omni
