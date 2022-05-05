#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "valence_potential.h"

namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateBondTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                   psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                   psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  return evaluateBondTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  return evaluateBondTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                   cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                   nullptr, ecard, EvaluateForce::NO,
                                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> &vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateBondTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateBondTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                    psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                    psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard,
                            eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                            eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                    cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                    nullptr, ecard, EvaluateForce::NO,
                                                    system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          const int system_index) {
  return evaluateAngleTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          const int system_index) {
  return evaluateAngleTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                       psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                       psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                               system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                               system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                       cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                       nullptr, ecard, EvaluateForce::NO,
                                                       system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              const int system_index) {
  return evaluateDihedralTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              const int system_index) {
  return evaluateDihedralTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, const EvaluateForce eval_force,
                                const int system_index) {
  return evaluateUreyBradleyTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                          psw.umat, psw.invu, psw.unit_cell,
                                                          psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                          eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {
  return evaluateUreyBradleyTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {
  return evaluateUreyBradleyTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                                  eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                          cfr.umat, cfr.invu, cfr.unit_cell,
                                                          nullptr, nullptr, nullptr, ecard,
                                                          EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateUreyBradleyTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateUreyBradleyTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, const EvaluateForce eval_force,
                                   const int system_index) {
  return evaluateCharmmImproperTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                             psw.umat, psw.invu, psw.unit_cell,
                                                             psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                             eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {
  return evaluateCharmmImproperTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {
  return evaluateCharmmImproperTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                             cfr.umat, cfr.invu, cfr.unit_cell,
                                                             nullptr, nullptr, nullptr, ecard,
                                                             EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph *ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                   psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                   psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                   cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                   nullptr, ecard, EvaluateForce::NO,
                                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateCmapTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateCmapTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateAttenuated14Terms<double, double, double>(vk, nbk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                           psw.umat, psw.invu, psw.unit_cell,
                                                           psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                           eval_elec_force, eval_vdw_force,
                                                           system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateAttenuated14Terms(ag.getDoublePrecisionValenceKit(),
                                   ag.getDoublePrecisionNonbondedKit(), ps->data(), ecard,
                                   eval_elec_force, eval_vdw_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateAttenuated14Terms(ag->getDoublePrecisionValenceKit(),
                                   ag->getDoublePrecisionNonbondedKit(), ps->data(), ecard,
                                   eval_elec_force, eval_vdw_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameReader cfr, ScoreCard *ecard,
                                  const int system_index) {
  return evaluateAttenuated14Terms<double, double, double>(vk, nbk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                           cfr.umat, cfr.invu, cfr.unit_cell,
                                                           nullptr, nullptr, nullptr, ecard,
                                                           EvaluateForce::NO, EvaluateForce::NO,
                                                           system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                  const int system_index) {
  return evaluateAttenuated14Terms(vk, nbk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, const int system_index) {
  return evaluateAttenuated14Terms(ag.getDoublePrecisionValenceKit(),
                                   ag.getDoublePrecisionNonbondedKit(), cf.data(), ecard,
                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph *ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, const int system_index) {
  return evaluateAttenuated14Terms(ag->getDoublePrecisionValenceKit(),
                                   ag->getDoublePrecisionNonbondedKit(), cf.data(), ecard,
                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatusDpReader rar, const double* xcrd,
                          const double* ycrd, const double* zcrd, const double* umat,
                          const double* invu, const UnitCellType unit_cell, double* xfrc,
                          double* yfrc, double* zfrc, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  double rest_energy = 0.0;
  llint rest_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all restraint terms
  for (int i = 0; i < rar.nposn; i++) {
    const double contrib =
      evalPosnRestraint(rar.rposn_atoms[i], rar.time_dependence, step_number,
                        rar.rposn_init_step[i], rar.rposn_finl_step[i],
                        Vec2<double>(rar.rposn_init_xy[i]), Vec2<double>(rar.rposn_finl_xy[i]),
                        rar.rposn_init_z[i], rar.rposn_finl_z[i],
                        Vec2<double>(rar.rposn_init_keq[i]), Vec2<double>(rar.rposn_finl_keq[i]),
                        Vec4<double>(rar.rposn_init_r[i]), Vec4<double>(rar.rposn_finl_r[i]), xcrd,
                        ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.nbond; pos++) {
    const double contrib =
      evalBondRestraint(rar.rbond_i_atoms[pos], rar.rbond_j_atoms[pos], rar.time_dependence,
                        step_number, rar.rbond_init_step[pos], rar.rbond_finl_step[pos],
                        Vec2<double>(rar.rbond_init_keq[pos]),
                        Vec2<double>(rar.rbond_finl_keq[pos]), Vec4<double>(rar.rbond_init_r[pos]),
                        Vec4<double>(rar.rbond_finl_r[pos]), xcrd, ycrd, zcrd, umat, invu,
                        unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.nangl; pos++) {
    const double contrib =
      evalAnglRestraint(rar.rangl_i_atoms[pos], rar.rangl_j_atoms[pos], rar.rangl_k_atoms[pos],
                        rar.time_dependence, step_number, rar.rangl_init_step[pos],
                        rar.rangl_finl_step[pos], Vec2<double>(rar.rangl_init_keq[pos]),
                        Vec2<double>(rar.rangl_finl_keq[pos]), Vec4<double>(rar.rangl_init_r[pos]),
                        Vec4<double>(rar.rangl_finl_r[pos]), xcrd, ycrd, zcrd, umat, invu,
                        unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.ndihe; pos++) {
    const double contrib =
      evalDiheRestraint(rar.rdihe_i_atoms[pos], rar.rdihe_j_atoms[pos], rar.rdihe_k_atoms[pos],
                        rar.rdihe_l_atoms[pos], rar.time_dependence, step_number,
                        rar.rdihe_init_step[pos], rar.rdihe_finl_step[pos],
                        Vec2<double>(rar.rdihe_init_keq[pos]),
                        Vec2<double>(rar.rdihe_finl_keq[pos]),
                        Vec4<double>(rar.rdihe_init_r[pos]), Vec4<double>(rar.rdihe_finl_r[pos]),
                        xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::RESTRAINT, rest_acc, system_index);
  
  // Return the double-precision energy sum, if of interest
  return rest_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatusDpReader rar, PhaseSpaceWriter psw,
                          ScoreCard *ecard, const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  return evaluateRestraints(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                            psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force, system_index,
                            step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus &ra, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  PhaseSpaceWriter psw = ps->data();
  return evaluateRestraints(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                            psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force,
                            system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus *ra, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  PhaseSpaceWriter psw = ps->data();
  return evaluateRestraints(ra->dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                            psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force,
                            system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatusDpReader rar, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  return evaluateRestraints(rar, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                            nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index,
                            step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatusDpReader rar, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  const CoordinateFrameReader cfr(cfw);
  return evaluateRestraints(rar, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                            nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index,
                            step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus &ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  return evaluateRestraints(ra.dpData(), cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                            cfr.unit_cell, nullptr, nullptr, nullptr, ecard, EvaluateForce::NO,
                            system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus *ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  return evaluateRestraints(ra->dpData(), cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                            cfr.unit_cell, nullptr, nullptr, nullptr, ecard, EvaluateForce::NO,
                            system_index, step_number);
}

} // namespace energy
} // namespace omni
