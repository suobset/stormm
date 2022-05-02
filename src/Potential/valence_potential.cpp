#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "valence_potential.h"

namespace omni {
namespace energy {

using topology::TorsionKind;
using topology::UnitCellType;

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                         const double* zcrd, const double* umat, const double* invu,
                         const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                         ScoreCard *ecard, const EvaluateForce eval_force,
                         const int system_index) {

  // Use two accumulators: one, a standard double-precision accumulator and the other a
  // fixed-precision long long integer with discretization at the global energy scaling factor
  // (see Constants/fixed_precision.h).
  double bond_energy = 0.0;
  llint bond_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int param_idx = vk.bond_param_idx[pos];
    const double du =
      evalHarmonicStretch<double, double, double>(vk.bond_i_atoms[pos], vk.bond_j_atoms[pos],
                                                  vk.bond_keq[param_idx],
                                                  fabs(vk.bond_leq[param_idx]),
                                                  xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                  yfrc, zfrc, eval_force);
    bond_energy += du;
    bond_acc += static_cast<llint>(llround(du * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::BOND, bond_acc, system_index);

  // Return the double-precision bond energy sum, if of interest
  return bond_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateBondTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                           psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force, system_index);
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
  return evaluateBondTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                           nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> &vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                          const double* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, const EvaluateForce eval_force,
                          const int system_index) {

  // As in the bond interaction computation above, use two accumulators.  The full double-precision
  // result will be returned if there is interest in comparing to the fixed-precision accumulation.
  // Also as above, the virial is not computed in full double precision; that quantity may affect
  // the way the system changes volume, but because its accumulation is spread over all of the
  // system's interactions and it in turn affects the motion of all particles collectively, it is
  // not of as much interest as the error in the fixed-precision net force on a single particle.
  double angl_energy = 0.0;
  llint angl_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all angle bending terms.
  for (int pos = 0; pos < vk.nangl; pos++) {
    const int param_idx = vk.angl_param_idx[pos];
    const double du =
      evalHarmonicBend<double, double, double>(vk.angl_i_atoms[pos], vk.angl_j_atoms[pos],
                                               vk.angl_k_atoms[pos], vk.angl_keq[param_idx],
                                               vk.angl_theta[param_idx], xcrd, ycrd, zcrd, umat,
                                               invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    angl_energy += du;
    angl_acc += static_cast<llint>(llround(du * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::ANGLE, angl_acc, system_index);

  // Return the double-precision angle energy sum, if of interest
  return angl_energy;
}


//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                            psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force, system_index);
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
  return evaluateAngleTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                            nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                              const double* zcrd, const double* umat, const double* invu,
                              const UnitCellType unit_cell, double* xfrc, double* yfrc,
                              double* zfrc, ScoreCard *ecard, const EvaluateForce eval_force,
                              const int system_index) {
  double2 dihe_energy = { 0.0, 0.0 };
  llint proper_acc   = 0LL;
  llint improper_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all dihedral terms.
  for (int pos = 0; pos < vk.ndihe; pos++) {
    const int param_idx = vk.dihe_param_idx[pos];
    const double du =
      evalDihedralTwist<double, double, double>(vk.dihe_i_atoms[pos], vk.dihe_j_atoms[pos],
                                                vk.dihe_k_atoms[pos], vk.dihe_l_atoms[pos],
                                                vk.dihe_amp[param_idx], vk.dihe_phi[param_idx],
                                                vk.dihe_freq[param_idx], DihedralStyle::COSINE,
                                                xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                yfrc, zfrc, eval_force);

    // Contribute the result to the correct pile: proper or improper
    const TorsionKind kind = static_cast<TorsionKind>(vk.dihe_modifiers[pos].w);
    switch (kind) {
    case TorsionKind::PROPER:
    case TorsionKind::PROPER_NO_14:
      dihe_energy.x += du;
      proper_acc += static_cast<llint>(llround(du * nrg_scale_factor));
      break;
    case TorsionKind::IMPROPER:
    case TorsionKind::IMPROPER_NO_14:
      dihe_energy.y += du;
      improper_acc += static_cast<llint>(llround(du * nrg_scale_factor));
      break;
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::PROPER_DIHEDRAL, proper_acc, system_index);
  ecard->contribute(StateVariable::IMPROPER_DIHEDRAL, improper_acc, system_index);

  // Return the double-precision dihedral energy sum, if of interest
  return dihe_energy;
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                               psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force, system_index);
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
  return evaluateDihedralTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                               nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const double* xcrd,
                                const double* ycrd, const double* zcrd, const double* umat,
                                const double* invu, const UnitCellType unit_cell, double* xfrc,
                                double* yfrc, double* zfrc, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {

  // Accumulate the results in two numerical precision models by looping over all terms
  double ubrd_energy = 0.0;
  llint ubrd_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nubrd; pos++) {
    const int param_idx = vk.ubrd_param_idx[pos];
    const double du =
      evalHarmonicStretch<double, double, double>(vk.ubrd_i_atoms[pos], vk.ubrd_k_atoms[pos],
                                                  vk.ubrd_keq[param_idx], vk.ubrd_leq[param_idx],
                                                  xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                  yfrc, zfrc, eval_force);
    ubrd_energy += du;
    ubrd_acc += static_cast<llint>(llround(du * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::UREY_BRADLEY, ubrd_acc, system_index);

  // Return the double-precision Urey-Bradley energy sum, if of interest
  return ubrd_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, const EvaluateForce eval_force,
                                const int system_index) {
  return evaluateUreyBradleyTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                  psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force,
                                  system_index);
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
  return evaluateUreyBradleyTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                                  cfr.unit_cell, nullptr, nullptr, nullptr, ecard,
                                  EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const double* xcrd,
                                   const double* ycrd, const double* zcrd, const double* umat,
                                   const double* invu, const UnitCellType unit_cell, double* xfrc,
                                   double* yfrc, double* zfrc, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {
  double cimp_energy = 0.0;
  llint cimp_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all CHARMM improper terms.
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  for (int pos = 0; pos < vk.ncimp; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int param_idx = vk.cimp_param_idx[pos];
    const double du = evalDihedralTwist(vk.cimp_i_atoms[pos], vk.cimp_j_atoms[pos],
                                        vk.cimp_k_atoms[pos], vk.cimp_l_atoms[pos],
                                        vk.cimp_keq[param_idx], vk.cimp_phi[param_idx],
                                        1.0, DihedralStyle::HARMONIC, xcrd, ycrd, zcrd, umat,
                                        invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    cimp_energy += du;
    cimp_acc += static_cast<llint>(llround(du * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::CHARMM_IMPROPER, cimp_acc, system_index);

  // Return the double-precision CHARMM improper energy sum, if of interest
  return cimp_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, const EvaluateForce eval_force,
                                   const int system_index) {
  return evaluateCharmmImproperTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard,
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
  return evaluateCharmmImproperTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                                     cfr.unit_cell, nullptr, nullptr, nullptr, ecard,
                                     EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const double* xcrd, const double* ycrd,
                         const double* zcrd, const double* umat, const double* invu,
                         const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                         ScoreCard *ecard, const EvaluateForce eval_force,
                         const int system_index) {
  double cmap_energy = 0.0;
  llint cmap_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all CMAP (here, meaning coupled dihedral) terms.
  for (int pos = 0; pos < vk.ncmap; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int i_atom = vk.cmap_i_atoms[pos];
    const int j_atom = vk.cmap_j_atoms[pos];
    const int k_atom = vk.cmap_k_atoms[pos];
    const int l_atom = vk.cmap_l_atoms[pos];
    const int m_atom = vk.cmap_m_atoms[pos];
    const int surf_idx = vk.cmap_surf_idx[pos];
    const int surf_dim = vk.cmap_dim[surf_idx];
    const double contrib = evalCmap(vk.cmap_patches, vk.cmap_patch_bounds, vk.cmap_surf_idx[pos],
                                    vk.cmap_dim[surf_idx], i_atom, j_atom, k_atom, l_atom, m_atom,
                                    xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                                    eval_force);
    cmap_energy += contrib;
    cmap_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::CMAP, cmap_acc, system_index);

  // Return the double-precision CMAP energy sum, if of interest
  return cmap_energy;
}


//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                           psw.xfrc, psw.yfrc, psw.zfrc, ecard, eval_force, system_index);
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
  return evaluateCmapTerms(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell,
                           nullptr, nullptr, nullptr, ecard, EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms(ag.getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms(ag->getDoublePrecisionValenceKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const double* xcrd, const double* ycrd, const double* zcrd,
                                  const double* umat, const double* invu,
                                  const UnitCellType unit_cell, double* xfrc, double* yfrc,
                                  double* zfrc, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {
  double vdw_energy = 0.0;
  llint vdw_acc = 0LL;
  double ele_energy = 0.0;
  llint ele_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();
  
  // Accumulate results in both electrostatic and van-der Waals (here, Lennard-Jones) energies by
  // looping over all 1:4 interactions.  The two 1:4 "non-bonded" energy components will be
  // combined into a tuple to return the single-threaded double-precision results.  Not only are
  // there two energies to accumulate, there are two possible sources of 1:4 interactions.  Start
  // with dihedrals that control 1:4 pairs between their I and L atoms.
  for (int pos = 0; pos < vk.ndihe; pos++) {

    // The zero index points to a pair interaction with zero electrostatic and zero Lennard-Jones
    // strength (full attenuation, a complete exclusion).  This is most often the case in dihedrals
    // that are part of a larger cosine series affecting the same atoms (only one of the dihedrals
    // will carry the responsibility of computing the I and L atom 1:4 interaction).  Another case
    // in which it can happen is dihedrals on either side of a six-membered ring.
    const int attn_idx = vk.dihe14_param_idx[pos];
    if (attn_idx == 0) {
      continue;
    }
    const Vec2<double> uc =
      evaluateAttenuated14Pair(vk.dihe_i_atoms[pos], vk.dihe_l_atoms[pos], attn_idx,
                               nbk.coulomb_constant, nbk.charge, nbk.lj_idx, vk.attn14_elec,
                               vk.attn14_vdw, nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.n_lj_types,
                               xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                               eval_elec_force, eval_vdw_force);
    ele_energy += uc.x;
    vdw_energy += uc.y;
    ele_acc += static_cast<llint>(llround(uc.x * nrg_scale_factor));
    vdw_acc += static_cast<llint>(llround(uc.y * nrg_scale_factor));
  }

  // Evaluate additional, inferred 1:4 attenuated interactions.  These occur between virtual sites
  // V and other atoms or virtual sites that are 1:4 to the parent atoms of V.
  for (int pos = 0; pos < vk.ninfr14; pos++) {
    const int attn_idx = vk.infr14_param_idx[pos];
    if (attn_idx == 0) {
      continue;
    }
    const Vec2<double> uc =
      evaluateAttenuated14Pair(vk.infr14_i_atoms[pos], vk.infr14_l_atoms[pos], attn_idx,
                               nbk.coulomb_constant, nbk.charge, nbk.lj_idx, vk.attn14_elec,
                               vk.attn14_vdw, nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.n_lj_types,
                               xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                               eval_elec_force, eval_vdw_force);
    ele_energy += uc.x;
    vdw_energy += uc.y;
    ele_acc += static_cast<llint>(llround(uc.x * nrg_scale_factor));
    vdw_acc += static_cast<llint>(llround(uc.y * nrg_scale_factor));
  }

  // Contribute results
  ecard->contribute(StateVariable::ELECTROSTATIC_ONE_FOUR, ele_acc, system_index);
  ecard->contribute(StateVariable::VDW_ONE_FOUR, vdw_acc, system_index);

  // Return the double-precision energy sums, if of interest
  return { ele_energy, vdw_energy };
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {
  return evaluateAttenuated14Terms(vk, nbk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                   psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                   eval_elec_force, eval_vdw_force, system_index);
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
  return evaluateAttenuated14Terms(vk, nbk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                                   cfr.unit_cell, nullptr, nullptr, nullptr, ecard,
                                   EvaluateForce::NO, EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                  const int system_index) {
  return evaluateAttenuated14Terms(vk, nbk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, const CoordinateFrameReader &cfr,
                                  ScoreCard *ecard, const int system_index) {
  return evaluateAttenuated14Terms(ag.getDoublePrecisionValenceKit(),
                                   ag.getDoublePrecisionNonbondedKit(), cfr, ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                  ScoreCard *ecard, const int system_index) {
  return evaluateAttenuated14Terms(ag->getDoublePrecisionValenceKit(),
                                   ag->getDoublePrecisionNonbondedKit(), cfr, ecard, system_index);
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
