#include "Constants/symbol_values.h"
#include "Math/vector_ops.h"
#include "Math/matrix_ops.h"
#include "Restraints/bounded_restraint.h"
#include "Structure/local_arrangement.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "valence_potential.h"

namespace omni {
namespace energy {

using math::crossProduct;
using math::matrixMultiply;
using math::matrixVectorMultiply;
using restraints::computeRestraintMixture;
using structure::imageCoordinates;
using structure::ImagingMethod;
using symbols::pi;
using symbols::twopi;
using symbols::inverse_twopi;
using topology::TorsionKind;
using topology::UnitCellType;

//-------------------------------------------------------------------------------------------------
double evalHarmonicStretch(const int i_atom, const int j_atom, const double stiffness,
                           const double equilibrium, const double* xcrd, const double* ycrd,
                           const double* zcrd, const double* umat, const double* invu,
                           const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                           const EvaluateForce eval_force) {
  double dx = xcrd[j_atom] - xcrd[i_atom];
  double dy = ycrd[j_atom] - ycrd[i_atom];
  double dz = zcrd[j_atom] - zcrd[i_atom];
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
  const double dl = dr - equilibrium;

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const double fmag = 2.0 * stiffness * dl / dr;
    const double fmag_dx = fmag * dx;
    const double fmag_dy = fmag * dy;
    const double fmag_dz = fmag * dz;
    xfrc[i_atom] += fmag_dx;
    yfrc[i_atom] += fmag_dy;
    zfrc[i_atom] += fmag_dz;
    xfrc[j_atom] -= fmag_dx;
    yfrc[j_atom] -= fmag_dy;
    zfrc[j_atom] -= fmag_dz;
  }
  return stiffness * dl * dl;
}
  
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
    const double du = evalHarmonicStretch(vk.bond_i_atoms[pos], vk.bond_j_atoms[pos],
                                          vk.bond_keq[param_idx], fabs(vk.bond_leq[param_idx]),
                                          xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                          zfrc, eval_force);
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
double evalHarmonicBend(const int i_atom, const int j_atom, const int k_atom,
                        const double stiffness, const double equilibrium, const double* xcrd,
                        const double* ycrd, const double* zcrd, const double* umat,
                        const double* invu, const UnitCellType unit_cell, double* xfrc,
                        double* yfrc, double* zfrc, const EvaluateForce eval_force) {

  // Compute displacements
  double ba[3], bc[3];
  ba[0] = xcrd[i_atom] - xcrd[j_atom];
  ba[1] = ycrd[i_atom] - ycrd[j_atom];
  ba[2] = zcrd[i_atom] - zcrd[j_atom];
  bc[0] = xcrd[k_atom] - xcrd[j_atom];
  bc[1] = ycrd[k_atom] - ycrd[j_atom];
  bc[2] = zcrd[k_atom] - zcrd[j_atom];
  imageCoordinates(&ba[0], &ba[1], &ba[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // On to the angle force computation
  const double mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
  const double mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
  const double invbabc = 1.0 / sqrt(mgba * mgbc);
  double costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  const double theta = acos(costheta);
  const double dtheta = theta - equilibrium;

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const double dA = -2.0 * stiffness * dtheta / sqrt(1.0 - (costheta * costheta));
    const double sqba = dA / mgba;
    const double sqbc = dA / mgbc;
    const double mbabc = dA * invbabc;
    double adf[3], cdf[3];
    for (int i = 0; i < 3; i++) {
      adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
      cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
    }
    xfrc[i_atom] -= adf[0];
    yfrc[i_atom] -= adf[1];
    zfrc[i_atom] -= adf[2];
    xfrc[j_atom] += adf[0] + cdf[0];
    yfrc[j_atom] += adf[1] + cdf[1];
    zfrc[j_atom] += adf[2] + cdf[2];
    xfrc[k_atom] -= cdf[0];
    yfrc[k_atom] -= cdf[1];
    zfrc[k_atom] -= cdf[2];
  }
  return stiffness * dtheta * dtheta;
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
    const double du = evalHarmonicBend(vk.angl_i_atoms[pos], vk.angl_j_atoms[pos],
                                       vk.angl_k_atoms[pos], vk.angl_keq[param_idx],
                                       vk.angl_theta[param_idx], xcrd, ycrd, zcrd, umat, invu,
                                       unit_cell, xfrc, yfrc, zfrc, eval_force);
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
double evalCosineTwist(const int i_atom, const int j_atom, const int k_atom, const int l_atom,
                       const double amplitude, const double phase_angle, double frequency,
                       const double* xcrd, const double* ycrd, const double* zcrd,
                       const double* umat, const double* invu, const UnitCellType unit_cell,
                       double* xfrc, double* yfrc, double* zfrc, const EvaluateForce eval_force) {

  // Compute displacements
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  ab[0] = xcrd[j_atom] - xcrd[i_atom];
  ab[1] = ycrd[j_atom] - ycrd[i_atom];
  ab[2] = zcrd[j_atom] - zcrd[i_atom];
  bc[0] = xcrd[k_atom] - xcrd[j_atom];
  bc[1] = ycrd[k_atom] - ycrd[j_atom];
  bc[2] = zcrd[k_atom] - zcrd[j_atom];
  cd[0] = xcrd[l_atom] - xcrd[k_atom];
  cd[1] = ycrd[l_atom] - ycrd[k_atom];
  cd[2] = zcrd[l_atom] - zcrd[k_atom];
  imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // Compute cross products and then the angle between the planes
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                   (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  const double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                            -acos(costheta);
  const double sangle = (frequency * theta) - phase_angle;

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {
    const double fr = amplitude * frequency * sin(sangle);
    const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    const double invab = 1.0 / mgab;
    const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
    const double invbc = 1.0 / mgbc;
    const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    const double invcd = 1.0 / mgcd;
    const double cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const double isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
                          fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
    const double cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const double isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
                          fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
    const double invabc = invab * invbc;
    const double invbcd = invbc * invcd;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
    }

    // Transform the rotational derivatives to cartesian coordinates
    const double fa = -invab * isinb2;
    const double fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
    const double fb2 = cosc * invbc * isinc2;
    const double fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
    const double fc2 = cosb * invbc * isinb2;
    const double fd = -invcd * isinc2;
    xfrc[i_atom] += crabbc[0] * fa;
    xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
    xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
    xfrc[l_atom] -= fd * crbccd[0];
    yfrc[i_atom] += crabbc[1] * fa;
    yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
    yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
    yfrc[l_atom] -= fd * crbccd[1];
    zfrc[i_atom] += crabbc[2] * fa;
    zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
    zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
    zfrc[l_atom] -= fd * crbccd[2];
  }
  return amplitude * (1.0 + cos(sangle));
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
    const TorsionKind kind = static_cast<TorsionKind>(vk.dihe_modifiers[pos].w);
    const double du = evalCosineTwist(vk.dihe_i_atoms[pos], vk.dihe_j_atoms[pos],
                                      vk.dihe_k_atoms[pos], vk.dihe_l_atoms[pos],
                                      vk.dihe_amp[param_idx], vk.dihe_phi[param_idx],
                                      vk.dihe_freq[param_idx], xcrd, ycrd, zcrd, umat, invu,
                                      unit_cell, xfrc, yfrc, zfrc, eval_force);

    // Contribute the result to the correct pile: proper or improper
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
    const double du = evalHarmonicStretch(vk.ubrd_i_atoms[pos], vk.ubrd_k_atoms[pos],
                                          vk.ubrd_keq[param_idx], vk.ubrd_leq[param_idx],
                                          xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                          zfrc, eval_force);
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
    const int i_atom = vk.cimp_i_atoms[pos];
    const int j_atom = vk.cimp_j_atoms[pos];
    const int k_atom = vk.cimp_k_atoms[pos];
    const int l_atom = vk.cimp_l_atoms[pos];
    const int param_idx = vk.cimp_param_idx[pos];
    const double keq = vk.cimp_keq[param_idx];
    const double phi = vk.cimp_phi[param_idx];

    // Compute displacements
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
    imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
    imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
    imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

    // Compute cross products and then the angle between the planes
    crossProduct(ab, bc, crabbc);
    crossProduct(bc, cd, crbccd);
    double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
    costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
    crossProduct(crabbc, crbccd, scr);
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    const double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                              -acos(costheta);
    const double dtheta = theta - phi;

    // Contribute the result to the correct pile: proper or improper
    const double contrib = keq * dtheta * dtheta;
    cimp_energy += contrib;
    cimp_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));

    // Compute forces
    if (eval_force == EvaluateForce::YES) {
      const double fr = -2.0 * keq * dtheta;
      const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      const double invab = 1.0 / mgab;
      const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      const double invbc = 1.0 / mgbc;
      const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      const double invcd = 1.0 / mgcd;
      const double cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
      const double isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
                            fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
      const double cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
      const double isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
                            fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
      const double invabc = invab * invbc;
      const double invbcd = invbc * invcd;
      for (int i = 0; i < 3; i++) {
        crabbc[i] *= invabc;
        crbccd[i] *= invbcd;
      }

      // Transform the rotational derivatives to cartesian coordinates
      const double fa = -invab * isinb2;
      const double fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
      const double fb2 = cosc * invbc * isinc2;
      const double fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
      const double fc2 = cosb * invbc * isinb2;
      const double fd = -invcd * isinc2;
      xfrc[i_atom] += crabbc[0] * fa;
      xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
      xfrc[l_atom] -= fd * crbccd[0];
      yfrc[i_atom] += crabbc[1] * fa;
      yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
      yfrc[l_atom] -= fd * crbccd[1];
      zfrc[i_atom] += crabbc[2] * fa;
      zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
      zfrc[l_atom] -= fd * crbccd[2];
    }
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
double evalCmap(const double* cmap_patches, const int* cmap_patch_bounds, const int surf_idx,
                const int surf_dim, const int i_atom, const int j_atom, const int k_atom,
                const int l_atom, const int m_atom, const double* xcrd, const double* ycrd,
                const double* zcrd, const double* umat, const double* invu,
                const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                const EvaluateForce eval_force) {

  // Compute displacements
  std::vector<double> acoef(16);
  double ab[3], bc[3], cd[3], de[3], crabbc[3], crbccd[3], crcdde[3], scr_phi[3], scr_psi[3];
  double phi_progression[4], psi_progression[4], acoef_psi[4];
  for (int i = 0; i < 4; i++) {
    acoef_psi[i] = 0.0;
  }
  ab[0] = xcrd[j_atom] - xcrd[i_atom];
  ab[1] = ycrd[j_atom] - ycrd[i_atom];
  ab[2] = zcrd[j_atom] - zcrd[i_atom];
  bc[0] = xcrd[k_atom] - xcrd[j_atom];
  bc[1] = ycrd[k_atom] - ycrd[j_atom];
  bc[2] = zcrd[k_atom] - zcrd[j_atom];
  cd[0] = xcrd[l_atom] - xcrd[k_atom];
  cd[1] = ycrd[l_atom] - ycrd[k_atom];
  cd[2] = zcrd[l_atom] - zcrd[k_atom];
  de[0] = xcrd[m_atom] - xcrd[l_atom];
  de[1] = ycrd[m_atom] - ycrd[l_atom];
  de[2] = zcrd[m_atom] - zcrd[l_atom];
  imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&de[0], &de[1], &de[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
    
  // Compute the first dihedral
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  double cos_phi = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  cos_phi /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                  (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  crossProduct(crabbc, crbccd, scr_phi);
  cos_phi = (cos_phi < -1.0) ? -1.0 : (cos_phi > 1.0) ? 1.0 : cos_phi;
  double phi = (scr_phi[0]*bc[0] + scr_phi[1]*bc[1] + scr_phi[2]*bc[2] > 0.0) ?  acos(cos_phi) :
                                                                                -acos(cos_phi);
  phi += pi;
  phi = (phi < 0.0) ? phi + twopi : (phi >= twopi) ? phi - twopi : phi;

  // Compute the second dihedral
  crossProduct(cd, de, crcdde);
  double cos_psi = crbccd[0]*crcdde[0] + crbccd[1]*crcdde[1] + crbccd[2]*crcdde[2];
  cos_psi /= sqrt((crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]) *
                  (crcdde[0]*crcdde[0] + crcdde[1]*crcdde[1] + crcdde[2]*crcdde[2]));
  crossProduct(crbccd, crcdde, scr_psi);
  cos_psi = (cos_psi < -1.0) ? -1.0 : (cos_psi > 1.0) ? 1.0 : cos_psi;
  double psi = (scr_psi[0]*bc[0] + scr_psi[1]*bc[1] + scr_psi[2]*bc[2] > 0.0) ?  acos(cos_psi) :
                                                                                -acos(cos_psi);
  psi += pi;
  psi = (psi < 0.0) ? psi + twopi : (psi >= twopi) ? psi - twopi : psi;

  // Compute the patch index (idx_phi, idx_psi) of the CMAP
  const double dsurf_dim = static_cast<double>(surf_dim);
  const double phi_grid = phi * dsurf_dim * inverse_twopi;
  const double psi_grid = psi * dsurf_dim * inverse_twopi;
  const int idx_phi = phi_grid;
  const int idx_psi = psi_grid;
  const double phifrac = phi_grid - idx_phi;
  const double psifrac = psi_grid - idx_psi;

  // Draw in the matrix of spline values and derivatives
  const int patch_idx = cmap_patch_bounds[surf_idx] + (((idx_psi * surf_dim) + idx_phi) * 16);
  for (int i = 0; i < 16; i++) {
    acoef[i] = cmap_patches[patch_idx + i];
  }

  // Perform the matrix multiplications to obtain the bicubic spline coefficients
  phi_progression[0] = 1.0;
  psi_progression[0] = 1.0;
  for (int i = 1; i < 4; i++) {
    phi_progression[i] = phi_progression[i - 1] * phifrac;
    psi_progression[i] = psi_progression[i - 1] * psifrac;
  }
  matrixVectorMultiply(acoef.data(), psi_progression, acoef_psi, 4, 4, 1.0, 1.0, 0.0);
  const double contrib = (phi_progression[0] * acoef_psi[0]) +
                         (phi_progression[1] * acoef_psi[1]) + 
                         (phi_progression[2] * acoef_psi[2]) + 
                         (phi_progression[3] * acoef_psi[3]);

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {

    // The derivatives along phi and psi follow from the energy expression, evaluation of the
    // guiding matrix equation for the bicubic spline potential:
    //
    //                                                    [ a00 a01 a02 a03     [ 1
    // Energy = [ 1 + phifrac + phifrac^2 + phifrac^3 ] *   a10 a11 a12 a13  *    psifrac
    //                                                      a20 a21 a22 a23       psifrac^2
    //                                                      a30 a31 a32 a33 ]     psifrac^3 ]
    double dphi = (((3.0 * acoef[15] * phifrac) + (2.0 * acoef[14])) * phifrac) + acoef[13];
    dphi *= psifrac;
    dphi +=       (((3.0 * acoef[11] * phifrac) + (2.0 * acoef[10])) * phifrac) + acoef[ 9];
    dphi *= psifrac;
    dphi +=       (((3.0 * acoef[ 7] * phifrac) + (2.0 * acoef[ 6])) * phifrac) + acoef[ 5];
    dphi *= psifrac;
    dphi +=       (((3.0 * acoef[ 3] * phifrac) + (2.0 * acoef[ 2])) * phifrac) + acoef[ 1];
    double dpsi = (((3.0 * acoef[15] * psifrac) + (2.0 * acoef[11])) * psifrac) + acoef[ 7];
    dpsi *= phifrac;
    dpsi +=       (((3.0 * acoef[14] * psifrac) + (2.0 * acoef[10])) * psifrac) + acoef[ 6];
    dpsi *= phifrac;
    dpsi +=       (((3.0 * acoef[13] * psifrac) + (2.0 * acoef[ 9])) * psifrac) + acoef[ 5];
    dpsi *= phifrac;
    dpsi +=       (((3.0 * acoef[12] * psifrac) + (2.0 * acoef[ 8])) * psifrac) + acoef[ 4];
    dphi *= dsurf_dim / twopi;
    dpsi *= dsurf_dim / twopi;

    // With the derivative in hand, evaluate the transformation of coordinates for either the
    // phi or psi dihedrals.
    const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    const double invab = 1.0 / mgab;
    const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
    const double invbc = 1.0 / mgbc;
    const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    const double invcd = 1.0 / mgcd;
    const double mgde = sqrt(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
    const double invde = 1.0 / mgde;
    const double invabc = invab * invbc;
    const double invbcd = invbc * invcd;
    const double invcde = invcd * invde;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
      crcdde[i] *= invcde;
    }

    // Feed the gradient, negative of the derivative, into the functions below
    dphi *= -1.0;
    dpsi *= -1.0;

    // Phi accumulation: transform the rotational derivatives to cartesian coordinates
    const double phi_cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const double phi_isinb2 = (phi_cosb * phi_cosb < asymptotic_to_one_lf) ?
                              dphi / (1.0 - (phi_cosb * phi_cosb)) :
                              dphi * inverse_one_minus_asymptote_lf;
    const double phi_cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const double phi_isinc2 = (phi_cosc * phi_cosc < asymptotic_to_one_lf) ?
                              dphi / (1.0 - (phi_cosc * phi_cosc)) :
                              dphi * inverse_one_minus_asymptote_lf;
    const double phi_fa = -invab * phi_isinb2;
    const double phi_fb1 = (mgbc - (mgab * phi_cosb)) * invabc * phi_isinb2;
    const double phi_fb2 = phi_cosc * invbc * phi_isinc2;
    const double phi_fc1 = (mgbc - (mgcd * phi_cosc)) * invbcd * phi_isinc2;
    const double phi_fc2 = phi_cosb * invbc * phi_isinb2;
    const double phi_fd = -invcd * phi_isinc2;

    // Psi accumulation: transform the rotational derivatives to cartesian coordinates
    const double psi_cosb = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const double psi_isinb2 = (psi_cosb * psi_cosb < asymptotic_to_one_lf) ?
                              dpsi / (1.0 - (psi_cosb * psi_cosb)) :
                              dpsi * inverse_one_minus_asymptote_lf;
    const double psi_cosc = -(cd[0]*de[0] + cd[1]*de[1] + cd[2]*de[2]) * invcd * invde;
    const double psi_isinc2 = (psi_cosc * psi_cosc < asymptotic_to_one_lf) ?
                              dpsi / (1.0 - (psi_cosc * psi_cosc)) : 
                              dpsi * inverse_one_minus_asymptote_lf;
    const double psi_fa = -invbc * psi_isinb2;
    const double psi_fb1 = (mgcd - (mgbc * psi_cosb)) * invbcd * psi_isinb2;
    const double psi_fb2 = psi_cosc * invcd * psi_isinc2;
    const double psi_fc1 = (mgcd - (mgde * psi_cosc)) * invcde * psi_isinc2;
    const double psi_fc2 = psi_cosb * invcd * psi_isinb2;
    const double psi_fd = -invde * psi_isinc2;

    // Accumulate the phi dihedral forces
    xfrc[i_atom] += crabbc[0] * phi_fa;
    yfrc[i_atom] += crabbc[1] * phi_fa;
    zfrc[i_atom] += crabbc[2] * phi_fa;
    xfrc[j_atom] += (phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0]);
    yfrc[j_atom] += (phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1]);
    zfrc[j_atom] += (phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2]);
    xfrc[k_atom] += (phi_fc2 * crabbc[0]) - (phi_fc1 * crbccd[0]);
    yfrc[k_atom] += (phi_fc2 * crabbc[1]) - (phi_fc1 * crbccd[1]);
    zfrc[k_atom] += (phi_fc2 * crabbc[2]) - (phi_fc1 * crbccd[2]);
    xfrc[l_atom] -= phi_fd * crbccd[0];
    yfrc[l_atom] -= phi_fd * crbccd[1];
    zfrc[l_atom] -= phi_fd * crbccd[2];

    // Accumulate the psi dihedral forces
    xfrc[j_atom] += crbccd[0] * psi_fa;
    yfrc[j_atom] += crbccd[1] * psi_fa;
    zfrc[j_atom] += crbccd[2] * psi_fa;
    xfrc[k_atom] += (psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0]);
    yfrc[k_atom] += (psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1]);
    zfrc[k_atom] += (psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2]);
    xfrc[l_atom] += (psi_fc2 * crbccd[0]) - (psi_fc1 * crcdde[0]);
    yfrc[l_atom] += (psi_fc2 * crbccd[1]) - (psi_fc1 * crcdde[1]);
    zfrc[l_atom] += (psi_fc2 * crbccd[2]) - (psi_fc1 * crcdde[2]);
    xfrc[m_atom] -= psi_fd * crcdde[0];
    yfrc[m_atom] -= psi_fd * crcdde[1];
    zfrc[m_atom] -= psi_fd * crcdde[2];
  }

  // Return the double-precision energy result.  It will be converted to a fixed-precision
  // quantity later.
  return contrib;
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
double2 evaluateAttenuated14Pair(const int i_atom, const int l_atom, const int attn_idx,
                                 const double coulomb_constant, const double* charges,
                                 const int* lj_param_idx, const double* attn14_elec_factors,
                                 const double* attn14_vdw_factors, const double* lja_14_coeff,
                                 const double* ljb_14_coeff, const int n_lj_types,
                                 const double* xcrd, const double* ycrd, const double* zcrd,
                                 const double* umat, const double* invu,
                                 const UnitCellType unit_cell, double* xfrc, double* yfrc,
                                 double* zfrc, const EvaluateForce eval_elec_force,
                                 const EvaluateForce eval_vdw_force) {
  const int ilj_t = lj_param_idx[i_atom];
  const int jlj_t = lj_param_idx[l_atom];
  double dx = xcrd[l_atom] - xcrd[i_atom];
  double dy = ycrd[l_atom] - ycrd[i_atom];
  double dz = zcrd[l_atom] - zcrd[i_atom];
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const double invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
  const double invr = sqrt(invr2);
  const double invr4 = invr2 * invr2;
  const double ele_scale = attn14_elec_factors[attn_idx];
  const double vdw_scale = attn14_vdw_factors[attn_idx];
  const double qiqj = (coulomb_constant * charges[i_atom] * charges[l_atom]) / ele_scale;
  const double lja = lja_14_coeff[(ilj_t * n_lj_types) + jlj_t] / vdw_scale;
  const double ljb = ljb_14_coeff[(ilj_t * n_lj_types) + jlj_t] / vdw_scale;
  const double ele_contrib = qiqj * invr;
  const double vdw_contrib = (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);

  // Evaluate the force, if requested
  if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
    double fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) : 0.0;
    if (eval_vdw_force == EvaluateForce::YES) {
      fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
    }
    const double fmag_dx = fmag * dx;
    const double fmag_dy = fmag * dy;
    const double fmag_dz = fmag * dz;
    xfrc[i_atom] += fmag_dx;
    yfrc[i_atom] += fmag_dy;
    zfrc[i_atom] += fmag_dz;
    xfrc[l_atom] -= fmag_dx;
    yfrc[l_atom] -= fmag_dy;
    zfrc[l_atom] -= fmag_dz;
  }

  return { ele_contrib, vdw_contrib };
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
    const double2 uc = evaluateAttenuated14Pair(vk.dihe_i_atoms[pos], vk.dihe_l_atoms[pos],
                                                attn_idx, nbk.coulomb_constant, nbk.charge,
                                                nbk.lj_idx, vk.attn14_elec, vk.attn14_vdw,
                                                nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.n_lj_types,
                                                xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                yfrc, zfrc, eval_elec_force, eval_vdw_force);
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
    const double2 uc = evaluateAttenuated14Pair(vk.infr14_i_atoms[pos], vk.infr14_l_atoms[pos],
                                                attn_idx, nbk.coulomb_constant, nbk.charge,
                                                nbk.lj_idx, vk.attn14_elec, vk.attn14_vdw,
                                                nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.n_lj_types,
                                                xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                yfrc, zfrc, eval_elec_force, eval_vdw_force);
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
double3 restraintDelta(const double2 init_k, const double2 final_k, const double4 init_r,
                       const double4 final_r, const double2 mixwt, const double dr) {
  const double r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const double r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const double r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const double r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const double k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const double k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  double dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return { keq, dl, du };
}

//-------------------------------------------------------------------------------------------------
double evalPosnRestraint(const int p_atom, const bool time_dependence, const int step_number,
                         const int kr_param_idx, const int xyz_param_idx, const int* init_step,
                         const int* finl_step, const double2* init_xy, const double2* finl_xy,
                         const double* init_z, const double* finl_z, const double2* init_keq,
                         const double2* finl_keq, const double4* init_r, const double4* finl_r,
                         const double* xcrd, const double* ycrd, const double* zcrd,
                         const double* umat, const double* invu, const UnitCellType unit_cell,
                         double* xfrc, double* yfrc, double* zfrc,
                         const EvaluateForce eval_force) {

  // Determine the weight to give to each endpoint of the restraint
  double2 mixwt = { 1.0, 0.0 };
  if (time_dependence) {
    mixwt = computeRestraintMixture(step_number, init_step[kr_param_idx], finl_step[kr_param_idx]);
  }
  double dx = xcrd[p_atom] - ((mixwt.x * init_xy[xyz_param_idx].x) +
                              (mixwt.y * finl_xy[xyz_param_idx].x));
  double dy = ycrd[p_atom] - ((mixwt.x * init_xy[xyz_param_idx].y) +
                              (mixwt.y * finl_xy[xyz_param_idx].y));
  double dz = zcrd[p_atom] - ((mixwt.x * init_z[xyz_param_idx]) +
                              (mixwt.y * finl_z[xyz_param_idx]));
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
  const double3 rst_eval = restraintDelta(init_keq[kr_param_idx], finl_keq[kr_param_idx],
                                          init_r[kr_param_idx], finl_r[kr_param_idx], mixwt, dr);

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    if (dr < constants::tiny) {

      // The case of positional restraints has a wrinkle when particles are already at their
      // exact target locations.  The force will likely be zero anyway, but it's possible to
      // define a positional restraint that forces a particle to be some finite distance away
      // from the target point, which would imply a non-zero force when the particle is at
      // the target location.  The proper direction of that force is an arbitrary thing in
      // such a case, so subdivide it among all three dimensions.
      const double fmag = 2.0 * rst_eval.x * rst_eval.y / sqrt(3.0);
      xfrc[p_atom] -= fmag;
      yfrc[p_atom] -= fmag;
      zfrc[p_atom] -= fmag;
    }
    else {
      const double fmag = 2.0 * rst_eval.x * rst_eval.y / dr;
      xfrc[p_atom] -= fmag * dx;
      yfrc[p_atom] -= fmag * dy;
      zfrc[p_atom] -= fmag * dz;
    }
  }
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
double evalBondRestraint(const int i_atom, const int j_atom, const bool time_dependence,
                         const int step_number, const int param_idx, const int* init_step,
                         const int* finl_step, const double2* init_keq, const double2* finl_keq,
                         const double4* init_r, const double4* finl_r, const double* xcrd,
                         const double* ycrd, const double* zcrd, const double* umat,
                         const double* invu, const UnitCellType unit_cell, double* xfrc,
                         double* yfrc, double* zfrc, const EvaluateForce eval_force) {
  double2 mixwt = { 1.0, 0.0 };
  if (time_dependence) {
    mixwt = computeRestraintMixture(step_number, init_step[param_idx], finl_step[param_idx]);
  }
  double dx = xcrd[j_atom] - xcrd[i_atom];
  double dy = ycrd[j_atom] - ycrd[i_atom];
  double dz = zcrd[j_atom] - zcrd[i_atom];
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
  const double3 rst_eval = restraintDelta(init_keq[param_idx], finl_keq[param_idx],
                                          init_r[param_idx], finl_r[param_idx], mixwt, dr);
  if (eval_force == EvaluateForce::YES) {
    const double fmag = 2.0 * rst_eval.x * rst_eval.y / dr;
    const double fmag_dx = fmag * dx;
    const double fmag_dy = fmag * dy;
    const double fmag_dz = fmag * dz;
    xfrc[i_atom] += fmag_dx;
    yfrc[i_atom] += fmag_dy;
    zfrc[i_atom] += fmag_dz;
    xfrc[j_atom] -= fmag_dx;
    yfrc[j_atom] -= fmag_dy;
    zfrc[j_atom] -= fmag_dz;
  }    
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
double evalAnglRestraint(const int i_atom, const int j_atom, const int k_atom,
                         const bool time_dependence, const int step_number, const int param_idx,
                         const int* init_step, const int* finl_step, const double2* init_keq,
                         const double2* finl_keq, const double4* init_r, const double4* finl_r,
                         const double* xcrd, const double* ycrd, const double* zcrd,
                         const double* umat, const double* invu, const UnitCellType unit_cell,
                         double* xfrc, double* yfrc, double* zfrc,
                         const EvaluateForce eval_force) {
  double ba[3], bc[3];
  ba[0] = xcrd[i_atom] - xcrd[j_atom];
  ba[1] = ycrd[i_atom] - ycrd[j_atom];
  ba[2] = zcrd[i_atom] - zcrd[j_atom];
  bc[0] = xcrd[k_atom] - xcrd[j_atom];
  bc[1] = ycrd[k_atom] - ycrd[j_atom];
  bc[2] = zcrd[k_atom] - zcrd[j_atom];
  imageCoordinates(&ba[0], &ba[1], &ba[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // On to the angle force computation
  const double mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
  const double mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
  const double invbabc = 1.0 / sqrt(mgba * mgbc);
  double costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  const double theta = acos(costheta);
  const double2 mixwt = computeRestraintMixture(step_number, init_step[param_idx],
                                                finl_step[param_idx]);
  const double3 rst_eval = restraintDelta(init_keq[param_idx], finl_keq[param_idx],
                                          init_r[param_idx], finl_r[param_idx], mixwt, theta);

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const double dA = -2.0 * rst_eval.x * rst_eval.y / sqrt(1.0 - (costheta * costheta));
    const double sqba = dA / mgba;
    const double sqbc = dA / mgbc;
    const double mbabc = dA * invbabc;
    double adf[3], cdf[3];
    for (int i = 0; i < 3; i++) {
      adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
      cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
    }
    xfrc[i_atom] -= adf[0];
    yfrc[i_atom] -= adf[1];
    zfrc[i_atom] -= adf[2];
    xfrc[j_atom] += adf[0] + cdf[0];
    yfrc[j_atom] += adf[1] + cdf[1];
    zfrc[j_atom] += adf[2] + cdf[2];
    xfrc[k_atom] -= cdf[0];
    yfrc[k_atom] -= cdf[1];
    zfrc[k_atom] -= cdf[2];
  }    
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
double evalDiheRestraint(const int i_atom, const int j_atom, const int k_atom, const int l_atom,
                         const bool time_dependence, const int step_number, const int param_idx,
                         const int* init_step, const int* finl_step, const double2* init_keq,
                         const double2* finl_keq, const double4* init_r, const double4* finl_r,
                         const double* xcrd, const double* ycrd, const double* zcrd,
                         const double* umat, const double* invu, const UnitCellType unit_cell,
                         double* xfrc, double* yfrc, double* zfrc,
                         const EvaluateForce eval_force) {
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  ab[0] = xcrd[j_atom] - xcrd[i_atom];
  ab[1] = ycrd[j_atom] - ycrd[i_atom];
  ab[2] = zcrd[j_atom] - zcrd[i_atom];
  bc[0] = xcrd[k_atom] - xcrd[j_atom];
  bc[1] = ycrd[k_atom] - ycrd[j_atom];
  bc[2] = zcrd[k_atom] - zcrd[j_atom];
  cd[0] = xcrd[l_atom] - xcrd[k_atom];
  cd[1] = ycrd[l_atom] - ycrd[k_atom];
  cd[2] = zcrd[l_atom] - zcrd[k_atom];
  imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // Compute cross products and then the angle between the planes
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);    
  double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                   (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                      -acos(costheta);
  const double2 mixwt = computeRestraintMixture(step_number, init_step[param_idx],
                                                finl_step[param_idx]);

  // As part of the setup, the restraint has been arranged so that r1, r2, r3, and r4 are
  // monotonically increasing and span at most two pi radians.  The center of this arrangement
  // may not be at zero, but will be within the range [-pi, pi).  Image the angle to align with
  // the center of the restraint displacements r2 and r3.
  const double midpoint = 0.5 * (mixwt.x * (init_r[param_idx].y + init_r[param_idx].z) +
                                 mixwt.y * (finl_r[param_idx].y + finl_r[param_idx].z));
  double midpoint_delta = imageValue(theta - midpoint, twopi, ImagingMethod::MINIMUM_IMAGE);
  theta += midpoint_delta - (theta - midpoint);
  const double3 rst_eval = restraintDelta(init_keq[param_idx], finl_keq[param_idx],
                                          init_r[param_idx], finl_r[param_idx], mixwt, theta);
  
  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const double fr = -2.0 * rst_eval.x * rst_eval.y;
    const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    const double invab = 1.0 / mgab;
    const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
    const double invbc = 1.0 / mgbc;
    const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    const double invcd = 1.0 / mgcd;
    const double cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const double isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
                          fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
    const double cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const double isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
                          fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
    const double invabc = invab * invbc;
    const double invbcd = invbc * invcd;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
    }

    // Transform the rotational derivatives to Cartesian coordinates
    const double fa = -invab * isinb2;
    const double fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
    const double fb2 = cosc * invbc * isinc2;
    const double fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
    const double fc2 = cosb * invbc * isinb2;
    const double fd = -invcd * isinc2;
    xfrc[i_atom] += crabbc[0] * fa;
    xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
    xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
    xfrc[l_atom] -= fd * crbccd[0];
    yfrc[i_atom] += crabbc[1] * fa;
    yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
    yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
    yfrc[l_atom] -= fd * crbccd[1];
    zfrc[i_atom] += crabbc[2] * fa;
    zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
    zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
    zfrc[l_atom] -= fd * crbccd[2];
  }
  return rst_eval.z;  
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
    const double contrib = evalPosnRestraint(rar.rposn_atoms[i], rar.time_dependence, step_number,
                                             i, i, rar.rposn_init_step, rar.rposn_finl_step,
                                             rar.rposn_init_xy, rar.rposn_finl_xy,
                                             rar.rposn_init_z, rar.rposn_finl_z,
                                             rar.rposn_init_keq, rar.rposn_finl_keq,
                                             rar.rposn_init_r, rar.rposn_finl_r, xcrd, ycrd, zcrd,
                                             umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.nbond; pos++) {
    const double contrib = evalBondRestraint(rar.rbond_i_atoms[pos], rar.rbond_j_atoms[pos],
                                             rar.time_dependence, step_number, pos,
                                             rar.rbond_init_step, rar.rbond_finl_step,
                                             rar.rbond_init_keq, rar.rbond_finl_keq,
                                             rar.rbond_init_r, rar.rbond_finl_r, xcrd, ycrd, zcrd,
                                             umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.nangl; pos++) {
    const double contrib = evalAnglRestraint(rar.rangl_i_atoms[pos], rar.rangl_j_atoms[pos],
                                             rar.rangl_k_atoms[pos], rar.time_dependence,
                                             step_number, pos, rar.rangl_init_step,
                                             rar.rangl_finl_step, rar.rangl_init_keq,
                                             rar.rangl_finl_keq, rar.rangl_init_r,
                                             rar.rangl_finl_r, xcrd, ycrd, zcrd, umat, invu,
                                             unit_cell, xfrc, yfrc, zfrc, eval_force);
    rest_energy += contrib;
    rest_acc += static_cast<llint>(llround(contrib * nrg_scale_factor));
  }
  for (int pos = 0; pos < rar.ndihe; pos++) {
    const double contrib = evalDiheRestraint(rar.rdihe_i_atoms[pos], rar.rdihe_j_atoms[pos],
                                             rar.rdihe_k_atoms[pos], rar.rdihe_l_atoms[pos],
                                             rar.time_dependence, step_number, pos,
                                             rar.rdihe_init_step, rar.rdihe_finl_step,
                                             rar.rdihe_init_keq, rar.rdihe_finl_keq,
                                             rar.rdihe_init_r, rar.rdihe_finl_r, xcrd, ycrd, zcrd,
                                             umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force);
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
