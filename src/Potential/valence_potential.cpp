#include "Constants/symbol_values.h"
#include "Math/vector_ops.h"
#include "Math/matrix_ops.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "valence_potential.h"

namespace omni {
namespace energy {

using math::crossProduct;
using math::matrixMultiply;
using math::matrixVectorMultiply;
using symbols::pi;
using symbols::twopi;
using symbols::inverse_twopi;
using topology::NonbondedKit;
using topology::TorsionKind;
using topology::ValenceKit;
using trajectory::PhaseSpaceWriter;

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();

  // Use two accumulators: one, a standard double-precision accumulator and the other a
  // fixed-precision long long integer with discretization at the global energy scaling factor
  // (see Constants/fixed_precision.h).
  double bond_energy = 0.0;
  llint bond_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int i_atom = vk.bond_i_atoms[pos];
    const int j_atom = vk.bond_j_atoms[pos];
    const int param_idx = vk.bond_param_idx[pos];
    const double keq = vk.bond_keq[param_idx];
    const double leq = fabs(vk.bond_leq[param_idx]);
    const double dx = psw.xcrd[j_atom] - psw.xcrd[i_atom];
    const double dy = psw.ycrd[j_atom] - psw.ycrd[i_atom];
    const double dz = psw.zcrd[j_atom] - psw.zcrd[i_atom];
    const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    const double dl = dr - leq;
    const double du = keq * dl * dl;
    bond_energy += du;
    bond_acc += static_cast<llint>(round(du * nrg_scale_factor));

    // Compute forces
    if (eval_force == EvaluateForce::YES) {
      const double fmag = 2.0 * keq * dl / dr;
      const double fmag_dx = fmag * dx;
      const double fmag_dy = fmag * dy;
      const double fmag_dz = fmag * dz;
      psw.xfrc[i_atom] += fmag_dx;
      psw.yfrc[i_atom] += fmag_dy;
      psw.zfrc[i_atom] += fmag_dz;
      psw.xfrc[j_atom] -= fmag_dx;
      psw.yfrc[j_atom] -= fmag_dy;
      psw.zfrc[j_atom] -= fmag_dz;
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::BOND, bond_acc, system_index);

  // Return the double-precision bond energy sum, if of interest
  return bond_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();

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
  double ac[3], ba[3], bc[3];
  for (int pos = 0; pos < vk.nangl; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int i_atom = vk.angl_i_atoms[pos];
    const int j_atom = vk.angl_j_atoms[pos];
    const int k_atom = vk.angl_k_atoms[pos];
    const int param_idx = vk.angl_param_idx[pos];
    const double keq = vk.angl_keq[param_idx];
    const double theta0 = vk.angl_theta[param_idx];

    // Compute displacements
    ba[0] = psw.xcrd[i_atom] - psw.xcrd[j_atom];
    ba[1] = psw.ycrd[i_atom] - psw.ycrd[j_atom];
    ba[2] = psw.zcrd[i_atom] - psw.zcrd[j_atom];
    bc[0] = psw.xcrd[k_atom] - psw.xcrd[j_atom];
    bc[1] = psw.ycrd[k_atom] - psw.ycrd[j_atom];
    bc[2] = psw.zcrd[k_atom] - psw.zcrd[j_atom];
    ac[0] = psw.xcrd[k_atom] - psw.xcrd[i_atom];
    ac[1] = psw.ycrd[k_atom] - psw.ycrd[i_atom];
    ac[2] = psw.zcrd[k_atom] - psw.zcrd[i_atom];

    // On to the angle force computation
    const double mgba = ba[0]*ba[0] + ba[1]*ba[1] + ba[2]*ba[2];
    const double mgbc = bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2];
    const double invbabc = 1.0 / sqrt(mgba * mgbc);
    double costheta = (ba[0]*bc[0] + ba[1]*bc[1] + ba[2]*bc[2]) * invbabc;
    costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
    const double theta = acos(costheta);
    const double dtheta = theta - theta0;
    const double du = keq * dtheta * dtheta;
    angl_energy += du;
    angl_acc += static_cast<llint>(round(du * nrg_scale_factor));

    // Compute forces
    if (eval_force == EvaluateForce::YES) {
      const double dA = -2.0 * keq * dtheta / sqrt(1.0 - costheta * costheta);
      const double sqba = dA/mgba;
      const double sqbc = dA/mgbc;
      const double mbabc = dA * invbabc;
      double adf[3], cdf[3];
      for (int i = 0; i < 3; i++) {
        adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
        cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
      }
      psw.xfrc[i_atom] -= adf[0];
      psw.yfrc[i_atom] -= adf[1];
      psw.zfrc[i_atom] -= adf[2];
      psw.xfrc[j_atom] += adf[0] + cdf[0];
      psw.yfrc[j_atom] += adf[1] + cdf[1];
      psw.zfrc[j_atom] += adf[2] + cdf[2];
      psw.xfrc[k_atom] -= cdf[0];
      psw.yfrc[k_atom] -= cdf[1];
      psw.zfrc[k_atom] -= cdf[2];
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::ANGLE, angl_acc, system_index);

  // Return the double-precision angle energy sum, if of interest
  return angl_energy;
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();
  double2 dihe_energy = { 0.0, 0.0 };
  llint proper_acc   = 0LL;
  llint improper_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all dihedral terms.
  double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  for (int pos = 0; pos < vk.ndihe; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int i_atom = vk.dihe_i_atoms[pos];
    const int j_atom = vk.dihe_j_atoms[pos];
    const int k_atom = vk.dihe_k_atoms[pos];
    const int l_atom = vk.dihe_l_atoms[pos];
    const int param_idx = vk.dihe_param_idx[pos];
    const double ampl = vk.dihe_amp[param_idx];
    const double freq = vk.dihe_freq[param_idx];
    const double phi  = vk.dihe_phi[param_idx];
    const TorsionKind kind = static_cast<TorsionKind>(vk.dihe_modifiers[pos].w);

    // Compute displacements
    ab[0] = psw.xcrd[j_atom] - psw.xcrd[i_atom];
    ab[1] = psw.ycrd[j_atom] - psw.ycrd[i_atom];
    ab[2] = psw.zcrd[j_atom] - psw.zcrd[i_atom];
    bc[0] = psw.xcrd[k_atom] - psw.xcrd[j_atom];
    bc[1] = psw.ycrd[k_atom] - psw.ycrd[j_atom];
    bc[2] = psw.zcrd[k_atom] - psw.zcrd[j_atom];
    cd[0] = psw.xcrd[l_atom] - psw.xcrd[k_atom];
    cd[1] = psw.ycrd[l_atom] - psw.ycrd[k_atom];
    cd[2] = psw.zcrd[l_atom] - psw.zcrd[k_atom];

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
    const double sangle = freq * theta - phi;

    // Contribute the result to the correct pile: proper or improper
    const double contrib = ampl * (1.0 + cos(sangle));
    switch (kind) {
    case TorsionKind::PROPER:
    case TorsionKind::PROPER_NO_14:
      dihe_energy.x += contrib;
      proper_acc += static_cast<llint>(round(contrib * nrg_scale_factor));
      break;
    case TorsionKind::IMPROPER:
    case TorsionKind::IMPROPER_NO_14:
      dihe_energy.y += ampl * (1.0 + cos(sangle));
      improper_acc += static_cast<llint>(round(contrib * nrg_scale_factor));
      break;
    }

    // Compute forces, if requested
    if (eval_force == EvaluateForce::YES) {
      const double fr = ampl * freq * sin(sangle);
      const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      const double invab = 1.0/mgab;
      const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      const double invbc = 1.0/mgbc;
      const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      const double invcd = 1.0/mgcd;
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
      psw.xfrc[i_atom] += crabbc[0] * fa;
      psw.xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      psw.xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
      psw.xfrc[l_atom] -= fd * crbccd[0];
      psw.yfrc[i_atom] += crabbc[1] * fa;
      psw.yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      psw.yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
      psw.yfrc[l_atom] -= fd * crbccd[1];
      psw.zfrc[i_atom] += crabbc[2] * fa;
      psw.zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      psw.zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
      psw.zfrc[l_atom] -= fd * crbccd[2];
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::PROPER_DIHEDRAL, proper_acc, system_index);
  ecard->contribute(StateVariable::IMPROPER_DIHEDRAL, improper_acc, system_index);

  // Return the double-precision dihedral energy sum, if of interest
  return dihe_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();

  // Accumulate the results in two numerical precision models by looping over all terms
  double ubrd_energy = 0.0;
  llint ubrd_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nubrd; pos++) {
    const int i_atom = vk.ubrd_i_atoms[pos];
    const int k_atom = vk.ubrd_k_atoms[pos];
    const int param_idx = vk.ubrd_param_idx[pos];
    const double keq = vk.ubrd_keq[param_idx];
    const double leq = fabs(vk.ubrd_leq[param_idx]);
    const double dx = psw.xcrd[k_atom] - psw.xcrd[i_atom];
    const double dy = psw.ycrd[k_atom] - psw.ycrd[i_atom];
    const double dz = psw.zcrd[k_atom] - psw.zcrd[i_atom];
    const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
    const double dl = dr - leq;
    const double du = keq * dl * dl;
    ubrd_energy += du;
    ubrd_acc += static_cast<llint>(round(du * nrg_scale_factor));

    // Compute forces
    if (eval_force == EvaluateForce::YES) {
      const double fmag = 2.0 * keq * dl / dr;
      const double fmag_dx = fmag * dx;
      const double fmag_dy = fmag * dy;
      const double fmag_dz = fmag * dz;
      psw.xfrc[i_atom] += fmag_dx;
      psw.yfrc[i_atom] += fmag_dy;
      psw.zfrc[i_atom] += fmag_dz;
      psw.xfrc[k_atom] -= fmag_dx;
      psw.yfrc[k_atom] -= fmag_dy;
      psw.zfrc[k_atom] -= fmag_dz;
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::UREY_BRADLEY, ubrd_acc, system_index);

  // Return the double-precision Urey-Bradley energy sum, if of interest
  return ubrd_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();
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
    ab[0] = psw.xcrd[j_atom] - psw.xcrd[i_atom];
    ab[1] = psw.ycrd[j_atom] - psw.ycrd[i_atom];
    ab[2] = psw.zcrd[j_atom] - psw.zcrd[i_atom];
    bc[0] = psw.xcrd[k_atom] - psw.xcrd[j_atom];
    bc[1] = psw.ycrd[k_atom] - psw.ycrd[j_atom];
    bc[2] = psw.zcrd[k_atom] - psw.zcrd[j_atom];
    cd[0] = psw.xcrd[l_atom] - psw.xcrd[k_atom];
    cd[1] = psw.ycrd[l_atom] - psw.ycrd[k_atom];
    cd[2] = psw.zcrd[l_atom] - psw.zcrd[k_atom];

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
    cimp_acc += static_cast<llint>(round(contrib * nrg_scale_factor));

    // Compute forces
    if (eval_force == EvaluateForce::YES) {
      const double fr = -2.0 * keq * dtheta;
      const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      const double invab = 1.0/mgab;
      const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      const double invbc = 1.0/mgbc;
      const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      const double invcd = 1.0/mgcd;
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
      psw.xfrc[i_atom] += crabbc[0] * fa;
      psw.xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      psw.xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
      psw.xfrc[l_atom] -= fd * crbccd[0];
      psw.yfrc[i_atom] += crabbc[1] * fa;
      psw.yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      psw.yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
      psw.yfrc[l_atom] -= fd * crbccd[1];
      psw.zfrc[i_atom] += crabbc[2] * fa;
      psw.zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      psw.zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
      psw.zfrc[l_atom] -= fd * crbccd[2];
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::CHARMM_IMPROPER, cimp_acc, system_index);

  // Return the double-precision CHARMM improper energy sum, if of interest
  return cimp_energy;
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  PhaseSpaceWriter psw = ps->data();
  double cmap_energy = 0.0;
  llint cmap_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();
  std::vector<double> acoef(16);

  // Accumulate results by looping over all CMAP (here, meaning coupled dihedral) terms.
  double ab[3], bc[3], cd[3], de[3], crabbc[3], crbccd[3], crcdde[3], scr_phi[3], scr_psi[3];
  double phi_progression[4], psi_progression[4], acoef_psi[4];
  for (int i = 0; i < 4; i++) {
    acoef_psi[i] = 0.0;
  }
  for (int pos = 0; pos < vk.ncmap; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int i_atom = vk.cmap_i_atoms[pos];
    const int j_atom = vk.cmap_j_atoms[pos];
    const int k_atom = vk.cmap_k_atoms[pos];
    const int l_atom = vk.cmap_l_atoms[pos];
    const int m_atom = vk.cmap_m_atoms[pos];
    const int surf_idx = vk.cmap_surf_idx[pos];
    const int surf_dim = vk.cmap_dim[surf_idx];

    // Compute displacements
    ab[0] = psw.xcrd[j_atom] - psw.xcrd[i_atom];
    ab[1] = psw.ycrd[j_atom] - psw.ycrd[i_atom];
    ab[2] = psw.zcrd[j_atom] - psw.zcrd[i_atom];
    bc[0] = psw.xcrd[k_atom] - psw.xcrd[j_atom];
    bc[1] = psw.ycrd[k_atom] - psw.ycrd[j_atom];
    bc[2] = psw.zcrd[k_atom] - psw.zcrd[j_atom];
    cd[0] = psw.xcrd[l_atom] - psw.xcrd[k_atom];
    cd[1] = psw.ycrd[l_atom] - psw.ycrd[k_atom];
    cd[2] = psw.zcrd[l_atom] - psw.zcrd[k_atom];
    de[0] = psw.xcrd[m_atom] - psw.xcrd[l_atom];
    de[1] = psw.ycrd[m_atom] - psw.ycrd[l_atom];
    de[2] = psw.zcrd[m_atom] - psw.zcrd[l_atom];

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
    const int patch_idx = vk.cmap_patch_bounds[surf_idx] + (((idx_psi * surf_dim) + idx_phi) * 16);
    for (int i = 0; i < 16; i++) {
      acoef[i] = vk.cmap_patches[patch_idx + i];
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
    cmap_energy += contrib;
    cmap_acc += static_cast<llint>(round(contrib * nrg_scale_factor));

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
      const double invab = 1.0/mgab;
      const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      const double invbc = 1.0/mgbc;
      const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      const double invcd = 1.0/mgcd;
      const double mgde = sqrt(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
      const double invde = 1.0/mgde;
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
      psw.xfrc[i_atom] += crabbc[0] * phi_fa;
      psw.yfrc[i_atom] += crabbc[1] * phi_fa;
      psw.zfrc[i_atom] += crabbc[2] * phi_fa;
      psw.xfrc[j_atom] += (phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0]);
      psw.yfrc[j_atom] += (phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1]);
      psw.zfrc[j_atom] += (phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2]);
      psw.xfrc[k_atom] += (phi_fc2 * crabbc[0]) - (phi_fc1 * crbccd[0]);
      psw.yfrc[k_atom] += (phi_fc2 * crabbc[1]) - (phi_fc1 * crbccd[1]);
      psw.zfrc[k_atom] += (phi_fc2 * crabbc[2]) - (phi_fc1 * crbccd[2]);
      psw.xfrc[l_atom] -= phi_fd * crbccd[0];
      psw.yfrc[l_atom] -= phi_fd * crbccd[1];
      psw.zfrc[l_atom] -= phi_fd * crbccd[2];

      // Accumulate the psi dihedral forces
      psw.xfrc[j_atom] += crbccd[0] * psi_fa;
      psw.yfrc[j_atom] += crbccd[1] * psi_fa;
      psw.zfrc[j_atom] += crbccd[2] * psi_fa;
      psw.xfrc[k_atom] += (psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0]);
      psw.yfrc[k_atom] += (psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1]);
      psw.zfrc[k_atom] += (psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2]);
      psw.xfrc[l_atom] += (psi_fc2 * crbccd[0]) - (psi_fc1 * crcdde[0]);
      psw.yfrc[l_atom] += (psi_fc2 * crbccd[1]) - (psi_fc1 * crcdde[1]);
      psw.zfrc[l_atom] += (psi_fc2 * crbccd[2]) - (psi_fc1 * crcdde[2]);
      psw.xfrc[m_atom] -= psi_fd * crcdde[0];
      psw.yfrc[m_atom] -= psi_fd * crcdde[1];
      psw.zfrc[m_atom] -= psi_fd * crcdde[2];
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::CMAP, cmap_acc, system_index);

  // Return the double-precision CMAP energy sum, if of interest
  return cmap_energy;
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  PhaseSpaceWriter psw = ps->data();
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
    const int i_atom = vk.dihe_i_atoms[pos];
    const int j_atom = vk.dihe_l_atoms[pos];
    const int ilj_t = nbk.lj_idx[i_atom];
    const int jlj_t = nbk.lj_idx[j_atom];
    const double dx = psw.xcrd[j_atom] - psw.xcrd[i_atom];
    const double dy = psw.ycrd[j_atom] - psw.ycrd[i_atom];
    const double dz = psw.zcrd[j_atom] - psw.zcrd[i_atom];
    const double invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
    const double invr = sqrt(invr2);
    const double invr4 = invr2 * invr2;
    const double ele_scale = vk.attn14_elec[attn_idx];
    const double vdw_scale = vk.attn14_vdw[attn_idx];
    const double qiqj = (nbk.coulomb_constant * nbk.charge[i_atom] * nbk.charge[j_atom]) /
                        ele_scale;
    const double lja = nbk.lja_14_coeff[(ilj_t * nbk.n_lj_types) + jlj_t] / vdw_scale;
    const double ljb = nbk.ljb_14_coeff[(ilj_t * nbk.n_lj_types) + jlj_t] / vdw_scale;
    const double ele_contrib = qiqj * invr;
    const double vdw_contrib = (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);
    ele_energy += ele_contrib;
    vdw_energy += vdw_contrib;
    ele_acc += static_cast<llint>(round(ele_contrib * nrg_scale_factor));
    vdw_acc += static_cast<llint>(round(vdw_contrib * nrg_scale_factor));

    // Evaluate the force, if requested
    if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
      double fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) : 0.0;
      if (eval_vdw_force == EvaluateForce::YES) {
        fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
      }
      const double fmag_dx = fmag * dx;
      const double fmag_dy = fmag * dy;
      const double fmag_dz = fmag * dz;
      psw.xfrc[i_atom] += fmag_dx;
      psw.yfrc[i_atom] += fmag_dy;
      psw.zfrc[i_atom] += fmag_dz;
      psw.xfrc[j_atom] -= fmag_dx;
      psw.yfrc[j_atom] -= fmag_dy;
      psw.zfrc[j_atom] -= fmag_dz;
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::ELECTROSTATIC_ONE_FOUR, ele_acc, system_index);
  ecard->contribute(StateVariable::VDW_ONE_FOUR, vdw_acc, system_index);

  // Return the double-precision energy sums, if of interest
  return { ele_energy, vdw_energy };
}

} // namespace energy
} // namespace omni
