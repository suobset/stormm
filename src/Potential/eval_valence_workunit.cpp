#include "Math/vector_ops.h"
#include "eval_valence_workunit.h"

namespace omni {
namespace energy {

using math::crossProduct;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
  
//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const double* xcrd, const double* ycrd, const double* zcrd,
                          const double* umat, const double* invu, const UnitCellType unit_cell,
                          double* xfrc, double* yfrc, double* zfrc, ScoreCard *ecard,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {  
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();
  for (size_t vidx = 0LLU; vidx < vwu_list.size(); vidx++) {

    // Get all tasks from the work unit
    const std::vector<int> task_counts = vwu_list[vidx].getTaskCounts();

    // Make local arrays of the coordinates and forces.  Import coordinates and extant forces
    // from the global arrays.
    const int max_atoms = vwu_list[vidx].getMaxAtoms();
    std::vector<double> sh_xcrd(max_atoms), sh_ycrd(max_atoms), sh_zcrd(max_atoms);
    std::vector<double> sh_xfrc(max_atoms), sh_yfrc(max_atoms), sh_zfrc(max_atoms);
    const int n_imp_atoms = vwu_list[vidx].getImportedAtomCount();
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      sh_xcrd[i] = xcrd[atom_idx];
      sh_ycrd[i] = ycrd[atom_idx];
      sh_zcrd[i] = zcrd[atom_idx];
      sh_xfrc[i] = xfrc[atom_idx];
      sh_yfrc[i] = yfrc[atom_idx];
      sh_zfrc[i] = zfrc[atom_idx];
    }
    
    // Evaluate bonds
    if (activity == VwuTask::BOND || activity == VwuTask::CBND || activity == VwuTask::ALL_TASKS) {
      llint bond_acc = 0LL;
      const int nbond = task_counts[static_cast<int>(VwuTask::CBND)];
      for (int pos = 0; pos < nbond; pos++) { 
        const uint2 tinsr = vwu_list[vidx].getCompositeBondInstruction(pos);
        const int i_atom = (tinsr.x & 0x3ff);
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);

        // Skip Urey-Bradley interactions in this loop
        if (activity == VwuTask::BOND && ((tinsr.x >> 20) & 0x1)) {
          continue;
        }
        const int param_idx = tinsr.y;
        const double keq = vk.bond_keq[param_idx];
        const double leq = fabs(vk.bond_leq[param_idx]);
        const double dx = sh_xcrd[j_atom] - sh_xcrd[i_atom];
        const double dy = sh_ycrd[j_atom] - sh_ycrd[i_atom];
        const double dz = sh_zcrd[j_atom] - sh_zcrd[i_atom];
        const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const double dl = dr - leq;
        const double du = keq * dl * dl;
        bond_acc += static_cast<llint>(round(du * nrg_scale_factor));

        // Compute forces
        if (eval_force == EvaluateForce::YES) {
          const double fmag = 2.0 * keq * dl / dr;
          const double fmag_dx = fmag * dx;
          const double fmag_dy = fmag * dy;
          const double fmag_dz = fmag * dz;
          sh_xfrc[i_atom] += fmag_dx;
          sh_yfrc[i_atom] += fmag_dy;
          sh_zfrc[i_atom] += fmag_dz;
          sh_xfrc[j_atom] -= fmag_dx;
          sh_yfrc[j_atom] -= fmag_dy;
          sh_zfrc[j_atom] -= fmag_dz;
        }
      }

      // Contribute results
      ecard->contribute(StateVariable::BOND, bond_acc, sysid);
    }
    
    // Evaluate harmonic bond angles
    if (activity == VwuTask::ANGL || activity == VwuTask::ALL_TASKS) {
      llint angl_acc = 0LL;
      const int nangl = task_counts[static_cast<int>(VwuTask::ANGL)];
      for (int pos = 0; pos < nangl; pos++) {
        const uint2 tinsr = vwu_list[vidx].getAngleInstruction(pos);
        const int i_atom = (tinsr.x & 0x3ff);
        const int j_atom = ((tinsr.x >> 10) & 0x3ff);
        const int k_atom = ((tinsr.x >> 20) & 0x3ff);
        const int param_idx = tinsr.y;
        const double keq = vk.angl_keq[param_idx];
        const double theta0 = vk.angl_theta[param_idx];

        // Compute displacements
        double ba[3], bc[3];
        ba[0] = sh_xcrd[i_atom] - sh_xcrd[j_atom];
        ba[1] = sh_ycrd[i_atom] - sh_ycrd[j_atom];
        ba[2] = sh_zcrd[i_atom] - sh_zcrd[j_atom];
        bc[0] = sh_xcrd[k_atom] - sh_xcrd[j_atom];
        bc[1] = sh_ycrd[k_atom] - sh_ycrd[j_atom];
        bc[2] = sh_zcrd[k_atom] - sh_zcrd[j_atom];

        // On to the angle force computation
        const double mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
        const double mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
        const double invbabc = 1.0 / sqrt(mgba * mgbc);
        double costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
        costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
        const double theta = acos(costheta);
        const double dtheta = theta - theta0;
        const double du = keq * dtheta * dtheta;
        angl_acc += static_cast<llint>(round(du * nrg_scale_factor));

        // Compute forces
        if (eval_force == EvaluateForce::YES) {
          const double dA = -2.0 * keq * dtheta / sqrt(1.0 - costheta * costheta);
          const double sqba = dA / mgba;
          const double sqbc = dA / mgbc;
          const double mbabc = dA * invbabc;
          double adf[3], cdf[3];
          for (int i = 0; i < 3; i++) {
            adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
            cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
          }
          sh_xfrc[i_atom] -= adf[0];
          sh_yfrc[i_atom] -= adf[1];
          sh_zfrc[i_atom] -= adf[2];
          sh_xfrc[j_atom] += adf[0] + cdf[0];
          sh_yfrc[j_atom] += adf[1] + cdf[1];
          sh_zfrc[j_atom] += adf[2] + cdf[2];
          sh_xfrc[k_atom] -= cdf[0];
          sh_yfrc[k_atom] -= cdf[1];
          sh_zfrc[k_atom] -= cdf[2];
        }
      }

      // Contribute results
      ecard->contribute(StateVariable::ANGLE, angl_acc, sysid);
    }
    
    // Add accumulated forces back to the global arrays (this is not done by all GPU kernels, as
    // in some cases the ValenceWorkUnits also move atoms and then leave the global force arrays
    // zero'ed for atoms they are responsible for moving).
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      xfrc[atom_idx] += sh_xfrc[i];
      yfrc[atom_idx] += sh_yfrc[i];
      zfrc[atom_idx] += sh_zfrc[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {
  PhaseSpaceWriter psw = ps->data();
  evalValenceWorkUnits(ag->getDoublePrecisionValenceKit(), ag->getDoublePrecisionVirtualSiteKit(),
                       psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                       psw.xfrc, psw.yfrc, psw.zfrc, ecard, sysid, vwu_list, eval_force, activity);
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, ScoreCard *ecard,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          const VwuTask activity) {
  PhaseSpaceReader psr = ps.data();
  evalValenceWorkUnits(ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionVirtualSiteKit(),
                       psr.xcrd, psr.ycrd, psr.zcrd, psr.umat, psr.invu, psr.unit_cell,
                       nullptr, nullptr, nullptr, ecard, sysid, vwu_list, EvaluateForce::NO,
                       activity);
}

} // namespace energy
} // namespace omni
