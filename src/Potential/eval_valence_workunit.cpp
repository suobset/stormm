#include "eval_valence_workunit.h"

namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const double* xcrd, const dobule* ycrd, const double* zcrd,
                          const double* umat, const double* invu, const UnitCellType unit_cell,
                          double* xfrc, double* yfrc, double* zfrc, ScoreCard *sc,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {  
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
    if (activity == VwuTask::BOND || activitiy == VwuTask::ALL_TASKS) {
      for (int pos = 0; pos < task_counts[VwuTask::BOND]; pos++) {

      }
    }
    
    // Evaluate harmonic bond angles
    if (activity == VwuTask::ANGL || activitiy == VwuTask::ALL_TASKS) {
      for (int pos = 0; pos < task_counts[VwuTask::ANGL]; pos++) {        
      }
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
void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *sc, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity) {
  PhaseSpaceWriter psw = ps.data();
  evalValenceWorkUnits(ag->getDoublePrecisionValenceKit(), ag->getDoublePrecisionVirtualSiteKit(),
                       psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                       psw.xfrc, psw.yfrc, psw.zfrc, sc, sysid, vwu_list, eval_force, activity);
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, ScoreCard *sc,
                          const int sysid, const std::vector<ValenceWorkUnit> &vwu_list,
                          const VwuTask activity) {
  PhaseSpaceReader psr = ps.data();
  evalValenceWorkUnits(ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionVirtualSiteKit(),
                       psr.xcrd, psr.ycrd, psr.zcrd, psr.umat, psr.invu, psr.unit_cell,
                       nullptr, nullptr, nullptr, sc, sysid, vwu_list, EvaluateForce::NO,
                       activity);
}

} // namespace energy
} // namespace omni
