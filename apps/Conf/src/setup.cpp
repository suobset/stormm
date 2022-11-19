#include <cmath>
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Math/tickcounter.h"
#include "../../../src/MoleculeFormat/molecule_parsing.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Topology/topology_util.h"
#include "../../../src/Trajectory/coordinate_series.h"
#include "setup.h"

// CHECK
#include "../../../src/DataTypes/stormm_vector_types.h"
// END CHECK

namespace conf_app {
namespace setup {

using stormm::chemistry::ConformationEdit;
using stormm::chemistry::MapRotatableGroups;
using stormm::constants::warp_size_int;
using stormm::math::roundUp;
using stormm::math::TickCounter;
using stormm::namelist::ConformerControls;
using stormm::structure::maskFromSdfDataItem;
using stormm::structure::rmsd;
using stormm::structure::RmsdMethod;
using stormm::structure::rotateAboutBond;
using stormm::structure::flipChiralCenter;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::isBonded;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::CoordinateSeries;
using stormm::trajectory::CoordinateSeriesWriter;

//-------------------------------------------------------------------------------------------------
AtomMask getCoreMask(const ConformerControls &conf_input, const MdlMol &sdf_example,
                     const PhaseSpace &ps, const AtomGraph *ag, const ChemicalFeatures &chemfe,
                     const ExceptionResponse policy) {
  AtomMask result(ag);
  bool core_found = false;
  if (conf_input.getCoreDataItemName().size() > 0LLU && sdf_example.getAtomCount() > 0) {
    result.addAtoms(maskFromSdfDataItem(conf_input.getCoreDataItemName(), sdf_example, ag, chemfe,
                                        sdf_example.exportCoordinateFrame(), policy),
                    sdf_example.exportCoordinateFrame(), chemfe);
    core_found = (result.getMaskedAtomCount() > 0);
  }
  if (core_found == false && conf_input.getCoreAtomMask().size() > 0LLU) {

    // If present, a core atom mask common to all systems will fill in for situations in which a
    // list of core atoms is not defined in an SD file entry.
    result.addAtoms(conf_input.getCoreAtomMask(), CoordinateFrame(ps), chemfe);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomMask> setGenerativeConditions(const UserSettings &ui, SystemCache *sc,
                                              const std::vector<MdlMol> &sdf_recovery,
                                              StopWatch *tm) {

  // Establish a category for this function and stash time up to this point in the "miscellaneous"
  // category.
  const int cat_no = tm->addCategory("Set generative conditions");
  tm->assignTime(0);
  
  // Add implicit solvent conditions.
  std::vector<AtomGraph*> all_topologies = sc->getTopologyPointer();
  const size_t ntop = all_topoligies.size();
  const ImplicitSolventModel igb = ui.getSolventNamelistInfo().getImplicitSolventModel();
  if (igb != ImplicitSolventModel::NONE) {
    for (size_t i = 0LLU; i < ntop; i++) {
      all_topologies[i]->setImplicitSolventModel(igb);
    }
  }

  // Create the core masks for all topologies.  Build a list of chemical features for each unique
  // topology in the systems cache.
  const ConformerControls& conf_input = ui.getConformerNamelistInfo();  
  std::vector<AtomMask> result;
  result.reserve(ntop);
  for (size_t i = 0LLU; i < ntop; i++) {
    AtomGraph* iag_ptr = sc->getTopologyPointer(i);
    const int example_system_idx = sc.getCoordinateExample(i);
    result.push_back(getCoreMask(conf_input, sdf_recovery[example_system_idx],
                                 sc->getCoordinateReference(example_system_idx),
                                 iag_ptr, sc->getFeaturesReference(i), ui.getExceptionBehavior()));
    iag_ptr->modifyAtomMobility(result.back(), MobilitySetting::OFF);

    // With the topology's core atoms now frozen, compute the rotatable bond groups.
    sc->getFeaturesPointer(i)->
  }

  // Record the time needed for this procedure.
  tm->assignTime(cat_no);
  return result;
}
  
//-------------------------------------------------------------------------------------------------
bool permutationsAreLinked(const std::vector<IsomerPlan> &isomerizers, const int permi,
                           const int permj, const NonbondedKit<double> &nbk) {
  const int root_i = isomerizers[permi].getRootAtom();
  const int pivt_i = isomerizers[permi].getPivotAtom();
  const int root_j = isomerizers[permj].getRootAtom();
  const int pivt_j = isomerizers[permj].getPivotAtom();
  switch (isomerizers[permi].getMotion()) {
  case ConformationEdit::BOND_ROTATION:
  case ConformationEdit::CIS_TRANS_FLIP:

    // If two isomerizations share atoms or have their root and pivot atoms bonded to one
    // another, consider them coupled and add the product of the number of possible states for
    // each to the running sum.
    switch (isomerizers[permj].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      return (root_i == root_j || root_i == pivt_j || pivt_i == root_j || pivt_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j) ||
              isBonded(nbk, pivt_i, root_j) || isBonded(nbk, pivt_i, pivt_j));
    case ConformationEdit::CHIRAL_INVERSION:
      return (root_i == root_j || pivt_i == root_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j));
    }
    break;
  case ConformationEdit::CHIRAL_INVERSION:
    switch (isomerizers[permj].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      return (root_i == root_j || root_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j));
    case ConformationEdit::CHIRAL_INVERSION:
      return (root_i == root_j || isBonded(nbk, root_i, root_j));
    }
    break;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
double computeLocalPermutations(const std::vector<int> &limits,
                                const std::vector<IsomerPlan> &isomerizers, const AtomGraph *ag) {
  double result = 0.0;
  const int n_var = limits.size();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  for (int i = 0; i < n_var; i++) {
    result += limits[i];
    for (int j = i + 1; j < n_var; j++) {
      if (permutationsAreLinked(isomerizers, i, j, nbk)) {
        result += limits[i] * limits[j];
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void forgeConformation(CoordinateSeries<double> *cseries, const int fc,
                       const std::vector<IsomerPlan> &isomerizers, const TickCounter &ptrack,
                       const std::vector<int> &chiral_centers,
                       const std::vector<ChiralInversionProtocol> &chiral_center_plans,
                       const std::vector<int> &chiral_variable_indices,
                       const std::vector<IsomerPlan> &invertible_groups) {
  const int n_variables = ptrack.getStateCount();
  CoordinateSeriesWriter<double> cseries_w   = cseries->data();
  const std::vector<int>& permutation_watch  = ptrack.getSettings();
  const std::vector<int>& permutation_limits = ptrack.getStateLimits();
  for (int i = 0; i < n_variables; i++) {
    switch (isomerizers[i].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      if (permutation_watch[i] > 0) {
        const double rval = static_cast<double>(permutation_watch[i]) /
                            static_cast<double>(permutation_limits[i]) * stormm::symbols::twopi;
        rotateAboutBond<double, double>(cseries_w, fc, isomerizers[i].getRootAtom(),
                                        isomerizers[i].getPivotAtom(),
                                        isomerizers[i].getMovingAtoms(), rval);
      }
      break;
    case ConformationEdit::CHIRAL_INVERSION:
      if (permutation_watch[i] == 1) {
        switch (isomerizers[i].getChiralPlan()) {
        case ChiralInversionProtocol::ROTATE:
        case ChiralInversionProtocol::REFLECT:
          flipChiralCenter<double, double>(cseries_w, fc, chiral_variable_indices[i],
                                           chiral_centers, chiral_center_plans, invertible_groups);
          break;
        case ChiralInversionProtocol::DO_NOT_INVERT:
          break;
        }
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc, 
                                     const std::vector<MdlMol> &sdf_recovery,
                                     Xoshiro256ppGenerator *xrs, StopWatch *tm) {
  const ConformerControls& conf_input = ui.getConformerNamelistInfo();

  // Count the expanded number of systems
  const int ntop = sc.getTopologyCount();
  const int nsys = sc.getSystemCount();
  const std::vector<const AtomGraph*> system_topologies = sc.getSystemTopologyPointer();
  const std::vector<const PhaseSpace*> system_coords    = sc.getCoordinatePointer();
  const std::vector<const AtomGraph*> unique_topologies = sc.getTopologyPointer();
  std::vector<ChemicalFeatures> chemfe_list;
  std::vector<AtomMask> core_masks;
  chemfe_list.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    const int example_system_idx = sc.getCoordinateExample(i);
    const ChemicalFeatures tmp_chemfe(sc.getTopologyPointer(i),
                                      sc.getCoordinateReference(example_system_idx),
                                      MapRotatableGroups::NO);
    core_masks.push_back(getCoreMask(conf_input, sdf_recovery[example_system_idx],
                                     sc.getCoordinateReference(example_system_idx),
                                     sc.getTopologyPointer(i), tmp_chemfe,
                                     ui.getExceptionBehavior()));
    chemfe_list.emplace_back(sc.getTopologyPointer(i),
                                sc.getCoordinateReference(example_system_idx),
                                MapRotatableGroups::NO)
  }
  tm->assignTime(2);

  // CHECK
  for (int i = 0; i < ntop; i++) {
    const AtomGraph *iag_ptr = sc.getTopologyPointer(i);
    printf("Topology: %s\n", iag_ptr->getFileName().c_str());
    const std::vector<int> masked_atoms = core_masks[i].getMaskedAtomList();
    for (int j = 0; j < core_masks[i].getMaskedAtomCount(); j++) {
      const stormm::char4 atom_name = iag_ptr->getAtomName(masked_atoms[j]);
      printf("  %c%c%c%c [ %2d ]\n", atom_name.x, atom_name.y, atom_name.z, atom_name.w,
             masked_atoms[j]);
    }
  }
  // END CHECK
  
  // Create lists of PhaseSpace objects and topology pointers to show how to model each of them
  std::vector<PhaseSpace> ps_list; 
  std::vector<AtomGraph*> ag_list;
  int conf_counter = 0;

  // Loop over all systems, grouping those with the same topology into a coherent group of
  // proto-conformers for coarse-grained sampling of rotatable bonds and chiral centers.
  const int nbond_rotations = conf_input.getRotationSampleCount();
  int nneighborhood = 0;
  for (int i = 0; i < ntop; i++) {

    // Obtain the topology pointer for possible use later
    const AtomGraph *iag_ptr = sc.getTopologyPointer(i);
    
    // Create a list of ways to change the conformation
    const std::vector<IsomerPlan> rotatable_groups  = chemfe_list[i].getRotatableBondGroups();
    const std::vector<IsomerPlan> cis_trans_groups  =
      chemfe_list[i].getCisTransIsomerizationGroups();
    const std::vector<IsomerPlan> invertible_groups = chemfe_list[i].getChiralInversionGroups();
    std::vector<IsomerPlan> isomerizers;
    const int nrot_bond = rotatable_groups.size();
    const int nctx_bond = cis_trans_groups.size();
    const int nchiral   = invertible_groups.size();
    const int n_variables = nrot_bond + nctx_bond + nchiral;
    isomerizers.reserve(n_variables);
    std::vector<int> permutation_states;
    permutation_states.reserve(n_variables);
    for (int j = 0; j < nrot_bond; j++) {
      isomerizers.emplace_back(rotatable_groups[j]);
      permutation_states.push_back(nbond_rotations);
    }
    for (int j = 0; j < nctx_bond; j++) {
      isomerizers.emplace_back(cis_trans_groups[j]);
      permutation_states.push_back(2);
    }
    for (int j = 0; j < nchiral; j++) {
      isomerizers.emplace_back(invertible_groups[j]);
      switch (invertible_groups[j].getChiralPlan()) {
      case ChiralInversionProtocol::ROTATE:
      case ChiralInversionProtocol::REFLECT:
        permutation_states.push_back(2);
        break;
      case ChiralInversionProtocol::DO_NOT_INVERT:
        permutation_states.push_back(1);
        break;
      }
    }
    TickCounter ptrack(permutation_states);
    
    // This reference, while const, will allow the program to see into the permutation tracking
    // TickCounter object as the state updates with each iteration.
    const std::vector<int> permutation_limits = ptrack.getStateLimits();
    const double ln_choices = ptrack.getLogPermutationCount();
    SamplingStrategy strat;
    int conformations_per_case;
    if (ln_choices < log(conf_input.getSystemTrialCount())) {
      conformations_per_case = ptrack.getExactPermutationCount();
      strat = SamplingStrategy::FULL;
    }
    else if (computeLocalPermutations(permutation_limits, isomerizers, iag_ptr) <
             static_cast<double>(conf_input.getSystemTrialCount())) {
      conformations_per_case = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::LIMITED;
    }
    else {
      conformations_per_case = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::SPARSE;
    }

    // Create a series to hold the conformers resulting from the coarse-grained search.  Populate
    // the series with copies of each case matching this topology, in sufficient quantities to
    // transform the original conformation into every permutation of its isomeric features.
    CoordinateSeries<double> cseries(sc.getTopologyPointer(i)->getAtomCount(),
                                     conformations_per_case);
    const std::vector<int> top_cases = sc.getTopologicalCases(i);
    const int ncases = top_cases.size();
    cseries.resize(0);
    int fc = 0;
    for (int j = 0; j < ncases; j++) {
      fc += conformations_per_case;
      cseries.resize(fc, sc.getCoordinateReference(top_cases[j]));
    }
    
    // Get the vectors of chiral centers and inversion plans for all centers in this molecule
    const std::vector<int> chiral_centers = chemfe_list[i].getChiralCenters();
    const std::vector<ChiralInversionProtocol> chiral_center_plans =
      chemfe_list[i].getChiralInversionMethods();
    std::vector<int> chiral_variable_indices(n_variables, -1);
    int ncen_found = 0;
    for (int j = 0; j < n_variables; j++) {
      switch (isomerizers[j].getMotion()) {
      case ConformationEdit::BOND_ROTATION:
      case ConformationEdit::CIS_TRANS_FLIP:
        break;
      case ConformationEdit::CHIRAL_INVERSION:
        chiral_variable_indices[j] = ncen_found;
        ncen_found++;
      }
    }
    
    // Manipulate this coordinate series using bond rotations as well as chiral inversions.
    // Compute the RMSD matrix and determine a set of diverse conformations.  Cull the results
    // to eliminate ring stabs or severe clashes between tertiary or quaternary atoms.
    fc = 0;
    CoordinateSeriesWriter<double> cseries_w = cseries.data();
    switch (strat) {
    case SamplingStrategy::FULL:
      for (int j = 0; j < ncases; j++) {
        ptrack.reset();
        for (int k = 0; k < conformations_per_case; k++) {
          forgeConformation(&cseries, fc, isomerizers, ptrack, chiral_centers, chiral_center_plans,
                            chiral_variable_indices, invertible_groups);
          fc++;
          ptrack.advance();
        }
        break;
      }
      break;
    case SamplingStrategy::LIMITED:
      {
        // In limited sampling, coupled isomerizing variables will be explored in concert.  In
        // sparse sampling, each variable will be explored in at least one case, provided that
        // there are enough samples of the various permutations to fill.
        const NonbondedKit<double> nbk = iag_ptr->getDoublePrecisionNonbondedKit();
        const int conf_goal = ncases * conformations_per_case;
        while (fc < conf_goal) {
          int k = 0;
          while (k < n_variables && fc < conf_goal) {
            int m = k + 1;
            while (m < n_variables && fc < conf_goal) {
              if (permutationsAreLinked(isomerizers, k, m, nbk)) {
                int pos_k = 0;
                while (pos_k < permutation_limits[k] && fc < conf_goal) {
                  int pos_m = 0;
                  while (pos_m < permutation_limits[m] && fc < conf_goal) {
                    ptrack.randomize(xrs);
                    ptrack.set(pos_k, k);
                    ptrack.set(pos_m, m);
                    forgeConformation(&cseries, fc, isomerizers, ptrack, chiral_centers,
                                      chiral_center_plans, chiral_variable_indices,
                                      invertible_groups);
                    fc++;
                    pos_m++;
                  }
                  pos_k++;
                }
              }
              m++;
            }
            k++;
          }
        }
      }
      break;
    case SamplingStrategy::SPARSE:
      {
        const NonbondedKit<double> nbk = iag_ptr->getDoublePrecisionNonbondedKit();
        const int conf_goal = ncases * conformations_per_case;
        while (fc < conf_goal) {
          int k = 0;
          while (k < n_variables && fc < conf_goal) {
            int m = 0;
            while (m < permutation_limits[k] && fc < conf_goal) {
              ptrack.randomize(xrs);
              ptrack.set(m, k);
              forgeConformation(&cseries, fc, isomerizers, ptrack, chiral_centers,
                                chiral_center_plans, chiral_variable_indices,
                                invertible_groups);
              fc++;
              m++;
            }
            k++;
          }
        }
      }
      break;
    }
    for (int j = 0; j < ncases * conformations_per_case; j++) {
      ps_list.push_back(cseries.exportPhaseSpace(j));
      ag_list.push_back(const_cast<AtomGraph*>(unique_topologies[i]));
    }
  }
  return PhaseSpaceSynthesis(ps_list, ag_list);
}

} // namespace setup
} // namespace conf_app
