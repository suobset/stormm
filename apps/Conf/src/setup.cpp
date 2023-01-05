#include <cmath>
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Math/summation.h"
#include "../../../src/MoleculeFormat/molecule_parsing.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Topology/topology_util.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::ConformationEdit;
using stormm::chemistry::MapRotatableGroups;
using stormm::constants::warp_size_int;
#ifndef STORMM_USE_HPC
using stormm::data_types::double2;
#endif
using stormm::diskutil::getBaseName;
using stormm::errors::rtErr;
using stormm::errors::rtWarn;
using stormm::math::prefixSumInPlace;
using stormm::math::PrefixSumType;
using stormm::math::roundUp;
using stormm::math::sum;
using stormm::math::loadScalarStateValues;
using stormm::namelist::ConformerControls;
using stormm::parse::char4ToString;
using stormm::structure::maskFromSdfDataItem;
using stormm::structure::rotateAboutBond;
using stormm::structure::flipChiralCenter;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::ImplicitSolventModel;
using stormm::topology::MobilitySetting;
using stormm::topology::isBonded;

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
void setGenerativeConditions(const UserSettings &ui, SystemCache *sc,
                             const std::vector<MdlMol> &sdf_recovery, StopWatch *tm) {

  // Establish a category for this function and stash time up to this point in the "miscellaneous"
  // category.
  const int cat_no = tm->addCategory("Set generative conditions");
  tm->assignTime(0);
  
  // Add implicit solvent conditions.
  std::vector<AtomGraph*> all_topologies = sc->getTopologyPointer();
  const size_t ntop = all_topologies.size();
  const ImplicitSolventModel igb = ui.getSolventNamelistInfo().getImplicitSolventModel();
  if (igb != ImplicitSolventModel::NONE) {
    for (size_t i = 0; i < ntop; i++) {
      all_topologies[i]->setImplicitSolventModel(igb);
    }
  }

  // Create the core masks for all topologies.  Build a list of chemical features for each unique
  // topology in the systems cache.
  const ConformerControls& conf_input = ui.getConformerNamelistInfo();  
  const int nsys = sc->getSystemCount();
  for (int i = 0; i < nsys; i++) {
    AtomGraph* iag_ptr = sc->getSystemTopologyPointer(i);
    const AtomMask imask = getCoreMask(conf_input, sdf_recovery[i], sc->getCoordinateReference(i),
                                       iag_ptr, sc->getFeaturesReference(i),
                                       ui.getExceptionBehavior());
    const std::vector<int> icore_atoms = imask.getMaskedAtomList();
    const size_t n_core_atoms = icore_atoms.size();
    std::vector<bool> icore_mobility(icore_atoms.size());
    for (size_t i = 0LLU; i < n_core_atoms; i++) {
      icore_mobility[i] = iag_ptr->getAtomMobility(icore_atoms[i]);
    }
    iag_ptr->modifyAtomMobility(icore_atoms, MobilitySetting::OFF);

    // With the topology's core atoms now frozen, compute the rotatable bond groups.
    sc->getFeaturesPointer(i)->findRotatableBondGroups(tm);

    // Restore the core atoms' original mobility settings (henceforth, they will be restrained to
    // their original positions using harmonic restraints).
    for (size_t i = 0; i < n_core_atoms; i++) {
      const MobilitySetting imb = (icore_mobility[i]) ? MobilitySetting::ON : MobilitySetting::OFF;
      iag_ptr->modifyAtomMobility(icore_atoms[i], imb);
    }
  }

  // Record the time needed for this procedure.
  tm->assignTime(cat_no);
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
CoordinateFrame forgeConformation(const CoordinateFrame &cf, const AtomGraph *ag,
                                  const StaticExclusionMask &exclusion_mask,
                                  const std::vector<IsomerPlan> &isomerizers,
                                  TickCounter<double> *ptrack, Xoshiro256ppGenerator *xrs,
                                  const std::vector<int> &chiral_centers,
                                  const std::vector<ChiralInversionProtocol> &chiral_center_plans,
                                  const std::vector<int> &chiral_variable_indices,
                                  const std::vector<IsomerPlan> &invertible_groups,
                                  const int max_seeding_attempts, const double self_clash_ratio) {
  const int n_variables = ptrack->getVariableCount();
  CoordinateFrame result = cf;
  int iter = 0;
  bool clash_detected = false;
  bool ptrack_altered = false;
  const std::vector<int> original_settings = ptrack->getSettings();
  do {
    const std::vector<int>& permutation_watch  = ptrack->getSettings();
    const std::vector<int>& permutation_limits = ptrack->getStateLimits();
    for (int i = 0; i < n_variables; i++) {
      switch (isomerizers[i].getMotion()) {
      case ConformationEdit::BOND_ROTATION:
      case ConformationEdit::CIS_TRANS_FLIP:
        if (permutation_watch[i] > 0) {
          const double rval = static_cast<double>(permutation_watch[i]) /
                              static_cast<double>(permutation_limits[i]) * stormm::symbols::twopi;
          rotateAboutBond(&result, isomerizers[i].getRootAtom(), isomerizers[i].getPivotAtom(),
                          isomerizers[i].getMovingAtoms(), rval);
        }
        break;
      case ConformationEdit::CHIRAL_INVERSION:
        if (permutation_watch[i] == 1) {
          switch (isomerizers[i].getChiralPlan()) {
          case ChiralInversionProtocol::ROTATE:
          case ChiralInversionProtocol::REFLECT:
            flipChiralCenter(&result, chiral_variable_indices[i], chiral_centers,
                             chiral_center_plans, invertible_groups);
            break;
          case ChiralInversionProtocol::DO_NOT_INVERT:
            break;
          }
        }
        break;
      }
    }

    // Randomize the configuration if the conformer contains a clash
#if 0
    if (detectVanDerWaalsClash(cf, ag, exclusion_mask, self_clash_ratio)) {
      clash_detected = true;
      ptrack->randomize(xrs);
      ptrack_altered = true;
    }
#endif
    iter++;
  } while (clash_detected && iter < max_seeding_attempts);

  // Reset the permutation tracker if necessary
  if (ptrack_altered) {
    ptrack->set(original_settings);
  }
  return (clash_detected) ? cf : result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> calculateReplicaCounts(const ConformerControls &conf_input,
                                        const SystemCache &sc, StopWatch *tm) {
  const double ln_rb = log(conf_input.getRotationSampleCount());
  const double ln_two = log(2.0);
  const double ln_max_states = log(conf_input.getSystemTrialCount());

  // Apply the topology-based counts to each system
  const int nsys = sc.getSystemCount();
  std::vector<int> result(nsys);
  for (int i = 0; i < nsys; i++) {
    const ChemicalFeatures& ichemfe = sc.getFeaturesReference(i);
    double ln_count = (ln_rb * static_cast<double>(ichemfe.getRotatableBondCount())) +
                      (ln_two * static_cast<double>(ichemfe.getCisTransBondCount()));
    const std::vector<ChiralInversionProtocol> chiral_ops = ichemfe.getChiralInversionMethods();
    const size_t nchir = chiral_ops.size();
    for (size_t j = 0; j < nchir; j++) {
      switch (chiral_ops[j]) {
      case ChiralInversionProtocol::ROTATE:
      case ChiralInversionProtocol::REFLECT:
        ln_count += ln_two;
        break;
      case ChiralInversionProtocol::DO_NOT_INVERT:
        break;
      }
    }
    if (ln_count < ln_max_states) {
      result[i] = ceil(exp(ln_count));
    }
    else {
      result[i] = conf_input.getSystemTrialCount();
    }
  }
  tm->assignTime(tm_coordinate_expansion);
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis buildReplicaWorkspace(const std::vector<int> &replica_counts,
                                          const SystemCache &sc, StopWatch *tm) {
  const int n_proto_sys = replica_counts.size();
  if (n_proto_sys != sc.getSystemCount()) {
    rtErr("Replica counts were provide for " + std::to_string(replica_counts.size()) + ", but the "
          "systems cache contains " + std::to_string(sc.getSystemCount()) + " systems.",
          "buildReplicaWorkspace");
  }
  const int n_total_rep = sum<int>(replica_counts);
  const std::vector<AtomGraph*> topology_pointers = sc.getSystemTopologyPointer();
  const std::vector<PhaseSpace>& proto_coords = sc.getCoordinateReference();
  std::vector<int> index_key(n_total_rep);
  int rep_con = 0;
  for (int i = 0; i < n_proto_sys; i++) {
    const int nrep = replica_counts[i];
    for (int j = 0; j < nrep; j++) {
      index_key[rep_con] = i;
      rep_con++;
    }
  }
  tm->assignTime(tm_coordinate_expansion);
  return PhaseSpaceSynthesis(proto_coords, topology_pointers, index_key);
}

//-------------------------------------------------------------------------------------------------
void expandConformers(PhaseSpaceSynthesis *poly_ps, const std::vector<int> &replica_counts,
                      const UserSettings &ui, const SystemCache &sc,
                      const std::vector<StaticExclusionMask> &exclusion_masks,
                      Xoshiro256ppGenerator *xrs, StopWatch *tm) {

  // Extract some user input directives
  const ConformerControls& conf_input = ui.getConformerNamelistInfo();
  const int nbond_rotations      = conf_input.getRotationSampleCount();
  const int max_seeding_attempts = conf_input.getMaxSeedingAttempts();
  const double self_clash_ratio  = conf_input.getSelfClashRatio();
  
  // Count the expanded number of systems
  std::vector<int> replica_bounds = replica_counts;
  replica_bounds.push_back(0);
  prefixSumInPlace(&replica_bounds, PrefixSumType::EXCLUSIVE);
  const int n_proto_sys = replica_counts.size();
  if (n_proto_sys != sc.getSystemCount()) {
    rtErr("Replica counts were provide for " + std::to_string(replica_counts.size()) + ", but the "
          "systems cache contains " + std::to_string(sc.getSystemCount()) + " systems.",
          "expandConformers");
  }
  for (int i = 0; i < n_proto_sys; i++) {
    
    // Obtain the topology pointer for possible use later
    const AtomGraph *iag_ptr = sc.getSystemTopologyPointer(i);
    const ChemicalFeatures& ichemfe = sc.getFeaturesReference(i);
    
    // Create a list of ways to change the conformation
    const std::vector<IsomerPlan> rotatable_groups  = ichemfe.getRotatableBondGroups();
    const std::vector<IsomerPlan> cis_trans_groups  = ichemfe.getCisTransIsomerizationGroups();
    const std::vector<IsomerPlan> invertible_groups = ichemfe.getChiralInversionGroups();
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
    TickCounter<double> ptrack(permutation_states);
    loadScalarStateValues(&ptrack,
                          std::vector<double2>(permutation_states.size(), { 0.0, stormm::twopi }));
    
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
    if (conformations_per_case != replica_bounds[i + 1] - replica_bounds[i]) {
      rtWarn("The number of conformations for system " + std::to_string(i) + ", based on " +
             getBaseName(iag_ptr->getFileName()) + ", was computed to be " +
             std::to_string(conformations_per_case) + " whereas " +
             std::to_string(replica_bounds[i + 1] - replica_bounds[i]) + " were expected.  The "
             "number of conformations will be adjusted to meet the number planned.",
             "expandConformers");
      conformations_per_case = replica_bounds[i + 1] - replica_bounds[i];
    }

    // Clone the original coordinates
    CoordinateFrame proto_cfi(sc.getCoordinateReference(i));

    // Get the vectors of chiral centers and inversion plans for all centers in this molecule
    const std::vector<int> chiral_centers = ichemfe.getChiralCenters();
    const std::vector<ChiralInversionProtocol> chiral_center_plans =
      ichemfe.getChiralInversionMethods();
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
    switch (strat) {
    case SamplingStrategy::FULL:
      for (int k = 0; k < conformations_per_case; k++) {
        const CoordinateFrame seed_cfi = forgeConformation(proto_cfi, iag_ptr, exclusion_masks[i],
                                                           isomerizers, &ptrack, xrs,
                                                           chiral_centers, chiral_center_plans,
                                                           chiral_variable_indices,
                                                           invertible_groups, max_seeding_attempts,
                                                           self_clash_ratio);
        poly_ps->import(seed_cfi, replica_bounds[i] + k);
        ptrack.advance();
      }
      break;
    case SamplingStrategy::LIMITED:
      {
        // In limited sampling, coupled isomerizing variables will be explored in concert.  In
        // sparse sampling, each variable will be explored in at least one case, provided that
        // there are enough samples of the various permutations to fill.
        const NonbondedKit<double> nbk = iag_ptr->getDoublePrecisionNonbondedKit();
        int fc = 0;
        while (fc < conformations_per_case) {
          int k = 0;
          while (k < n_variables && fc < conformations_per_case) {
            int m = k + 1;
            while (m < n_variables && fc < conformations_per_case) {
              if (permutationsAreLinked(isomerizers, k, m, nbk)) {
                int pos_k = 0;
                while (pos_k < permutation_limits[k] && fc < conformations_per_case) {
                  int pos_m = 0;
                  while (pos_m < permutation_limits[m] && fc < conformations_per_case) {
                    ptrack.randomize(xrs);
                    ptrack.set(pos_k, k);
                    ptrack.set(pos_m, m);
                    const CoordinateFrame seed_cfi =
                      forgeConformation(proto_cfi, iag_ptr, exclusion_masks[i], isomerizers,
                                        &ptrack, xrs, chiral_centers, chiral_center_plans,
                                        chiral_variable_indices, invertible_groups,
                                        max_seeding_attempts, self_clash_ratio);
                    poly_ps->import(seed_cfi, replica_bounds[i] + fc);
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
        int fc = 0;
        while (fc < conformations_per_case) {
          int k = 0;
          while (k < n_variables && fc < conformations_per_case) {
            int m = 0;
            while (m < permutation_limits[k] && fc < conformations_per_case) {
              ptrack.randomize(xrs);
              ptrack.set(m, k);
              const CoordinateFrame seed_cfi =
                forgeConformation(proto_cfi, iag_ptr, exclusion_masks[i], isomerizers, &ptrack,
                                  xrs, chiral_centers, chiral_center_plans,
                                  chiral_variable_indices, invertible_groups,
                                  max_seeding_attempts, self_clash_ratio);
              poly_ps->import(seed_cfi, replica_bounds[i] + fc);
              fc++;
              m++;
            }
            k++;
          }
        }
      }
      break;
    }
  }
  tm->assignTime(2);
}

} // namespace setup
} // namespace conf_app
