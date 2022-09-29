#include <cmath>
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Math/tickcounter.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Topology/topology_util.h"
#include "../../../src/Trajectory/coordinate_series.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::ChiralInversionProtocol;
using stormm::chemistry::ConformationEdit;
using stormm::chemistry::MapRotatableGroups;
using stormm::chemistry::IsomerPlan;
using stormm::constants::warp_size_int;
using stormm::math::roundUp;
using stormm::math::TickCounter;
using stormm::namelist::ConformerControls;
using stormm::structure::rmsd;
using stormm::structure::RmsdMethod;
using stormm::structure::rotateAboutBond;
using stormm::structure::flipChiralCenter;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::isBonded;
using stormm::topology::NonbondedKit;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::CoordinateSeries;
using stormm::trajectory::CoordinateSeriesWriter;
  
//-------------------------------------------------------------------------------------------------
double computeLocalPermutations(const std::vector<int> &limits,
                                const std::vector<IsomerPlan> &isomerizers, const AtomGraph *ag) {
  double result = 0.0;
  const int n_var = limits.size();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  std::vector<double> log_limits(n_var);
  for (int i = 0; i < n_var; i++) {
    log_limits[i] = log(static_cast<double>(limits[i]));
  }
  for (int i = 0; i < n_var; i++) {
    result += log_limits[i];
    const int root_i = isomerizers[i].getRootAtom();
    const int pivt_i = isomerizers[i].getPivotAtom();
    switch (isomerizers[i].getMotion()) {
    case ConformationEdit::BOND_ROTATION:
    case ConformationEdit::CIS_TRANS_FLIP:
      for (int j = i + 1; j < n_var; j++) {
        const int root_j = isomerizers[j].getRootAtom();
        const int pivt_j = isomerizers[j].getPivotAtom();

        // If two isomerizations share atoms or have their root and pivot atoms bonded to one
        // another, consider them coupled and add the product of the number of possible states for
        // each to the running sum.
        switch (isomerizers[j].getMotion()) {
        case ConformationEdit::BOND_ROTATION:
        case ConformationEdit::CIS_TRANS_FLIP:
          if (root_i == root_j || root_i == pivt_j || pivt_i == root_j || pivt_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j) ||
              isBonded(nbk, pivt_i, root_j) || isBonded(nbk, pivt_i, pivt_j)) {
            result += log_limits[i] + log_limits[j];
          }
          break;
        case ConformationEdit::CHIRAL_INVERSION:
          if (root_i == root_j || pivt_i == root_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j)) {
            result += log_limits[i] + log_limits[j];
          }
          break;
        }
      }
      break;
    case ConformationEdit::CHIRAL_INVERSION:
      for (int j = i + 1; j < n_var; j++) {
        const int root_j = isomerizers[j].getRootAtom();
        const int pivt_j = isomerizers[j].getPivotAtom();
        switch (isomerizers[j].getMotion()) {
        case ConformationEdit::BOND_ROTATION:
        case ConformationEdit::CIS_TRANS_FLIP:
          if (root_i == root_j || root_i == pivt_j ||
              isBonded(nbk, root_i, root_j) || isBonded(nbk, root_i, pivt_j)) {
            result += log_limits[i] + log_limits[j];
          }
          break;
        case ConformationEdit::CHIRAL_INVERSION:
          if (root_i == root_j || isBonded(nbk, root_i, root_j)) {
            result += log_limits[i] + log_limits[j];
          }
          break;
        }
      }
      break;
    }
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     Xoshiro256ppGenerator *xrs, StopWatch *tm) {
  const ConformerControls conf_input = ui.getConformerNamelistInfo();

  // Count the expanded number of systems
  const int ntop = sc.getTopologyCount();
  const int nsys = sc.getSystemCount();
  const std::vector<const AtomGraph*> system_topologies = sc.getSystemTopologyPointer();
  const std::vector<const PhaseSpace*> system_coords    = sc.getCoordinatePointer();
  const std::vector<const AtomGraph*> unique_topologies = sc.getTopologyPointer();
  std::vector<ChemicalFeatures> chemfe_list;
  chemfe_list.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    const int example_system_idx = sc.getCoordinateExample(i);
    chemfe_list.emplace_back(sc.getTopologyPointer(i),
                             sc.getCoordinateReference(example_system_idx),
                             MapRotatableGroups::YES);
  }
  tm->assignTime(2);

  // Create lists of PhaseSpace objects and topology pointers to show how to model each of them
  std::vector<PhaseSpace> ps_list; 
  std::vector<AtomGraph*> ag_list;
  int conf_counter = 0;

  // Loop over all systems, grouping those with the same topology into a coherent group of
  // proto-conformers for coarse-grained sampling of rotatable bonds and chiral centers.
  const int nbond_rotations = conf_input.getRotationSampleCount();
  int nneighborhood = 0;
  for (int i = 0; i < ntop; i++) {

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
    const std::vector<int>& permutation_watch  = ptrack.getSettings();
    const std::vector<int> permutation_limits = ptrack.getStateLimits();
    const double ln_choices = ptrack.getLogPermutationCount();
    SamplingStrategy strat;
    int conformations_per_case;
    if (ln_choices < log(conf_input.getSystemTrialCount())) {
      conformations_per_case = ptrack.getExactPermutationCount();
      strat = SamplingStrategy::FULL;
    }
    else if (computeLocalPermutations(permutation_limits, isomerizers, sc.getTopologyPointer(i)) <
             log(static_cast<double>(conf_input.getSystemTrialCount()))) {
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
    for (int j = 0; j < ncases; j++) {
      ptrack.reset();
      switch (strat) {
      case SamplingStrategy::FULL:
        for (int k = 0; k < conformations_per_case; k++) {
          for (int m = 0; m < n_variables; m++) {
            switch (isomerizers[m].getMotion()) {
            case ConformationEdit::BOND_ROTATION:
            case ConformationEdit::CIS_TRANS_FLIP:
              if (permutation_watch[m] > 0) {
                const double rval = static_cast<double>(permutation_watch[m]) /
                                    static_cast<double>(permutation_limits[m]) *
                                    stormm::symbols::twopi;
                rotateAboutBond<double, double>(cseries_w, fc, isomerizers[m].getRootAtom(),
                                                isomerizers[m].getPivotAtom(),
                                                isomerizers[m].getMovingAtoms(), rval);
              }
              break;
            case ConformationEdit::CHIRAL_INVERSION:
              if (permutation_watch[m] == 1) {
                switch (isomerizers[m].getChiralPlan()) {
                case ChiralInversionProtocol::ROTATE:
                case ChiralInversionProtocol::REFLECT:
                  flipChiralCenter<double, double>(cseries_w, fc, chiral_variable_indices[m],
                                                   chiral_centers, chiral_center_plans,
                                                   invertible_groups);
                  break;
                case ChiralInversionProtocol::DO_NOT_INVERT:
                  break;
                }
              }
              break;
            }
          }
          fc++;
          ptrack.advance();
        }
        break;
      case SamplingStrategy::LIMITED:
      case SamplingStrategy::SPARSE:
        {
          // In limited sampling, coupled isomerizing variables will be explored in concert.  In
          // sparse sampling, each variable will be explored in at least one case, provided that
          // there are enough samples of the various permutations to fill.  In either case, the
          // rest of the available samples will be filled with randomized permutations.
          if (strat == SamplingStrategy::LIMITED) {
            
          }
          else {
            
          }
          for (int k = 0; k < conformations_per_case; k++) {
            for (int m = 0; m < n_variables; m++) {
              switch (isomerizers[m].getMotion()) {
              case ConformationEdit::BOND_ROTATION:
              case ConformationEdit::CIS_TRANS_FLIP:
                if (permutation_watch[m] > 0) {
                  const double rval = static_cast<double>(permutation_watch[m]) /
                                      static_cast<double>(permutation_limits[m]) *
                                      stormm::symbols::twopi;
                  rotateAboutBond<double, double>(cseries_w, fc, isomerizers[m].getRootAtom(),
                                                  isomerizers[m].getPivotAtom(),
                                                  isomerizers[m].getMovingAtoms(), rval);
                }
                break;
              case ConformationEdit::CHIRAL_INVERSION:
                if (permutation_watch[m] == 1) {
                  switch (isomerizers[m].getChiralPlan()) {
                  case ChiralInversionProtocol::ROTATE:
                  case ChiralInversionProtocol::REFLECT:
                    flipChiralCenter<double, double>(cseries_w, fc, chiral_variable_indices[m],
                                                     chiral_centers, chiral_center_plans,
                                                     invertible_groups);
                    break;
                  case ChiralInversionProtocol::DO_NOT_INVERT:
                    break;
                  }
                }
                break;
              }
            }
            fc++;
            ptrack.randomize(xrs);
          }
        }
        break;
      }
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
