#include <cmath>
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Structure/rmsd.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
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
using stormm::namelist::ConformerControls;
using stormm::structure::rmsd;
using stormm::structure::RmsdMethod;
using stormm::structure::rotateAboutBond;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::CoordinateSeries;
using stormm::trajectory::CoordinateSeriesWriter;

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     StopWatch *tm) {
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
    const double ln_choices = ptrack.getLogPermutationCount();
    SamplingStrategy strat;
    int conformations_per_case;
    if (ln_choices < log(conf_input.getSystemTrialCount())) {
      conformations_per_case = ptrack.getExactPermutationCount();
      strat = SamplingStrategy::FULL;
    }
    else if (ln_choices < log(static_cast<double>(conf_input.getSystemTrialCount()) * 32.0)) {
      conformations_per_case = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::LIMITED;
    }
    else {
      conformations_per_case = conf_input.getSystemTrialCount();
      strat = SamplingStrategy::SPARSE;
    }

    // This reference, while const, will allow the program to see into the permutation tracking
    // TickCounter object as the state updates with each iteration.
    const std::vector<int>& permutation_watch  = ptrack.getSettings();
    const std::vector<int> permutation_limits = ptrack.getStateLimits();

    // Create a series to hold the conformers resulting from the coarse-grained search.  Populate
    // the series with copies of each case matching this topology, in sufficient quantities to
    // transform the original conformation into every permutation of its isomeric features.
    CoordinateSeries<float> cseries(sc.getTopologyPointer(i)->getAtomCount(),
                                    conformations_per_case);
    const std::vector<int> top_cases = sc.getTopologicalCases(i);
    cseries.resize(0);
    int fc = 0;
    for (int j = 0; j < ncases; j++) {
      fc += conformations_per_case;
      cseries.resize(fc, sc.getCoordinateReference(top_cases[j]));
    }

    // Get the vectors of chiral centers and inversion plans for all centers in this molecule
    const std::vector<int> chiral_centers = chemfe_list[i].getChiralCenters();
    const std::vector<ChiralInversionPlan> chiral_center_plans =
      chemfe_list[i].getChiralInversionMethods();
    
    // Manipulate this coordinate series using bond rotations as well as chiral inversions.
    // Compute the RMSD matrix and determine a set of diverse conformations.  Cull the results
    // to eliminate ring stabs or severe clashes between tertiary or quaternary atoms.
    fc = 0;
    CoordinateSeriesWriter<float> cseries_w = cseries.data();
    for (int j = 0; j < ncases; j++) {
      ptrack.reset();
      switch (strat) {
      case SamplingStrategy::FULL:
        for (int k = 0; k < conformations_per_case; k++) {
          for (int m = 0; m < n_variables; m++) {
            switch (isomerizers[m].getMotion()) {
            case BOND_ROTATION:
            case CIS_TRANS_FLIP:
              {
                const double rval = static_cast<double>(permutation_watch[m]) /
                                    static_cast<double>(permutation_limits[m]) *
                                    stormm::symbols::twopi;
                rotateAboutBond(cseries_w, fc, isomerizers[m].getRootAtom(),
                                isomerizers[m].getPivotAtom(), isomerizers[m].getMovingAtoms(),
                                rval);
              }
              break;s
            case CHIRAL_INVERSION:
              {
                switch (isomerizers[m].getChiralPlan()) {
                case ChiralInversionProtocol::ROTATE:
                case ChiralInversionProtocol::REFLECT:
                  flipChiralCenter(cseries_w, fc, isomerizers[m].getRootAtom(),
                                   chiral_centers, chiral_center_plans, invertible_groups);
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
        for (int k = 0; k < conformations_per_case; k++) {
          
        }
        break;
      }
    }
    ps_list.push_back(
  }
  
  for (int i = 0; i < nsys; i++) {
    const int top_idx = sc.getSystemTopologyIndex(i);
    const int nrot_bond = std::min(chemfe_list[top_idx].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int ncopy = nrot_bond * conf_input.getRotationSampleCount();
    for (int j = 0; j < ncopy; j++) {
      ag_list[conf_counter] = const_cast<AtomGraph*>(sc.getTopologyPointer(i));
      ps_list.push_back(PhaseSpace(sc.getCoordinateReference(i)));
      conf_counter++;
    }
  }
  return PhaseSpaceSynthesis(ps_list, ag_list);
}

} // namespace setup
} // namespace conf_app
