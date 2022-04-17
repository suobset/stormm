#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Trajectory/coordinate_series.h"
#include "user_settings.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using user_input::ConformerControls;
using omni::chemistry::ChemicalFeatures;
using omni::chemistry::ChiralInversionProtocol;
using omni::chemistry::MapRotatableGroups;
using omni::topology::AtomGraph;
using omni::trajectory::PhaseSpace;
using omni::trajectory::CoordinateSeries;
  
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     StopWatch *tm) {
  const ConformerControls conf_input = ui.getConformerNamelistInfo();

  // Count the expanded number of systems
  const int ntop = sc.getTopologyCount();
  std::vector<ChemicalFeatures> chemfe_list;
  chemfe_list.reserve(ntop);
  for (int i = 0; i < ntop; i++) {
    const int example_system_idx = sc.getCoordinateExample(i);
    chemfe_list.emplace_back(sc.getTopologyPointer(example_system_idx),
                             sc.getCoordinateReference(example_system_idx),
                             MapRotatableGroups::YES);
  }

  // Loop over all systems, grouping those with the same topology into a coherent group of
  // proto-conformers for coarse-grained sampling of rotatable bonds and chiral centers.
  const int nsys = sc.getSystemCount();
  int nconformer = 0;
  for (int i = 0; i < ntop; i++) {
    const int nrot_bond = std::min(chemfe_list[i].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int ncases = sc.getTopologyCaseCount(i);
    const int bond_rotation_reps = nrot_bond * conf_input.getRotationSampleCount();
    int nproto_conf = ncases * bond_rotation_reps;
    int nchiral_samples = 1;
    if (conf_input.sampleChirality()) {
      const std::vector<ChiralInversionProtocol> chiral_protocols;
      const int nchiral = chemfe_list[i].getChiralCenterCount();
      for (int i = 0; i < nchiral; i++) {
        switch (chiral_protocols[i]) {
        case ChiralInversionProtocol::ROTATE:
        case ChiralInversionProtocol::REFLECT:
          nproto_conf *= 2;
          nchiral_samples *= 2;
          break;
        case ChiralInversionProtocol::DO_NOT_INVERT:
          break;
        }
        if (nproto_conf > 32768) {
          break;
        }
      }
    }
    CoordinateSeries<float> cseries(sc.getCoordinateReference(i), nproto_conf);
    const std::vector<int> top_cases = sc.getTopologicalCases(i);
    int pos = 0;
    for (int j = 0; j < ncases; j++) {
      for (int k = 0; k < bond_rotation_reps * nchiral_samples; k++) {
        cseries.importCoordinateSet(sc.getCoordinateReference(top_cases[j]), pos);
        pos++;
      }
    }

    // Manipulate this coordinate series using bond rotations as well as chiral inversions.
    // Compute the RMSD matrix and determine a set of diverse conformations.  Cull the results
    // to eliminate ring stabs or severe clashes between tetiary or quaternary atoms.
    pos = 0;
    for (int j = 0; j < ncases; j++) {
      for (int k = 0; k < bond_rotation_reps * nchiral_samples; k++) {
        
        pos++;
      }
    }
  }
  tm->assignTime(2);
  
  // Create lists of PhaseSpace objects and topology pointers to show how to model each of them
  std::vector<PhaseSpace> ps_list;
  ps_list.reserve(nconformer);
  std::vector<AtomGraph*> ag_list(nconformer);
  int conf_counter = 0;
  for (int i = 0; i < nsys; i++) {
    const int top_idx = sc.getTopologyIndex(i);
    const int nrot_bond = std::min(chemfe_list[top_idx].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    const int ncopy = pow(nrot_bond, conf_input.getRotationSampleCount());
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
