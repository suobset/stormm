#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Topology/atomgraph.h"
#include "user_settings.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using user_input::ConformerControls;
using omni::chemistry::ChemicalFeatures;
using omni::chemistry::MapRotatableGroups;
using omni::topology::AtomGraph;
using omni::trajectory::PhaseSpace;
  
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
  const int nsys = sc.getSystemCount();  
  int nconformer = 0;
  for (int i = 0; i < nsys; i++) {
    const int top_idx = sc.getTopologyIndex(i);
    const int nrot_bond = std::min(chemfe_list[top_idx].getRotatableBondCount(),
                                   conf_input.getRotatableBondLimit());
    nconformer += pow(nrot_bond, conf_input.getRotationSampleCount());
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
