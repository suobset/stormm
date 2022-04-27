#include <vector>
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"

using omni::chemistry::MapRotatableGroups;
using omni::energy::ScoreCard;
using omni::mm::minimize;
using omni::namelist::AppName;
using omni::namelist::MinimizeControls;
using omni::namelist::UserSettings;
using omni::restraints::RestraintApparatus;
using omni::synthesis::PhaseSpaceSynthesis;
using omni::synthesis::SystemCache;
using omni::testing::StopWatch;
using omni::topology::UnitCellType;
using omni::trajectory::PhaseSpace;

// CHECK
using omni::restraints::RestraintApparatus;
// END CHECK

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch timer("Timings for dynamics.omni");
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::DYNAMICS);
  
  // Read topologies and coordinate files.  Assemble critical deatils about each system.
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::NO,
                 &timer);

  // Perform minimizations as requested.
  const int system_count = sc.getSystemCount();
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  const int mini_timings = timer.addCategory("Minimization");
  if (ui.getMinimizePresence()) {
    std::vector<ScoreCard> all_mme;
    all_mme.reserve(system_count);
    for (int i = 0; i < system_count; i++) {
      PhaseSpace *ps = sc.getCoordinatePointer(i);

      // CHECK
      RestraintApparatus ra(sc.getSystemTopologyPointer(i));
      // END CHECK
    
      switch(ps->getUnitCellType()) {
      case UnitCellType::NONE:
        all_mme.emplace_back(minimize(ps, sc.getSystemTopologyReference(i), ra,
                                      sc.getSystemStaticMaskReference(i), mincon));
        timer.assignTime(mini_timings);
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
  }

  // Print restart files from energy minimization
  if (mincon.getCheckpointProduction()) {
    for (int i = 0; i < system_count; i++) {
      const PhaseSpace ps = sc.getCoordinateReference(i);
      ps.exportToFile(sc.getSystemCheckpointName(i));
    }
  }

  // Collate the topologies in preparation to operate on the new coordinate population
  timer.assignTime(0);
  timer.printResults();
  
  return 0;
}
