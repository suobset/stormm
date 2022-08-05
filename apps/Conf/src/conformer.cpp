#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "setup.h"

using stormm::chemistry::MapRotatableGroups;
using stormm::namelist::AppName;
using stormm::namelist::UserSettings;
using stormm::testing::StopWatch;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SystemCache;
using conf_app::setup::expandConformers;

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch master_timer("Master timings for conformer.omni");
  master_timer.addCategory("Input parsing");
  master_timer.addCategory("Chemical feature detection");
  master_timer.addCategory("Coordinate expansion");
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::CONFORMER);
  
  // Read topologies and coordinate files
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::YES,
                 &master_timer);
  master_timer.assignTime(1);

  // Parse the rotatable bonds, cis-trans isomers, and chiral centers in each system to prepare
  // a much larger list of all the possible conformers that each system might be able to access.
  PhaseSpaceSynthesis conformer_population = expandConformers(ui, sc, &master_timer);
  master_timer.assignTime(3);

  // Collate the topologies in preparation to operate on the new coordinate population
  
  master_timer.printResults();
  
  // CHECK
  printf("Back in main.\n");
  // END CHECK
  
  return 0;
}
