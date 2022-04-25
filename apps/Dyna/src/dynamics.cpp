#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"

using omni::chemistry::MapRotatableGroups;
using omni::namelist::AppName;
using omni::namelist::UserSettings;
using omni::testing::StopWatch;
using omni::synthesis::PhaseSpaceSynthesis;
using omni::synthesis::SystemCache;

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch timer("Timings for dynamics.omni");
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::DYNAMICS);
  
  // Read topologies and coordinate files
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::NO,
                 &timer);

  // Parse the rotatable bonds, cis-trans isomers, and chiral centers in each system to prepare
  // a much larger list of all the possible conformers that each system might be able to access.

  // Collate the topologies in preparation to operate on the new coordinate population
  timer.printResults();
  
  return 0;
}
