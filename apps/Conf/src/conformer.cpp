#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Random/random.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "setup.h"

using stormm::random::Xoshiro256ppGenerator;
using stormm::testing::StopWatch;
using conf_app::setup::expandConformers;

using namespace stormm::chemistry;
using namespace stormm::diskutil;
using namespace stormm::namelist;
using namespace stormm::synthesis;

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch master_timer("Master timings for conformer.omni");
  master_timer.addCategory("Input parsing");
  master_timer.addCategory("Chemical feature detection");
  master_timer.addCategory("Coordinate expansion");
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::CONFORMER);

  // Boot up the master random number generator
  const RandomControls rngcon = ui.getRandomNamelistInfo();
  Xoshiro256ppGenerator xrs(rngcon.getRandomSeed(), rngcon.getWarmupCycleCount());
  
  // Read topologies and coordinate files
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::YES,
                 &master_timer);
  master_timer.assignTime(1);

  // Parse the rotatable bonds, cis-trans isomers, and chiral centers in each system to prepare
  // a much larger list of all the possible conformers that each system might be able to access.
  PhaseSpaceSynthesis conformer_population = expandConformers(ui, sc, &xrs, &master_timer);
  master_timer.assignTime(3);

  // CHECK
  std::vector<bool> is_printed(conformer_population.getSystemCount(), false);
  const int nconf = conformer_population.getSystemCount();
  for (int i = 0; i < nconf; i++) {
    if (is_printed[i]) {
      continue;
    }
    const AtomGraph* iag_ptr = conformer_population.getSystemTopologyPointer(i);
    std::vector<int> like_confs(1, i);
    for (int j = i + 1; j < nconf; j++) {
      const AtomGraph* jag_ptr = conformer_population.getSystemTopologyPointer(j);
      if (iag_ptr == jag_ptr) {
        like_confs.push_back(j);
        is_printed[j] = true;
      }
    }
    is_printed[i] = true;
    const std::string fname = substituteNameExtension("bloom_" +
                                                      getBaseName(iag_ptr->getFileName()), "crd");
    conformer_population.printTrajectory(like_confs, fname, 0.0, CoordinateFileKind::AMBER_CRD,
                                         PrintSituation::OVERWRITE);
  }
  // END CHECK
  
  // Collate the topologies in preparation to operate on the new coordinate population
  
  master_timer.printResults();

  return 0;
}
