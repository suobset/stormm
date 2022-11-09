#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Random/random.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "setup.h"

using stormm::random::Xoshiro256ppGenerator;
using stormm::testing::StopWatch;
using conf_app::setup::expandConformers;

using namespace stormm::chemistry;
using namespace stormm::diskutil;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::synthesis;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// Display a general help message for this program.
//-------------------------------------------------------------------------------------------------
void displayGeneralHelpMessage() {
  const std::string base_msg =
    terminalFormat("A program for probing conformations of chemical structures and returning "
                   "those with the lowest energy.\n\n", "conformer"); 
  const std::string nml_msg =
    terminalFormat("Applicable namelists (re-run with one of these terms as the command-line "
                   "argument, IN QUOTES, i.e. \"&files\" or '&files', for further details):\n"
                   "  - &files\n  - &restraint\n  - &conformer\n  - &minimize\n  - &report\n",
                   nullptr, nullptr, 0, 2, 2);
  printf("%s", base_msg.c_str());
  printf("%s", nml_msg.c_str());
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Check for a help message
  if (detectHelpSignal(argc, argv)) {
    displayGeneralHelpMessage();
    return 0;
  }
  if (displayNamelistHelp(argc, argv, { "&files", "&conformer", "&restraint", "&solvent",
                                        "&minimize", "&report" })) {
    return 0;
  }

  // Wall time tracking
  StopWatch master_timer("Master timings for conformer.omni");
  master_timer.addCategory("Input parsing");
  master_timer.addCategory("Chemical feature detection");
  master_timer.addCategory("Coordinate expansion");
  master_timer.addCategory("Energy minimization");

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::CONFORMER);

  // Boot up the master random number generator
  const RandomControls rngcon = ui.getRandomNamelistInfo();
  Xoshiro256ppGenerator xrs(rngcon.getRandomSeed(), rngcon.getWarmupCycleCount());
  
  // Read topologies and coordinate files
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::YES,
                 ui.getPrintingPolicy(), &master_timer);
  master_timer.assignTime(1);

  // Parse the rotatable bonds, cis-trans isomers, and chiral centers in each system to prepare
  // a much larger list of all the possible conformers that each system might be able to access.
  PhaseSpaceSynthesis conformer_population = expandConformers(ui, sc, &xrs, &master_timer);
  std::vector<AtomGraph*> unique_topologies = conformer_population.getUniqueTopologies();
  const int n_unique_topologies = unique_topologies.size();
  const ImplicitSolventModel igb = ui.getSolventNamelistInfo().getImplicitSolventModel();
  if (igb != ImplicitSolventModel::NONE) {
    for (int i = 0; i < n_unique_topologies; i++) {
      unique_topologies[i]->setImplicitSolventModel(igb);
    }
  }
  const std::vector<int> conformer_topology_indices =
    conformer_population.getUniqueTopologyIndices();
  AtomGraphSynthesis conf_poly_ag(unique_topologies, conformer_topology_indices);
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

  // Loop over all systems and perform energy minimizations
  std::vector<StaticExclusionMask> conformer_masks;
  conformer_masks.reserve(unique_topologies.size());
  for (int i = 0; i < unique_topologies.size(); i++) {
    conformer_masks.emplace_back(unique_topologies[i]);
  }
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  for (int i = 0; i < nconf; i++) {
    PhaseSpace psi = conformer_population.exportSystem(i);
    AtomGraph *agi = conf_poly_ag.getSystemTopologyPointer(i);
    const ScoreCard emin = minimize(&psi, agi, conformer_masks[conformer_topology_indices[i]],
                                    mincon);
  }
  master_timer.assignTime(4);

  // Print timings results
  master_timer.printResults();

  return 0;
}
