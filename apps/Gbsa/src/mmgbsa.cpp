#include "../../../src/copyright.h"
#include "../../../src/Accelerator/hpc_config.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Namelists/command_line_parser.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/stopwatch.h"

using stormm::testing::StopWatch;
using namespace stormm::constants;
using namespace stormm::display;
using namespace stormm::errors;
using namespace stormm::namelist;
using namespace stormm::review;
using namespace stormm::synthesis;

// Define some parameters for MM-GBSA calculations
constexpr int default_sphere_points = 1000;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch master_timer("Master timings for mmgbsa.stormm");
  master_timer.addCategory("Input parsing");
  master_timer.addCategory("Work unit building");

  // Parse the command line
  CommandLineParser clip("mmgbsa.stormm", "A program for calculating MMGBSA energies of ligands "
                         "based on a .");
  clip.addStandardApplicationInputs({ "-i", "-O", "-o", "-except" });
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-nsph", NamelistType::INTEGER, std::to_string(default_sphere_points));
  t_nml->addHelp("-nsph", "Set the number of sphere points that will be used to estimate surface "
                 "area");
  t_nml->addKeyword("-rec_p", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-rec_p", "File path to the receptor topology file.  If this structure is not "
                 "provided on the command line, the information will be sought from a -sys "
                 "keyword will be sought from the &files namelist of the input file with the "
                 "label \"receptor\".");
  t_nml->addKeyword("-rec_c", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-rec_c", "File path to the receptor input coordinates file.  If this structure "
                 "is not provided on the command line, the information will be sought from a -sys "
                 "keyword will be sought from the &files namelist of the input file with the "
                 "label \"receptor\".");
  t_nml->addKeyword("-lig_p", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-lig_p", "File path to the ligand topology file");
  t_nml->addKeyword("-lig_c", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-lig_c", "File path to the ligand input coordinates file.  This may contain "
                 "many frames with different poses of the ligand, but must exist in the same "
                 "'coordinate frame as the receptor.");
  const std::vector<std::string> my_namelists = { "&files", "&minimize", "&restraint", "&solvent",
                                                  "&random", "&report", "&precision" };
  clip.addControlBlocks(my_namelists);
  if (displayNamelistHelp(argc, argv, my_namelists) && clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv); 

  // Take in the user input from the input file.
  const UserSettings ui(clip, { "-pe", "-ce", "-rg" });
  const ExceptionResponse policy = translateExceptionResponse(t_nml->getStringValue("-except"));
  const PrintSituation prnt_protocol = t_nml->getBoolValue("-O") ? PrintSituation::OVERWRITE :
                                                                   PrintSituation::OPEN_NEW;

  // Search for a receptor among the various structures.  This will make a copy of the &files
  // namelist in the UserSettings class object, but the copy will be modifiable in the context of
  // this program.  Check for input errors.
  FilesControls ficon = ui.getFilesNamelistInfo();
  const bool cli_top = (t_nml->getKeywordStatus("-rec_p") == InputStatus::USER_SPECIFIED);
  const bool cli_crd = (t_nml->getKeywordStatus("-rec_c") == InputStatus::USER_SPECIFIED);
  if (cli_top && cli_crd) {
    MoleculeSystem recmol(t_nml->getStringValue("-rec_p"), t_nml->getStringValue("-rec_c"), "", "",
                          "receptor", 0, 1, 1, CoordinateFileKind::UNKNOWN,
                          CoordinateFileKind::UNKNOWN, CoordinateFileKind::UNKNOWN);
    ficon.addSystem(recmol);
  }
  else if ((cli_top && cli_crd == false) || (cli_top == false && cli_crd)) {
    rtErr("A receptor " +
          ((cli_top) ? std::string("topology") : std::string("input coordinates")) + " file was "
          "provided but the corresponding " +
          ((cli_top) ? std::string("input_coordinates") : std::string("topology")) + " was not.",
          "mmgbsa");
  }
  
  // Load all of the systems
  SystemCache sysche(ficon, policy, MapRotatableGroups::NO, prnt_protocol, &master_timer);

  // CHECK
  printf("There are %d systems in the cache with %d label groups.\n", sysche.getSystemCount(),
         sysche.getLabelCount());
  for (int i = 0; i < sysche.getSystemCount(); i++) {
    const AtomGraph& ag_i = sysche.getSystemTopology(i);
    printf("  %s\n", ag_i.getFileName().c_str());
  }
  // END CHECK
  
  // Find the receptor among all of the different structures in the &files namelist (including
  // command line edits).
  const std::vector<int> rec_cache_idx = sysche.getMatchingSystemIndices("receptor");
  const std::vector<int> lig_cache_idx = sysche.getMatchingSystemIndices("ligand");

  // CHECK
  printf("The receptor is system index %d in the cache (%zu systems under this label).\n",
         rec_cache_idx[0], rec_cache_idx.size());
  printf("The ligands are system indices [ ");
  for (size_t i = 0; i < lig_cache_idx.size(); i++) {
    printf("%d ", lig_cache_idx[i]);
  }
  printf("] in the cache (%zu systems under this label).\n", lig_cache_idx.size());
  // END CHECK

  // Combine the topologies of each ligand with that of the receptor.
  AtomGraph ag_complex(sysche.getSystemTopology(rec_cache_idx[0]),
                       sysche.getSystemTopology(lig_cache_idx[0]));
  
  return 0;
}
