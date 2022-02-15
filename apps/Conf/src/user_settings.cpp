#include <cstdlib>
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/polynumeric.h"
#include "../../../src/Parsing/textfile.h"
#include "../../../src/Reporting/error_format.h"
#include "user_settings.h"

namespace conf_app {
namespace user_input {

using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtErr;
using omni::errors::rtWarn;
using omni::namelist::translateImplicitSolventModel;
using omni::namelist::NamelistEmulator;
using omni::namelist::minimizeInput;
using omni::namelist::randomInput;
using omni::namelist::solventInput;
using omni::parse::NumberFormat;
using omni::parse::TextFile;
using omni::parse::TextOrigin;
using omni::parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
    topology_file_name{}, coordinate_file_name{}, frame_start{0}, frame_end{0}, replica_count{0},
    coordinate_kind{CoordinateFileKind::UNKNOWN}
{}

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem(const std::string &topology_file_in,
                               const std::string &coordinate_file_in, const int frame_start_in,
                               const int frame_end_in, const int replica_count_in,
                               const CoordinateFileKind coordinate_kind_in) :
    topology_file_name{topology_file_in},
    coordinate_file_name{coordinate_file_in},
    frame_start{frame_start_in},
    frame_end{(replica_count_in > 1) ? frame_start : frame_end_in},
    replica_count{replica_count_in},
    coordinate_kind{coordinate_kind_in}
{}

//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const int argc, const char* argv[]) :
    input_file{std::string("cgen.in")},
    topology_file_names{},
    coordinate_file_names{},
    systems{},
    common_core_mask{std::string("")},
    report_file{std::string("cgen.out")},
    conf_file_base{std::string("conf_")},
    conf_file_ext{std::string(".crd")},
    all_free_trajectory_frames{false},
    system_count{0},
    igseed{omni::namelist::default_random_seed},
    random_stream_count{omni::namelist::default_random_streams},
    minimization_steps{omni::namelist::default_minimize_maxcyc},
    steepest_descent_steps{omni::namelist::default_minimize_ncyc},
    initial_move_length{omni::namelist::default_minimize_dx0},
    convergence_target{omni::namelist::default_minimize_drms},
    igb{translateImplicitSolventModel(omni::namelist::default_solvent_igb)},
    born_radii_cutoff{omni::namelist::default_solvent_rgbmax},
    dielectric{omni::namelist::default_solvent_extdiel},
    salt_concentration{omni::namelist::default_solvent_saltcon},
    topology_cache{},
    initial_coordinates_cache{},
    topology_indices{}
{
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile  = false;
  bool cli_igseed   = false;
  int n_cli_topfile = 0;
  int n_cli_crdfile = 0;
  bool cli_report   = false;
  bool cli_confbase = false;
  bool cli_confext  = false;
  for (int i = 1; i < argc; i++) {
    if (i < argc - 1 && strcmp(argv[i], "-i") == 0) {
      input_file = std::string(argv[i + 1]);
      cli_inpfile = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-igseed") == 0) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER) == false) {
        rtErr("The random seed " + std::string(argv[i + 1]) + " is not valid.", "UserSettings");
      }
      igseed = atoi(argv[i + 1]);
      cli_igseed = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-p") == 0) {
      std::string tmp_top(argv[i + 1]);
      topology_file_names.push_back(tmp_top);
      n_cli_topfile++;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-c") == 0) {
      std::string tmp_crd(argv[i + 1]);
      coordinate_file_names.push_back(tmp_crd);
      n_cli_crdfile++;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-o") == 0) {
      report_file = std::string(argv[i + 1]);
      cli_report = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-xbase") == 0) {
      conf_file_base = std::string(argv[i + 1]);
      cli_confbase = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-xext") == 0) {
      conf_file_ext = std::string(argv[i + 1]);
      cli_confext = true;
      i++;
    }
    else {
      rtErr("Command line argument " + std::string(argv[i]) + " was not recognized.",
            "UserSettings");
    }
  }

  // Process the input file
  if (getDrivePathType(input_file) != DrivePathType::FILE) {
    const std::string descriptor = (cli_inpfile) ? std::string("user specified") :
                                                   std::string("default");
    rtErr("The " + descriptor + " input file " + input_file + " was not found or could not be "
          "read.", "UserSettings");
  }
  TextFile inp_tf(input_file, TextOrigin::DISK, "Input deck for conformer.omni", "UserSettings");
  int start_line = 0;
  const NamelistEmulator solvent_info  = solventInput(inp_tf, &start_line);
  start_line = 0;
  const NamelistEmulator minimize_info = minimizeInput(inp_tf, &start_line);
  start_line = 0;
  const NamelistEmulator random_info   = randomInput(inp_tf, &start_line);

}

//-------------------------------------------------------------------------------------------------
int UserSettings::getSystemCount() const {
  return system_count;
}

} // namespace user_input
} // namespace conf_app
