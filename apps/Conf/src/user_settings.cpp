#include <cstdlib>
#include <string>
#include <vector>
#include "../../../src/Constants/behavior.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/polynumeric.h"
#include "../../../src/Parsing/textfile.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "user_settings.h"

namespace conf_app {
namespace user_input {

using omni::constants::ExceptionResponse;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtErr;
using omni::errors::rtWarn;
using omni::parse::NumberFormat;
using omni::parse::TextFile;
using omni::parse::TextOrigin;
using omni::parse::verifyNumberFormat;
using omni::trajectory::detectCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const int argc, const char* argv[]) :
    input_file{std::string(default_conformer_input_file)},
    common_core_mask{std::string("")},
    file_io_input{}, line_min_input{}, solvent_input{}, prng_input{},
    topology_cache{},
    initial_coordinates_cache{},
    topology_indices{}
{
  // Local variables to store command line arguments
  int cval_igseed = 0;
  std::string cval_report_file;
  std::string cval_conf_file_name;
  std::vector<std::string> cval_topology_file_names;
  std::vector<std::string> cval_coordinate_file_names;
  
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile  = false;
  bool cli_igseed   = false;
  int n_cli_topfile = 0;
  int n_cli_crdfile = 0;
  bool cli_report   = false;
  bool cli_confname = false;
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
      cval_igseed = atoi(argv[i + 1]);
      cli_igseed = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-p") == 0) {
      std::string tmp_top(argv[i + 1]);
      cval_topology_file_names.push_back(tmp_top);
      n_cli_topfile++;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-c") == 0) {
      std::string tmp_crd(argv[i + 1]);
      cval_coordinate_file_names.push_back(tmp_crd);
      n_cli_crdfile++;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-o") == 0) {
      cval_report_file = std::string(argv[i + 1]);
      cli_report = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-xname") == 0) {
      cval_conf_file_name = std::string(argv[i + 1]);
      cli_confname = true;
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
  file_io_input = FilesControls(inp_tf, &start_line, ExceptionResponse::WARN);
  start_line = 0;
  line_min_input = MinimizeControls(inp_tf, &start_line, ExceptionResponse::WARN);
  start_line = 0;
  solvent_input = SolventControls(inp_tf, &start_line, ExceptionResponse::WARN);
  start_line = 0;
  prng_input = RandomControls(inp_tf, &start_line, ExceptionResponse::WARN);

  // Read all free topologies and free coordinate sets, then determine which free topology
  // matches which free coordinate set.
  const int n_free_top = file_io_input.getFreeTopologyCount();
  topology_cache.reserve(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    topology_cache.push_back(AtomGraph(file_io_input.getFreeTopologyName(i)));
  }
  const int n_free_crd = file_io_input.getFreeCoordinatesCount();
  for (int i = 0; i < n_free_crd; i++) {
    const std::string crd_name = file_io_input.getFreeCoordinateName(i);
    const CoordinateFileKind kind = detectCoordinateFileKind(crd_name);
  }
  
  // Filter the unique topologies and match them to systems
}

} // namespace user_input
} // namespace conf_app
