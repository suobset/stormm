#include <cstdlib>
#include <string>
#include <vector>
#include "FileManagement/file_listing.h"
#include "Namelists/namelist_emulator.h"
#include "Namelists/namelist_element.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "user_settings.h"

namespace omni {
namespace namelist {

using constants::CaseSensitivity;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using errors::rtErr;
using errors::rtWarn;
using namelist::NamelistElement;
using namelist::NamelistType;
using parse::NumberFormat;
using parse::TextOrigin;
using parse::verifyNumberFormat;
using parse::WrapTextSearch;

//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const int argc, const char* argv[], const AppName prog_set) :
    policy{ExceptionResponse::DIE},
    input_file{std::string(default_conformer_input_file)},
    file_io_input{}, line_min_input{}, solvent_input{}, prng_input{}
{
  // Local variables to store command line arguments
  int cval_igseed = 0;
  std::string cval_report_file, cval_traj_file_name;
  std::vector<std::string> cval_topology_file_names;
  std::vector<std::string> cval_coordinate_file_names;
  
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile  = false;
  bool cli_igseed   = false;
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
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-c") == 0) {
      std::string tmp_crd(argv[i + 1]);
      cval_coordinate_file_names.push_back(tmp_crd);
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-o") == 0) {
      cval_report_file = std::string(argv[i + 1]);
      cli_report = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-xname") == 0) {
      cval_traj_file_name = std::string(argv[i + 1]);
      cli_confname = true;
      i++;
    }
    else if (strcmp(argv[i], "-warn") == 0) {
      policy = ExceptionResponse::WARN;
    }
    else if (strcmp(argv[i], "-silent") == 0) {
      policy = ExceptionResponse::SILENT;
    }
    else {
      rtErr("Command line argument " + std::string(argv[i]) + " was not recognized.",
            "UserSettings");
    }
  }

  // Process the input file.  Take only the first instance of each namelist, as found by searching
  // from the beginning.
  if (getDrivePathType(input_file) != DrivePathType::FILE) {
    const std::string descriptor = (cli_inpfile) ? std::string("user specified") :
                                                   std::string("default");
    rtErr("The " + descriptor + " input file " + input_file + " was not found or could not be "
          "read.", "UserSettings");
  }
  TextFile inp_tf(input_file, TextOrigin::DISK, "Input deck for OMNI executable", "UserSettings");
  int start_line = 0;
  file_io_input = FilesControls(inp_tf, &start_line, policy);
  start_line = 0;
  line_min_input = MinimizeControls(inp_tf, &start_line, policy);
  start_line = 0;
  solvent_input = SolventControls(inp_tf, &start_line, policy);
  start_line = 0;
  prng_input = RandomControls(inp_tf, &start_line, policy);
  switch (prog_set) {
  case AppName::CONFORMER:
    start_line = 0;
    conf_input = ConformerControls(inp_tf, &start_line, policy);
    break;
  case AppName::DYNAMICS:
    break;
  }

  // Superimpose, or contribute, command line directives
  if (cli_igseed) {
    prng_input.setRandomSeed(cval_igseed);    
  }
  if (cli_report) {
    file_io_input.setReportFileName(cval_report_file);
  }
  if (cli_confname) {
    file_io_input.setGeneralTrajectoryFileName(cval_traj_file_name);
  }
  if (cval_topology_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_topology_file_names.size(); i++) {
      file_io_input.addFreeTopologyName(cval_topology_file_names[i]);
    }
  }
  if (cval_coordinate_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_coordinate_file_names.size(); i++) {
      file_io_input.addFreeCoordinateName(cval_coordinate_file_names[i]);
    }
  }
}

//-------------------------------------------------------------------------------------------------
ExceptionResponse UserSettings::getExceptionBehavior() const {
  return policy;
}
  
//-------------------------------------------------------------------------------------------------
std::string UserSettings::getInputFileName() const {
  return input_file;
}

//-------------------------------------------------------------------------------------------------
FilesControls UserSettings::getFilesNamelistInfo() const {
  return file_io_input;
}

//-------------------------------------------------------------------------------------------------
MinimizeControls UserSettings::getMinimizeNamelistInfo() const {
  return line_min_input;
}

//-------------------------------------------------------------------------------------------------
SolventControls UserSettings::getSolventNamelistInfo() const {
  return solvent_input;
}

//-------------------------------------------------------------------------------------------------
RandomControls UserSettings::getRandomNamelistInfo() const {
  return prng_input;
}

//-------------------------------------------------------------------------------------------------
ConformerControls UserSettings::getConformerNamelistInfo() const {
  return conf_input;
}

} // namespace namelist
} // namespace omni
