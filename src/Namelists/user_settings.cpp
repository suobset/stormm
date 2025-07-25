#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "Namelists/namelist_emulator.h"
#include "Namelists/namelist_element.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "Trajectory/trajectory_enumerators.h"
#include "user_settings.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using constants::translateExceptionResponse;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using errors::rtErr;
using errors::rtWarn;
using parse::NumberFormat;
using parse::TextOrigin;
using parse::verifyNumberFormat;
using parse::WrapTextSearch;
using trajectory::getEnumerationName;
using trajectory::translateCoordinateFileKind;
  
//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const CommandLineParser &clip,
                           const std::vector<std::string> &sys_reqs) :
    policy{ExceptionResponse::DIE}, print_policy{default_file_writing_directive},
    has_files_nml{false}, has_minimize_nml{false}, has_solvent_nml{false}, has_random_nml{false},
    has_precision_nml{false}, has_conformer_nml{false}, has_receptor_nml{false},
    has_pppm_nml{false}, has_dynamics_nml{false}, has_remd_nml{false}, has_ffmorph_nml{false},
    has_emulator_nml{false}, has_report_nml{false}, restraint_nml_count{0},
    input_file{std::string(default_conformer_input_file)},
    file_io_input{}, line_min_input{}, solvent_input{}, prng_input{}, conf_input{},
    receptor_input{}, pppm_input{}, dyna_input{}, remd_input{}, ffmod_input{}, emul_input{},
    diagnostic_input{}, rstr_inputs{}
{
  // Local variables to store command line arguments
  int cval_igseed = 0;
  std::string cval_report_file, cval_traj_file_name, cval_input_transcript_file;
  std::vector<std::string> cval_topology_file_names;
  std::vector<std::string> cval_coordinate_file_names;
  
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile          = false;
  bool cli_igseed           = false;
  bool cli_report           = false;
  bool cli_input_transcript = false;
  bool cli_trajname         = false;
  CoordinateFileKind c_kind = default_filecon_inpcrd_type;
  CoordinateFileKind x_kind = default_filecon_outcrd_type;
  CoordinateFileKind r_kind = default_filecon_chkcrd_type;
  const NamelistEmulator *t_nml = clip.getNamelistPointer();
  const InputStatus user_spec = InputStatus::USER_SPECIFIED;
  if (t_nml->hasKeyword("-i") && t_nml->getKeywordStatus("-i") == user_spec) {
    input_file = t_nml->getStringValue("-i");
    cli_inpfile = true;
  }
  if (t_nml->hasKeyword("-ig_seed") &&
      t_nml->getKeywordStatus("-ig_seed") == user_spec) {
    cval_igseed = t_nml->getIntValue("-ig_seed");
    cli_igseed = true;
  }
  if (t_nml->hasKeyword("-p") && t_nml->getKeywordStatus("-p") == user_spec) {
    cval_topology_file_names = t_nml->getAllStringValues("-p");
  }
  if (t_nml->hasKeyword("-c") && t_nml->getKeywordStatus("-c") == user_spec) {
    cval_coordinate_file_names = t_nml->getAllStringValues("-c");
  }
  if (t_nml->hasKeyword("-o") && t_nml->getKeywordStatus("-o") == user_spec) {
    cval_report_file = t_nml->getStringValue("-o");
    cli_report = true;
  }
  if (t_nml->hasKeyword("-t") && t_nml->getKeywordStatus("-t") == user_spec) {
    cval_input_transcript_file = t_nml->getStringValue("-t");
    cli_input_transcript = true;
  }
  if (t_nml->hasKeyword("-x") && t_nml->getKeywordStatus("-x") == user_spec) {
    cval_traj_file_name = t_nml->getStringValue("-x");
    cli_trajname = true;
  }
  if (t_nml->hasKeyword("-c_kind") && t_nml->getKeywordStatus("-c_kind") == user_spec) {
    c_kind = translateCoordinateFileKind(t_nml->getStringValue("-c_kind"));
  }
  if (t_nml->hasKeyword("-x_kind") && t_nml->getKeywordStatus("-x_kind") == user_spec) {
    x_kind = translateCoordinateFileKind(t_nml->getStringValue("-x_kind"));
  }
  if (t_nml->hasKeyword("-r_kind") && t_nml->getKeywordStatus("-r_kind") == user_spec) {
    r_kind = translateCoordinateFileKind(t_nml->getStringValue("-r_kind"));
  }
  if (t_nml->hasKeyword("-except") && t_nml->getKeywordStatus("-except") == user_spec) {
    policy = translateExceptionResponse(t_nml->getStringValue("-except"));
  }
  if (t_nml->hasKeyword("-O") && t_nml->getBoolValue("-O")) {
    print_policy = PrintSituation::OVERWRITE;
  }
  
  // Process the input file.  Take only the first instance of each namelist, as found by searching
  // from the beginning.
  if (getDrivePathType(input_file) != DrivePathType::FILE) {
    const std::string descriptor = (cli_inpfile) ? std::string("user specified") :
                                                   std::string("default");
    rtErr("The " + descriptor + " input file " + input_file + " was not found or could not be "
          "read.", "UserSettings");
  }
  TextFile inp_tf(input_file, TextOrigin::DISK, "Input deck for STORMM executable",
                  "UserSettings");

  std::vector<std::string> alternatives = {
    "coordinate_input_format",      getEnumerationName(c_kind),
    "coordinate_output_format",     getEnumerationName(x_kind),
    "coordinate_checkpoint_format", getEnumerationName(r_kind)
  };
  int start_line = 0;
  file_io_input = FilesControls(inp_tf, &start_line, &has_files_nml, policy, WrapTextSearch::NO,
                                alternatives, sys_reqs);
  start_line = 0;
  line_min_input = MinimizeControls(inp_tf, &start_line, &has_minimize_nml, policy);
  start_line = 0;
  solvent_input = SolventControls(inp_tf, &start_line, &has_solvent_nml, policy);
  start_line = 0;
  prng_input = RandomControls(inp_tf, &start_line, &has_random_nml, policy);
  start_line = 0;
  prec_input = PrecisionControls(inp_tf, &start_line, &has_precision_nml, policy);
  start_line = 0;
  conf_input = ConformerControls(inp_tf, &start_line, &has_conformer_nml, policy);
  start_line = 0;
  receptor_input = ReceptorControls(inp_tf, &start_line, &has_receptor_nml, policy);
  start_line = 0;
  pppm_input = PPPMControls(inp_tf, &start_line, &has_pppm_nml, policy);
  start_line = 0;
  dyna_input = DynamicsControls(inp_tf, &start_line, &has_dynamics_nml, policy);
  start_line = 0;
  remd_input = RemdControls(inp_tf, &start_line, &has_remd_nml, policy);
  start_line = 0;
  ffmod_input = FFMorphControls(inp_tf, &start_line, &has_ffmorph_nml, policy);
  start_line = 0;
  emul_input = EmulatorControls(inp_tf, &start_line, &has_emulator_nml, policy);
  start_line = 0;
  diagnostic_input = ReportControls(inp_tf, &start_line, &has_report_nml, policy);
  start_line = 0;
  while (start_line < inp_tf.getLineCount()) {
    bool restraint_nml_found = false;
    RestraintControls tmp_rstr_input = RestraintControls(inp_tf, &start_line, &restraint_nml_found,
                                                         policy);
    if (restraint_nml_found) {
      restraint_nml_count += 1;
      rstr_inputs.push_back(tmp_rstr_input);
    }
  }
  
  // Superimpose, or contribute, command line directives
  if (cli_igseed) {
    prng_input.setRandomSeed(cval_igseed);    
  }
  if (cli_report) {
    file_io_input.setReportFileName(cval_report_file);
  }
  if (cli_trajname) {
    file_io_input.setGeneralTrajectoryFileName(cval_traj_file_name);
  }
  if (cli_input_transcript) {
    file_io_input.setInputTranscriptFileName(cval_input_transcript_file);
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
const std::string& UserSettings::getInputFileName() const {
  return input_file;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getFilesPresence() const {
  return has_files_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getMinimizePresence() const {
  return has_minimize_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getSolventPresence() const {
  return has_solvent_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getRandomPresence() const {
  return has_random_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getPrecisionPresence() const {
  return has_precision_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getConformerPresence() const {
  return has_conformer_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getPPPMPresence() const {
  return has_pppm_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getDynamicsPresence() const {
  return has_dynamics_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getRemdPresence() const {
  return has_remd_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getFFMorphPresence() const {
  return has_ffmorph_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getEmulatorPresence() const {
  return has_emulator_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getReportPresence() const {
  return has_report_nml;
}

//-------------------------------------------------------------------------------------------------
const FilesControls& UserSettings::getFilesNamelistInfo() const {
  return file_io_input;
}

//-------------------------------------------------------------------------------------------------
const MinimizeControls& UserSettings::getMinimizeNamelistInfo() const {
  return line_min_input;
}

//-------------------------------------------------------------------------------------------------
const SolventControls& UserSettings::getSolventNamelistInfo() const {
  return solvent_input;
}

//-------------------------------------------------------------------------------------------------
const RandomControls& UserSettings::getRandomNamelistInfo() const {
  return prng_input;
}

//-------------------------------------------------------------------------------------------------
const PrecisionControls& UserSettings::getPrecisionNamelistInfo() const {
  return prec_input;
}

//-------------------------------------------------------------------------------------------------
const ConformerControls& UserSettings::getConformerNamelistInfo() const {
  return conf_input;
}

//-------------------------------------------------------------------------------------------------
const ReceptorControls& UserSettings::getReceptorNamelistInfo() const {
  return receptor_input;
}

//-------------------------------------------------------------------------------------------------
const PPPMControls& UserSettings::getPPPMNamelistInfo() const {
  return pppm_input;
}

//-------------------------------------------------------------------------------------------------
const DynamicsControls& UserSettings::getDynamicsNamelistInfo() const {
  return dyna_input;
}

//-------------------------------------------------------------------------------------------------
const RemdControls& UserSettings::getRemdNamelistInfo() const {
  return remd_input;
}

//-------------------------------------------------------------------------------------------------
const FFMorphControls& UserSettings::getFFMorphNamelistInfo() const {
  return ffmod_input;
}

//-------------------------------------------------------------------------------------------------
const EmulatorControls& UserSettings::getEmulatorNamelistInfo() const {
  return emul_input;
}

//-------------------------------------------------------------------------------------------------
const ReportControls& UserSettings::getReportNamelistInfo() const {
  return diagnostic_input;
}

//-------------------------------------------------------------------------------------------------
const std::vector<RestraintControls>& UserSettings::getRestraintNamelistInfo() const {
  return rstr_inputs;
}

//-------------------------------------------------------------------------------------------------
const RestraintControls& UserSettings::getRestraintNamelistInfo(const int index) const {
  if (index < 0 || index >= restraint_nml_count) {
    rtErr("The input contained " + std::to_string(restraint_nml_count) + " &restraint namelists.  "
          "Index " + std::to_string(index) + " is invalid.", "UserSettings",
          "getRestraintNamelistInfo");
  }
  return rstr_inputs[index];
}

//-------------------------------------------------------------------------------------------------
PrintSituation UserSettings::getPrintingPolicy() const {
  return print_policy;
}

} // namespace namelist
} // namespace stormm
