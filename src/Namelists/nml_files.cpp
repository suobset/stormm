#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "Trajectory/trajectory_enumerators.h"
#include "nml_files.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using diskutil::listDirectory;
using diskutil::listFilesInPath;
using diskutil::SearchStyle;
using parse::findStringInVector;
using parse::strcmpCased;
using parse::strncmpCased;
using trajectory::getEnumerationName;
using trajectory::translateCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
    topology_file_name{}, coordinate_file_name{}, coordinate_output_name{}, checkpoint_name{},
    label{}, frame_start{0}, frame_end{0}, replica_count{0},
    coordinate_kind{CoordinateFileKind::UNKNOWN},
    trajectory_kind{CoordinateFileKind::AMBER_CRD},
    checkpoint_kind{CoordinateFileKind::AMBER_ASCII_RST}
{}

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem(const std::string &topology_file_in,
                               const std::string &coordinate_file_in,
                               const std::string &trajectory_file_in,
                               const std::string &checkpoint_file_in, const std::string &label_in,
                               const int frame_start_in, const int frame_end_in,
                               const int replica_count_in,
                               const CoordinateFileKind coordinate_kind_in,
                               const CoordinateFileKind trajectory_kind_in,
                               const CoordinateFileKind checkpoint_kind_in) :
    topology_file_name{topology_file_in},
    coordinate_file_name{coordinate_file_in},
    coordinate_output_name{trajectory_file_in},
    checkpoint_name{checkpoint_file_in},
    label{label_in},
    frame_start{frame_start_in},
    frame_end{(replica_count_in > 1) ? frame_start : frame_end_in},
    replica_count{replica_count_in},
    coordinate_kind{coordinate_kind_in},
    trajectory_kind{trajectory_kind_in},
    checkpoint_kind{checkpoint_kind_in}
{}

//-------------------------------------------------------------------------------------------------
std::string MoleculeSystem::getTopologyFileName() const {
  return topology_file_name;
}

//-------------------------------------------------------------------------------------------------
std::string MoleculeSystem::getInputCoordinateFileName() const {
  return coordinate_file_name;
}

//-------------------------------------------------------------------------------------------------
std::string MoleculeSystem::getTrajectoryFileName() const {
  return coordinate_output_name;
}

//-------------------------------------------------------------------------------------------------
std::string MoleculeSystem::getCheckpointFileName() const {
  return checkpoint_name;
}

//-------------------------------------------------------------------------------------------------
std::string MoleculeSystem::getLabel() const {
  return label;
}

//-------------------------------------------------------------------------------------------------
int MoleculeSystem::getStartingFrame() const {
  return frame_start;
}

//-------------------------------------------------------------------------------------------------
int MoleculeSystem::getFinalFrame() const {
  return frame_end;
}

//-------------------------------------------------------------------------------------------------
int MoleculeSystem::getTotalFrames() const {
  return frame_end - frame_start + 1;
}

//-------------------------------------------------------------------------------------------------
int MoleculeSystem::getReplicaCount() const {
  return replica_count;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind MoleculeSystem::getInputCoordinateFileKind() const {
  return coordinate_kind;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind MoleculeSystem::getTrajectoryFileKind() const {
  return trajectory_kind;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind MoleculeSystem::getCheckpointFileKind() const {
  return checkpoint_kind;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setTopologyFileName(const std::string &file_name) {
  topology_file_name = file_name;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setInputCoordinateFileName(const std::string &file_name) {
  coordinate_file_name = file_name;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setTrajectoryFileName(const std::string &file_name) {
  coordinate_output_name = file_name;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setCheckpointFileName(const std::string &file_name) {
  checkpoint_name = file_name;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setStartingFrame(const int frame_number) {
  frame_start = frame_number;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setFinalFrame(const int frame_number) {
  frame_end = frame_number;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setReplicaCount(const int count) {
  replica_count = count;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setInputCoordinateFileKind(const std::string &kind) {
  coordinate_kind = translateCoordinateFileKind(kind);
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setInputCoordinateFileKind(const CoordinateFileKind kind) {
  coordinate_kind = kind;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setTrajectoryFileKind(const std::string &kind) {
  trajectory_kind = translateCoordinateFileKind(kind);
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setTrajectoryFileKind(const CoordinateFileKind kind) {
  trajectory_kind = kind;
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setCheckpointFileKind(const std::string &kind) {
  checkpoint_kind = translateCoordinateFileKind(kind);
}

//-------------------------------------------------------------------------------------------------
void MoleculeSystem::setCheckpointFileKind(const CoordinateFileKind kind) {
  checkpoint_kind = kind;
}

//-------------------------------------------------------------------------------------------------
bool MoleculeSystem::validateTopologyFile() const {
  return (getDrivePathType(topology_file_name) == DrivePathType::FILE);
}

//-------------------------------------------------------------------------------------------------
bool MoleculeSystem::validateInputCoordinateFile() const {
  return (getDrivePathType(topology_file_name) == DrivePathType::FILE);
}

//-------------------------------------------------------------------------------------------------
FilesControls::FilesControls(const ExceptionResponse policy_in,
                             const WrapTextSearch wrap) :
    policy{policy_in}, structure_count{0}, free_topology_count{0}, free_coordinate_count{0},
    system_count{0},
    all_free_frames{default_filecon_read_all_free},
    coordinate_input_format{default_filecon_inpcrd_type},
    coordinate_output_format{default_filecon_outcrd_type},
    coordinate_checkpoint_format{default_filecon_chkcrd_type},
    topology_file_names{}, coordinate_file_names{}, systems{},
    report_file{std::string(default_filecon_report_name)},
    coordinate_output_name{std::string(default_filecon_trajectory_name)},
    checkpoint_name{std::string(default_filecon_checkpoint_name)}
{}

//-------------------------------------------------------------------------------------------------
FilesControls::FilesControls(const TextFile &tf, int *start_line,
                             const ExceptionResponse policy_in, const WrapTextSearch wrap,
                             const std::vector<std::string> &alternatives,
                             const std::vector<std::string> &sys_requirements) :
    FilesControls(policy_in)
{
  // Set some alternative defaults.  This takes a vector of strings, the even-numbered strings
  // being names of actual member variables and the odd-numbered strings being the new defaults to
  // apply.  Different applications will then be able to call the constructor with different
  // default settings.
  const int n_alt = alternatives.size();
  for (int i = 0; i < n_alt; i++) {
    if (i < n_alt - 1) {
      if (alternatives[i] == std::string("coordinate_input_format")) {
        coordinate_input_format = translateCoordinateFileKind(alternatives[i + 1]);
        i++;
      }
      else if (alternatives[i] == std::string("coordinate_output_format")) {
        coordinate_output_format = translateCoordinateFileKind(alternatives[i + 1]);
        i++;
      }
      else if (alternatives[i] == std::string("coordinate_checkpoint_format")) {
        coordinate_checkpoint_format = translateCoordinateFileKind(alternatives[i + 1]);
        i++;
      }
      else if (alternatives[i] == std::string("report_file")) {
        report_file = alternatives[i + 1];
        i++;
      }
      else if (alternatives[i] == std::string("coordinate_output_name")) {
        coordinate_output_name = alternatives[i + 1];
        i++;
      }
      else if (alternatives[i] == std::string("checkpoint_name")) {
        checkpoint_name = alternatives[i + 1];
        i++;
      }
      else if (alternatives[i] == std::string("all_free_frames")) {
        all_free_frames = (alternatives[i + 1] == std::string("true"));
        i++;
      }
    }
  }
  
  // Get the systems requirements in order
  bool top_name_required = false;
  bool crd_name_required = false;
  bool trj_name_required = false;
  bool rst_name_required = false;
  bool top_must_exist = false;
  bool crd_must_exist = false;
  bool trj_must_exist = false;
  bool rst_must_exist = false;
  bool top_is_bogus = false;
  bool crd_is_bogus = false;
  bool trj_is_bogus = false;
  bool rst_is_bogus = false;
  for (size_t i = 0; i < sys_requirements.size(); i++) {
    top_name_required = (top_name_required ||
                         (sys_requirements[i] == "-p" || sys_requirements[i] == "-pe"));
    crd_name_required = (crd_name_required ||
                         (sys_requirements[i] == "-c" || sys_requirements[i] == "-ce"));
    trj_name_required = (trj_name_required ||
                         (sys_requirements[i] == "-x" || sys_requirements[i] == "-xe"));
    rst_name_required = (rst_name_required ||
                         (sys_requirements[i] == "-r" || sys_requirements[i] == "-re"));
    top_must_exist = (top_must_exist || sys_requirements[i] == "-pe");
    crd_must_exist = (crd_must_exist || sys_requirements[i] == "-ce");
    trj_must_exist = (trj_must_exist || sys_requirements[i] == "-xe");
    rst_must_exist = (rst_must_exist || sys_requirements[i] == "-re");
    top_is_bogus = (top_is_bogus || sys_requirements[i] == "-pg");
    crd_is_bogus = (crd_is_bogus || sys_requirements[i] == "-cg");
    trj_is_bogus = (trj_is_bogus || sys_requirements[i] == "-xg");
    rst_is_bogus = (rst_is_bogus || sys_requirements[i] == "-rg");
  }
  const SubkeyRequirement short_req = SubkeyRequirement::REQUIRED;
  const SubkeyRequirement short_opt = SubkeyRequirement::OPTIONAL;
  const SubkeyRequirement short_bog = SubkeyRequirement::BOGUS;
  const std::vector<SubkeyRequirement> sys_keyword_reqs = {
    (top_name_required) ? short_req : (top_is_bogus) ? short_bog : short_opt,
    (crd_name_required) ? short_req : (crd_is_bogus) ? short_bog : short_opt,
    (trj_name_required) ? short_req : (trj_is_bogus) ? short_bog : short_opt,
    (rst_name_required) ? short_req : (rst_is_bogus) ? short_bog : short_opt,
    short_opt
  };
  NamelistEmulator t_nml = filesInput(tf, start_line, sys_keyword_reqs, policy, wrap,
                                      coordinate_input_format, coordinate_output_format,
                                      coordinate_checkpoint_format);
  const InputStatus stt_missing = InputStatus::MISSING;
  const int nsys = t_nml.getKeywordEntries("-sys");
  for (int i = 0; i < nsys; i++) {
    bool complete = true;
    const std::string top_name = (t_nml.getKeywordStatus("-sys", "-p", i) != stt_missing) ?
                                 t_nml.getStringValue("-sys", "-p", i) : std::string("");
    const std::string crd_name = (t_nml.getKeywordStatus("-sys", "-c", i) != stt_missing) ?
                                 t_nml.getStringValue("-sys", "-c", i) : std::string("");
    const std::string trj_name = (t_nml.getKeywordStatus("-sys", "-x", i) != stt_missing) ?
                                 t_nml.getStringValue("-sys", "-x", i) : std::string("");
    const std::string rst_name = (t_nml.getKeywordStatus("-sys", "-r", i) != stt_missing) ?
                                 t_nml.getStringValue("-sys", "-r", i) : std::string("");
    const std::string sys_label = (t_nml.getKeywordStatus("-sys", "-label", i) != stt_missing) ?
                                  t_nml.getStringValue("-sys", "-label", i) : std::string("");
    std::string missing_elements("");
    if (top_name.size() == 0LLU && top_name_required) {
      missing_elements += "-p";
      complete = false;
    }
    if (crd_name.size() == 0LLU && crd_name_required) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -c" : "-c";
      complete = false;
    }
    if (trj_name.size() == 0LLU && trj_name_required) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -x" : "-x";
      complete = false;
    }
    if (rst_name.size() == 0LLU && rst_name_required) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -r" : "-r";
      complete = false;
    }
    if (complete == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Instance " + std::to_string(i) + " of the \"-sys\" keyword is missing elements: " +
              missing_elements + ".", "FilesControls");
      case ExceptionResponse::WARN:
        rtWarn("Instance " + std::to_string(i) + " of the \"-sys\" keyword is missing elements: " +
               missing_elements + ".  These must be supplied in order for this system to be "
               "considered.", "FilesControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    if (complete == false) {
      continue;
    }
    const CoordinateFileKind c_kind =
      (t_nml.getKeywordStatus("-sys", "c_kind", i) == InputStatus::MISSING) ?
      coordinate_input_format :
      translateCoordinateFileKind(t_nml.getStringValue("-sys", "c_kind", i));
    const CoordinateFileKind x_kind =
      (t_nml.getKeywordStatus("-sys", "x_kind", i) == InputStatus::MISSING) ?
      coordinate_input_format :
      translateCoordinateFileKind(t_nml.getStringValue("-sys", "x_kind", i));
    const CoordinateFileKind r_kind =
      (t_nml.getKeywordStatus("-sys", "r_kind", i) == InputStatus::MISSING) ?
      coordinate_input_format :
      translateCoordinateFileKind(t_nml.getStringValue("-sys", "r_kind", i));
    
    // Check for the existence of critical files
    const std::vector<std::string> sys_components = { "-p", "-c", "-x", "-r" };
    const std::vector<bool> existence_checks = { top_must_exist, crd_must_exist, trj_must_exist,
                                                 rst_must_exist };
    const std::vector<std::string> component_names = { "topology", "starting coordinates",
                                                       "trajectory", "checkpoint" };
    for (size_t j = 0; j < sys_components.size(); j++) {
      if (existence_checks[j] &&
          getDrivePathType(t_nml.getStringValue("-sys", sys_components[j], i)) !=
          DrivePathType::FILE) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("System " + std::to_string(i + 1) + " names a non-existent " + component_names[j] +
                " file " + t_nml.getStringValue("-sys", sys_components[j], i) + ".",
                "FilesControls");
        case ExceptionResponse::WARN:
          rtWarn("System " + std::to_string(i + 1) + " names a non-existent " +
                 component_names[j] + " file " +
                 t_nml.getStringValue("-sys", sys_components[j], i) + " and will be omitted.",
                 "FilesControls");
          complete = false;
          break;
        case ExceptionResponse::SILENT:
          complete = false;
          break;
        }
      }
    }
    if (complete) {
      if (strcmpCased(sys_label, "all", CaseSensitivity::NO) ||
          strcmpCased(sys_label, "all_possible", CaseSensitivity::NO) ||
          strncmpCased(sys_label, "pairedsystem", CaseSensitivity::NO, 12)) {
        rtErr("The system built from topology " + top_name + " has a label '" + sys_label +
              "', which collides with a reserved value for labelling systems.", "MoleculeSystem");
      }
      systems.push_back(MoleculeSystem(top_name, crd_name, trj_name, rst_name, sys_label,
                                       t_nml.getIntValue("-sys", "frame_start", i),
                                       t_nml.getIntValue("-sys", "frame_end", i),
                                       t_nml.getIntValue("-sys", "-n", i),
                                       c_kind, x_kind, r_kind));
    }
  }
  system_count = systems.size();
  
  // Get free topologies
  const int n_free_top = t_nml.getKeywordEntries("-p") *
                         (t_nml.getKeywordStatus("-p") != InputStatus::MISSING);
  for (int i = 0; i < n_free_top; i++) {
    const std::string fipath = t_nml.getStringValue("-p", i);
    const DrivePathType fitype = getDrivePathType(fipath);
    switch (fitype) {
    case DrivePathType::FILE:
      topology_file_names.push_back(fipath);
      break;
    case DrivePathType::DIRECTORY:
      {
        const std::vector<std::string> allfi = listDirectory(fipath);
        topology_file_names.insert(topology_file_names.end(), allfi.begin(), allfi.end());
      }
    case DrivePathType::REGEXP:
      {
        const std::vector<std::string> allfi = listFilesInPath(fipath, SearchStyle::RECURSIVE);
        topology_file_names.insert(topology_file_names.end(), allfi.begin(), allfi.end());        
      }
    }
  }
  free_topology_count = topology_file_names.size();

  // Get free coordinates
  const int n_free_crd = t_nml.getKeywordEntries("-c") *
                         (t_nml.getKeywordStatus("-c") != InputStatus::MISSING);
  for (int i = 0; i < n_free_crd; i++) {
    const std::string fipath = t_nml.getStringValue("-c", i);
    const DrivePathType fitype = getDrivePathType(fipath);
    switch (fitype) {
    case DrivePathType::FILE:
      coordinate_file_names.push_back(fipath);
      break;
    case DrivePathType::DIRECTORY:
      {
        const std::vector<std::string> allfi = listDirectory(fipath);
        coordinate_file_names.insert(coordinate_file_names.end(), allfi.begin(), allfi.end());
      }
    case DrivePathType::REGEXP:
      {
        const std::vector<std::string> allfi = listFilesInPath(fipath, SearchStyle::RECURSIVE);
        coordinate_file_names.insert(coordinate_file_names.end(), allfi.begin(), allfi.end());
      }
    }
  }
  free_coordinate_count = coordinate_file_names.size();  
}

//-------------------------------------------------------------------------------------------------
int FilesControls::getStructureCount() const {
  return structure_count;
}

//-------------------------------------------------------------------------------------------------
int FilesControls::getFreeTopologyCount() const {
  return free_topology_count;
}

//-------------------------------------------------------------------------------------------------
int FilesControls::getFreeCoordinatesCount() const {
  return free_coordinate_count;
}

//-------------------------------------------------------------------------------------------------
int FilesControls::getSystemDefinitionCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
bool FilesControls::readAllFreeFrames() const {
  return all_free_frames;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind FilesControls::getOutputCoordinateFormat() const {
  return coordinate_output_format;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind FilesControls::getCheckpointFormat() const {
  return coordinate_checkpoint_format;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getFreeTopologyName(const int index) const {
  return topology_file_names[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> FilesControls::getFreeTopologyNames() const {
  return topology_file_names;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getFreeCoordinateName(const int index) const {
  return coordinate_file_names[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> FilesControls::getFreeCoordinateNames() const {
  return coordinate_file_names;
}

//-------------------------------------------------------------------------------------------------
MoleculeSystem FilesControls::getSystem(int index) const {
  return systems[index];
}


//-------------------------------------------------------------------------------------------------
std::string FilesControls::getReportFile() const {
  return report_file;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getTrajectoryFileName() const {
  return coordinate_output_name;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getCheckpointFileName() const {
  return checkpoint_name;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getWarningFileName() const {
  return warning_file_name;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setAllFreeFrameReading(const bool active) {
  all_free_frames = active;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setOutputCoordinateFormat(const std::string &traj_kind) {
  coordinate_output_format = translateCoordinateFileKind(traj_kind);
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setOutputCoordinateFormat(const CoordinateFileKind traj_kind) {
  coordinate_output_format = traj_kind;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setCheckpointFormat(const std::string &chk_kind) {
  coordinate_checkpoint_format = translateCoordinateFileKind(chk_kind);
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setCheckpointFormat(const CoordinateFileKind chk_kind) {
  coordinate_checkpoint_format = chk_kind;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::addFreeTopologyName(const std::string &file_name) {
  if (findStringInVector(topology_file_names, file_name) == topology_file_names.size()) {
    topology_file_names.push_back(file_name);
    free_topology_count++;
  }
}

//-------------------------------------------------------------------------------------------------
void FilesControls::removeFreeTopologyName(const int starting_index, const int stretch) {
  if (starting_index < 0 || starting_index >= static_cast<int>(topology_file_names.size())) {
    rtErr("Unable to remove elements starting at index " + std::to_string(starting_index) +
          "from a list of " + std::to_string(topology_file_names.size()) + ".", "FilesControls",
          "removeFreeTopologyName");
  }
  if (starting_index + stretch >= static_cast<int>(topology_file_names.size())) {
    rtErr("Unable to remove " + std::to_string(stretch) + " elements starting at index " +
          std::to_string(starting_index) + "from a list of " +
          std::to_string(topology_file_names.size()) + ".", "FilesControls",
          "removeFreeTopologyName");
  }
  topology_file_names.erase(topology_file_names.begin() + starting_index,
                            topology_file_names.begin() + starting_index + stretch);
}

//-------------------------------------------------------------------------------------------------
void FilesControls::removeFreeTopologyName(const std::string &fname) {
  const size_t findex = findStringInVector(topology_file_names, fname);
  if (findex < topology_file_names.size()) {
    topology_file_names.erase(topology_file_names.begin() + findex);
  }
}

//-------------------------------------------------------------------------------------------------
void FilesControls::addFreeCoordinateName(const std::string &file_name) {
  if (findStringInVector(coordinate_file_names, file_name) == coordinate_file_names.size()) {
    coordinate_file_names.push_back(file_name);
    free_coordinate_count++;
    structure_count++;
  }
}

//-------------------------------------------------------------------------------------------------
void FilesControls::removeFreeCoordinateName(const int starting_index, const int stretch) {
  if (starting_index < 0 || starting_index >= static_cast<int>(coordinate_file_names.size())) {
    rtErr("Unable to remove elements starting at index " + std::to_string(starting_index) +
          "from a list of " + std::to_string(coordinate_file_names.size()) + ".", "FilesControls",
          "removeFreeCoordinateName");
  }
  if (starting_index + stretch >= static_cast<int>(coordinate_file_names.size())) {
    rtErr("Unable to remove " + std::to_string(stretch) + " elements starting at index " +
          std::to_string(starting_index) + "from a list of " +
          std::to_string(coordinate_file_names.size()) + ".", "FilesControls",
          "removeFreeCoordinateName");
  }
  coordinate_file_names.erase(coordinate_file_names.begin() + starting_index,
                              coordinate_file_names.begin() + starting_index + stretch);
}

//-------------------------------------------------------------------------------------------------
void FilesControls::removeFreeCoordinateName(const std::string &fname) {
  const size_t findex = findStringInVector(coordinate_file_names, fname);
  if (findex < coordinate_file_names.size()) {
    coordinate_file_names.erase(coordinate_file_names.begin() + findex);
  }
}

//-------------------------------------------------------------------------------------------------
void FilesControls::addSystem(const MoleculeSystem &new_mol) {
  systems.push_back(new_mol);
  system_count++;
  structure_count += new_mol.getTotalFrames();
}

//-------------------------------------------------------------------------------------------------
void FilesControls::removeSystem(const int starting_index, const int stretch) {
  if (starting_index < 0 || starting_index >= static_cast<int>(systems.size())) {
    rtErr("Unable to remove elements starting at index " + std::to_string(starting_index) +
          "from a list of " + std::to_string(systems.size()) + ".", "FilesControls",
          "removeSystem");
  }
  if (starting_index + stretch >= static_cast<int>(systems.size())) {
    rtErr("Unable to remove " + std::to_string(stretch) + " elements starting at index " +
          std::to_string(starting_index) + "from a list of " +
          std::to_string(systems.size()) + ".", "FilesControls", "removeSystem");
  }
  systems.erase(systems.begin() + starting_index, systems.begin() + starting_index + stretch);
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setReportFileName(const std::string &file_name) {
  report_file = file_name;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setGeneralTrajectoryFileName(const std::string &proto_name) {
  coordinate_output_name = proto_name;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setGeneralCheckpointFileName(const std::string &proto_name) {
  checkpoint_name = proto_name;
}

//-------------------------------------------------------------------------------------------------
void FilesControls::setWarningFileName(const std::string &file_name) {
  warning_file_name = file_name;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator filesInput(const TextFile &tf, int *start_line,
                            const std::vector<SubkeyRequirement> &sys_keyword_reqs,
                            const ExceptionResponse policy, const WrapTextSearch wrap,
                            const CoordinateFileKind crd_input_format,
                            const CoordinateFileKind crd_output_format,
                            const CoordinateFileKind crd_chkpt_format) {
  NamelistEmulator t_nml("files", CaseSensitivity::AUTOMATIC, policy, "Collects file names for "
                         "STORMM programs, offloading work that would otherwise require "
                         "command-line arguments.");
  t_nml.addKeyword(NamelistElement("-p", NamelistType::STRING,
                                   std::string(default_filecon_topology_name),
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("-c", NamelistType::STRING,
                                   std::string(default_filecon_coordinate_name),
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  const std::string sys_help("Expression for a complete system, linking a topology file "
                             "explicitly to a starting coordinates file, with the option of that "
                             "coordinates file being a trajectory with more than one frame.  This "
                             "keyword provides a means to read more than one frame from a "
                             "trajectory starting coordinates file, if the frame_end subkey is "
                             "given and greater than frame_start.  All starting coordinates will "
                             "be paired to the same topology object.  Like several other "
                             "specifiers in this namelist, this keyword is repeatable.");
  const std::vector<std::string> sys_keys_help = {
    "Topology file", "Starting coordinates file", "Output trajectory file", "Checkpoint file",
    "System label", "Starting frame (if the coordinates are a trajectory)", "Ending frame (if the "
    "coordinates are a trajectory).  If unspecified, only the starting frame will be read.  "
    "Otherwise, distinct systems will be made for the given topology and every frame between "
    "frame_start and frame_end.", "Type of coordinates file to expect (if unspecified, the type "
    "will be detected automatically)", "Type of trajectory file to write", "Type of checkpoint "
    "(restart) coordinates file to write" };
  
  // Prepare the requirements list for the -sys keyword
  t_nml.addKeyword(NamelistElement("-sys", { "-p", "-c", "-x", "-r", "-label", "frame_start",
                                             "frame_end", "-n", "c_kind", "x_kind", "r_kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::INTEGER,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), "0", "0", "1",
                                     getEnumerationName(crd_input_format),
                                     getEnumerationName(crd_output_format),
                                     getEnumerationName(crd_chkpt_format) },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sys_help,
                                   sys_keys_help, sys_keyword_reqs));
  t_nml.addKeyword(NamelistElement("-o", NamelistType::STRING,
                                   std::string(default_filecon_report_name)));
  t_nml.addKeyword(NamelistElement("-x", NamelistType::STRING,
                                   std::string(default_filecon_trajectory_name)));
  t_nml.addKeyword(NamelistElement("-r", NamelistType::STRING,
                                   std::string(default_filecon_checkpoint_name)));
  t_nml.addKeyword(NamelistElement("-wrn", NamelistType::STRING,
                                   std::string(default_filecon_warnings_name)));
  t_nml.addHelp("-p", "System topology file.  Repeatable for multiple systems.  Also accepts "
                "regular expressions.");
  t_nml.addHelp("-c", "Input coordinates file.  Repeatable for multiple systems.  Also accepts "
                "regular expressions.  Input coordinate files will be matched to topology files "
                "using atom counts and sanity of valence parameters, if free -c and -p parameters "
                "are provided.  Otherwise, use the system keyword and its subkeys to tie specific "
                "sets of starting coordinates to each topology.");
  t_nml.addHelp("-o", "Output diagnostics file, equivalent to mdout from Amber's sander program.  "
                "Reports for all systems will be included in this file.  Output can become "
                "voluminous for multiple systems if all dump their details into one file.  Use "
                "the \"outfmt\" keyword with setting \"INDIVIDUAL\" to obtain separate files for "
                "each system along with a \".master\" output file providing details of the entire "
                "run.");
  t_nml.addHelp("-x", "Trajectory output file (base name) for each system.  The actual name of "
                "each output file will be \"yyy_(sysID).zzz\", where \"yyy\" is any part of the "
                "-x string value preceding the final dot [.], \"_(sysID)\" is based on the name "
                "of the initial coordinates file, perhaps with a frame number appended, and "
                "\"zzz\" is an extension obtained from all content of the -x string following the "
                "final dot [.], or nothing if there is no dot.");
  t_nml.addHelp("-r", "Checkpoint (coordinate and velocity restart) file (base name) for each "
                "system.  As in the case of the trajectory output file specification, this is a "
                "fallback for free topology / coordinate pairs or systems with no specified "
                "restart file name.");
  t_nml.addHelp("-wrn", "Warnings reported for the run, collecting results from all systems.");
  
  // There is expected to be one unique &files namelist in a given input file.  Seek it out by
  // wrapping back to the beginning of the input file if necessary.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount());
  
  return t_nml;
}
  
} // namespace namelist
} // namespace stormm
