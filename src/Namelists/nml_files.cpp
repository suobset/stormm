#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "Trajectory/trajectory_enumerators.h"
#include "nml_files.h"

namespace omni {
namespace namelist {

using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using parse::findStringInVector;
using parse::WrapTextSearch;
using trajectory::translateCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem() :
    topology_file_name{}, coordinate_file_name{}, coordinate_output_name{}, checkpoint_name{},
    frame_start{0}, frame_end{0}, replica_count{0}, coordinate_kind{CoordinateFileKind::UNKNOWN},
    trajectory_kind{CoordinateFileKind::AMBER_CRD},
    checkpoint_kind{CoordinateFileKind::AMBER_ASCII_RST}
{}

//-------------------------------------------------------------------------------------------------
MoleculeSystem::MoleculeSystem(const std::string &topology_file_in,
                               const std::string &coordinate_file_in,
                               const std::string &trajectory_file_in,
                               const std::string &checkpoint_file_in, const int frame_start_in,
                               const int frame_end_in, const int replica_count_in,
                               const CoordinateFileKind coordinate_kind_in,
                               const CoordinateFileKind trajectory_kind_in,
                               const CoordinateFileKind checkpoint_kind_in) :
    topology_file_name{topology_file_in},
    coordinate_file_name{coordinate_file_in},
    coordinate_output_name{trajectory_file_in},
    checkpoint_name{checkpoint_file_in},
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
FilesControls::FilesControls(const ExceptionResponse policy_in) :
    policy{policy_in}, structure_count{0}, free_topology_count{0}, free_coordinate_count{0},
    system_count{0},
    all_free_frames{default_filecon_read_all_free},
    coordinate_output_format{translateCoordinateFileKind(default_filecon_outcrd_type)},
    coordinate_checkpoint_format{translateCoordinateFileKind(default_filecon_chkcrd_type)},
    topology_file_names{}, coordinate_file_names{}, systems{},
    report_file{std::string(default_filecon_report_name)},
    coordinate_output_name{std::string(default_filecon_trajectory_name)},
    checkpoint_name{std::string(default_filecon_checkpoint_name)}
{}

//-------------------------------------------------------------------------------------------------
FilesControls::FilesControls(const TextFile &tf, int *start_line,
                             const ExceptionResponse policy_in,
                             const std::vector<std::string> &alternatives) :
    FilesControls(policy_in)
{
  // Set some alternative defaults.  This takes a vector of strings, the even-numbered strings
  // being names of actual member variables and the odd-numbered strings being the new defaults to
  // apply.  Different applications will then be able to call the constructor with different
  // default settings.
  const int n_alt = alternatives.size();
  for (int i = 0; i < n_alt; i++) {
    if (i < n_alt - 1) {
      if (alternatives[i] == std::string("coordinate_output_format")) {
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
  NamelistEmulator t_nml = filesInput(tf, start_line, policy);
  const int nsys = t_nml.getKeywordEntries("-sys") *
                   (t_nml.getKeywordStatus("-sys") != InputStatus::MISSING);
  for (int i = 0; i < nsys; i++) {
    bool complete = true;
    const std::string top_name = t_nml.getStringValue("-sys", "-p", i);
    const std::string crd_name = t_nml.getStringValue("-sys", "-c", i);
    const std::string trj_name = t_nml.getStringValue("-sys", "-x", i);
    const std::string rst_name = t_nml.getStringValue("-sys", "-r", i);
    std::string missing_elements("");
    if (top_name.size() == 0LLU) {
      missing_elements += "-p";
    }
    if (crd_name.size() == 0LLU) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -c" : "-c";
    }
    if (trj_name.size() == 0LLU) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -x" : "-x";
    }
    if (rst_name.size() == 0LLU) {
      missing_elements += (missing_elements.size() > 0LLU) ? ", -r" : "-r";
    }
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
    CoordinateFileKind c_kind, x_kind, r_kind;
    c_kind = translateCoordinateFileKind(t_nml.getStringValue("-sys", "-c_kind", i));
    x_kind = translateCoordinateFileKind(t_nml.getStringValue("-sys", "-x_kind", i));
    r_kind = translateCoordinateFileKind(t_nml.getStringValue("-sys", "-r_kind", i));
    systems.push_back(MoleculeSystem(t_nml.getStringValue("-sys", "-p", i),
                                     t_nml.getStringValue("-sys", "-c", i),
                                     t_nml.getStringValue("-sys", "-x", i),
                                     t_nml.getStringValue("-sys", "-r", i),
                                     t_nml.getIntValue("-sys", "frame_start", i),
                                     t_nml.getIntValue("-sys", "frame_end", i),
                                     t_nml.getIntValue("-sys", "-n", i), c_kind, x_kind, r_kind));
  }
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
void FilesControls::addFreeCoordinateName(const std::string &file_name) {
  if (findStringInVector(coordinate_file_names, file_name) == coordinate_file_names.size()) {
    coordinate_file_names.push_back(file_name);
    free_coordinate_count++;
    structure_count++;
  }
}

//-------------------------------------------------------------------------------------------------
void FilesControls::addSystem(const MoleculeSystem &new_mol) {
  systems.push_back(new_mol);
  system_count++;
  structure_count += new_mol.getTotalFrames();
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
NamelistEmulator filesInput(const TextFile &tf, int *start_line, const ExceptionResponse policy) {
  NamelistEmulator t_nml("files", CaseSensitivity::AUTOMATIC, policy, "Collects file names for "
                         "OMNI programs, offloading work that would otherwise require "
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
    "Starting frame (if the coordinates are a trajectory)", "Ending frame (if the coordinates are "
    "a trajectory).  If unspecified, only the starting frame will be read.  Otherwise, distinct "
    "systems will be made for the given topology and every frame between frame_start and "
    "frame_end.", "Type of coordinates file to expect (if unspecified, the type will be detected "
    "automatically)", "Type of trajectory file to write", "Type of checkpoint (restart) "
    "coordinates file to write"
  };
  t_nml.addKeyword(NamelistElement("-sys", { "-p", "-c", "-x", "-r", "frame_start", "frame_end",
                                             "-n", "c_kind", "x_kind", "r_kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::INTEGER, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), "0", "0", "1",
                                     std::string(default_filecon_inpcrd_type),
                                     std::string(default_filecon_outcrd_type),
                                     std::string(default_filecon_chkcrd_type) },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sys_help,
                                   sys_keys_help));
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
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount());
  
  return t_nml;
}
  
} // namespace namelist
} // namespace omni
