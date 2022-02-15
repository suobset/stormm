#include "nml_files.h"

namespace omni {
namespace namelist {

using parse::WrapTextSearch;

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
FilesControls::FilesControls() :
    structure_count{0}, free_topology_count{0}, free_coordinate_count{0}, system_count{0},
    coordinate_output_format{CoordinateFileKind::AMBER_CRD},
    coordinate_checkpoint_format{CoordinateFileKind::AMBER_ASCII_RST},
    topology_file_names{}, coordinate_file_names{}, systems{},
    report_file{std::string(default_filecon_report_name)},
    coordinate_output_base{std::string(default_filecon_trajectory_base)},
    coordinate_output_ext{std::string(default_filecon_trajectory_ext)},
    checkpoint_base{std::string(default_filecon_checkpoint_base)},
    checkpoint_ext{std::string(default_filecon_checkpoint_ext)}
{

}

//-------------------------------------------------------------------------------------------------
FilesControls::FilesControls(const TextFile &tf, int *start_line, const ExceptionResponse policy,
                             const std::string &report_name_in, const std::string &coord_base_in,
                             const std::string &coord_ext_in, const std::string &chk_base_in,
                             const std::string &chk_ext_in) :
    structure_count{0}, free_topology_count{0}, free_coordinate_count{0}, system_count{0},
    coordinate_output_format{CoordinateFileKind::AMBER_CRD},
    coordinate_checkpoint_format{CoordinateFileKind::AMBER_ASCII_RST},
    topology_file_names{}, coordinate_file_names{}, systems{},
    report_file{report_name_in},
    coordinate_output_base{coord_base_in},
    coordinate_output_ext{coord_ext_in},
    checkpoint_base{chk_base_in},
    checkpoint_ext{chk_ext_in}
{
  NamelistEmulator t_nml = filesInput(tf, start_line, policy);

  
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
CoordinateFileKind FilesControls::getOutputCoordinateFormat() const {
  return coordinate_output_format;
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind FilesControls::getOutputCheckpointFormat() const {
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
std::string FilesControls::getTrajectoryFileBase() const {
  return coordinate_output_base;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getTrajectoryFileExtension() const {
  return coordinate_output_ext;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getCheckpointFileBase() const {
  return checkpoint_base;
}

//-------------------------------------------------------------------------------------------------
std::string FilesControls::getCheckpointFileExtension() const {
  return checkpoint_ext;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator filesInput(const TextFile &tf, int *start_line, const ExceptionResponse policy) {
  NamelistEmulator t_nml("files", CaseSensitivity::AUTOMATIC, policy, "Collects file names for "
                         "OMNI programs, offloading work that would otherwise require "
                         "command-line arguments.");
  t_nml.addKeyword(NamelistElement("-p", NamelistType::STRING, "prmtop", DefaultIsObligatory::NO,
                                   InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("-c", NamelistType::STRING, "inpcrd", DefaultIsObligatory::NO,
                                   InputRepeats::YES));
  const std::string sys_help("Expression for a complete system, linking a topology file "
                             "explicitly to a starting coordinates file, with the option of that "
                             "coordinates file being a trajectory with more than one frame.  This "
                             "keyword provides a means to read more than one frame from a "
                             "trajectory starting coordinates file, if the frame_end subkey is "
                             "given and greater than frame_start.  All starting coordinates will "
                             "be paired to the same topology object.  Like several other "
                             "specifiers in this namelist, this keyword is repeatable.");
  const std::vector<std::string> sys_keys_help = {
    "Topology file", "Starting coordinates file", "Starting frame (if the coordinates are a "
    "trajectory)", "Ending frame (if the coordinates are a trajectory).  If unspecified, only the "
    "starting frame will be read.  Otherwise, distinct systems will be made for the given "
    "topology and every frame between frame_start and frame_end.", "Type of coordinates file (if "
    "unspecified, the type will be detected automatically."
  };
  t_nml.addKeyword(NamelistElement("-sys", { "-p", "-c", "frame_start", "frame_end", "-n",
                                             "kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::STRING },
                                   { "prmtop", "inpcrd", "0", "0", "1", "AMBER_ASCII_RST" },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sys_help,
                                   sys_keys_help));
  t_nml.addKeyword(NamelistElement("-o", NamelistType::STRING, "mdout"));
  t_nml.addKeyword(NamelistElement("-x", NamelistType::STRING, "mdcrd"));
  t_nml.addKeyword(NamelistElement("-warn", NamelistType::STRING, "warnings"));
  t_nml.addKeyword(NamelistElement("-error", NamelistType::STRING, "errors"));
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
                "-x string value preceding the final dot [.]. \"_(sysID)\" is the number of the "
                "system in some internal list (key printed in the master output file) for free "
                "coordinate inputs matched to free topologies, or an identifier of the system "
                "for couple coordinate and topology pairs specified with the -sys command.  In "
                "both cases \"zzz\" is an extension obtained from all content of the -x string "
                "following the final dot [.], or nothing if there is no dot.");
  t_nml.addHelp("-warn", "Warnings reported for the run, collecting results from all systems.");
  t_nml.addHelp("-error", "Errors reported for the run, collecting results from all systems.");

  // There is expected to be one unique &files namelist in a given input file.  Seek it out by
  // wrapping back to the beginning of the input file if necessary.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount());
  
  return t_nml;
}
  
} // namespace namelist
} // namespace omni
