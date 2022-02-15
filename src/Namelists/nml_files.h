// -*-c++-*-
#ifndef OMNI_NML_RST_H
#define OMNI_NML_RST_H

#include <string>
#include <vector>
#include "input.h"
#include "namelist_emulator.h"
#include "Parsing/textfile.h"
#include "Trajectory/trajectory_enumerators.h"

namespace omni {
namespace namelist {

using trajectory::CoordinateFileKind;

/// \brief Default file names and extensions.  The default input file name varies according to the
///        particular application and is therefore not defined within these libraries.
/// \{
constexpr char default_filecon_topology_name[] = "prmtop";
constexpr char default_filecon_coordinate_name[] = "inpcrd";
constexpr char default_filecon_report_name[] = "mdout";
constexpr char default_filecon_trajectory_base[] = "mdcrd";
constexpr char default_filecon_trajectory_ext[] = ".crd";
constexpr char default_filecon_checkpoint_base[] = "mdrst";
constexpr char default_filecon_checkpoint_ext[] = ".rst";
/// \}

/// \brief Object to encapsulate a system, a coupled set of coordinates and a single topology.
///        This can be used to enforce conf.omni to read certain coordinates in the context of a
///        particular topology, even when other viable topologies might be available.
struct MoleculeSystem {

  /// \brief The constructor can make a blank system or automatically fill values.
  ///
  /// Overloaded:
  ///   - Create a blank object
  ///   - Create an object with a complete list of files and read settings
  /// \{
  MoleculeSystem();
  MoleculeSystem(const std::string &topology_file_in, const std::string &coordinate_file_in,
                 int frame_start_in, int frame_end_in, int replica_count_in,
                 CoordinateFileKind coordinate_kind_in);
  /// \}

private:
  std::string topology_file_name;      ///< Topology file describing the molecular system
  std::string coordinate_file_name;    ///< Coordinate file (may be a whole trajectory)
  int frame_start;                     ///< Strating frame of the coordinates file to read
  int frame_end;                       ///< Final frame of the coordinates file to read (the
                                       ///<   program will read frames over the range
                                       ///<   [ start, end ], so [ 0, 0 ] gets the first frame)
  int replica_count;                   ///< Number of times to replicate this system
  CoordinateFileKind coordinate_kind;  ///< Kind of coordinate file to be expected
};

/// \brief Distill the results of file identification, producing clean lists of free topologies,
///        free coordinate files, and linked topology / coordinate systems.  Record output and
///        trajectory files, if available.
struct FilesControls {

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param policy_in       Requested error handling behavior
  /// \param tf              Input file translated into RAM
  /// \param start_line      Line of the input file to begin searching for the &solvent namelist
  /// \param report_name_in  New default name of the report file (for application-specific naming
  ///                        conventions)
  /// \param coord_base_in   New default base name of trajectory or other coordinate files (for
  ///                        application-specific naming conventions)
  /// \param coord_ext_in    New default extension of trajectory or other coordinate files (for
  ///                        application-specific naming conventions)
  /// \{
  FilesControls();
  FilesControls(const TextFile &tf, int *start_line,
                ExceptionResponse policy = ExceptionResponse::DIE,
                const std::string &report_name_in = std::string(default_filecon_report_name),
                const std::string &coord_base_in = std::string(default_filecon_trajectory_base),
                const std::string &coord_ext_in = std::string(default_filecon_trajectory_ext),
                const std::string &chk_base_in = std::string(default_filecon_checkpoint_base),
                const std::string &chk_ext_in = std::string(default_filecon_checkpoint_ext));
  /// \}

  /// \brief Get the structure count, based on the number of free coordinate files as well as the
  ///        number of systems, with frame counts and replicas therein.  If one of the systems
  ///        does not have the requisite number of frames, this is an error and may be trapped.
  ///        The structure count reported here is therefore a maximum that can be expected.
  int getStructureCount() const;

  /// \brief Get the free topology count.
  int getFreeTopologyCount() const;

  /// \brief Get the free coordinate file count.
  int getFreeCoordinatesCount() const;

  /// \brief Get the number of system specifications made with the -sys keyword.
  int getSystemDefinitionCount() const;

  /// \brief Get the coordinate (trajectory) file output format
  CoordinateFileKind getOutputCoordinateFormat() const;
  
  /// \brief Get the coordinate (checkpoint) file output format
  CoordinateFileKind getOutputCheckpointFormat() const;
  
  /// \brief Get one or more free topology names.
  ///
  /// Overloaded:
  ///   - Get the free topology named at one index of the array
  ///   - Get all free topology names in the array
  ///
  /// \param index  The index to query
  /// \{
  std::string getFreeTopologyName(int index) const;
  std::vector<std::string> getFreeTopologyNames() const;
  /// \}

  /// \brief Get one or more free coordinate file names.
  ///
  /// Overloaded:
  ///   - Get the free coordinate file named at one index of the array
  ///   - Get all free coordinate file names in the array
  ///
  /// \param index  The index to query
  /// \{
  std::string getFreeCoordinateName(int index) const;
  std::vector<std::string> getFreeCoordinateNames() const;
  /// \}

  /// \brief Get a molecule system from this object's array.
  ///
  /// \param index  The index to query
  MoleculeSystem getSystem(int index) const;

  /// \brief Get the name of the report file (equivalent to mdout in sander or pmemd)
  std::string getReportFile() const;

  /// \brief Get the base name of trajectory files to write
  std::string getTrajectoryFileBase() const;

  /// \brief Get the file extension of trajectory files to write
  std::string getTrajectoryFileExtension() const;

  /// \brief Get the base name of (coordinate) checkpoint files to write
  std::string getCheckpointFileBase() const;

  /// \brief Get the file extension of (coordinate) checkpoint files to write
  std::string getCheckpointFileExtension() const;

private:

  // Counts of critical data
  int structure_count;        ///< Total number of initial structures, the sum of all free
                              ///<   coordinate files (for which only the first frame will be read,
                              ///<   if multiple frames are present) and all frames arising from
                              ///<   system definitions (which can include trajectory files and
                              ///<   span multiple frames)
  int free_topology_count;    ///< The number of free topologies, which will be matched to free
                              ///<   coordinate files
  int free_coordinate_count;  ///< The number of free coordinate files, which will be paired to
                              ///<   free trajectories
  int system_count;           ///< The number of system keyword specifications

  /// Format of the coordinate output files
  CoordinateFileKind coordinate_output_format;

  /// Format of the coordinate checkpoint files
  CoordinateFileKind coordinate_checkpoint_format;

  /// List of free topologies.  These names can indicate specific files, directories containing
  /// files, or regular expressions.  All such files will be evaluated as possible topologies,
  /// and read as individual objects if they appear promising.
  std::vector<std::string> topology_file_names;

  /// List of free coordinate files or trajectories.  These names can indicate specific files,
  /// directories containing files, or regular expressions.  All such files will be evaluated as
  /// possible sources of coordinates.  If trajectories are supplied, the default behavior will be
  /// to take only the first frame from each trajectory, unless the flag all_free_trajectory_frames
  /// is set to true.
  std::vector<std::string> coordinate_file_names;

  /// List of specific systems.  This will be initialized with user input, if there are any
  /// -sys arguments supplied in the &files namelist.  Afterwards, the list will be extended by
  /// adding any viable topology and coordinate pairs that can be detected from the two lists of
  /// independent topologies and coordinate file names.
  std::vector<MoleculeSystem> systems;

  /// Name of the output file.  This is akin to sander's mdout but much more involved as it
  /// spans all systems.
  std::string report_file;

  /// Conformation output base name.  Each molecule will become a set of conformers, expressed in
  /// a trajectory file with this as the base name.  The typical process will then be to use the
  /// system topology's root name (the file path base name, with any subdirectories clipped from
  /// the front and any extensions after a dot [.] clipped form the back) as the next part, and
  /// finally an extension set with the next control variable.
  std::string coordinate_output_base;

  /// Conformation output extension.  Each molecule's set of conformations will get a separate
  /// file, the name determined from a variety of user-specified factors, with this extension.
  std::string coordinate_output_ext;

  /// Checkpoint file base name, i.e. the restart file in Amber's sander or pmemd
  std::string checkpoint_base;

  /// Checkpoint file extension, i.e. the restart file in Amber's sander or pmemd
  std::string checkpoint_ext;
};
  
/// \brief Produce a namelist for specifying basic input and output files, which can take the place
///        of a great deal of command line input in the Amber pmemd and sander programs.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this will
///                    wrap back to the beginning of the file in search of a unique &files
///                    namelist)
/// \param policy      Reaction to exceptions encountered during namelist reading
NamelistEmulator filesInput(const TextFile &tf, int *start_line,
                            ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace namelist
} // namespace omni

#endif
