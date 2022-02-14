// -*-c++-*-
#ifndef OMNI_USER_SETTINGS_H
#define OMNI_USER_SETTINGS_H

#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"

namespace conf_app {
namespace user_input {

using omni::topology::AtomGraph;
using omni::topology::ImplicitSolventModel;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFileKind;

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
  std::string topology_file_name;
  std::string coordinate_file_name;
  int frame_start;
  int frame_end;
  int replica_count;
  CoordinateFileKind coordinate_kind;
};

/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program
  ///
  /// \param argc  Number of command-line variables
  /// \param argv  List of command line argument strings
  UserSettings(int argc, const char* argv[]);

  /// \brief Get the number of systems
  int getSystemCount() const;
  
private:

  /// Name of the original input file
  std::string input_file;
  
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

  /// Atom mask specifying common core atoms.  This will be applied to all systems to generate
  /// some subset of the atoms, which will then be subject to additional restraints during
  /// energy minimizations.
  std::string common_core_mask;

  /// Name of the output file.  This is akin to sander's mdout but much more involved as it
  /// spans all systems.
  std::string report_file;

  // Scalars and user-specified run settings
  bool all_free_trajectory_frames;  ///< Flag to have all frames, rather than just the first, read
                                    ///<   from free trajectories
  int system_count;                 ///< Total number of systems, after processing user-defined
                                    ///<   systems and then the free topology and coordinate lists

  // Energy minimization and target planning
  int igseed;                       ///< Random number seed
  int random_stream_count;          ///< Number of distinct random streams, the kth stream being
                                    ///<   based on the kth jump forward from the initial seed in
                                    ///<   the xoshiro256++ generator.
  int minimization_steps;           ///< Maximum number of minimization steps to attempt for any
                                    ///<   individual system (equivalent to maxcyc in sander)
  int steepest_descent_steps;       ///< Number of steepest descent steps to perform in each
                                    ///<   system prior to beginning conjugate gradient moves
                                    ///<   (equivalent to ncyc in sander)
  double initial_move_length;       ///< Magnitude of the initial displacement along the gradient
                                    ///<   vector.  The size of subsequent moves will grow or
                                    ///<   shrink based on the history of success in previous
                                    ///<   optimizations.
  double convergence_criterion;     ///< Convergence criterion for the minimization, in kcal/mol.

  // Implicit solvent model
  ImplicitSolventModel igb;         ///< The type of Generalized Born or other implicit solvent
                                    ///<   model to use
  double born_radii_cutoff;         ///< Cutoff for pairs of atoms to participate in Generalized
                                    ///<   Born radii (equivalent to rgbmax in sander)
  double dielectric;                ///< Sets the dielectric constant for all molecules (applicable
                                    ///<   to implicit solvent situations only)
  double salt_concentration;        ///< Sets the salt concentration for all molecules (applicable
                                    ///<   to implicit solvent situations only)

  /// An array of all topologies to be read by the system: all free topologies and all topologies
  /// read as part of a MoleculeSystem.
  std::vector<AtomGraph> topology_cache;

  /// An array of all coordinate sets to be read by the system: all free coordinate sets and all
  /// coordinates read as part of a MoleculeSystem.
  std::vector<CoordinateFrame> initial_coordinates_cache;

  /// The vector of all topology indices guiding each simulation.  This may contain repeats, if
  /// the various MoleculeSystem objects contain the same topology, but the list will be reduced
  /// when composing the synthesis objects.
  std::vector<int> topology_indices;
};

} // namespace user_input
} // namespace conf_app

#endif
