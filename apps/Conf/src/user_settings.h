// -*-c++-*-
#ifndef OMNI_USER_SETTINGS_H
#define OMNI_USER_SETTINGS_H

#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"

namespace conf_app {
namespace user_input {

using omni::namelist::FilesControls;
using omni::namelist::MinimizeControls;
using omni::namelist::RandomControls;
using omni::namelist::SolventControls;
using omni::topology::AtomGraph;
using omni::topology::ImplicitSolventModel;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFileKind;

/// \brief Default input settings for conformer.omni
/// \{
constexpr char default_conformer_input_file[] = "cgen.in";
/// \}

/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program
  ///
  /// \param argc  Number of command-line variables
  /// \param argv  List of command line argument strings
  UserSettings(int argc, const char* argv[]);

private:

  /// Name of the original input file
  std::string input_file;
  
  /// Atom mask specifying common core atoms.  This will be applied to all systems to generate
  /// some subset of the atoms, which will then be subject to additional restraints during
  /// energy minimizations.
  std::string common_core_mask;

  // Control parameters: these structs each encapsulate their own namelist from the input file.
  FilesControls file_io_input;      ///< All input and output file names, save for the command file
  MinimizeControls line_min_input;  ///< Line minimization directives
  SolventControls solvent_input;    ///< Implicit solvent specifications
  RandomControls prng_input;        ///< Random number generator specifications

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

  /// \brief Create vectors of topologies, starting coordinates, and indices of which topology
  ///        controlsthe motion of each coordinate set.
  void gatherUniqueTopologies();
};

} // namespace user_input
} // namespace conf_app

#endif
