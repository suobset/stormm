// -*-c++-*-
#ifndef OMNI_USER_SETTINGS_H
#define OMNI_USER_SETTINGS_H

#include "../../../src/Constants/behavior.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Parsing/textfile.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "nml_conformer.h"

namespace conf_app {
namespace user_input {

using omni::constants::ExceptionResponse;
using omni::namelist::FilesControls;
using omni::namelist::MinimizeControls;
using omni::namelist::RandomControls;
using omni::namelist::SolventControls;
using omni::parse::TextFile;
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

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program.
  ///
  /// \param argc  Number of command-line variables
  /// \param argv  List of command line argument strings
  UserSettings(int argc, const char* argv[]);

  /// \brief Get the policy (for passing to other operations that may trigger an exception).
  ExceptionResponse getExceptionBehavior() const;

  /// \brief Get the name of the input file.
  std::string getInputFileName() const;

  /// \brief Get the block of information associated with the &files namelist.
  FilesControls getFilesNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &minimize namelist.
  MinimizeControls getMinimizeNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &solvent namelist.
  SolventControls getSolventNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &random namelist.
  RandomControls getRandomNamelistInfo() const;
  
private:

  /// Action in the event of bad input
  ExceptionResponse policy;
  
  /// Name of the original input file
  std::string input_file;
  
  // Control parameters: these structs each encapsulate their own namelist from the input file.
  FilesControls file_io_input;      ///< All input and output file names, save for the command file
  MinimizeControls line_min_input;  ///< Line minimization directives
  SolventControls solvent_input;    ///< Implicit solvent specifications
  RandomControls prng_input;        ///< Random number generator specifications
  ConformerControls conf_input;     ///< Conformer generation instructions

  /// Atom mask specifying common core atoms.  This will be applied to all systems to generate
  /// some subset of the atoms, which will then be subject to additional restraints during
  /// energy minimizations.
  std::string common_core_mask;
};

} // namespace user_input
} // namespace conf_app

#endif
