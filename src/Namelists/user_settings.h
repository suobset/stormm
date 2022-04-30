// -*-c++-*-
#ifndef OMNI_USER_SETTINGS_H
#define OMNI_USER_SETTINGS_H

#include "Constants/behavior.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_random.h"
#include "Namelists/nml_solvent.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/trajectory_enumerators.h"
#include "nml_conformer.h"

namespace omni {
namespace namelist {

using constants::ExceptionResponse;
using namelist::DynamicsControls;
using namelist::FFMorphControls;
using namelist::FilesControls;
using namelist::MinimizeControls;
using namelist::RandomControls;
using namelist::SolventControls;
using parse::TextFile;
using topology::AtomGraph;
using topology::ImplicitSolventModel;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFileKind;

/// \brief Default input settings for conformer.omni
/// \{
constexpr char default_conformer_input_file[] = "cgen.in";
/// \}

/// \brief Enumerate the various apps using OMNI libraries, and which use this UserSettings object
///        as a way to collect common namelists and other command line input.  This enumerator
///        makes it possible to tailor the contents of UserSettings to a specific application, and
///        to tune the behavior in response to particular inputs.
enum class AppName {
  CONFORMER,  ///< The OMNI conformer generator
  DYNAMICS,   ///< The OMNI molecular dynamics program
  FFREFINE    ///< The OMNI force field refinement tool
};
  
/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

  /// \brief The constructor requires an input file, similar to mdin for the Amber sander program.
  ///
  /// \param argc  Number of command-line variables
  /// \param argv  List of command line argument strings
  UserSettings(int argc, const char* argv[], AppName prog_set);

  /// \brief Get the policy (for passing to other operations that may trigger an exception).
  ExceptionResponse getExceptionBehavior() const;

  /// \brief Get the name of the input file.
  std::string getInputFileName() const;
  
  /// \brief Detect whether a &minimize namelist was present
  bool getMinimizePresence() const;

  /// \brief Detect whether a &conformer namelist was present
  bool getConformerPresence() const;

  /// \brief Get the block of information associated with the &files namelist.
  FilesControls getFilesNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &minimize namelist.
  MinimizeControls getMinimizeNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &solvent namelist.
  SolventControls getSolventNamelistInfo() const;
  
  /// \brief Get the block of information associated with the &random namelist.
  RandomControls getRandomNamelistInfo() const;

  /// \brief Get the block of information associated with the &conformer namelist.
  ConformerControls getConformerNamelistInfo() const;

private:

  ExceptionResponse policy;   /// Action in the event of bad input
  bool has_minimize_nml;      /// Indicate the presence of a &minimize namelist in the input file
  bool has_solvent_nml;       /// Indicate the presence of a &solvent namelist in the input file
  bool has_random_nml;        /// Indicate the presence of a &random namelist in the input file
  bool has_conformer_nml;     /// Indicate the presence of a &conformer namelist in the input file
  bool has_dynamics_nml;      /// Indicate the presence of a &dynamics namelist in the input file
  bool has_ffmorph_nml;       /// Indicate the presence of an &ffmorph namelist in the input file
  
  /// Name of the original input file
  std::string input_file;
  
  // Control parameters: these structs each encapsulate their own namelist from the input file.
  FilesControls file_io_input;      ///< All input and output file names, save for the command file
  MinimizeControls line_min_input;  ///< Line minimization directives
  SolventControls solvent_input;    ///< Implicit solvent specifications
  RandomControls prng_input;        ///< Random number generator specifications
  ConformerControls conf_input;     ///< Conformer generation instructions
  DynamicsControls dyna_input;      ///< Molecular dynamics instructions
  FFMorphControls ffmod_input;      ///< Force field modification instructions
};

} // namespace namelist
} // namespace omni

#endif
