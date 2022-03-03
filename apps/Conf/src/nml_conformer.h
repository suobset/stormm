// -*-c++-*-
#ifndef OMNI_NML_CONFORMER_H
#define OMNI_NML_CONFORMER_H

#include "../../../src/Constants/behavior.h"
#include "../../../src/Namelists/namelist_emulator.h"
#include "../../../src/Namelists/namelist_element.h"
#include "../../../src/Parsing/textfile.h"

namespace conf_app {
namespace user_input {

using omni::constants::ExceptionResponse;
using omni::namelist::NamelistEmulator;
using omni::parse::TextFile;

/// \brief Default input settings for the &conformer namelist
/// \{
constexpr int default_conf_rotation_samples     = 3;
constexpr char default_conf_chirality[]         = "false";
constexpr char default_conf_cis_trans[]         = "false";
constexpr char default_conf_stop_hbonds[]       = "false";
constexpr int default_conf_running_states       = 16;
constexpr int default_conf_final_states         = 100;
constexpr double default_conf_rmsd_tolerance    = 0.5;
constexpr int default_conf_reshuffle_iterations = 0;
/// \}

/// \brief Object to encapsulate the data that can be extracted from the &conformer namelist.
struct ConformerControls {

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param policy_in   Requested error handling behavior
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \{
  ConformerControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  ConformerControls(const TextFile &tf, int *start_line,
                    ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

private:
  ExceptionResponse policy;       ///< Set the behavior when bad inputs are encountered.  DIE =
                                  ///<   abort program, WARN = warn the user, and likely reset to
                                  ///<   the default value if one is available, SILENT = do not
                                  ///<   warn the user, but also likely reset to the default value
                                  ///<   if one is available.
  std::string common_atom_mask;   ///< Mask of common core atoms, applied to all systems
  bool sample_chirality;          ///< Flag to have chiral enantiomers sampled
  bool sample_cis_trans;          ///< Flag to have cis-trans isomers sampled
  bool prevent_hbonds;            ///< Flag to apply restraints that will prevent hydrogen bonds
                                  ///<   from forming during energy minimizations
  int running_states;             ///< Number of states to try minimizing at one time
  int final_states;               ///< Number of final states to collect
  int rotation_samples;           ///< Number of times to sample about a rotatable bond
  double rmsd_tolerance;          ///< Minimum mass-weighted root-mean squared deviation between
                                  ///<   unique conformers
};

/// \brief Free function to read the &conformer namelist.  This works in analogous fashion to
///        various namelist-assoicated objects implemented in src/Namelists/, except that the
///        &conformer namelist is private to this application.  The UserSettings object is the
///        object associated in this manner to the &conformer namelist.  Unlike the objects linked
///        to general-purpose namelists in src/Namelists, the UserSettings object does not need
///        additional setter functions for keyword-associated content in a &confomer namelist, as
///        the namelist in an input file is the one way that this information will enter the
///        program.
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param policy      Response to bad inputs
NamelistEmulator conformerInput(const TextFile &tf, int *start_line, ExceptionResponse policy);
  
} // namespace user_input
} // namespace conf_app

#endif
