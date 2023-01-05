// -*-c++-*-
#ifndef STORMM_NML_CONFORMER_H
#define STORMM_NML_CONFORMER_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/namelist_emulator.h"
#include "Namelists/namelist_element.h"
#include "Parsing/textfile.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using parse::TextFile;
using parse::WrapTextSearch;

/// \brief Default input settings for the &conformer namelist
/// \{
constexpr int default_conf_rotation_samples     = 3;
constexpr int default_conf_max_rotatable_bonds  = 4;
constexpr int default_conf_max_seeding_attempts = 3;
constexpr int default_conf_max_system_trials    = 16384;
constexpr int default_conf_running_states       = 16;
constexpr int default_conf_final_states         = 100;
constexpr int default_conf_reshuffle_iterations = 0;
constexpr int active_states_limit               = 524288;
constexpr double default_conf_rmsd_tolerance    = 0.5;
constexpr double default_conf_core_restraint    = 16.0;
constexpr double default_conf_self_clash_ratio  = 0.5;
constexpr char default_conf_chirality[]         = "false";
constexpr char default_conf_cis_trans[]         = "false";
constexpr char default_conf_stop_hbonds[]       = "false";
/// \}

/// \brief Object to encapsulate the data that can be extracted from the &conformer namelist.
struct ConformerControls {

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param policy_in   Requested error handling behavior
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  ConformerControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                    WrapTextSearch wrap = WrapTextSearch::NO);
  ConformerControls(const TextFile &tf, int *start_line, bool *found_nml,
                    ExceptionResponse policy_in = ExceptionResponse::DIE,
                    WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief Get the core atom mask string.
  const std::string& getCoreAtomMask() const;

  /// \brief Get the name of the data item expected to contain core atoms.
  const std::string& getCoreDataItemName() const;

  /// \brief Get the core restraint r2 penalty (near-repulsive potential) stiffness.
  double getCoreRK2Value() const;

  /// \brief Get the core restraint r3 penalty (position-keeping potential) stiffness.
  double getCoreRK3Value() const;

  /// \brief Get the core restraint r2 displacement.
  double getCoreR2Value() const;

  /// \brief Get the core restraint r3 displacement.
  double getCoreR3Value() const;

  /// \brief Get an indicator of whether to sample chirality.
  bool sampleChirality() const;

  /// \brief Get an indicator of whether to sample cis- and trans- isomers.
  bool sampleCisTrans() const;

  /// \brief Get an indicator as to whether to apply restraints that will prevent hydrogen bond
  ///        formation in the resulting conformers.
  bool preventHydrogenBonding() const;

  /// \brief Get the total number of states to attempt minimizing at one time.  This will put a
  ///        limit on the expanded population of conformer systems that the program will attempt
  ///        to model and minimize on the GPU, which takes a significant amount of memory.  If
  ///        there are more systems, the program will expand each and attempt minimizations as
  ///        this limit (a safeguard against overrunning available resources) permits.
  int getRunningStateCount() const;

  /// \brief Get the number of final states to produce for each initial system.
  int getFinalStateCount() const;

  /// \brief Get the number of samples to apply to each explicitly sampled rotatable bond.
  int getRotationSampleCount() const;

  /// \brief Get the maximum number of rotatable bonds to sample.
  int getRotatableBondLimit() const;

  /// \brief Get the maximum number of conformer seeding attempts
  int getMaxSeedingAttempts() const;

  /// \brief Get the minimum ratio for self interactions below which they will be considered
  ///        van-der Waals clashes.
  double getSelfClashRatio() const;

  /// \brief Get the maximum number of minimizations to attempt with any one molecule.  Each
  ///        initial state provided by the user will be subject to this limit, so if the limit
  ///        is 5000 and one molecule has two initial states listed in the input deck, the total
  ///        number of conformations sampled will be no greater than 10000.
  int getSystemTrialCount() const;
  
  /// \brief Get the positional root mean squared deviation that will distinguish each reported
  ///        confomer.
  double getRMSDTolerance() const;

  /// \brief Validate the restraining potential that will define the common core.
  void validateCoreRestraint() const;
  
  /// \brief Validate input pertaining to chiral sampling .
  ///
  /// \param directive  The keyword setting for chirality sampling (must be 'true' or 'false',
  ///                   without case sensitivity)
  void validateSampleChirality(const std::string &directive) const;

  /// \brief Validate input pertaining to sampling cis- and trans- states of molecules.
  ///
  /// \param directive  The keyword setting for cis- and trans- sampling (must be 'true' or
  ///                   'false', without case sensitivity)
  void validateSampleCisTrans(const std::string &directive) const;

  /// \brief Validate input pertaining to hydrogen bonding prevention.
  ///
  /// \param directive  The keyword setting for cis- and trans- sampling (must be 'true' or
  ///                   'false', without case sensitivity)
  void validatePreventHBonds(const std::string &directive) const;  
  
  /// \brief Validate the replica counts and criteria for distinguishing unique conformers.
  void validateStateCounts();
  
private:
  ExceptionResponse policy;         ///< Set the behavior when bad inputs are encountered.  DIE =
                                    ///<   abort program, WARN = warn the user, and likely reset to
                                    ///<   the default value if one is available, SILENT = do not
                                    ///<   warn the user, but also likely reset to the default
                                    ///<   value if one is available.
  std::string core_atom_mask;       ///< Mask of common core atoms, applied to all systems
  std::string core_data_item;       ///< Name of a data item from a Biovia SD file indicating the
                                    ///<   core atom mask.  Various formats of the input in the
                                    ///<   SD file will be automatically detected and accepted.
  double core_rk2;                  ///< Core atom repulsive positional restraint stiffness
  double core_rk3;                  ///< Core atom attractive positional restraint stiffness
  double core_r2;                   ///< Core atom maximum repulsive potential range
  double core_r3;                   ///< Range for onset of core atom attractive potential
  std::string anchor_conformation;  ///< The other half of the common core mask, the thing to align
                                    ///<   core atoms against.  This anchor conformation will be
                                    ///<   docked or placed within the receptor of interest and
                                    ///<   must contain atoms that register in the common atom
                                    ///<   mask.  All conformations will be aligned against this
                                    ///<   reference in order to bias the search towards results
                                    ///<   that fit inside of the receptor.
  bool sample_chirality;            ///< Flag to have chiral enantiomers sampled
  bool sample_cis_trans;            ///< Flag to have cis-trans isomers sampled
  bool prevent_hbonds;              ///< Flag to apply restraints that will prevent hydrogen bonds
                                    ///<   from forming during energy minimizations
  int running_states;               ///< Number of states to try minimizing at one time
  int final_states;                 ///< Number of final states to collect
  int rotation_samples;             ///< Number of times to sample about a rotatable bond
  int rotatable_bond_limit;         ///< Maximum number of rotatable bonds to explicitly sample
  int max_seeding_attempts;         ///< Maximum number of attempts to make in seeding each
                                    ///<   conformer.  If the seeding fails, the conformer will be
                                    ///<   left in its original state from the user input.
  double self_clash_ratio;          ///< The minimum ratio of the distance between two particles to
                                    ///<   their pair sigma value.  Ratios below this level will be
                                    ///<   considered clashes.
  int system_trials;                ///< Maximum number of distinct minimizations to attempt with
                                    ///<   one molecule
  double rmsd_tolerance;            ///< Minimum mass-weighted root-mean squared deviation between
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
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator conformerInput(const TextFile &tf, int *start_line, bool *found,
                                ExceptionResponse policy = ExceptionResponse::DIE,
                                WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace namelist
} // namespace stormm

#endif
