// -*-c++-*-
#ifndef EMULATOR_NAMELIST_H
#define EMULATOR_NAMELIST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/input.h"
#include "Namelists/namelist_emulator.h"
#include "Namelists/namelist_element.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "Potential/nbemulator_util.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using energy::EmulationTarget;
using energy::ForceZero;
using energy::NBAtomSource;
using energy::StateVariable;
using parse::TextFile;
using parse::WrapTextSearch;

/// \brief Class to collect user inputs taken for fitting non-bonded pair potentials to reproduce
///        energies for a specific set of structures or complexes
class EmulatorControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  EmulatorControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  EmulatorControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  EmulatorControls(const EmulatorControls &original) = default;
  EmulatorControls(EmulatorControls &&original) = default;
  EmulatorControls& operator=(const EmulatorControls &original) = default;
  EmulatorControls& operator=(EmulatorControls &&original) = default;

  /// \brief Get the number of unique source definitions for pair potentials.
  int getSourceCount() const;

  /// \brief Get the definition of one of the atom sources.  The order in which sources are
  ///        dispensed from this object may be different than the order in which sources are listed
  ///        in the user input, as internal re-arrangement occurs to ensure that more specific
  ///        source definitions are placed ahead of more general, matter-encompassing ones.
  const NBAtomSource& getSource(int index) const;

  /// \brief Get the number of fitting targets.
  int getTargetCount() const;

  /// \brief Get a specific fitting target from the stored list.
  ///
  /// \param t_index  Index of the target of interest.  This will be checked against the size of
  ///                 the list, for validity.
  const EmulationTarget& getTarget(int t_index) const;

  /// \brief Get the number of atomic force balancing directives.
  int getBalanceCount() const;
  
  /// \brief Get a specific atom force balancing directive.
  ///
  /// \param b_index  Index of the balancing instruction of interest.  This will be checked against
  ///                 the size of the list, for validity.
  const ForceZero& getBalance(int b_index) const;

  /// \brief Get the number of basis functions.
  int getBasisFunctionCount() const;
  
  /// \brief Get the width of the support region for a given basis function.
  ///
  /// \param index  The index of the basis function of interest
  double getSupportWidth(int index) const;

  /// \brief Get the start of the support region for a given basis function.
  ///
  /// \param index  The index of the basis function of interest
  double getSupportStart(int index) const;

  /// \brief Get the regulation for various basis functions.
  ///
  /// Overloaded:
  ///   - Get the generic regulation used when no value was given for a particular basis function
  ///   - Get the specific regulation used for the indexed basis function of each pair potential
  ///
  /// \param index  The index of the basis function of interest
  /// \{
  double getRegularization() const;
  double getRegularization(int index) const;
  /// \}

  /// \brief Get the list of molecular mechanics terms which provide the context for the fitted
  ///        potentials.
  const std::vector<StateVariable>& getMMContext() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator getTranscript() const;

private:

  ExceptionResponse policy;  ///< The course of action to take in response to bad input 

  /// The list of atom source definitions.  This is pruned and ordered at the point of input so
  /// that stepping through the list will categorize atoms in the most specific terms first, then
  /// into broader categories if other classifications fail.
  std::vector<NBAtomSource> source_list;

  /// The list of targets for fitting will reference labels in a system cache.
  std::vector<EmulationTarget> target_list;

  /// The list of targets may alos include specific atoms for which forces should be balanced to
  /// the greatest degree possible.
  std::vector<ForceZero> balance_list;
  
  /// The list of molecular mechanics terms in which the emulation potentials are to function
  std::vector<StateVariable> mm_context;
  
  // The basis functions for each pair potential are identical in form, with different scaling
  // factors.
  double generic_zero_hold;            ///< General-purpose weight for regulating basis function
                                       ///<   scaling factors.  The overall weight will be scaled
                                       ///<   by the number of times any given basis function
                                       ///<   appears in the fitting problem.
  std::vector<double> support_widths;  ///< Width of the support for each basis function
  std::vector<double> support_starts;  ///< The start of the support range for each basis function
  std::vector<double> bf_zero_holds;   ///< Regulating weights for each basis function.  The
                                       ///<   overall weight will be scaled by the number of times
                                       ///<   the basis function appears in the fitting problem.
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};

/// \brief A free function to assemble a non-bonded potential atom source definition from
///        user-supplied data in a &namelist input struct
///
/// \param t_nml     The namelist bearing source structs.  The keyword for eahc atomic source must
///                  match the provided key_name.  Subkeys within the keyword must match an
///                  expected pattern.
/// \param src_idx   Index of the source to extract information from
/// \param key_name  Name of the keyword bearing atom source details
NBAtomSource extractNBAtomSourceInput(const NamelistEmulator &t_nml, const int src_idx,
                                      const std::string &key_name = std::string("source"));

/// \brief A free function to prepare the &emulator namelist, wrapping user commands that modulate
///        the workflow for generating non-bonded potentials that reproduce target energies when
///        applied to a set of known structures.
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator emulatorInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace emulation
} // namespace stormm

#endif
