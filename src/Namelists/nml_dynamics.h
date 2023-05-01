// -*-c++-*-
#ifndef STORMM_NML_DYNAMICS_H
#define STORMM_NML_DYNAMICS_H

#include "copyright.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using parse::WrapTextSearch;
  
/// \brief Default number of molecular dynamics cycles
constexpr int default_dynamics_nstlim = 100;

/// \brief Default time step for molecular dynamics, in femtoseconds
constexpr double default_dynamics_time_step = 1.0;

/// \brief Default tolerance for RATTLE bond constraints
constexpr double default_rattle_tolerance = 1.0e-6;

/// \brief The minimum molecular dynamics time step, in units of femtoseconds
constexpr double minimum_dynamics_time_step = 0.015625;

/// \brief The tightest possible RATTLE tolerance
constexpr double minimum_rattle_tolerance = 1.0e-9;

/// \brief Object to encapsulate molecular dynamics control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class DynamicsControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indicator that the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &dynamics namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  DynamicsControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  DynamicsControls(const DynamicsControls &original) = default;
  DynamicsControls(DynamicsControls &&original) = default;
  DynamicsControls& operator=(const DynamicsControls &original) = default;
  DynamicsControls& operator=(DynamicsControls &&original) = default;
  /// \}
  
  /// \brief Get the total number of dynamics steps
  int getStepCount() const;

  /// \brief Get the simulation time step
  double getTimeStep() const;

  /// \brief Get the RATTLE tolerance
  double getRattleTolerance() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the total number of dynamics steps
  ///
  /// \param nstlim_in  The number of steps to take
  void setStepCount(int nstlim_in);

  /// \brief Set the simulation time step
  ///
  /// \param time_step_in  The requested time step
  void setTimeStep(double time_step_in);

  /// \brief Set the simulation RATTLE tolerance
  ///
  /// \param tol_in  The requested tolerance
  void setRattleTolerance(double rattle_tolerance_in);
  
private:
  ExceptionResponse policy;     ///< Set the behavior when bad inputs are encountered.  DIE =
                                ///<   abort program, WARN = warn the user, and likely reset to
                                ///<   the default value if one is available, SILENT = do not
                                ///<   warn the user, but also likely reset to the default value
                                ///<   if one is available.
  int nstlim;                   ///< Total number of dynamics steps to perform (equivalent to
                                ///<   maxcyc in sander)
  double time_step;             ///< Time step to take after each force evaluation
  double rattle_tolerance;      ///< The tolerance to apply to bond constraint calculations

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the total number of steps
  void validateStepCount() const;

  /// \brief Validate the time step
  void validateTimeStep() const;

  /// \brief Validate the RATTLE tolerance
  void validateRattleTolerance() const;
};
  
/// \brief Produce a namelist for specifying molecular dynamics directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the energy
///        minimization input.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
/// \param found       Indicator that the namelist was present in the input file
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for an &dynamics namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
