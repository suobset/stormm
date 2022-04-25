// -*-c++-*-
#ifndef OMNI_NML_DYNAMICS_H
#define OMNI_NML_DYNAMICS_H

#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

/// \brief Default values for molecular dynamics
/// \{
constexpr int default_dynamics_nstlim = 100;
/// \}

/// \brief Object to encapsulate molecular dynamics control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using OMNI libraries.
class DynamicsControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indicator that the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \{
  DynamicsControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

  /// \brief Get the total number of dynamics steps
  int getStepCount() const;

  /// \brief Set the total number of dynamics steps
  ///
  /// \param nstlim_in  The number of steps to take
  void setStepCount(int nstlim_in);

private:
  ExceptionResponse policy;     ///< Set the behavior when bad inputs are encountered.  DIE =
                                ///<   abort program, WARN = warn the user, and likely reset to
                                ///<   the default value if one is available, SILENT = do not
                                ///<   warn the user, but also likely reset to the default value
                                ///<   if one is available.
  int nstlim;                   ///< Total number of dynamics steps to perform (equivalent to
                                ///<   maxcyc in sander)

  /// \brief Validate the total number of steps
  void validateStepCount() const;
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
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace namelist
} // namespace omni

#endif
