// -*-c++-*-
#ifndef OMNI_NML_MINIMIZE_H
#define OMNI_NML_MINIMIZE_H

#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

/// \brief Default values for energy minimization
/// \{
constexpr int default_minimize_maxcyc        = 200;
constexpr int default_minimize_ncyc          = 50;
constexpr int default_minimize_ntpr          = 50;
constexpr char default_minimize_checkpoint[] = "true";
constexpr double default_minimize_cut        = 8.0;
constexpr double default_minimize_dx0        = 0.01;
constexpr double default_minimize_drms       = 0.0001;
/// \}
  
/// \brief Object to encapsulate energy minimization control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using OMNI libraries.
class MinimizeControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &minimize namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \{
  MinimizeControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  MinimizeControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

  /// \brief Get the total number of minimization cycles.
  int getTotalCycles() const;

  /// \brief Get the number of steepest descent cycles.
  int getSteepestDescentCycles() const;

  /// \brief Get the diagnostic output printing frequency, akin to the major contribution to pmemd
  ///        and sander mdout files.
  int getDiagnosticPrintFrequency() const;

  /// \brief Get the directive on whether to produce a checkpoint file for the final state of each
  ///        energy minimization run.
  bool getCheckpointProduction() const;
  
  /// \brief Get the electrostatic cutoff.
  double getElectrostaticCutoff() const;
  
  /// \brief Get the Lennard-Jones cutoff.
  double getLennardJonesCutoff() const;
  
  /// \brief Get the initial step length.
  double getInitialStep() const;

  /// \brief Get the convergence criterion.
  double getConvergenceTarget() const;

  /// \brief Set the total number of minimization cycles.
  ///
  /// \param cycles_in  The requested number of minimization cycles
  void setTotalCycles(int cycles_in);
  
  /// \brief Set the number of steepest descent cycles.
  ///
  /// \param cycles_in  The requested number of steepest descent cycles
  void setSteepestDescentCycles(int cycles_in);

  /// \brief Set the diagnostic printing frequency.
  ///
  /// \param frequency_in  The chosen printing interval
  void setDiagnosticPrintFrequency(int frequency_in);

  /// \brief Set the checkpoint production flag.
  ///
  /// \param produce_in  Whether to produce a checkpoint at the end of the run
  void setCheckpointProduction(bool produce_in);
  
  /// \brief Set the electrostatic cutoff
  void setElectrostaticCutoff(double cutoff_in);
  
  /// \brief Set the Lennard-Jones cutoff
  void setLennardJonesCutoff(double cutoff_in);
  
  /// \brief Set the initial step length.
  ///
  /// \param step_size_in  The requested initial step length
  void setInitialStep(double step_size_in);

  /// \brief Set the convergence criterion, the target for the root mean squared value of all
  ///        gradients obtained after the minimization, in kcal/mol.
  ///
  /// \param target_in  The requested convergence target
  void setConvergenceTarget(double target_in);
  
private:
  ExceptionResponse policy;     ///< Set the behavior when bad inputs are encountered.  DIE =
                                ///<   abort program, WARN = warn the user, and likely reset to
                                ///<   the default value if one is available, SILENT = do not
                                ///<   warn the user, but also likely reset to the default value
                                ///<   if one is available.
  int total_cycles;             ///< Maximum number of minimization steps to attempt (equivalent
                                ///<   to maxcyc in sander)
  int steepest_descent_cycles;  ///< Number of steepest descent steps to perform prior to beginning
                                ///<   conjugate gradient moves (equivalent to ncyc in sander)
  int print_frequency;          ///< Print results at step 0 and, thereafter, after each interval
                                ///<   of this many line minimizations.  The default of 0
                                ///<   suppresses output except at the outset of the run.
  bool produce_checkpoint;      ///< Indicate that a checkpoint file should be produced at the end
                                ///<   of the energy minimization run (default TRUE), with the name
                                ///<   of the checkpoint file for each system found in the &files
                                ///<   namelist.
  double electrostatic_cutoff;  ///< Cutoff for (short-ranged) electrostatic interactions, or for
                                ///<   all gas-phase Coulombic electrostatics in a non-periodic
                                ///<   system.  Units of Angstroms (A).
  double lennard_jones_cutoff;  ///< Cutoff for van-der Waals interactions, in Angstroms (A).
  double initial_step;          ///< Magnitude of the initial displacement along the gradient
                                ///<   vector.  The size of subsequent moves will grow or shrink
                                ///<   based on the history of success in previous optimizations.
                                ///<   Units of Angstroms (A).
  double convergence_target;    ///< Convergence target for root mean squared value of all
                                ///<   gradients obtained after the minimization, in kcal/mol-A.
  
  /// \brief Validate the total number of minimization cycles
  void validateTotalCycles();

  /// \brief Validate the number of steepest descent cycles.  This does NOT ensure that this number
  ///        is less than the number of total cycles.  That behavior is enforced by only taking up
  ///        to the total number of steps, and taking a steepest descent approach if the step
  ///        number is less than the appropriate threshold.
  void validateSteepestDescentCycles();

  /// \brief Validate the diagnostic printing frequency
  void validatePrintFrequency();

  /// \brief Validate the checkpointing behavior.
  ///
  /// \param directive  The command given for whether to write a checkpoint file
  void validateCheckpointProduction(const std::string &directive) const;
  
  /// \brief Validate the electrostatic cutoff.
  void validateElectrostaticCutoff();

  /// \brief Validate the Lennard-Jones cutoff.
  void validateLennardJonesCutoff();
  
  /// \brief Validate the initial step size.
  void validateInitialStep();

  /// \brief Validate the convergence target value.
  void validateConvergenceTarget();
};
  
/// \brief Produce a namelist for specifying energy minimization directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the
///        molecular dynamics input, obviating the need for the imin setting found in the general
///        &cntrl namelist of sander and pmemd.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if needed,
///                    to find a &minimize namelist) 
/// \param found       Indicate that the namelist was found
/// \param policy      Reaction to exceptions encountered during namelist reading
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace namelist
} // namespace omni

#endif
