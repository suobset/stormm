#include "Constants/scaling.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "namelist_element.h"
#include "nml_minimize.h"

namespace omni {
namespace namelist {

using constants::CaseSensitivity;
using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
MinimizeControls::MinimizeControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    total_cycles{default_minimize_maxcyc},
    steepest_descent_cycles{default_minimize_maxcyc},
    initial_step{default_minimize_dx0},
    convergence_target{default_minimize_drms}
{}

//-------------------------------------------------------------------------------------------------
MinimizeControls::MinimizeControls(const TextFile &tf, int *start_line,
                                   const ExceptionResponse policy_in) :
    MinimizeControls(policy_in)
{
  NamelistEmulator t_nml = minimizeInput(tf, start_line, policy);
  total_cycles = t_nml.getIntValue("maxcyc");
  steepest_descent_cycles = t_nml.getIntValue("ncyc");
  initial_step = t_nml.getRealValue("dx0");
  convergence_target = t_nml.getRealValue("drms");

  // Validate input
  validateTotalCycles();
  validateSteepestDescentCycles();
  validateInitialStep();
  validateConvergenceTarget();
}

//-------------------------------------------------------------------------------------------------
int MinimizeControls::getTotalCycles() const {
  return total_cycles;
}

//-------------------------------------------------------------------------------------------------
int MinimizeControls::getSteepestDescentCycles() const {
  return steepest_descent_cycles;
}

//-------------------------------------------------------------------------------------------------
double MinimizeControls::getInitialStep() const {
  return initial_step;
}

//-------------------------------------------------------------------------------------------------
double MinimizeControls::getConvergenceTarget() const {
  return convergence_target;
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setTotalCycles(const int cycles_in) {
  total_cycles = cycles_in;
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setSteepestDescentCycles(const int cycles_in) {
  steepest_descent_cycles = cycles_in;
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setInitialStep(const double step_size_in) {
  initial_step = step_size_in;
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setConvergenceTarget(const double target_in) {
  convergence_target = target_in;
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateTotalCycles() {
  if (total_cycles < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The total cycle count of " + std::to_string(total_cycles) + " is invalid.",
            "MinimizeControls", "validateTotalCycles");
    case ExceptionResponse::WARN:
      rtWarn("The total cycle count of " + std::to_string(total_cycles) + " is invalid and will "
             "be reset to the default value of " + std::to_string(default_minimize_maxcyc) + ".",
             "MinimizeControls", "validateTotalCycles");      
      total_cycles = default_minimize_maxcyc;
      break;
    case ExceptionResponse::SILENT:
      total_cycles = default_minimize_maxcyc;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateSteepestDescentCycles() {
  if (total_cycles < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The cycle count of " + std::to_string(steepest_descent_cycles) + " is invalid for "
            "the steepest descent phase of minimization.", "MinimizeControls",
            "validateTotalCycles");
    case ExceptionResponse::WARN:
      rtWarn("The cycle count of " + std::to_string(steepest_descent_cycles) + " is invalid for "
             "the steepest descent phase and will be reset to the default value of " +
             std::to_string(default_minimize_ncyc) + ".", "MinimizeControls",
             "validateTotalCycles");      
      steepest_descent_cycles = default_minimize_ncyc;
      break;
    case ExceptionResponse::SILENT:
      steepest_descent_cycles = default_minimize_ncyc;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateInitialStep() {
  if (convergence_target < constants::verytiny || convergence_target) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The initial step size must be nontrivial in order for the algorithm to begin moving "
            "particles.  " + realToString(initial_step, 4, 11, NumberFormat::SCIENTIFIC) +
            " is not a valid initial step.", "MinimizeControls", "validateConvergenceTarget");
    case ExceptionResponse::WARN:
      rtWarn("The initial step size must be nontrivial in order for the algorithm to begin moving "
             "particles.  " + realToString(initial_step, 4, 11, NumberFormat::SCIENTIFIC) +
             " is not a valid initial step and will be replaced by the default of " +
             realToString(default_minimize_dx0, 4, 11, NumberFormat::SCIENTIFIC) + ".",
             "MinimizeControls", "validateConvergenceTarget");
      initial_step = default_minimize_dx0;
      break;
    case ExceptionResponse::SILENT:
      initial_step = default_minimize_dx0;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateConvergenceTarget() {
  if (convergence_target < constants::tiny) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A convergence target of " +
            realToString(convergence_target, 4, 11, NumberFormat::SCIENTIFIC) +
            " is not realistic.", "MinimizeControls", "validateConvergenceTarget");
    case ExceptionResponse::WARN:
      rtErr("A convergence target of " +
            realToString(convergence_target, 4, 11, NumberFormat::SCIENTIFIC) +
            " is not realistic and will be replaced by the default value of " +
            realToString(default_minimize_drms, 4, 11, NumberFormat::SCIENTIFIC) + ".",
            "MinimizeControls", "validateConvergenceTarget");
      convergence_target = default_minimize_drms;
      break;
    case ExceptionResponse::SILENT:
      convergence_target = default_minimize_drms;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line,
                               const ExceptionResponse policy) {
  NamelistEmulator t_nml("minimize", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to carry out energy minimization of a molecular system guided by an "
                         "energy surface based on a force field.");
  t_nml.addKeyword(NamelistElement("maxcyc", NamelistType::INTEGER,
                                   std::to_string(default_minimize_maxcyc)));
  t_nml.addKeyword(NamelistElement("ncyc", NamelistType::INTEGER,
                                   std::to_string(default_minimize_ncyc)));
  t_nml.addKeyword(NamelistElement("dx0", NamelistType::REAL,
                                   std::to_string(default_minimize_dx0)));
  t_nml.addKeyword(NamelistElement("drms", NamelistType::REAL,
                                   std::to_string(default_minimize_drms)));
  t_nml.addHelp("maxcyc", "Maximum number of line minimization cycles to pursue.");
  t_nml.addHelp("ncyc", "Number of steepest-descent optimization steps to perform, prior to the "
                "onset of conjugate gradient optimization for the remaining maxcyc - ncyc steps.");
  t_nml.addHelp("dx0", "Magnitude of the initial displacement along the gradient vector.  The "
                "size of subsequent moves will grow or shrink based on the history of success in "
                "previous optimizations.");
  t_nml.addHelp("drms", "Convergence criterion for the minimization, based on the root mean "
                "squared value of the Cartesian forces on all particles.  Units of kcal/mol-A.");  
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.  All calls to this function should
  // proceed in consecutive calls, to make use of the updates to start_line and avoid reading any
  // instance of this namelist twice or skipping instances of it in the search for some other
  // namelist.  An alternative is to keep an independent counter to track progress through the
  // input file in search for &rst namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::NO, tf.getLineCount());

  return t_nml;
}

} // namespace namelist
} // namespace omni
