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
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
MinimizeControls::MinimizeControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    total_cycles{default_minimize_maxcyc},
    steepest_descent_cycles{default_minimize_ncyc},
    print_frequency{default_minimize_ntpr},
    produce_checkpoint{true},
    electrostatic_cutoff{default_minimize_cut},
    lennard_jones_cutoff{default_minimize_cut},
    initial_step{default_minimize_dx0},
    convergence_target{default_minimize_drms}
{}

//-------------------------------------------------------------------------------------------------
MinimizeControls::MinimizeControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in) :
    MinimizeControls(policy_in)
{
  NamelistEmulator t_nml = minimizeInput(tf, start_line, found_nml, policy);
  total_cycles = t_nml.getIntValue("maxcyc");
  steepest_descent_cycles = t_nml.getIntValue("ncyc");
  print_frequency = t_nml.getIntValue("ntpr");
  const std::string chkp_behavior = t_nml.getStringValue("checkpoint");
  validateCheckpointProduction(chkp_behavior);
  produce_checkpoint = strcmpCased(chkp_behavior, "true");
  if (t_nml.getKeywordStatus("cut") != InputStatus::MISSING) {
    electrostatic_cutoff = t_nml.getRealValue("cut");
    lennard_jones_cutoff = t_nml.getRealValue("cut");
  }
  else {
    electrostatic_cutoff = t_nml.getRealValue("es_cutoff");
    lennard_jones_cutoff = t_nml.getRealValue("vdw_cutoff");
  }
  initial_step = t_nml.getRealValue("dx0");
  convergence_target = t_nml.getRealValue("drms");

  // Validate input
  validateTotalCycles();
  validateSteepestDescentCycles();
  validatePrintFrequency();
  validateElectrostaticCutoff();
  validateLennardJonesCutoff();
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
int MinimizeControls::getDiagnosticPrintFrequency() const {
  return print_frequency;
}

//-------------------------------------------------------------------------------------------------
bool MinimizeControls::getCheckpointProduction() const {
  return produce_checkpoint;
}
  
//-------------------------------------------------------------------------------------------------
double MinimizeControls::getElectrostaticCutoff() const {
  return electrostatic_cutoff;
}

//-------------------------------------------------------------------------------------------------
double MinimizeControls::getLennardJonesCutoff() const {
  return lennard_jones_cutoff;
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
  validateTotalCycles();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setSteepestDescentCycles(const int cycles_in) {
  steepest_descent_cycles = cycles_in;
  validateSteepestDescentCycles();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setDiagnosticPrintFrequency(const int frequency_in) {
  print_frequency = frequency_in;
  validatePrintFrequency();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setCheckpointProduction(const bool produce_in) {
  produce_checkpoint = produce_in;
}
  
//-------------------------------------------------------------------------------------------------
void MinimizeControls::setElectrostaticCutoff(const double cutoff_in) {
  electrostatic_cutoff = cutoff_in;
  validateElectrostaticCutoff();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setLennardJonesCutoff(const double cutoff_in) {
  lennard_jones_cutoff = cutoff_in;
  validateLennardJonesCutoff();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setInitialStep(const double step_size_in) {
  initial_step = step_size_in;
  validateInitialStep();
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::setConvergenceTarget(const double target_in) {
  convergence_target = target_in;
  validateConvergenceTarget();
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
  if (steepest_descent_cycles < 0) {
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
      steepest_descent_cycles = std::min(default_minimize_ncyc, total_cycles);
      break;
    case ExceptionResponse::SILENT:
      steepest_descent_cycles = std::min(default_minimize_ncyc, total_cycles);
      break;
    }
  }
  if (steepest_descent_cycles > total_cycles) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The cycle count of " + std::to_string(steepest_descent_cycles) + " is greater than "
            "the stated number of total cycles " + std::to_string(total_cycles) + ".",
            "MinimizeControls", "validateTotalCycles");
    case ExceptionResponse::WARN:
      rtErr("The cycle count of " + std::to_string(steepest_descent_cycles) + " is greater than "
            "the stated number of total cycles " + std::to_string(total_cycles) + " and will be "
            "reduced to the maximum overall cycle count.", "MinimizeControls",
            "validateTotalCycles");
      steepest_descent_cycles = total_cycles;
      break;
    case ExceptionResponse::SILENT:
      steepest_descent_cycles = total_cycles;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validatePrintFrequency() {
  if (print_frequency < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Diagnostics must be printed with a non-negative frequency, not " +
             std::to_string(print_frequency) + ".  The default of " +
             std::to_string(default_minimize_ntpr) + " will be restored, which suppresses most "
             "output.", "MinimizeControls", "validatePrintFrequency");
      print_frequency = default_minimize_ntpr;
      break;
    case ExceptionResponse::SILENT:
      print_frequency = default_minimize_ntpr;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateCheckpointProduction(const std::string &directive) const {
  if (strcmpCased(directive, "true", CaseSensitivity::NO) ||
      strcmpCased(directive, "false", CaseSensitivity::NO)) {
    return;
  }
  rtErr("Acceptable values for the checkpoint string are 'true' or 'false' (case-insensitive).  "
        "The value of " + directive + " is invalid.", "MinimizeControls",
        "validateCheckpointProduction");
}
    
//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateElectrostaticCutoff() {
  if (electrostatic_cutoff < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("No electrostatic interaction will be computed with a cutoff of " +
             std::to_string(electrostatic_cutoff) + ".", "MinimizeControls",
             "validateElectrostaticCutoff");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateLennardJonesCutoff() {
  if (electrostatic_cutoff < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("No van-der Waals interaction will be computed with a cutoff of " +
             std::to_string(lennard_jones_cutoff) + ".", "MinimizeControls",
             "validateLennardJonesCutoff");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MinimizeControls::validateInitialStep() {
  if (initial_step < constants::verytiny) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The initial step size must be nontrivial in order for the algorithm to begin moving "
            "particles.  " + realToString(initial_step, 11, 4, NumberFormat::SCIENTIFIC) +
            " is not a valid initial step.", "MinimizeControls", "validateInitialStep");
    case ExceptionResponse::WARN:
      rtWarn("The initial step size must be nontrivial in order for the algorithm to begin moving "
             "particles.  " + realToString(initial_step, 11, 4, NumberFormat::SCIENTIFIC) +
             " is not a valid initial step and will be replaced by the default of " +
             realToString(default_minimize_dx0, 11, 4, NumberFormat::SCIENTIFIC) + ".",
             "MinimizeControls", "validateInitialStep");
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
            realToString(convergence_target, 11, 4, NumberFormat::SCIENTIFIC) +
            " is not realistic.", "MinimizeControls", "validateConvergenceTarget");
    case ExceptionResponse::WARN:
      rtErr("A convergence target of " +
            realToString(convergence_target, 11, 4, NumberFormat::SCIENTIFIC) +
            " is not realistic and will be replaced by the default value of " +
            realToString(default_minimize_drms, 11, 4, NumberFormat::SCIENTIFIC) + ".",
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
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy) {
  NamelistEmulator t_nml("minimize", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to carry out energy minimization of a molecular system guided by an "
                         "energy surface based on a force field.");
  t_nml.addKeyword(NamelistElement("maxcyc", NamelistType::INTEGER,
                                   std::to_string(default_minimize_maxcyc)));
  t_nml.addKeyword(NamelistElement("ncyc", NamelistType::INTEGER,
                                   std::to_string(default_minimize_ncyc)));
  t_nml.addKeyword(NamelistElement("ntpr", NamelistType::INTEGER,
                                   std::to_string(default_minimize_ntpr)));
  t_nml.addKeyword(NamelistElement("checkpoint", NamelistType::STRING,
                                   std::string(default_minimize_checkpoint)));
  t_nml.addKeyword(NamelistElement("cut", NamelistType::REAL, std::string("MISSING")));
  t_nml.addKeyword(NamelistElement("es_cutoff", NamelistType::REAL,
                                   std::to_string(default_minimize_cut)));
  t_nml.addKeyword(NamelistElement("vdw_cutoff", NamelistType::REAL,
                                   std::to_string(default_minimize_cut)));
  t_nml.addKeyword(NamelistElement("dx0", NamelistType::REAL,
                                   std::to_string(default_minimize_dx0)));
  t_nml.addKeyword(NamelistElement("drms", NamelistType::REAL,
                                   std::to_string(default_minimize_drms)));
  t_nml.addHelp("maxcyc", "Maximum number of line minimization cycles to pursue.");
  t_nml.addHelp("ncyc", "Number of steepest-descent optimization steps to perform, prior to the "
                "onset of conjugate gradient optimization for the remaining maxcyc - ncyc steps.");
  t_nml.addHelp("ntpr", "Interval at which to report energy diagnostics for the minimization run, "
                "akin to the mdout results in Amber's sander and pmemd programs.  The default "
                "of " + std::to_string(default_minimize_ntpr) + " suppresses output except at the "
                "outset of the run.");
  t_nml.addHelp("cut", "Cutoff to apply to all short-ranged interactions.");
  t_nml.addHelp("es_cutoff", "Cutoff to apply to electrostatic (short-ranged) interactions.");
  t_nml.addHelp("vdw_cutoff", "Cutoff to apply to Lennard-Jones interactions.");
  t_nml.addHelp("dx0", "Magnitude of the initial displacement along the gradient vector.  The "
                "size of subsequent moves will grow or shrink based on the history of success in "
                "previous optimizations.");
  t_nml.addHelp("drms", "Convergence criterion for the minimization, based on the root mean "
                "squared value of the Cartesian forces on all particles.  Units of kcal/mol-A.");  
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount(),
                             found);

  return t_nml;
}

} // namespace namelist
} // namespace omni
