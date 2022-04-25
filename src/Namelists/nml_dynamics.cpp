#include "Reporting/error_format.h"
#include "namelist_element.h"
#include "nml_dynamics.h"

namespace omni {
namespace namelist {

//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const ExceptionResponse policy_in) :
  policy{policy_in},
  nstlim{default_dynamics_nstlim}
{}

//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in) :
    DynamicsControls(policy_in)
{
  NamelistEmulator t_nml = dynamicsInput(tf, start_line, found_nml, policy);
  nstlim = t_nml.getIntValue("nstlim");

  // Validate input
  validateStepCount();
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getStepCount() const {
  return nstlim;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setStepCount(const int nstlim_in) {
  nstlim = nstlim_in;
  validateStepCount();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateStepCount() const {
  if (nstlim < 0) {
    rtErr("A negative value for the number of dynamics steps is invalid.  This error may be the "
          "result of trying to supply too large a number of steps (greater than 2.1 billion, "
          "2^31, which overflows the signed integer format.  Use checkpoint files to carry out "
          "runs with very large numbers of total steps in multiple segments.", "DynamicsControls",
          "validateStepCount");
  }
}
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy) {
  NamelistEmulator t_nml("dynamics", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to propagate dynamics of a molecular system.");
  t_nml.addKeyword(NamelistElement("nstlim", NamelistType::INTEGER,
                                   std::to_string(default_dynamics_nstlim)));
  t_nml.addHelp("nstlim", "Number of dynamics steps to carry out");

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
