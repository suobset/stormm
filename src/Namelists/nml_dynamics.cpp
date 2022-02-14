#include "namelist_element.h"
#include "nml_dynamics.h"

namespace omni {
namespace namelist {

//-------------------------------------------------------------------------------------------------
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("dynamics", CaseSensitivity::AUTOMATIC, ExceptionResponse::DIE,
                         "Wraps directives needed to propagate dynamics of a molecular system.");
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
