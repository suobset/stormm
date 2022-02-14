#include "namelist_element.h"
#include "nml_minimize.h"

namespace omni {
namespace namelist {

//-------------------------------------------------------------------------------------------------
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("minimize", CaseSensitivity::AUTOMATIC, ExceptionResponse::DIE,
                         "Wraps directives needed to carry out energy minimization of a molecular "
                         "system guided by an energy surface based on a force field.");
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
