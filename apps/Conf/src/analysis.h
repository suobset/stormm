// -*-c++-*-
#ifndef CONFORMER_ANALYSIS_H
#define CONFORMER_ANALYSIS_H

#include "../../../src/Structure/local_arrangement.h"
#include "../../../src/Structure/rmsd.h"

namespace conf_app {
namespace analysis {
    
/// \brief Compute the angles made by each rotatable bond.  This will compute the angle based on
///        the highest Z-number atoms attached to either end of the bond (not counting the atoms at
///        the other end of the bond).  If there is a tie, the atom with the lower topological
///        index will be preferred.


} // namespace analysis
} // namespace conf_app

#endif
