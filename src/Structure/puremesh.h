// -*-c++-*-
#ifndef STORMM_PUREMESH_H
#define STORMM_PUREMESH_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
namespace structure {

/// \brief A workspace for constructing a pure potential mesh based on the frozen atoms of a
///        large molecule.  If the large molecule has nonrigid components, they will be excluded
///        from contributing to the grid.  In addition, any atoms up to 1:3 (connected by two bonds
///        of less) will also be excluded from the grid-based potential.  Computations on these
///        atoms will not be accurate off the grid, but since they are frozen the consequences are
///        mitigated.  Partially attenuated 1:4 interactions between the frozen atoms and mobile
///        atoms will be corrected post-hoc while computing other potentials.

  
} // namespace structure
} // namespace stormm

#endif
