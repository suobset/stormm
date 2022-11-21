// -*-c++-*-
#ifndef STORMM_RECEPTOR_GRID_H
#define STORMM_RECEPTOR_GRID_H

#include <string>
#include <vector>
#include "copyright.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

class ReceptorGrid {
public:

  /// \brief The constructor requires a topology and may take coordinates of the receptor as one
  ///        of the coordinate objects or a handful of arrays.
  
  
  
private:
  GridDetailLevel contents;
  int length_x;            ///< Number of grid points along the Cartesian X axis
  int length_y;            ///< Number of grid points along the Cartesian Y axis
  int length_z;            ///< Number of grid points along the Cartesian Z axis
  int max_explicit_atoms;  ///< Maximum number of atoms from the rigid components of the receptor
                           ///<   that will contribute explicit interactions to any one grid bin.
                           ///<   If this limit is exceeded, the furthest atoms' potentials will
                           ///<   be mapped to the grid corners to become part of the tricubic
                           ///<   spline function.
  double grid_spacing;     ///< Spacing between grid points in all directions
  double explicit_cutoff;  ///< Cutoff below which atoms near any point in a grid bin will be
                           ///<   represented explicitly.  This corresponds to the demarcation
                           ///<   described in the NONBONDED_ATOMIC case of the GridDetail
                           ///<   enumerator.
  
  const AtomGraph *ag_pointer;
};
 
} // namespace structure
} // namespace stormm

#endif
