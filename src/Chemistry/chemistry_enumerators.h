// -*-c++-*-
#ifndef STORMM_CHEMISTRY_ENUMERATORS_H
#define STORMM_CHEMISTRY_ENUMERATORS_H

#include "copyright.h"

namespace stormm {
namespace chemistry {

/// \brief Enumerate the chiral orientations of a center with four unique bonded groups
enum class ChiralOrientation {
  RECTUS,    ///< R- or D-chiral
  SINISTER,  ///< S- or L-chiral
  NONE       ///< No chirality, or no preference
};

/// \brief Enumerate the ways for inverting a chiral center
enum class ChiralInversionProtocol {
  ROTATE,        ///< A simple rotation of two arms about the bisector between them will invert
                 ///<   the chiral center without severe bond stretching.
  REFLECT,       ///< A reflection of the entire molecule across a plane will be necessary to
                 ///<   invert the chiral center in a way that does not break covalent bonds.
                 ///<   Because this inverts all centers in the molecule, those that can be
                 ///<   inverted individually by rotation must be reset to their original
                 ///<   orientations following the reflection.
  DO_NOT_INVERT  ///< Applied to the 2nd, 3rd, and subsequent chiral centers that require
                 ///<   reflection across a plane in order to invert them.  All such centers must
                 ///<   be inverted together, so there is no point in sampling orientations beyond
                 ///<   those of the first.
};
  
/// \brief Enumerate the options for mapping rotatable groups in the molecule:
enum class MapRotatableGroups {
  YES, NO
};

/// \brief Enumerate a series of actions for retrieving rotatable group information from a list
///        held within a ChemicalFeatures object.  Whichever criterion is chosen, the results
///        will be ordered in descending order by the other criterion.
enum class RotatorCriterion {
  COM_PROXIMITY,  ///< Proximity to the center of mass
  GROUP_SIZE      ///< Size of the group of rotating atoms (excludes the pivot and root atoms at
                  ///<   the ends of the rotatable bond)
};
  
} // namespace chemistry
} // namespace stormm

#endif

