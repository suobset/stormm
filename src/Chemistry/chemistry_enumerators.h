// -*-c++-*-
#ifndef OMNI_CHEMISTRY_ENUMERATORS_H
#define OMNI_CHEMISTRY_ENUMERATORS_H

namespace omni {
namespace chemistry {

/// \brief Enumerate the chiral orientations of a center with four unique bonded groups
enum class ChiralOrientation {
  RECTUS,    ///< R- or D-chiral
  SINISTER,  ///< S- or L-chiral
  NONE       ///< No chirality, or no preference
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
} // namespace omni

#endif

