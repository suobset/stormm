// -*-c++-*-
#include "copyright.h"
#include <string>

namespace stormm {
namespace restraints {

/// \brief An enumerator for the various restraint assortments that one can apply to a system.
enum class RestraintEnsemble {
  PREVENT_HBONDS,            ///< Penalize hydrogen bonds that may form between identified donor
                             ///<   and acceptor pairs.
  PRESERVE_HEAVY_DIHEDRALS,  ///< Apply dihedral restraints to dihedrals involving purely heavy
                             ///<   atoms.
  PRESERVE_DISTANCES         ///< Select heavy atoms throughout the molecule and apply distance
                             ///<   restraints to maintain the relative displacements.  This is
                             ///<   a distinct method of preserving molecular geometry, but works
                             ///<   towards the same goal as PRESERVE_HEAVY_DIHEDRALS
};

/// \brief Translate a user input value or other string code into the appropriate RestraintEnsemble
///        enumeration.
RestraintEnsemble translateRestraintEnsemble(const std::string &rst_group);

} // namespace restraints
} // namespace stormm
