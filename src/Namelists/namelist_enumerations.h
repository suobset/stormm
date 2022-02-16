// -*-c++-*-
#ifndef OMNI_NAMELIST_ENUMERATIONS_H
#define OMNI_NAMELIST_ENUMERATIONS_H

namespace omni {
namespace namelist {

/// \brief Enumeration to indicate whether a restraint encodes its atoms by atom masks (strings) or
///        by topological indices.
enum class RestraintAnchoring {
  ATOMMASK,  ///< Atom masks encode the connected atoms
  INDICES,   ///< Topological indices specify the restrained atoms
  MIXED,     ///< Both masks and topological indices specify the mask
  UNKNOWN    ///< Who knows?  This should never be possible to achieve.
};

} // namespace namelist
} // namespace omni

#endif
