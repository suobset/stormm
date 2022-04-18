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

/// \brief Enumeration to indicate the requirement of STRUCT subkeys
enum class SubkeyRequirement {
  OPTIONAL,  ///< The keyword is not required, but can be specified
  REQUIRED,  ///< The keyword is essential and must be supplied in order to create a valid STRUCT
  BOGUS      ///< The option should actually NOT be present--useful for situations in which a
             ///<   pre-compiled namelist is re-used but not intended to take arguments from
             ///<   certain STRUCT options
};

} // namespace namelist
} // namespace omni

#endif
