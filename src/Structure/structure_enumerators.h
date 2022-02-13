// -*-c++-*-
#ifndef OMNI_STRUCTURE_ENUMERATORS_H
#define OMNI_STRUCTURE_ENUMERATORS_H

namespace omni {
namespace structure {
  
/// \brief There are two typical ways of re-imaging particles.
enum class ImagingMethod {
  PRIMARY_UNIT_CELL,  ///< Place coordinates such that, in fractional space, they lie in the first
                      ///<   octant, all components within [0, 1).
  MINIMUM_IMAGE,      ///< Place coordinates (more typically, displacements) such that, in
                      ///<   fractional space, all components fall in the range [ -0.5, 0.5 ).
};

} // namespace structure
} // namespace omni

#endif
