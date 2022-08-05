// -*-c++-*-
#ifndef STORMM_STRUCTURE_ENUMERATORS_H
#define STORMM_STRUCTURE_ENUMERATORS_H

namespace stormm {
namespace structure {
  
/// \brief There are two typical ways of re-imaging particles.
enum class ImagingMethod {
  PRIMARY_UNIT_CELL,  ///< Place coordinates such that, in fractional space, they lie in the first
                      ///<   octant, all components within [0, 1).
  MINIMUM_IMAGE,      ///< Place coordinates (more typically, displacements) such that, in
                      ///<   fractional space, all components fall in the range [ -0.5, 0.5 ).
};

/// \brief List methods for computing (positional) Root Mean Squared Deviation (RMSD) between
///        two sets of coordinates
enum class RmsdMethod {
  ALIGN_MASS,     ///< Align the two structures, centering each by their respective centers of
                  ///<   mass, prior to computing mass-weighted positional RMSD
  ALIGN_GEOM,     ///< Align the two structures, centering each by their respective centers of
                  ///<   geometry, prior to computing positional RMSD with no mass weighting
  NO_ALIGN_MASS,  ///< Do not align the two structures prior to computing positional RMSD with
                  ///<   mass weighting
  NO_ALIGN_GEOM   ///< Do not align the two structures prior to computing positional RMSD without
                  ///<   mass weighting
};

} // namespace structure
} // namespace stormm

#endif
