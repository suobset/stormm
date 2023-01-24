// -*-c++-*-
#ifndef STORMM_STRUCTURE_ENUMERATORS_H
#define STORMM_STRUCTURE_ENUMERATORS_H

#include "copyright.h"

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
enum class RMSDMethod {
  ALIGN_MASS,     ///< Align the two structures, centering each by their respective centers of
                  ///<   mass, prior to computing mass-weighted positional RMSD
  ALIGN_GEOM,     ///< Align the two structures, centering each by their respective centers of
                  ///<   geometry, prior to computing positional RMSD with no mass weighting
  NO_ALIGN_MASS,  ///< Do not align the two structures prior to computing positional RMSD with
                  ///<   mass weighting
  NO_ALIGN_GEOM   ///< Do not align the two structures prior to computing positional RMSD without
                  ///<   mass weighting
};

/// \brief List the strategies by which RMSD can be computed.  Depending on the mass of equivalent
///        atom groups, it may be necessary to do a combinatorial search of all symmetry-related
///        atoms until enough of the molecule's mass is placed, then test the minutiae of smaller
///        symmetric groups.  On the other end of the spectrum, there may be no symmetric atom
///        groups, which would allow the RMSD to be computed in a single pass.  If no alignment is
///        necessary, then it is possible to sample any symmetry-related groups in order of their
///        dependence on one another, in any order.
enum class RMSDAlignmentProtocol {
  BUILD_CORE,  ///< Sample the possible arrangements of the largest symmetry-related groups until
               ///<   enough mass has been placed, then proceed to...
  ALIGN_CORE,  ///< Align the molecule's massive core (starting with asymmetric atoms) and
               ///<   afterwards test the various arrangements of each symmetry-related group in
               ///<   the order of their dependencies.
  ALIGN_ALL    ///< With no symmetry-related groups, perform a standard RMSD calculation according
               ///<   to the RMSDMethod indicated.
};

/// \brief There are two orders of RMSD calculation: all to one (reference structure), or all to
///        all (matrix).
enum class RMSDTask {
  REFERENCE,  ///< Compute the RMSD of all structures to a reference structure
  MATRIX      ///< Compute the RMSD of all structures following a particular topology to all others
              ///<   following that same topology.
};

/// \brief Virtual site standalone functiosn and kernels call into two categories.
enum class VirtualSiteActivity {
  PLACEMENT,       ///< Place virtual sites after motion of the underlying frame atoms
  TRANSMIT_FORCES  ///< Transmit forces accumulated on virtual sites to their frame atoms
};

/// \brief Enumerate the various levels of grid detail by which a receptor may be mapped.
enum class GridDetail {
  OCCLUSION,        ///< The grid contains only an occlusion mask, with one bit per grid point
                    ///<   stored in a manner like std::vector<bool>.
  NONBONDED_FIELD,  ///< The grid stores electrostatic, Generalized Born radial derivatives,
                    ///<   Generalized Born pairwise energy, and van-der Waals potentials with
                    ///<   tricibic interpolation.
  NONBONDED_ATOMIC  ///< The grid stores all of the potentials found in the NONBONDED_FIELD case,
                    ///<   but only from atom contributions of rigid components of the receptor
                    ///<   which lie far enough from the bin that tricubic interpolation is a sound
                    ///<   approximation.  Other rigid atoms will be represented in a list of
                    ///<   explicit particles with unique segments for each bin.
};

/// \brief The types of clashes that can occur.
enum class ClashKind {
  VAN_DER_WAALS,  ///< The distance between two particles is less than the minimum van-der Waals
                  ///<   (Lennard-Jones) sigma ratio.
  PURE_DISTANCE   ///< The distance between two particles is less than the minimum absolute
                  ///<   distance allowed before a clash is declared
};
  
} // namespace structure
} // namespace stormm

#endif
