// -*-c++-*-
#ifndef	OMNI_RESTRAINT_BUILDER_H
#define	OMNI_RESTRAINT_BUILDER_H

#include <vector>
#include "Chemistry/atommask.h"
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

using chemistry::AtomMask;
using topology::AtomGraph;
using trajectory::CoordinateFrameReader;

/// \brief Perform basic checks the validity of the topology and coordinate frame needed by the
///        following restraint builders.
///
/// \param ag                    Pointer to the topology for the system of interest
/// \param cframe                Coordinates of the system in its current state
/// \param reference_cframe      Reference coordinates of the system, if different from
/// \param mask                  Atom mask (must match the topology by some basic checks))
void restraintTopologyChecks(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                             const AtomMask &mask);
  
/// \brief Apply positional restraints to a topology based on an AtomMask, with general parameters
///        for the permittivity and stiffness.
///
/// \param ag                    Pointer to the topology for the system of interest
/// \param cframe                Coordinates of the system in its current state
/// \param reference_cframe      Reference coordinates of the system, if different from
/// \param mask                  Atom mask (must match the topology by some basic checks))
/// \param displacement_penalty  Harmonic stiffness constant for displacing particles away from
///                              the target locations.  This is k3 in the terminology of Amber NMR
///                              restraints.
/// \param displacement_onset    Point at which the restraint penalty begins for particles moving
///                              away from their target coordinates.  This is r3 in the Amber NMR
///                              restraint terminology and the default is 0.0.
/// \param displacement_plateau  Maximum displacment away from the target location beyond which
///                              the penalty goes linear and the restraining force is constant.
///                              This is r4 in the Amber NMR restraint terminology and the default
///                              behavior is to set this 5.0 Angstroms past displacement_onset.
/// \param proximity_penalty     Atoms that are not at least proximity_onset far from the target
///                              location will be subject to this harmonic penalty as their
///                              proximity to the target point increases (default 0.0, no
///                              proximity penalty is applied).  k2 in the Amber NMR nomenclature.
/// \param proximity_onset       The distance beneath which a penalty is applied for atoms being
///                              too close to the target location.  This is r2 in the Amber NMR
///                              nomenclature.
/// \param proximity_plateau     The distance beneath which the proximity penality increases
///                              linearly rather than quadratically and the penalty force
///                              flatlines.  This is r1 in the Amber NMR nomeclature.
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                          const CoordinateFrameReader &reference_cframe, const AtomMask &mask,
                          double displacement_penalty, double displacement_onset,
                          double displacement_plateau, double proximity_penalty = 0.0,
                          double proximity_onset = 0.0, double proximity_plateau = 0.0);

/// \brief Build restraints needed to maintain elements of the conformation not intended to change
///        their shape.  This will be accomplished by distance, angle, and dihedral restraints
///        between heavy atoms not intended to move.  
///
/// \param ag                      Pointer to the topology for the system of interest
/// \param crd                     Coordinates of the system in its current state
/// \param mask                    Mask of atoms to be held in place, relative to one another
/// \param penalty                 The stiffness of the harmonic penalty to apply.  This is
///                                identical for distance, angle, and dihedral restraints.
/// \param flat_bottom_half_width  The flat bottom of the restraint well extends this amount in
///                                either direction from the current, observed value of the
///                                distance, angle, or dihedral coordinate, in the appropriate
///                                unit system.
/// \param harmonic_maximum        The maximum value of the restraint in its harmonic range,
///                                defining the point at which the penalty goes linear.
std::vector<BoundedRestraint>
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                       const AtomMask &mask, double penalty, double flat_bottom_half_width,
                       double harmonic_maximum);

/// \brief Build a series of restraints to rotate about a particular bond.
///  
/// \param ag      System topology (this will be used to branch out from the rotatable bond
/// \param atom_i  The first atom of the bond about which to rotate
/// \param atom_j  The second atom of the bond about which to rotate

} // namespace restraints
} // namespace omni

#endif
