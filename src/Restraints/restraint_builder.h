// -*-c++-*-
#ifndef	OMNI_RESTRAINT_BUILDER_H
#define	OMNI_RESTRAINT_BUILDER_H

#include <vector>
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

/// \brief Apply positional restraints to a topology based on an AtomMask, with general parameters
///        for the permittivity and stiffness.
///
/// \param ag                    Topology for the system of interest
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
applyPositionalRestraints(const AtomGraph &ag, const CoordinateFrameReader &cframe,
                          const CoordinateFrameReader &reference_cframe, const AtomMask &mask,
                          double displacement_penalty, double displacement_onset,
                          double displacement_plateau, double proximity_penalty,
                          double proximity_onset, double proximity_plateau);

/// \brief Build restraints needed to invert chiral centers in a molecule.  Given a chiral center
///        detected by the chemical perception code (see chemical_features.h in the src/Chemistry/
///        directory), the strategy is to construct positional restraints that will keep two of the
///        arms (R1, R2) and the center (C) itself in place while swapping the position of the
///        other two arms with a 180-degree rotation about the bisector of R1--C--R2.  The swap
///        will be made based on what appears easiest to rotate, and will detect whether two of
///        the bonds connecting to a chiral center are part of some ring system (which would
///        designate them R1 and R2, not subject to rotation).  In the rare case that a
///        chiral center is part of two rings, the rotation will be chosen to swap the orientations
///        of the rings.  In the rare, rare case that a chiral center is an intersection of three
///        rings, there will be another center, chiral or not, which lies at a second intersection
///        of the rings, which will undergo a complementary rotation.
///
/// \param ag           Topology for the system of interest
/// \param crd          Coordinates of the system in its current state
/// \param chiral_mask  Mask of all atoms with detected chirality (it is not critical what the
///                     value of the chirality is, or the IUPAC priority of each of the chiral
///                     arms--all that will happen is that two of the arms will be selected and
///                     their places swapped)
std::vector<BoundedRestraint>
applyChiralityInversionRestraints(const AtomGraph &ag, const CoordinateFrameReader &cframe,
                                  const AtomMask &chiral_mask);

} // namespace restraints
} // namespace omni

#endif
