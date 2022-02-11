// -*-c++-*-
#ifndef	OMNI_RESTRAINT_BUILDER_H
#define	OMNI_RESTRAINT_BUILDER_H

#include <vector>
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
bool restraintTopologyChecks(const AtomGraph *ag, const CoordinateFrameReader &cframe,
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

/// \brief Build restraints needed to manage the inversion of chiral centers in a molecule.  The
///        strategy will be to maintain as many details of the original conformation as possible,
///        as expressed in the unique dihedral angles not involving a chiral center that is about
///        to undergo inversion by the function invertChirality() in the conformers.cpp library of
///        the src/Conformations/ folder.  The 
///
/// \param ag            Pointer to the topology for the system of interest
/// \param crd           Coordinates of the system in its current state
/// \param chiral_atoms  List of all atoms with detected chirality that are to be inverted (these
///                      atoms will help define what dihedrals to avoid constraining)
std::vector<BoundedRestraint>
applyChiralityInversionRestraints(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                                  const AtomMask &chiral_mask);

} // namespace restraints
} // namespace omni

#endif
