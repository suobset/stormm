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

/// \brief Unguarded struct to hold elemens of a flat bottom restraint.  Two such objects can hold
///        the initial and final states of a restraint, as well as the steps at which each state
///        takes effect.
struct FlatBottomPlan {

  /// \brief The default constructor sets stiffnesses to zero and the activation step to zero.
  ///        Additional constructors touch on the most frequently used combinations of parameters.
  FlatBottomPlan();
  FlatBottomPlan(double k_in, double r_in);
  FlatBottomPlan(double k_in, double r2_in, double r3_in);
  FlatBottomPlan(double k2_in, double k3_in, double r2_in, double r3_in);
  FlatBottomPlan(double k2_in, double k3_in, double r1_in, double r2_in, double r3_in,
                 double r4_in, int step_in);
  
  int activation_step;  // Step number at which this state of the restraint will take (full) effect
  double k2;            // Left-hand harmonic stiffness
  double k3;            // Right-hand harmonic stiffness
  double r1;            // Leftmost boundary between linear and harmonic behavior
  double r2;            // Leftmost boundary of the flat-bottom potential
  double r3;            // Rightmost boundary of the flat-bottom potential
  double r4;            // Rightmost boundard between harmonic and linear behavior
};
  
/// \brief Build a restraint to hold two atoms at a particular distance from one another.
///  
/// \param ag      System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i  The first atom of the angle
/// \param atom_j  The second atom of the angle
BoundedRestraint applyDistanceRestraint(const AtomGraph *ag, int atom_i, int atom_j,
                                        const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final = FlatBottomPlan());

/// \brief Build a restraint to hold three atoms in a particular angle arrangement.
///  
/// \param ag      System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i  The first atom of the angle
/// \param atom_j  The second atom of the angle
/// \param atom_k  The third atom of the angle
BoundedRestraint applyAngleRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                     const FlatBottomPlan fb_init,
                                     const FlatBottomPlan fb_final = FlatBottomPlan());
  
/// \brief Build a restraint to hold a dihedral towards a particular value.  This can rotate the
///        group about a particular bond, or secure a conformation in a position to which it has
///        already been rotated.
///  
/// \param ag      System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i  The first atom of the dihedral
/// \param atom_j  The second atom of the dihedral (first atom of the bond about which to rotate)
/// \param atom_k  The third atom of the dihedral (second atom of the bond about which to rotate)
/// \param atom_l  The fourth atom of the dihedral
BoundedRestraint applyDihedralRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                        int atom_l, const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final = FlatBottomPlan());

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
  
} // namespace restraints
} // namespace omni

#endif
