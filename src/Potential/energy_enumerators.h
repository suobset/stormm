// -*-c++-*-
#ifndef STORMM_ENERGY_ENUMERATORS_H
#define STORMM_ENERGY_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace energy {

/// \brief Enumerate the choices on whether to evaluate the energy... yes or no.  CPU functions
///        are obligated to evaluate the energy, but GPU functions may or may not due to register
///        pressure that it creates.
enum class EvaluateEnergy {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the force... yes or no.
enum class EvaluateForce {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the virial... yes or no.
enum class EvaluateVirial {
  NO, YES
};

/// \brief Enumerate the types of dihedral potentials that can be evaluated in the context of the
///        abstracted routine.
enum class DihedralStyle {
  COSINE, HARMONIC
};

/// \brief Enumerate all state variables that STORMM will track.  This must always increment by one
///        from one enumeration to the next, start at zero, and always end with ALL_STATES.
enum class StateVariable {
  BOND = 0,               ///< Harmonic bond stretching energy
  ANGLE,                  ///< Harmonic angle bending energy
  PROPER_DIHEDRAL,        ///< Proper dihedral energy
  IMPROPER_DIHEDRAL,      ///< Improper (plane-enforcing) dihedral energy
  UREY_BRADLEY,           ///< CHARMM Urey-Bradley angle stretching energy
  CHARMM_IMPROPER,        ///< CHARMM harmonic improper torsion energy
  CMAP,                   ///< Correction map (typically, coupled dihedral-dihedral) energy
  VDW,                    ///< van-der Waals (typically, Lennard Jones) energy from non-bonded
                          ///<   interactions
  VDW_ONE_FOUR,           ///< van-der Waals (typically, Lennard Jones) energy from 1-4 attenuated
                          ///<   interactions
  ELECTROSTATIC,          ///< Electrostatic energy from non-bonded interactions
  ELECTROSTATIC_ONE_FOUR, ///< Electrostatic energy from 1-4 attenuated interactions
  GENERALIZED_BORN,       ///< Generalized Born (implicit solvent) energy
  RESTRAINT,              ///< Energy due to flat-bottom bimodal harmonic potential restraints
  KINETIC,                ///< Energy due to particle motion
  PRESSURE,               ///< System pressure (only computed if virials are accumulated)
  VIRIAL_11,              ///< Virial tensor (1,1) element
  VIRIAL_12,              ///< Virial tensor (1,2) element
  VIRIAL_22,              ///< Virial tensor (2,2) element
  VIRIAL_13,              ///< Virial tensor (1,3) element
  VIRIAL_23,              ///< Virial tensor (2,3) element
  VIRIAL_33,              ///< Virial tensor (3,3) element
  VOLUME,                 ///< Unit cell volume (only relevant to periodic simulations)
  TEMPERATURE_ALL,        ///< Overall system temperature
  TEMPERATURE_PROTEIN,    ///< Temperature of atoms in the "protein" component of the system (this,
                          ///<   like the other temperature subcategories that collow, is an
                          ///<   arbitrary subset of the atoms defined by a special mask in the
                          ///<   control input)
  TEMPERATURE_LIGAND,     ///< Overall system temperature
  TEMPERATURE_SOLVENT,    ///< Overall system temperature
  DU_DLAMBDA,             ///< Derivative of the mixed potential energy function with respect to
                          ///<   the mixing parameter Lambda (relevant to thermodynamic
                          ///<   integration applications only)
  POTENTIAL_ENERGY,       ///< Sum of all potential energy contributions in the system
  TOTAL_ENERGY,           ///< Sum of all potential and kinetic energy components in the system
  ALL_STATES              ///< This must always be the final entry.  The number of tracked
                          ///<   quantities is equal to the value of this entry (ALL_STATES does
                          ///<   not define its own index in the subsequent tracking arrays).
};

/// \brief Enumerate the specific kinds of non-bonded potentials.
enum class NonbondedPotential {
  ELECTROSTATIC,  ///< Interactions between charges
  VAN_DER_WAALS,  ///< Dispersion iteractions between particles, likely a Lennard-Jones potential
  CLASH           ///< Zero (no clash) or one (clash) based on the van-der Waals radii and the
                  ///<   width of some probe sphere representing heavy atoms in the ligand
};

/// \brief List the methods for combining Lennard-Jones (or other van-der Waals potential)
///        parameters in pairs of interacting atoms.
enum class VdwCombiningRule {
  LORENTZ_BERTHELOT,  ///< Sigma parameters s_i and s_j are combined by taking (s_i + s_j) / 2.
                      ///<   Epsilon parameters are always combined by taking sqrt(E_i * E_j).
  GEOMETRIC,          ///< Sigma parameters s_i and s_j are combined by taking sqrt(s_i * s_j)
  NBFIX               ///< Pair-specific combinations, with some off-diagonal terms not conforming
                      ///<   to either of the other rules
};

/// \brief Functions and kernels can be configured to dampen the effects of clashes (at a minor
///        expense in computation and registers), or not.
enum class ClashResponse {
  NONE,    ///< Do not attempt to dampen the effects of clashes between particles
  FORGIVE  ///< Forgive clashes by not letting the perceived interparticle distance drop below
           ///<   some multiple of the pairwise Lennard-Jones sigma value, or if that is also zero
           ///<   then some arbitrary minimum value.
};

/// \brief List the different modes in which to extract energies from a ScoreCard object.
enum class EnergySample {
  TIME_SERIES,  ///< Return the entire time series of all values stored for a state variable
  FINAL,        ///< Return the current (most recent, or final) value of a state variable
  TIME_AVERAGE  ///< Return the average of the time series of a state variable
};
  
/// \brief Get a human-readable name for the enumerations detailed above.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(EvaluateForce input);
std::string getEnumerationName(EvaluateEnergy input);
std::string getEnumerationName(EvaluateVirial input);
std::string getEnumerationName(DihedralStyle input);
std::string getEnumerationName(StateVariable input);
std::string getEnumerationName(NonbondedPotential input);
std::string getEnumerationName(VdwCombiningRule input);
std::string getEnumerationName(ClashResponse input);
std::string getEnumerationName(EnergySample input);
/// \}

} // namespace energy
} // namespace stormm

#endif
