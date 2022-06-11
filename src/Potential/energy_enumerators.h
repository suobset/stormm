// -*-c++-*-
#ifndef OMNI_ENERGY_ENUMERATORS_H
#define OMNI_ENERGY_ENUMERATORS_H

namespace omni {
namespace energy {

/// \brief Enumerate the choices on whether to evaluate the force... yes or no.
enum class EvaluateForce {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the energy... yes or no.  CPU functions
///        are obligated to evaluate the energy, but GPU functions may or may not due to register
///        pressure that it creates.
enum class EvaluateEnergy {
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

/// \brief Enumerate all state variables that OMNI will track.  This must always increment by one
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
  ALL_STATES              ///< This must always be the final entry.  The number of tracked
                          ///<   quantities is equal to the value of this entry (ALL_STATES does
                          ///<   not define its own index in the subsequent tracking arrays).
};

} // namespace energy
} // namespace omni

#endif
