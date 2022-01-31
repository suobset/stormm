// -*-c++-*-
#ifndef OMNI_SYMBOLS_H
#define OMNI_SYMBOLS_H

namespace omni {
namespace symbols {

/// The number pi, to twenty decimal places (deeper than double precision)
constexpr double pi = 3.14159265358979323846;
constexpr double twopi = 6.28318530717958647692;
constexpr double inverse_twopi = 0.159154943091895345608;

/// The tetrahedral bond angle, in radians to twenty decimal places
constexpr double tetrahedral_angle = 1.91063323624901859610;

/// The "biological charge" conversion, taking atomic units of charge, Coulomb's constant, and
/// Angstroms into kilocalories per mole.  Amber and CHARMM disagree on this value.  This is the
/// energy, in kcal/mol, of a pair of charges with 1 a.u. charge being brought from infinity to
/// within 1 Angstrom of each other.
/// \{
constexpr double amber_ancient_bioq  = 332.0522172900000;
constexpr double charmm_gromacs_bioq = 332.0636974382250;
/// \}

/// \brief Avogadro's number
constexpr double avogadro_number = 6.02214076e+23;
  
/// \brief Boltzmann's constant in kcal/mol-K
constexpr double boltzmann_constant = (1.38064852e-23) / 4184.0 * avogadro_number;

/// \brief Hartree to kcal/mol conversion
constexpr double hartree_to_kcal = 627.509474;
  
} // namespace constants
} // namespace omni

#endif
