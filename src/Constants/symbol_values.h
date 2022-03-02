// -*-c++-*-
#ifndef OMNI_SYMBOLS_H
#define OMNI_SYMBOLS_H

namespace omni {
namespace symbols {

/// \brief Pi and related quantities, represented in values that can be stored exactly in
///        double-precision floating point numbers according to the IEEE_754 standard.
constexpr double pi = 3.141592653589793115997963468544185161590576171875;
constexpr double twopi = 2.0 * pi;
constexpr double inverse_pi = 1.0 / pi;
constexpr double inverse_twopi = 1.0 / twopi;

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
