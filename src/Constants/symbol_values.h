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
constexpr float pi_f = 3.1415927410125732421875f;
constexpr float twopi_f = 6.283185482025146484375f;
constexpr float inverse_pi_f = 0.3183098733425140380859375f;
constexpr float inverse_twopi_f = 0.15915493667125701904296875f;

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

/// \brief Values which approach one from below.  They are used in dihedral and similar
///        computations to detect when a value is nearing 1.0 and might generate a singularity
///        in some denominator, or otherwise become numerically ill-conditioned.  The first is
///        nearly the closest value to 1.0 that can be represented in an IEEE-754 format 32-bit
///        floating point number (the format can go a couple of bits' worth of precision closer,
///        but this is 1 part in about a million).  The second is essential for single-precision
///        dihedral computations.
/// \{
constexpr double asymptotic_to_one_lf = 0.99999904632568359375;
constexpr float  asymptotic_to_one_f  = (float)asymptotic_to_one_lf;
constexpr float  near_to_one_f        = 0.99993896484375f;
/// \}                                                                                              

/// \brief A value which captures 1 / (1 - asymptotic_to_one), to put a cap on the value of such    
///        fractions 1 / (1 - x) as x -> 1.                                                         
/// \{                                                                                              
constexpr double inverse_one_minus_asymptote_lf = 1048576.0;
constexpr float  inverse_one_minus_asymptote_f = (float)1048576.0;
/// \}                                  

} // namespace constants
} // namespace omni

#endif
