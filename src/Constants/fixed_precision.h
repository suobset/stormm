// -*-c++-*-
#ifndef OMNI_FIXED_PRECISION_H
#define OMNI_FIXED_PRECISION_H

namespace omni {
namespace numerics {

/// \brief The fixed-precision discretizations of global and local coordinate frame positions, in
///        parts per Angstrom (the internal unit of length).  
/// \{
constexpr double default_globalpos_scale_lf = 268435456.0;
constexpr float  default_globalpos_scale_f  = (float)default_globalpos_scale_lf;
constexpr int    default_globalpos_scale_bits = 28;
constexpr double default_inverse_globalpos_scale_lf = 1.0 / default_globalpos_scale_lf;
constexpr float  default_inverse_globalpos_scale_f  = (float)1.0 / default_globalpos_scale_f;
constexpr double default_localpos_scale_lf = 16777216.0;
constexpr float  default_localpos_scale_f  = (float)default_localpos_scale_lf;
constexpr int    default_localpos_scale_bits = 24;
constexpr double default_inverse_localpos_scale_lf = 1.0 / default_localpos_scale_lf;
constexpr float  default_inverse_localpos_scale_f  = (float)1.0 / default_localpos_scale_f;
/// \}

/// \brief Velocities are expressed in A / sqrt(418.4) fs, and as such the velocity scaling should
///        be high in order to preserve bits commensurate with the force and position
///        quantities.
/// \{
constexpr double default_velocity_scale_lf = 17179869184.0;
constexpr float  default_velocity_scale_f  = (float)default_velocity_scale_lf;
constexpr int    default_velocity_scale_bits = 34;
constexpr double default_inverse_velocity_scale_lf = 1.0 / default_velocity_scale_lf;
constexpr float  default_inverse_velocity_scale_f = (float)1.0 / default_velocity_scale_f;
/// \}
  
/// \brief Time is expressed in units of femtoseconds, and forces are discretized into increments
///        one part in 8,388,608 of one kcal/mol-A.  For 32-bit integer accumulation, this gives a
///        range of [ -256.0, +256.0 ) for each of three force components, which will suffice for
///        nearly all interactions.  The rare very large interaction will overflow a 32-bit "minor"
///        accumulator, but this will be detected and the excess will flow into a second "major"
///        accumulator, ensuring that there is no practical upper bound to the magnitudes of
///        forces that can be accumulated.
/// \{
constexpr double default_force_scale_lf = 8388608.0;
constexpr float  default_force_scale_f  = (float)default_force_scale_lf;
constexpr int    default_force_scale_bits = 23;
constexpr double default_inverse_force_scale_lf = 1.0 / default_force_scale_lf;
constexpr float  default_inverse_force_scale_f  = (float)default_inverse_force_scale_lf;
/// \}

/// \brief Energies are accumulated in units of kcal/mol, discretized into one part in 33554432
///        of one kcal/mol.  Summation funnels directly into 64-bit, long long int accumulators,
///        again with no practical upper bound on the numbers.
/// \{
constexpr double default_energy_scale_lf = 33554432.0;
constexpr float  default_energy_scale_f  = (float)default_energy_scale_lf;
constexpr int    default_energy_scale_bits = 25;
constexpr double default_inverse_energy_scale_lf = 1.0 / default_energy_scale_lf;
constexpr float  default_inverse_energy_scale_f  = (float)1.0 / default_energy_scale_f;
/// \}
  
} // namespace numerics
} // namespace omni

#endif
