// -*-c++-*-
#ifndef OMNI_FIXED_PRECISION_H
#define OMNI_FIXED_PRECISION_H

namespace omni {
namespace numerics {

/// \brief The fixed-precision discretizations of global and local coordinate frame positions, in
///        parts per Angstrom (the internal unit of length).  
/// \{
constexpr double global_position_scale_lf = 268435456.0;
constexpr float  global_position_scale_f  = (float)global_position_scale_lf;
constexpr int    global_position_scale_bits = 28;
constexpr double inverse_global_position_scale_lf = 1.0 / global_position_scale_lf;
constexpr float  inverse_global_position_scale_f  = (float)1.0 / global_position_scale_f;
constexpr double local_position_scale_lf = 16777216.0;
constexpr float  local_position_scale_f  = (float)local_position_scale_lf;
constexpr int    local_position_scale_bits = 24;
constexpr double inverse_local_position_scale_lf = 1.0 / local_position_scale_lf;
constexpr float  inverse_local_position_scale_f  = (float)1.0 / local_position_scale_f;
/// \}

/// \brief Time is expressed in units of femtoseconds, and forces are discretized into increments
///        one part in 8,388,608 of one kcal/mol-A.  For 32-bit integer accumulation, this gives a
///        range of [ -256.0, +256.0 ) for each of three force components, which will suffice for
///        nearly all interactions.  The rare very large interaction will overflow a 32-bit "minor"
///        accumulator, but this will be detected and the excess will flow into a second "major"
///        accumulator, ensuring that there is no practical upper bound to the magnitudes of
///        forces that can be accumulated.
/// \{
constexpr double global_force_scale_lf = 8388608.0;
constexpr float  global_force_scale_f  = (float)global_force_scale_lf;
constexpr int    global_force_scale_bits = 23;
constexpr double nonbond_force_scale_lf = 8388608.0;
constexpr float  nonbond_force_scale_f  = (float)nonbond_force_scale_lf;
constexpr int    nonbond_force_scale_bits = 23;
/// \}

/// \brief Energies are accumulated in units of kcal/mol, discretized into one part in 33554432
///        of one kcal/mol.  Summation funnels directly into 64-bit, long long int accumulators,
///        again with no practical upper bound on the numbers.
/// \{
constexpr double global_energy_scale_lf = 33554432.0;
constexpr float  global_energy_scale_f  = (float)global_energy_scale_lf;
constexpr int    global_energy_scale_bits = 25;
constexpr double inverse_global_energy_scale_lf = 1.0 / global_energy_scale_lf;
constexpr float  inverse_global_energy_scale_f  = (float)1.0 / global_energy_scale_f;
/// \}
  
} // namespace numerics
} // namespace omni

#endif
