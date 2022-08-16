// -*-c++-*-
#ifndef STORMM_FIXED_PRECISION_H
#define STORMM_FIXED_PRECISION_H

#include <cmath>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace numerics {

using constants::ExceptionResponse;
using constants::int_bit_count_int;
using constants::llint_bit_count_int;
using constants::PrecisionModel;

/// \brief The fixed-precision discretizations of global and local coordinate frame positions, in
///        parts per Angstrom (the internal unit of length).  The local position scaling is
///        intended to work in local coordinate frames, perhaps in an implicit solvent enviroment
///        but more likely with an explicit solvent neighbor list, to encode positions at lower
///        resolution with 32-bit numbers for non-bonded calculations.  The "nonoverflow" bit count
///        indicates a threshold beneath which the overflow counters will be assumed never to be
///        touched.  The "transform" bit count and scaling factors are fixed values to ensure that
///        large unit cells (up to +/- 4096 Angstroms) can be handled while making the most out of
///        the 64-bit integer format.
/// \{
constexpr double default_globalpos_scale_lf = 268435456.0;
constexpr float  default_globalpos_scale_f  = (float)default_globalpos_scale_lf;
constexpr int    default_globalpos_scale_bits = 28;
constexpr double default_inverse_globalpos_scale_lf = 1.0 / default_globalpos_scale_lf;
constexpr float  default_inverse_globalpos_scale_f  = (float)1.0 / default_globalpos_scale_f;
constexpr int    min_globalpos_scale_bits = 24;
constexpr int    max_globalpos_scale_bits = 72;
constexpr int    globalpos_scale_nonoverflow_bits = 32;
constexpr double default_localpos_scale_lf = 16777216.0;
constexpr float  default_localpos_scale_f  = (float)default_localpos_scale_lf;
constexpr int    default_localpos_scale_bits = 24;
constexpr double default_inverse_localpos_scale_lf = 1.0 / default_localpos_scale_lf;
constexpr float  default_inverse_localpos_scale_f  = (float)1.0 / default_localpos_scale_f;
constexpr int    min_localpos_scale_bits = 20;
constexpr int    max_localpos_scale_bits = 26;
constexpr double default_trajpos_scale_lf = 16384.0;
constexpr float  default_trajpos_scale_f  = (float)default_trajpos_scale_lf;
constexpr int    default_trajpos_scale_bits = 14;
constexpr double default_inverse_trajpos_scale_lf = 1.0 / default_trajpos_scale_lf;
constexpr float  default_inverse_trajpos_scale_f  = (float)1.0 / default_trajpos_scale_f;
constexpr int    min_trajpos_scale_bits = 10;
constexpr int    max_trajpos_scale_bits = 48;
constexpr int    transform_scale_bits = 51;
constexpr double transform_scale_lf   = 2251799813685248.0;
constexpr float  transform_scale_f    = (float)transform_scale_lf;
/// \}

/// \brief Velocities are expressed in A / sqrt(418.4) fs, and as such the velocity scaling should
///        be high in order to preserve bits commensurate with the force and position
///        quantities.  The "nonoverflow" bit count indicates a threshold beneath which the
///        overflow counters will be assumed never to be touched.
/// \{
constexpr double default_velocity_scale_lf = 17179869184.0;
constexpr float  default_velocity_scale_f  = (float)default_velocity_scale_lf;
constexpr int    default_velocity_scale_bits = 34;
constexpr double default_inverse_velocity_scale_lf = 1.0 / default_velocity_scale_lf;
constexpr float  default_inverse_velocity_scale_f = (float)1.0 / default_velocity_scale_f;
constexpr int    min_velocity_scale_bits = 28;
constexpr int    max_velocity_scale_bits = 72;
constexpr int    velocity_scale_nonoverflow_bits = 44;
/// \}
  
/// \brief Time is expressed in units of femtoseconds, and forces are discretized into increments
///        one part in 8,388,608 of one kcal/mol-A.  For 32-bit integer accumulation, this gives a
///        range of [ -256.0, +256.0 ) for each of three force components, which will suffice for
///        nearly all interactions.  The rare very large interaction will overflow a 32-bit "minor"
///        accumulator, but this will be detected and the excess will flow into a second "major"
///        accumulator, ensuring that there is no practical upper bound to the magnitudes of
///        forces that can be accumulated.  The "nonoverflow" bit count indicates a threshold
///        beneath which overflow counters will be assumed never to be touched.
/// \{
constexpr double default_force_scale_lf = 8388608.0;
constexpr float  default_force_scale_f  = (float)default_force_scale_lf;
constexpr int    default_force_scale_bits = 23;
constexpr double default_inverse_force_scale_lf = 1.0 / default_force_scale_lf;
constexpr float  default_inverse_force_scale_f  = (float)default_inverse_force_scale_lf;
constexpr int    min_force_scale_bits = 18;
constexpr int    max_force_scale_bits = 72;
constexpr int    force_scale_nonoverflow_bits = 40;
/// \}

/// \brief Energies are accumulated in units of kcal/mol, discretized into one part in 33554432
///        of one kcal/mol.  Summation funnels directly into 64-bit, long long int accumulators,
///        again with no practical upper bound on the numbers.
/// \{
constexpr double default_energy_scale_lf = 67108864.0;
constexpr float  default_energy_scale_f  = (float)default_energy_scale_lf;
constexpr int    default_energy_scale_bits = 26;
constexpr double default_inverse_energy_scale_lf = 1.0 / default_energy_scale_lf;
constexpr float  default_inverse_energy_scale_f  = (float)default_inverse_energy_scale_lf;
constexpr int    min_energy_scale_bits = 22;
constexpr int    max_energy_scale_bits = 40;
/// \}

/// \brief Charges are mapped to the mesh in atomic units--the forces on atoms are later multiplied
///        by Coulomb's constant, just like the electrostatic short-ranged non-bonded forces.
/// \{
constexpr double default_charge_mesh_scale_lf = 268435456.0;
constexpr float  default_charge_mesh_scale_f = (float)default_charge_mesh_scale_lf;
constexpr int    default_charge_mesh_scale_bits = 28;
constexpr double default_inverse_charge_mesh_scale_lf = 1.0 / default_charge_mesh_scale_lf;
constexpr float  default_inverse_charge_mesh_scale_f = (float)default_inverse_charge_mesh_scale_lf;
constexpr int    min_charge_mesh_scale_bits = 24;
constexpr int    max_charge_mesh_scale_bits = 48;
/// \}

} // namespace numerics
} // namespace stormm

#endif
