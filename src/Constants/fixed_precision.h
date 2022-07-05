// -*-c++-*-
#ifndef OMNI_FIXED_PRECISION_H
#define OMNI_FIXED_PRECISION_H

#include <cmath>
#include <string>
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"

namespace omni {
namespace numerics {

using constants::ExceptionResponse;
using constants::int_bit_count_int;
using constants::llint_bit_count_int;
  
/// \brief Enumerate different available precision models
enum class PrecisionLevel {
  SINGLE,       ///< Computations will take place in fp32 or int64 fixed-precision.  This may limit
                ///<   the upper bound on some fixed-precision bit settings.
  SINGLE_PLUS,  ///< Computations will take place in fp32 or int64 fixed-precision, with some
                ///<   enhancements depending on the situation.  This may limit the upper bound on
                ///<   some fixed-precision bit settings.
  DOUBLE        ///< Computations will take place in fp64 or int64 fixed-precision
};

/// \brief Enumerate the choices for carrying out fixed-precision accumulation
enum class ForceAccumulationMethod {
  SPLIT,     ///< Use split accumulation, stashing the low 32 bits in a locally cached int and the
             ///<   high 32 bits in a secondary accumulator probably located further away in main
             ///<   memory.  So long as most of the work happens in the low 32 bits, this reduces
             ///<   local memory demand and overall memory bandwidth by a factor of two, lowers
             ///<   GPU kernel register pressure on many architectures, and has shown 2.2 - 2.9x
             ///<   the speed of accumulating in int64.
  WHOLE,     ///< Sum fixed-precision numbers in int64 accumulators.  This is needed when the
             ///<   fixed-precision work cannot be mostly confined to the low 32 bits.
  AUTOMATIC  ///< Determine the accumulation method by looking at the number of fixed-precision
             ///<   bits after the decimal and making some assumptions about typical molecular
             ///<   mechanics forces.
};

/// \brief The fixed-precision discretizations of global and local coordinate frame positions, in
///        parts per Angstrom (the internal unit of length).  The local position scaling is
///        intended to work in local coordinate frames, perhaps in an implicit solvent enviroment
///        but more likely with an explicit solvent neighbor list, to encode positions at lower
///        resolution with 32-bit numbers for non-bonded calculations.
/// \{
constexpr double default_globalpos_scale_lf = 268435456.0;
constexpr float  default_globalpos_scale_f  = (float)default_globalpos_scale_lf;
constexpr int    default_globalpos_scale_bits = 28;
constexpr double default_inverse_globalpos_scale_lf = 1.0 / default_globalpos_scale_lf;
constexpr float  default_inverse_globalpos_scale_f  = (float)1.0 / default_globalpos_scale_f;
constexpr int    min_globalpos_scale_bits = 24;
constexpr int    max_globalpos_scale_bits = 48;
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
constexpr int    min_velocity_scale_bits = 28;
constexpr int    max_velocity_scale_bits = 52;
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
constexpr int    min_force_scale_bits = 18;
constexpr int    max_force_scale_bits = 48;
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

/// \brief The maximum contribution for signed integer accumulation
/// \{
constexpr llint max_llint_accumulation_ll = (1LL << (llint_bit_count_int - 1));
constexpr double max_llint_accumulation   = max_llint_accumulation_ll;
constexpr float max_llint_accumulation_f  = max_llint_accumulation;
constexpr llint max_int_accumulation_ll = (1LL << (int_bit_count_int - 1));
constexpr double max_int_accumulation   = max_int_accumulation_ll;
constexpr float max_int_accumulation_f  = max_int_accumulation;
/// \}

/// \brief Translate a string into a known precision level enumeration.
///
/// \param choice  The named precision level (will be checked for validity)
PrecisionLevel translatePrecisionLevel(const std::string &choice, ExceptionResponse policy);

/// \brief Get a descriptive string corresponding to each enumerated compute precision level.
///
/// \param plevel  The precision level to name
std::string getPrecisionLevelName(PrecisionLevel plevel);

/// \brief Translate a string specifying a force accumulation method into the numerical code.
///
/// \param method  The string to translate
ForceAccumulationMethod translateForceAccumulationMethod(const std::string &choice,
                                                         ExceptionResponse policy);

/// \brief Get a string for the name of a force accumulation method.
///
/// \param method  The method in question
std::string getForceAccumulationMethodName(ForceAccumulationMethod method);

/// \brief Produce an error message describing range violations in user choices for various
///        fixed-precision methods.
///
/// \param choice   The selected number of bits in the precision model
/// \param min_val  The minimum allowed number of bits
/// \param max_val  The maximum allowed number of bits
std::string fixedPrecisionRangeErrorMessage(int choice, int min_val, int max_val);

/// \brief Check user input regarding the global position scaling.  Positions are represented in
///        internal units of Angstroms.
///
/// \param choice  The fixed-precision bits for representing global positions of particles
void checkGlobalPositionBits(int choice);

/// \brief Check user input regarding the local position scaling.  Positions are represented in
///        internal units of Angstroms.
///
/// \param choice  The fixed-precision bits for representing local positions of particles
void checkLocalPositionBits(int choice);

/// \brief Check user input regarding the fixed-precision velocity scaling.  Velocities are
///        represented in units of Angstroms per sqrt(418.4) * femtoseconds.
///
/// \param choice  The fixed-precision bits for representing particle velocities
void checkVelocityBits(int choice);

/// \brief Check user input regarding the fixed-precision force accumulation.  Forces are
///        represented in units of kcal/mol - Angstroms.
///
/// \param choice  The fixed-precision bits for representing forces acting on particles
void checkForceBits(int choice);

/// \brief Check user input regarding the fixed-precision energy accumulation.  Energies are
///        represented in units of kcal/mol.
///
/// \param choice  The fixed-precision bits for representing energy contributions
void checkEnergyBits(int choice);

/// \brief Check user input regarding the fixed-precision charge denisty accumulation on the mesh.
///        Charge density per grid point is represented in atomic units.
///
/// \param choice  The fixed-precision bits for performing charge density accumulation
/// \param pmdoel  The fixed-precision model, which implies the accumulation range
void checkChargeMeshBits(int choice, PrecisionLevel pmodel);

} // namespace numerics
} // namespace omni

#endif
