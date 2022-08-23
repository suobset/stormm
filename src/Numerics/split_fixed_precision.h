// -*-c++-*-
#ifndef STORMM_SPLIT_FIXED_PRECISION_H
#define STORMM_SPLIT_FIXED_PRECISION_H

#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace numerics {

using constants::ExceptionResponse;
using constants::PrecisionModel;
using data_types::int95_t;
  
/// \brief Enumerate the choices for carrying out fixed-precision accumulation
enum class AccumulationMethod {
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

/// \brief The maximum contributions for signed integer accumulation.  There is no long long
///        integer form of the maximum long long integer accumulation as the number would wrap
///        the format to become the negative of its intended value.  The long long integer
///        representation of the maximum int accumulation is present to facilitate conversion to
///        a unified 64-bit integer, but even that is only reliable on the CPU (the NVIDIA CUDA
///        compiler seems to take type specifications of long long int as a suggestion, not a
///        command, and this may relate to thread counts and register pressure as the compiler
///        tries to optimize a kernel for given launch bounds).
/// \{
constexpr llint max_int_accumulation_ll = (1LL << (int_bit_count_int - 1));
constexpr double max_int_accumulation   = max_int_accumulation_ll;
constexpr float max_int_accumulation_f  = max_int_accumulation;
constexpr double max_llint_accumulation = max_int_accumulation * max_int_accumulation * 2.0;
constexpr float max_llint_accumulation_f  = max_llint_accumulation;
/// \}

/// \brief Translate a string specifying a force accumulation method into the numerical code.
///
/// \param method  The string to translate
AccumulationMethod translateAccumulationMethod(const std::string &choice,
                                                         ExceptionResponse policy);

/// \brief Get a string for the name of a force accumulation method.
///
/// \param method  The method in question
std::string getAccumulationMethodName(AccumulationMethod method);

/// \brief Determine the best accumulation method based on the precision level of the forces.
///
/// \param frc_bits  Number of bits stored after the decimal in fixed-precision force
///                  representations
AccumulationMethod chooseAccumulationMethod(int frc_bits);

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
void checkChargeMeshBits(int choice, PrecisionModel pmodel);

/// \brief Convert floating point numbers into fixed-precision representations with two integers.
///        This is similar to splitRealAccumulation below, but will set the values rather than
///        add new contributions.
///
/// Overloaded:
///   - Convert a single-precision floating point number into two 32-bit signed integers.
///   - Convert a double-precision floating point number into a 64-bit primary integer and a
///     32-bit secondary / overflow integer.
///   - Return the fixed-precision representation or assign it directly to the corresponding
///     locations in two appropriate arrays.
///
/// \param fval      Single-precision value to convert to fixed-precision
/// \param dval      Double-precision value to convert to fixed-precision
/// \param primary   The primary accumulator (the low 32 bits)
/// \param overflow  The secondary accumulator (the high 31 bits)
/// \{
int2 floatToInt63(const float fval);

void floatToInt63(const float fval, int *primary, int *overflow);

int95_t doubleToInt95(const double fval);

void doubleToInt95(const double fval, llint *primary, int *overflow);
/// \}

/// \brief Accumulate floating point numbers into fixed-precision representations with two
///        integers.
///
/// Overloaded:
///   - Convert a single-precision floating point number into two 32-bit signed integers.
///   - Convert a double-precision floating point number into a 64-bit primary integer and a
///     32-bit secondary / overflow integer.
///
/// \param fval      Single-precision value to convert to fixed-precision
/// \param dval      Double-precision value to convert to fixed-precision
/// \param primary   The primary accumulator (the low 32 bits)
/// \param overflow  The secondary accumulator (the high 31 bits)
/// \{
void splitAccumulation(const float fval, int *primary, int *overflow);

void splitAccumulation(const double fval, llint *primary, int *overflow);
/// \}

} // namespace numerics
} // namespace stormm

#endif

