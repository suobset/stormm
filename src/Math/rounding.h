// -*-c++-*-
#ifndef STORMM_ROUNDING_H
#define STORMM_ROUNDING_H

#include <cmath>
#include <vector>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

using symbols::pi;
using symbols::pi_f;
using symbols::near_to_one_f;
using symbols::near_to_one_lf;
  
/// \brief Round a number to the increment.  A templated function can deal with each integral
///        type.  Running this function with float or double will just produce the original number
///        plus the increment minus one, all divided by the increment.
///
/// \param jagged     The original number, a jagged cut that needs rounding
/// \param increment  The rounding increment
template <typename T> T roundUp(T jagged, T increment);

/// \brief Return a vector containing the unique prime factors found within a positive integer, out
///        of a list of primes.
///
/// \param number    The number to factorize
/// \param primes    The vector of prime numbers to attempt
/// \param n_primes  The number of unique primes to test as factors (for taking a subset of the
///                  available set, the default of zero indicating that all entries in primes
///                  should be tried)
std::vector<uint> primeFactors(ullint number, const std::vector<uint> &primes, int n_primes = 0);

/// \brief Return a vector containing the number of times each of a series of prime factors appear
///        in a given number.  Descriptions of input parameters follow from primeFactors() above.
std::vector<uint> primeFactorCounts(ullint number, const std::vector<uint> &primes,
                                    int n_primes = 0);

/// \brief Find some small prime factorizations of a given container size and a target memory
///        increment, then return the minimum number of containers which together take up memory
///        equal to a multiple of the increment.
///
/// \param element_size   The container size
/// \param increment      Minimum increment of memory (in bytes) to allocate
/// \param n_primes       [Optional] The number of prime factors to test (max 8)
ulint getSmallestLot(int element_size, int increment, int n_primes = 8);

/// \brief Calculate the nearest multiple of prime factors to some target, approaching from below
///        without going over or from above without going under.  This calculation can be repeated
///        in order to strike a balance between fine tuning the result and using the smallest
///        possible primes.
///
/// \param number    The number for which to develop a factor as close as possible to target
/// \param target    The number to approximate with prime factors of number
/// \param primes    The selection of primes to try in factorizing the number
/// \param approach  Indicate whether to approach the target from above or from below
/// \param n_primes  The extent of elements within the array primes to use (the default of zero
///                  triggers use of the whole array)
ullint nearestFactor(ullint number, ullint target, const std::vector<uint> &primes,
                     LimitApproach approach, int n_primes = 0);
  
/// \brief Find a sensible maximum memory size to accommodate a given number of elements.  The size
///        of each element may be supplied to make the thing work in terms of bytes, otherwise it
///        trims down the amount of padding by the number of existing elements.
///
/// \param length            The current length of the array
/// \param growth_increment  The minimum growth increment of the array
/// \param element_size      Size of each element, in bytes
size_t getPaddedMemorySize(size_t length, size_t growth_increment, size_t element_size = 1);

/// \brief Guard against the ill-conditioned region of the arccos function when its argument is
///        close to 1.0.  This will produce a much more accurate dihedral angle representation
///        when computing torsions and CMAPs, particularly in improper terms when the force
///        constants are larger and the angles tend to occupy the ill-conditioned regions.
///
/// \param costheta  Argument to the arccos function, naively computed as the dot product of
///                  crabbc and crbccd over the product of their magnitudes.
/// \param crabbc    Normal vector to the ABC (atoms I-J-K) plane
/// \param crbccd    Normal vector to the BCD (atoms J-K-L) plane
/// \param bc        B-C (atoms J-K) vector
/// \param scr       Cross product of crabbc and crbccd, pre-computed for prior use in the naive
///                  computation
template <typename Tcalc>
Tcalc angleVerification(const Tcalc costheta, const Tcalc* crabbc, const Tcalc* crbccd,
                        const Tcalc* bc, const Tcalc* scr);

} // namespace stmath
} // namespace stormm

#include "rounding.tpp"

#endif
