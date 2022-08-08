#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Constants/scaling.h"
#include "rounding.h"

namespace stormm {
namespace math {

using constants::mega;

//-------------------------------------------------------------------------------------------------
std::vector<ulint> primeFactors(const ulint number, const std::vector<ulint> &primes,
                                const int n_primes) {
  std::vector<ulint> factors(n_primes, 0);
  int residual = number;
  for (int i = 0; i < n_primes; i++) {
    while ((residual / primes[i]) * primes[i] == residual) {
      factors[i] += 1;
      residual /= primes[i];
    }
  }
  return factors;
}

//-------------------------------------------------------------------------------------------------
ulint getSmallestLot(const int element_size, const int increment, const int n_primes) {

  // Small primes
  std::vector<ulint> primes = { 2, 3, 5, 7, 11, 13, 17, 19 };
  const int useable_primes = (n_primes > 8) ? 8 : n_primes;

  // Naive factorization of each number
  std::vector<ulint> efac = primeFactors(element_size, primes, useable_primes);
  std::vector<ulint> ifac = primeFactors(increment, primes, useable_primes);

  // Find the lowest common multiple
  ulint common_multiple = 1;
  for (int i = 0; i < useable_primes; i++) {
    for (int j = 0; j < std::min(efac[i], ifac[i]); j++) {
      common_multiple *= primes[i];
    }
  }

  // Divide the increment by the common multiple, then multiply the element size by the result
  return increment / common_multiple;
}

//-------------------------------------------------------------------------------------------------
size_t getPaddedMemorySize(const size_t length, const size_t growth_increment,
                           const size_t element_size) {

  // The actual growth increment
  const size_t incr = growth_increment * element_size;

  // Handle the case of zero
  if (length == 0) {
    return incr;
  }
  else if (length < 4 * incr) {
    return roundUp(length, incr);
  }
  else if (length < 8 * incr) {
    return roundUp(length, 2 * incr);
  }
  else if (length < 4 * constants::mega) {
    return roundUp(2 * length, incr);
  }
  else if (length < 16 * constants::mega) {
    return roundUp((3 * length) / 2, incr);
  }
  else if (length < 64 * constants::mega) {
    return roundUp((5 * length) / 4, incr);
  }
  else {
    return roundUp((9 * length) / 8, incr);
  }
  __builtin_unreachable();
}

} // namespace math
} // namespace stormm
