#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Constants/scaling.h"
#include "rounding.h"
#include "tickcounter.h"

namespace stormm {
namespace stmath {

using constants::mega;

//-------------------------------------------------------------------------------------------------
std::vector<uint> primeFactors(const ullint number, const std::vector<uint> &primes,
                               const int n_primes) {
  std::vector<uint> factors;
  ullint residual = number;
  const int actual_prime_count = (n_primes == 0) ? primes.size() : n_primes;
  for (int i = 0; i < actual_prime_count; i++) {
    const ullint ul_pi = primes[i];
    int nfac = 0;
    while ((residual / ul_pi) * ul_pi == residual) {
      nfac++;
      residual /= ul_pi;
    }
    if (nfac > 0) {
      factors.push_back(primes[i]);
    }
  }
  return factors;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> primeFactorCounts(const ullint number, const std::vector<uint> &primes,
                                    const int n_primes) {
  std::vector<uint> factors((n_primes <= 0 || n_primes > primes.size()) ? primes.size() : n_primes,
                            0);
  ullint residual = number;
  const int actual_prime_count = (n_primes == 0) ? factors.size() : n_primes;
  for (int i = 0; i < actual_prime_count; i++) {
    const ullint ul_pi = primes[i];
    while ((residual / ul_pi) * ul_pi == residual) {
      factors[i] += 1;
      residual /= ul_pi;
    }
  }
  return factors;
}

//-------------------------------------------------------------------------------------------------
ulint getSmallestLot(const int element_size, const int increment, const int n_primes) {

  // Small primes
  std::vector<uint> primes = { 2, 3, 5, 7, 11, 13, 17, 19 };
  const int useable_primes = (n_primes > 8) ? 8 : n_primes;

  // Naive factorization of each number
  std::vector<uint> efac = primeFactorCounts(element_size, primes, useable_primes);
  std::vector<uint> ifac = primeFactorCounts(increment, primes, useable_primes);

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
ullint nearestFactor(const ullint number, const ullint target, const std::vector<uint> &primes,
                     const LimitApproach approach, const int n_primes) {
  const std::vector<uint> factors = primeFactorCounts(number, primes, n_primes);
  const int actual_pc = (n_primes <= 0 || n_primes >= primes.size()) ? primes.size() : n_primes;
  std::vector<int> maximum_prime_inclusions(actual_pc);
  int unique_prime_factors = 0;
  for (int i = 0; i < actual_pc; i++) {
    uint ilog_calc;
    switch (approach) {
    case LimitApproach::BELOW:
      ilog_calc = floor(log2(target) / log2(primes[i]));
      break;
    case LimitApproach::ABOVE:
      ilog_calc = ceil(log2(target) / log2(primes[i]));
      break;
    }
    maximum_prime_inclusions[i] = std::min(factors[i], ilog_calc);
    unique_prime_factors += (maximum_prime_inclusions[i] > 0);
  }
  if (unique_prime_factors > 0) {
    std::vector<std::vector<int>> pf_settings(unique_prime_factors);
    unique_prime_factors = 0;
    for (int i = 0; i < actual_pc; i++) {
      if (maximum_prime_inclusions[i] > 0) {
        std::vector<int> tmp_powers(maximum_prime_inclusions[i] + 1);
        int k = 1;
        for (int j = 0; j <= maximum_prime_inclusions[i]; j++) {
          tmp_powers[j] = k;
          k *= primes[i];
        }
        pf_settings[unique_prime_factors] = tmp_powers;
        unique_prime_factors++;
      }
    }
    
    // The problem of finding the best combination of prime factors to get as close as possible
    // to the target number of entries per line is hard, perhaps NP-hard, but the combinatorics
    // is not excessively large for a small target number.
    TickCounter<int> dials(pf_settings);
    const std::vector<int>& dial_settings = dials.getSettings();
    int ncombo = (dials.getLogPermutations() > 3.0) ? 1000 : dials.getExactPermutationCount();
    int best_approx = 1;
    for (int i = 0; i < ncombo; i++) {
      int trial_approx = 1;
      for (int j = 0; j < unique_prime_factors; j++) {
        trial_approx *= dials.getState(j);
      }
      switch (approach) {
      case LimitApproach::BELOW:
        if (trial_approx <= target && target - trial_approx < target - best_approx) {
          best_approx = trial_approx;
        }
        break;
      case LimitApproach::ABOVE:
        if (trial_approx >= target && trial_approx - target < best_approx - target) {
          best_approx = trial_approx;
        }
        break;
      }
      dials.advance();
    }
    return best_approx;
  }
  else {
    return number;
  }
  __builtin_unreachable();
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

} // namespace stmath
} // namespace stormm
