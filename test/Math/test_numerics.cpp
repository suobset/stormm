#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/file_snapshot.h"

using omni::data_types::llint;
using omni::data_types::float2;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::parse::NumberFormat;
using omni::parse::polyNumericVector;
using omni::random::Ran2Generator;
using omni::symbols::pi;
using namespace omni::math;
using namespace omni::numerics;
using namespace omni::testing;

int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Long long int <=> float2 conversion");

  // Section 2
  section("Double <=> float2 conversion");

  // Section 3
  section("Compariosn of precision formats for unit cell operations");

  // Make a series of prime numbers
  std::vector<llint> primes(1, 2);
  int p = 3;
  int g = 3;
  while (g < 10000) {
    const int q = g / p;
    if (g - (p * q) == 0) {
      g += 2;
      p = 3;
      continue;
    }
    if (p > q) {
      primes.push_back(g);
      g += 2;
      p = 3;
    }
    else {
      p += 2;
    }
  }

  // Make a series of integers covering a large portion of the long long int range
  const int n_pts = 16277216;
  const int half_pts = n_pts / 2;
  std::vector<llint> exact_test_range(n_pts);
  llint acc = 1LL;
  llint inc = 2LL;
  int pcon = 0;
  const int n_primes = primes.size();
  for (int i = half_pts; i < n_pts; i++) {
    exact_test_range[i] = acc;
    exact_test_range[n_pts - i - 1] = -acc + 1LL;
    acc += inc;
    inc = (inc > 134217728LL) ? primes[pcon] : inc + primes[pcon];
    pcon++;
    if (pcon == n_primes) {
      pcon = 0;
    }
  }

  // Check the condition of the number series
  const llint etr_expected_min = -545794316914904LL;  
  const llint etr_expected_max =  545794316914905LL;
  const bool exact_passes = (minValue(exact_test_range) == etr_expected_min &&
                             maxValue(exact_test_range) == etr_expected_max);
  if (exact_passes == false) {
    rtWarn("A series of long long integers needed for many tests was not constructed properly.  "
           "Its contents should span a range [" + std::to_string(etr_expected_min) + ", " +
           std::to_string(etr_expected_max) + "], but instead it spans " +
           std::to_string(minValue(exact_test_range)) + " to " +
           std::to_string(maxValue(exact_test_range)) + ".  Dependent tests will be skipped.",
           "test_numerics");
  }
  const TestPriority do_exact_tests = (exact_passes) ? TestPriority::CRITICAL :
                                                       TestPriority::ABORT;
  
  // Compute for a slightly larger series of numbers
  std::vector<llint> large_test_range(n_pts);
  acc = 1LL;
  inc = 2LL;
  pcon = 0;
  for (int i = half_pts; i < n_pts; i++) {
    large_test_range[i] = acc;
    large_test_range[n_pts - i - 1] = -acc + 1LL;
    acc += inc;
    inc = (inc > 260000000LL) ? primes[pcon] : inc + primes[pcon];
    pcon++;
    if (pcon == n_primes) {
      pcon = 0;
    }
  }
  const llint ltr_expected_max =  1056755308863635LL;
  const llint ltr_expected_min = -1056755308863634LL;
  const bool large_passes = (minValue(large_test_range) == ltr_expected_min &&
                             maxValue(large_test_range) == ltr_expected_max);
  if (large_passes == false) {
    rtWarn("A series of long long integers needed for many tests was not constructed properly.  "
           "Its contents should span a range [" + std::to_string(ltr_expected_min) + ", " +
           std::to_string(ltr_expected_max) + "], but instead it spans " +
           std::to_string(minValue(large_test_range)) + " to " +
           std::to_string(maxValue(large_test_range)) + ".  Dependent tests will be skipped.",
           "test_numerics");
  }
  const TestPriority do_large_tests = (exact_passes) ? TestPriority::CRITICAL :
                                                       TestPriority::ABORT;

  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
