#include "../../src/Constants/behavior.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/file_snapshot.h"

using omni::constants::PrecisionModel;
using omni::constants::tiny;
using omni::data_types::llint;
using omni::data_types::float2;
using omni::data_types::int95_t;
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

//-------------------------------------------------------------------------------------------------
// Test the split accumulation method over a segment of the number line.
//
// Arguments:
//   llim:        Low limit for sampling
//   hlim:        High limit for sampling
//   incr:        Incrementation value
//   scale_bits:  The number of bits after the decimal
//-------------------------------------------------------------------------------------------------
void testSplitAccumulation(const double llim, const double hlim, const double incr,
                           const int scale_bits, const PrecisionModel lvl) {
  const double scale_factor = pow(2.0, scale_bits);
  const double inv_scale_factor = 1.0 / scale_factor;
  const double inv_overflow_factor = max_llint_accumulation / scale_factor;
  const float scale_factorf = scale_factor;
  int n_basic_fail = 0;
  for (double r = llim; r < hlim; r += incr) {
    const double scaled_r = r * scale_factor;
    switch (lvl) {
    case PrecisionModel::SINGLE:
      {
        const llint dp_result = scaled_r;
        float fr = r;
        const float scaled_fr = fr * scale_factor;
        const llint fp_result = scaled_fr;
        int overflow, workbin;
        splitRealConversion(scaled_fr, &workbin, &overflow);
        const llint lloverflow = overflow;
        const llint llworkbin  = workbin;
        const llint fp_reconst = (lloverflow * max_int_accumulation_ll) + llworkbin;
        if (fp_reconst != fp_result) {
          n_basic_fail++;
        }
      }
      break;
    case PrecisionModel::DOUBLE:
      {
        llint workbin;
        int overflow;
        splitRealConversion(r * scale_factor, &workbin, &overflow);
        double dp_reconst = static_cast<double>(workbin) * inv_scale_factor;
        dp_reconst += static_cast<double>(overflow) * inv_overflow_factor;
        if (fabs(dp_reconst - r) > tiny) {
          n_basic_fail++;
        }
      }
      break;
    }
  }
  const double ellim = llim * 0.125;
  const double ehlim = hlim * 0.125;
  const double eincr = incr * 0.125;
  int n_inc_fail = 0;
  for (double r = ellim; r < ehlim; r += eincr) {
    switch (lvl) {
    case PrecisionModel::SINGLE:
      {
        llint dp_result = 0LL;
        llint fp_result = 0LL;
        int overflow = 0;
        int workbin = 0;
        for (int i = 0; i < 9; i++) {
          float fr = r;
          double scaled_r = r * scale_factor;
          float scaled_fr = fr * scale_factorf;
          splitRealAccumulation(scaled_fr, &workbin, &overflow);
          dp_result += static_cast<llint>(scaled_r);
          fp_result += static_cast<llint>(scaled_fr);
        }
        const llint lloverflow = overflow;
        const llint llworkbin  = workbin;
        const llint fp_reconst = (lloverflow * max_int_accumulation_ll) + llworkbin;
        if (fp_reconst != fp_result) {
          n_inc_fail++;
        }
      }
      break;
    case PrecisionModel::DOUBLE:
      {
        double dp_tracker = 0.0;
        int overflow = 0;
        llint workbin = 0LL;
        for (int i = 0; i < 9; i++) {
          double scaled_r = r * scale_factor;
          splitRealAccumulation(scaled_r, &workbin, &overflow);
          dp_tracker += r;
        }
        double dp_reconst = static_cast<double>(workbin) * inv_scale_factor;
        dp_reconst += static_cast<double>(overflow) * inv_overflow_factor;
        if (fabs(dp_reconst - dp_tracker) > tiny) {
          n_inc_fail++;
        }
      }
      break;
    }
  }
  check(n_basic_fail == 0, "A total of " + std::to_string(n_basic_fail) + " failures were "
        "recorded converting the range " + realToString(llim, 8, 4) + " : " +
        realToString(hlim, 8, 4) + " to split integer fixed-precision values when sampled with "
        "a " + realToString(incr, 9, 2, NumberFormat::SCIENTIFIC) + " increment.  Precision "
        "model: " + getPrecisionModelName(lvl) + ".");
  check(n_inc_fail == 0, "A total of " + std::to_string(n_inc_fail) + " failures were recorded "
        "converting the range " + realToString(ellim, 8, 4) + " : " + realToString(ehlim, 8, 4) +
        " to split integer fixed-precision values by incrementing 9x when sampled with a " +
        realToString(incr, 9, 2, NumberFormat::SCIENTIFIC) + " increment.  Precision model: " +
        getPrecisionModelName(lvl) + ".");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Long long int <=> float2 conversion");

  // Section 2
  section("Real to int2, 63-bit conversion");

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
  section(1);
  const llint etr_expected_min = -545794316914904LL;  
  const llint etr_expected_max =  545794316914905LL;
  const bool exact_passes = (minValue(exact_test_range) == etr_expected_min &&
                             maxValue(exact_test_range) == etr_expected_max);
  check(exact_passes, "A series of long long integers needed for many tests was not constructed "
        "properly.  Its contents should span a range [" + std::to_string(etr_expected_min) + ", " +
        std::to_string(etr_expected_max) + "], but instead it spans " +
        std::to_string(minValue(exact_test_range)) + " to " +
        std::to_string(maxValue(exact_test_range)) + ".  Dependent tests will be skipped.");
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
  check(large_passes, "A series of long long integers needed for many tests was not constructed "
        "properly.  Its contents should span a range [" + std::to_string(ltr_expected_min) + ", " +
        std::to_string(ltr_expected_max) + "], but instead it spans " +
        std::to_string(minValue(large_test_range)) + " to " +
        std::to_string(maxValue(large_test_range)) + ".  Dependent tests will be skipped.");
  const TestPriority do_large_tests = (exact_passes) ? TestPriority::CRITICAL :
                                                       TestPriority::ABORT;


  // Test direct conversion
  section(2);
  testSplitAccumulation(-48.5, -47.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation(-32.5, -31.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( -0.5,   0.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 31.5,  32.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 47.5,  48.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation(-64.5, -63.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 63.5,  64.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation(-48.5, -47.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation(-32.5, -31.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( -0.5,   0.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 31.5,  32.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 47.5,  48.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation(-64.5, -63.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 63.5,  64.5, 1.0e-5, 58, PrecisionModel::DOUBLE);

  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
