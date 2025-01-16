// -*-c++-*-
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

extern "C" {
#include <pocketfft.h>
}
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/summation.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Random/random.h"
#include "../../src/Random/hpc_random.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Tests the real FFT (RFFT) functionality by generating random data, performing forward and
// backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   maxlen:   The maximum length of the real FFT to test (several values in the range will be
//             tested);
//   epsilon:  Any given point in the FFT grid must return to its oiginal value after the forward
//             FFT, reverse FFT, and normalization cycle.  This is the tolerance for judging
//             whether that condition is met.
//-------------------------------------------------------------------------------------------------
void testRealFFT(Xoshiro256ppGenerator *rng, const int maxlen = 8192,
                 const double epsilon = 5.0e-10) {

  // Perform FFT test for different lengths of the data
  std::vector<int> test_lengths;
  int tl = 1;
  while (tl < maxlen) {
    test_lengths.push_back(tl);
    if (tl < 8) {
      tl++;
    }
    else if (tl < maxlen / 8) {
      tl = maxlen / 8;
    }
    else {
      tl += maxlen / 8;
    }
  }
  for (size_t l_idx = 0; l_idx < test_lengths.size(); l_idx++) {
    std::vector<double> grid_values  = gaussianRand(rng, test_lengths[l_idx], 1.0);
    std::vector<double> original_grid_values = grid_values;

    // Create and apply the forward and backward real FFT (RFFT) plan
    rfft_plan plan = make_rfft_plan(test_lengths[l_idx]);

    // Forward FFT
    rfft_forward(plan, grid_values.data(), 1.0);
    
    // Backward FFT (normalize by length)  
    rfft_backward(plan, grid_values.data(), 1.0 / static_cast<double>(test_lengths[l_idx]));

    // Clean up the FFT plan
    destroy_rfft_plan(plan);

    // Calculate error between transformed data and original data
    check(grid_values, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
          "The root mean-squared error in real FFT transformation of dimension " +
          std::to_string(test_lengths[l_idx]) + " exceeds the permissible threshold.");
  }
}

//-------------------------------------------------------------------------------------------------
// Tests the complex FFT (CFFT) functionality by generating random complex data, performing forward
// and backward FFT operations, and comparing the result with the original complex data.
// Descriptions of input arguments follow from testRealFFT(), above.
//-------------------------------------------------------------------------------------------------
void testComplexFFT(Xoshiro256ppGenerator* rng, const int maxlen = 8192,
                    const double epsilon = 5.0e-10) {

  // Perform FFT test for different lengths of the data
  std::vector<int> test_lengths;
  int tl = 1;
  while (tl < maxlen) {
    test_lengths.push_back(tl);
    if (tl < 4) {
      tl++;
    }
    else if (tl < maxlen / 8) {
      tl = maxlen / 8;
    }
    else {
      tl += maxlen / 8;
    }
  }
  for (size_t l_idx = 0; l_idx < test_lengths.size(); l_idx++) {
    std::vector<double> grid_values  = gaussianRand(rng, 2 * test_lengths[l_idx], 1.0);
    std::vector<double> original_grid_values = grid_values;

    // Create and apply the forward and backward real FFT (RFFT) plan
    cfft_plan plan = make_cfft_plan(test_lengths[l_idx]);

    // Forward FFT
    cfft_forward(plan, grid_values.data(), 1.0);
    
    // Backward FFT (normalize by length)
    cfft_backward(plan, grid_values.data(), 1.0 / static_cast<double>(test_lengths[l_idx]));

    // Clean up the FFT plan
    destroy_cfft_plan(plan);

    // Calculate error between transformed data and original data
    check(grid_values, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
          "The root mean-squared error in complex FFT transformation of dimension " +
          std::to_string(test_lengths[l_idx]) + " exceeds the permissible threshold.");
  }
}

//-------------------------------------------------------------------------------------------------
// Test the real FFT (RFFT) functionality in a three-dimensional array by generating random data,
// performing forward and backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   dims:     Dimensions of the transform to set up and perform
//   epsilon:  Tolerance for deviations in the input and output arguments
//-------------------------------------------------------------------------------------------------
void testReal3DFFT(Xoshiro256ppGenerator *rng, const std::vector<int> &dims,
                   const double epsilon = 5.0e-10) {
  if (dims.size() != 3) {
    rtErr("A three-dimensional array is expected (" + std::to_string(dims.size()) +
          " dimensions were provided).", "testReal3DFFT");
  }
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) {
    rtErr("Positive dimensions for the grid must be specified (" + std::to_string(dims[0]) + ", " +
          std::to_string(dims[1]) + ", and " + std::to_string(dims[2]) + " were provided).",
          "testReal3DFFT");
  }
  std::vector<double> grid_values = gaussianRand(rng, dims[0] * dims[1] * dims[2], 1.0);
  std::vector<double> original_grid_values = grid_values;

  // Calculate error between transformed data and original data
  check(grid_values, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
        "The root mean-squared error in complex FFT transformation of dimension " +
        std::to_string(dims[0]) + " x " + std::to_string(dims[1]) + " x " +
        std::to_string(dims[2]) + " exceeds the permissible threshold.");
}

//-------------------------------------------------------------------------------------------------
// Main function to run the FFT tests (both real and complex) and report results.
//
// Arguments:
//   argc: The number of command-line arguments.
//   argv: The array of command-line arguments (unused in this case).
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Initialize STORMM test environment
  TestEnvironment oe(argc, argv);
  Xoshiro256ppGenerator rng(7183529);  // Using STORMM's random number generator

  // Section 1
  section("Running a real FFT test in PocketFFT");

  // Section 2
  section("Running a complex FFT in PocketFFT");

  // Section 3
  section("Test real-to-(half) complex 3D FFTs in PocketFFT");
  
  section(1);
  testRealFFT(&rng);
  section(2);
  testComplexFFT(&rng);
  section(3);
  testReal3DFFT(&rng, { 8, 8, 8 });
  
  // Print results
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
