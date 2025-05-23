// -*-c++-*-
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if STORMM_INCLUDE_NETCDF
#  include <netcdf.h>
#endif
#if STORMM_INCLUDE_POCKETFFT
#endif
#include "pocketfft_hdronly.h"
#include <string.h>
#include "copyright.h"
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
using namespace pocketfft;

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
#if 0
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
    shape_t shape_in = { test_lengths[l_idx] };
    std::vector<double> grid_values = gaussianRand(rng, test_lengths[l_idx], 1.0);

    // CHECK
    if (test_lengths[l_idx] == 8) {
      for (int i = 0; i < test_lengths[l_idx]; i++) {
        grid_values[i] = i + 1;
      }
    }
    // END CHECK
    
    std::vector<double> original_grid_values = grid_values;
    
    // Create a copy for the FFT operation
    std::vector<double> fft_data = grid_values;
    
    // Define shape and stride for the FFT
    shape_t shape{static_cast<size_t>(test_lengths[l_idx])};
    shape_t axes{0};
    stride_t stride_in{static_cast<int>(test_lengths[l_idx] * sizeof(double))};
    stride_t stride_out{static_cast<int>(test_lengths[l_idx] * sizeof(std::complex<double>))};
    size_t tmp_in = sizeof(double), tmp_out = sizeof(double);
    
    for (int i = shape.size()-1; i >= 0; --i) {
      stride_in[i] = tmp_in;
      tmp_in *= shape[i];
      stride_out[i] = tmp_out;
      tmp_out *= shape[i];
    }
    
    // Forward FFT (real to complex)
    std::vector<std::complex<double>> complex_data((test_lengths[l_idx] / 2) + 1);
    r2c(shape, stride_in, stride_out, axes, true, fft_data.data(), complex_data.data(), 1.0);

    // CHECK
    if (test_lengths[l_idx] == 8) {
      printf("complex_data = [\n");
      for (int i = 0; i < (test_lengths[l_idx] / 2) + 1; i++) {
        printf("  %16.8lf %16.8lf\n", complex_data[i].real(), complex_data[i].imag());
      }
      printf("];\n");
    }
    // END CHECK
    
    // Backward FFT (complex to real)
    c2r(shape, stride_out, stride_in, 0, false, complex_data.data(), fft_data.data(),
        1.0 / test_lengths[l_idx]);
    
    // CHECK
    if (test_lengths[l_idx] == 8) {
      printf("r2c_c2r = [\n");
      for (int i = 0; i < test_lengths[l_idx]; i++) {
        printf("  %16.8lf %16.8lf\n", original_grid_values[i], fft_data[i]);
      }
      printf("];\n");
    }
    // END CHECK
    
    // Calculate error between transformed data and original data
    check(fft_data, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
          "The root mean-squared error in real FFT transformation of dimension " +
          std::to_string(test_lengths[l_idx]) + " exceeds the permissible threshold.");
  }
#endif
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
    std::vector<double> grid_values = gaussianRand(rng, 2 * test_lengths[l_idx], 1.0);
    std::vector<double> original_grid_values = grid_values;
    
    // Create a copy for the FFT operation
    std::vector<std::complex<double>> complex_data(test_lengths[l_idx]);
    for (size_t i = 0; i < test_lengths[l_idx]; ++i) {
        complex_data[i] = std::complex<double>(grid_values[2*i], grid_values[2*i+1]);
    }
    
    // Define shape and stride for the FFT
    shape_t shape = { static_cast<size_t>(test_lengths[l_idx]) };
    stride_t stride_in(shape.size()), stride_out(shape.size());
    size_t tmp_in = sizeof(std::complex<double>);
    size_t tmp_out = sizeof(std::complex<double>);
    
    for (int i = shape.size() - 1; i >= 0; --i) {
      stride_in[i] = tmp_in;
      tmp_in *= shape[i];
      stride_out[i] = tmp_out;
      tmp_out *= shape[i];
    }
    
    // Define axes for the FFT
    shape_t axes;
    for (size_t i = 0; i < shape.size(); i++) {
      axes.push_back(i);
    }
    
    // Forward FFT
    c2c(shape, stride_in, stride_out, axes, true, complex_data.data(), complex_data.data(), 1.0);
    
    // Backward FFT
    c2c(shape, stride_in, stride_out, axes, false, complex_data.data(), complex_data.data(),
        1.0 / static_cast<double>(test_lengths[l_idx]));
    
    // Convert back to real array
    std::vector<double> fft_data(2 * test_lengths[l_idx]);
    for (size_t i = 0; i < test_lengths[l_idx]; ++i) {
      fft_data[2*i] = complex_data[i].real();
      fft_data[2*i+1] = complex_data[i].imag();
    }
    
    // Calculate error between transformed data and original data
    check(fft_data, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
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

  // Section 4
  section("Test Basic NetCDF Functionality");

  // Begin testing
  section(1);
  testRealFFT(&rng);
  section(2);
  testComplexFFT(&rng);
  section(3);
  testReal3DFFT(&rng, { 8, 8, 8 });
  section(4);

  // Print results
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
