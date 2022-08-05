#include <algorithm>
#include "../../src/Constants/symbol_values.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/rounding.h"
#include "../../src/Math/sorting.h"
#include "../../src/Math/summation.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/file_snapshot.h"

using stormm::double3;
using stormm::ulint;
using stormm::ullint;
using stormm::ullint2;
using stormm::constants::tiny;
using stormm::card::Hybrid;
using stormm::card::HybridTargetLevel;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::parse::TextFile;
using stormm::parse::polyNumericVector;
using stormm::random::Ran2Generator;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::RandomNumberMill;
using stormm::random::Xoshiro256ppGenerator;
using stormm::symbols::pi;
using namespace stormm::math;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Produce a series of random numbers from a generator based on an initial state and some
// directions about how to cycle the generator after each sample.
//-------------------------------------------------------------------------------------------------
template <typename Tgen, typename Tstate>
std::vector<double> generateExpectedSeries(Tgen *xrs, const Tstate init_state, const int gen_idx,
                                           const int ngen, const int depth, const int maxcyc) {

  // Initialize the pseudo-random number generator
  xrs->setState(init_state);
  for (int i = 0; i < gen_idx; i++) {
    xrs->longJump();
  }
  
  // Create a buffer to represent the series that should be held within the generator series
  std::vector<double> buffer(depth);
  for (int i = 0; i < depth; i++) {
    buffer[i] = xrs->gaussianRandomNumber();
  }

  // Allocate the result
  std::vector<double> result(maxcyc * depth);
  
  // Compute the refresh schedule
  const int rstride = (ngen + depth - 1) / depth;
  const int rsched  = gen_idx / rstride;
  for (int i = 0; i < maxcyc * depth; i++) {

    // Pull a value from the table
    const int table_idx = i - ((i / depth) * depth);
    result[i] = buffer[table_idx];
    
    // Refresh the table as appropriate
    if (table_idx == rsched) {
      for (int j = 0; j < depth; j++) {
        buffer[j] = xrs->gaussianRandomNumber();
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Vector processing capabilities");

  // Section 2
  section("Random number generation");

  // Section 3
  section("Rounding and prime factorization");

  // Section 4
  section("Lightweight matrix math from the hand-coded routines");
  
  // Check vector processing capabilities
  section(1);
  const std::vector<double> dv_i = { 0.1, 0.5, 0.9, 0.7, 0.8 };
  check(mean(dv_i), RelationalOperator::EQUAL, Approx(0.6).margin(1.0e-12), "The mean value of a "
        "simple real number vector is incorrect.");
  const std::vector<int> iv_i = { 1, 5, 9, 7, 38 };
  check(mean(iv_i), RelationalOperator::EQUAL, 12, "The mean value of a "
        "simple integer vector is incorrect.");
  const std::vector<double> dv_i_long = { -0.3, -0.6, 0.1, 0.5, 0.9, 0.7, 0.8 };
  CHECK_THROWS(maxAbsoluteDifference(dv_i_long, dv_i), "Vectors of different lengths were "
               "compared to obtain a maximum difference.");
  const std::vector<double> dv_ii = { 2.1, 3.5, -9.9, 4.7, 7.8 };
  check(maxAbsoluteDifference(dv_i, dv_ii), RelationalOperator::EQUAL, Approx(10.8).margin(1.0e-8),
        "Incorrect evaluation of the maximum absolute difference between two short vectors.");
  const std::vector<int> iv_ii = { -1, 15, -9, 72, -3 };
  check(maxAbsoluteDifference(iv_i, iv_ii), RelationalOperator::EQUAL, Approx(65).margin(1.0e-10),
        "Incorrect evaluation of the maximum absolute difference between two short vectors.");
  const std::vector<int> uv_i = { 92, 47, 81, 36, 29 };
  const std::vector<int> uv_ii = { 82, 37, 101, 19, 147 };
  check(maxAbsoluteDifference(uv_i, uv_ii), RelationalOperator::EQUAL,
        maxAbsoluteDifference(uv_ii, uv_i), "Non-reflexive evaluation of the maximum absolute "
        "difference between two vectors of unsigned integers.");
  check(maxRelativeDifference(dv_ii, dv_i), RelationalOperator::EQUAL, Approx(20.0).margin(1.0e-8),
        "Incorrect evaluation of maximum relative error between two short vectors.");
  check(maxRelativeDifference(uv_ii, uv_i), RelationalOperator::EQUAL,
        Approx(4.06896552).margin(1.0e-7), "Incorrect evaluation of maximum relative error "
        "between two short unsigned integer vectors.");
  check(mean(uv_i), RelationalOperator::EQUAL, 57.0, "Incorrect evaluation of the mean of a "
        "vector of unsigned integers.");
  const std::vector<double> cp_a = {  1.4, 5.8, -8.5 };
  const std::vector<double> cp_b = { -3.4, 2.7, 10.1 };
  std::vector<double> cp_c(3);
  crossProduct(cp_a, cp_b, &cp_c);
  const std::vector<double> cp_ab = { 81.53, 14.76, 23.50 };
  check(cp_c, RelationalOperator::EQUAL, cp_ab, "Vector cross product does not function as "
        "expected when provided Standard Template Library vectors.");
  const std::vector<float> cp_af = {  1.4, 5.8, -8.5 };
  const std::vector<float> cp_bf = { -3.4, 2.7, 10.1 };
  std::vector<float> cp_cf(3);
  crossProduct(cp_af.data(), cp_bf.data(), cp_cf.data());
  check(cp_cf, RelationalOperator::EQUAL, Approx(cp_ab).margin(1.0e-5), "Vector cross product "
        "does not function as expected when provided C-style pointers.");
  const double3 d3_a = { cp_a[0],  cp_a[1],  cp_a[2] };
  const double3 d3_b = { cp_b[0],  cp_b[1],  cp_b[2] };
  const double3 d3_c = crossProduct(d3_a, d3_b);
  check(d3_c.x == Approx(cp_ab[0]).margin(tiny) && d3_c.y == Approx(cp_ab[1]).margin(tiny) &&
        d3_c.z == Approx(cp_ab[2]).margin(tiny), "Vector cross product does not function as "
        "expected when provided double-precision HPC tuples.");
  const std::vector<int> i_searchable_asc = {  0,  1,  2,  3,  4,  5, 16, 17, 28, 99 };
  const std::vector<int> i_searchable_des = {  9,  8,  6,  5,  4,  2,  1, -3, -4, -6 };
  const std::vector<int> i_searchable_non = {  9, 18, 46,  5, -4,  2,  1, -3, 15, -8 };
  std::vector<int> search_order = { 0, 4, 8, 2, 6, 1, 5, 3, 7, 9 };
  const int nsearch = i_searchable_asc.size();
  std::vector<int> asc_result(nsearch), des_result(nsearch), non_result(nsearch);
  for (int i = 0; i < nsearch; i++) {
    asc_result[i] = locateValue(i_searchable_asc, i_searchable_asc[search_order[i]],
                                DataOrder::ASCENDING);
    des_result[i] = locateValue(i_searchable_des, i_searchable_des[search_order[i]],
                                DataOrder::DESCENDING);
    non_result[i] = locateValue(i_searchable_non, i_searchable_non[search_order[i]],
                                DataOrder::NONE);
  }
  check(asc_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "ascending order yielded incorrect results.");
  check(des_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "descending order yielded incorrect results.");
  check(non_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "no discernible order yielded incorrect results.");
  for (int i = 0; i < nsearch - 1; i++) {
    asc_result[i] = locateValue(i_searchable_asc.data(), i_searchable_asc[search_order[i]],
                                nsearch - 1, DataOrder::ASCENDING);
    des_result[i] = locateValue(i_searchable_des.data(), i_searchable_des[search_order[i]],
                                nsearch - 1, DataOrder::DESCENDING);
    non_result[i] = locateValue(i_searchable_non.data(), i_searchable_non[search_order[i]],
                                nsearch - 1, DataOrder::NONE);
  }
  asc_result.resize(nsearch - 1);
  des_result.resize(nsearch - 1);
  non_result.resize(nsearch - 1);
  search_order.resize(nsearch - 1);
  check(asc_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) + ", arranged in ascending order, yielded incorrect results.");
  check(des_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) +
        ", arranged in descending order, yielded incorrect results.");
  check(non_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) +
        ", arranged in no discernible order, yielded incorrect results.");
  const std::vector<double> d_searchable_asc = { 0.1, 0.7, 2.3, 6.8, 7.1, 7.2, 9.3, 9.5, 9.7 };
  const std::vector<double> d_searchable_des = { 6.1, 5.7, 4.3, 3.8, 3.1, 2.2, 1.3, 0.5, 0.2 };
  const std::vector<double> d_searchable_non = { 5.1, 1.7, 4.3, 9.8, 3.1, 0.2, 0.3, 7.5, 1.8 };
  check(locateValue(Approx(2.35).margin(0.06), d_searchable_asc, DataOrder::ASCENDING),
        RelationalOperator::EQUAL, 2, "Binary search in ascending order with an approximate "
        "target value that should hit its mark from above fails.");
  check(locateValue(Approx(7.05).margin(0.06), d_searchable_asc, DataOrder::ASCENDING),
        RelationalOperator::EQUAL, 4, "Binary search in ascending order with an approximate "
        "target value that should hit its mark from below fails.");
  check(locateValue(Approx(0.55).margin(0.06), d_searchable_des, DataOrder::DESCENDING),
        RelationalOperator::EQUAL, 7, "Binary search in descending order with an approximate "
        "target value that should hit its mark from above fails.");
  check(locateValue(Approx(3.75).margin(0.06), d_searchable_des, DataOrder::DESCENDING),
        RelationalOperator::EQUAL, 3, "Binary search in descending order with an approximate "
        "target value that should hit its mark from below fails.");
  
  // Verify that the internal random number generation is consistent with expectations.  Create
  // three generators, the first two initialized in the same way (which thus should track one
  // another precisely) and the third initialized to a different value (which should decorrelate
  // from the first two immediately thanks to some priming that is done in the initialization).
  section(2);
  check(oe.getRandomSeed(), RelationalOperator::EQUAL, 827493, "Pseudo-random number seed is not "
        "set to the expected value.  While permitted as part of the baseline test environment "
        "initialization, changing the random seed will, in this case, result in several bad "
        "results in later tests.");
  const int n_pts = 2000;
  Ran2Generator prng_a(oe.getRandomSeed());
  Ran2Generator prng_b(oe.getRandomSeed());
  Ran2Generator prng_c(71277);
  std::vector<double> result_a(n_pts, 0.0);
  std::vector<double> result_b(n_pts, 0.0);
  std::vector<double> result_c(n_pts, 0.0);
  for (int i = 0; i < n_pts; i++) {
    result_a[i] = prng_a.uniformRandomNumber();
    result_b[i] = prng_b.uniformRandomNumber();
    result_c[i] = prng_c.gaussianRandomNumber();
  }
  check(result_a, RelationalOperator::EQUAL,
        Approx(result_b, ComparisonType::MEAN_UNSIGNED_ERROR).margin(1.0e-12), "Differences were "
        "detected between random number series created starting from the same seed.");
  check(mean(result_a), RelationalOperator::EQUAL, Approx(0.49674889).margin(1.0e-7), "The mean "
	"of a set of " + std::to_string(n_pts) + " random numbers is incorrect.");
  check(mean(result_c), RelationalOperator::EQUAL, Approx(-0.01931485).margin(1.0e-7), "The mean "
        "value of a normal distribution of random numbers is incorrect.");

  // Additional checks, using the file reference system
  const std::string randoms_snp = oe.getStormmHomePath() + osSeparator() + "test" + osSeparator() +
                                  "Math" + osSeparator() + "randoms.m";
  TestPriority snp_found = (getDrivePathType(randoms_snp) == DrivePathType::FILE) ?
                           TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_found == TestPriority::ABORT && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot file " + randoms_snp + " was not found.  Check the $STORMM_SOURCE environment "
           "variable and make sure it indicates the root source directory where src/ and test/ "
           "can be found.  Some subsequent tests will be skipped.", "test_math");
  }
  snapshot(oe.getStormmHomePath() + osSeparator() + "test" + osSeparator() + "Math" + osSeparator() +
           "randoms.m", polyNumericVector(result_c), "rngvec", 1.0e-4, "Series of random numbers "
           "created by the ran2 method does not conform to expectations.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE, snp_found);
  
  // Check scrambled linear random number generators
  Xoroshiro128pGenerator xrs128p(798031);
  Xoshiro256ppGenerator xrs256pp(901835);
  Xoroshiro128pGenerator xrs128pj(798031);
  Xoshiro256ppGenerator xrs256ppj(901835);
  xrs128pj.jump();
  xrs256ppj.jump();
  const int nrand_trial = 16;
  std::vector<double> xrs128p_result(nrand_trial);
  std::vector<double> xrs128p_jump_result(nrand_trial);
  std::vector<double> xrs256pp_result(nrand_trial);
  std::vector<double> xrs256pp_jump_result(nrand_trial);
  for (int i = 0; i < nrand_trial; i++) {
    xrs128p_result[i] = xrs128p.uniformRandomNumber();
    xrs256pp_result[i] = xrs256pp.uniformRandomNumber();
    xrs128pj.uniformRandomNumber();
    xrs256ppj.uniformRandomNumber();
  }
  xrs128p.jump();
  xrs256pp.jump();
  const std::vector<double> ans_128p = { 0.7453648164, 0.4049923254, 0.8584963726, 0.9833355535,
                                         0.0066062865, 0.6311114017, 0.9820136114, 0.2413733841,
                                         0.2753459418, 0.4993040685, 0.0806123499, 0.3691566725,
                                         0.4001401073, 0.0590209187, 0.5804605659, 0.5293153466 };
  const std::vector<double> ans_256pp = { 0.9570185700, 0.9443435021, 0.5241518529, 0.1428242166,
                                          0.5687186981, 0.9839369524, 0.9751010737, 0.8695307120,
                                          0.4293373290, 0.1555431764, 0.0005355056, 0.0820197480,
                                          0.8509827801, 0.7310430980, 0.7770864125, 0.9021266541 };
  check(xrs128p_result, RelationalOperator::EQUAL, Approx(ans_128p).margin(1.0e-8),
        "Random numbers generated by the xoroshift128+ method do not meet expectations.");
  check(xrs256pp_result, RelationalOperator::EQUAL, Approx(ans_256pp).margin(1.0e-8),
        "Random numbers generated by the xoshift256++ method do not meet expectations.");
  for (int i = 0; i < nrand_trial; i++) {
    xrs128p_result[i] = xrs128p.uniformRandomNumber();
    xrs128p_jump_result[i] = xrs128pj.uniformRandomNumber();
    xrs256pp_result[i] = xrs256pp.uniformRandomNumber();
    xrs256pp_jump_result[i] = xrs256ppj.uniformRandomNumber();
  }
  check(xrs128p_result, RelationalOperator::EQUAL, xrs128p_jump_result, "Two xoroshiro128+ "
        "generators do not re-synchronize after different combinations of random bit string "
        "generation and a jump.");
  check(xrs256pp_result, RelationalOperator::EQUAL, xrs256pp_jump_result, "Two xoroshiro256++ "
        "generators do not re-synchronize after different combinations of random bit string "
        "generation and a jump.");
  Xoshiro256ppGenerator xrs256ppb({ 0x7743a154e17a5e9bLLU, 0x7823a1cd9453899bLLU,
                                    0x976589eefbb1c7f5LLU, 0x702cf168260fa29eLLU });
  const std::vector<ullint> orbit = { 0xd5c766557e6e16e4LLU, 0xf8eb8f8747a8cc67LLU,
                                      0xfc18365710a653eeLLU, 0xc698a193593f232LLU,
                                      0xa44ddeaac93b1be7LLU, 0x678dcd0e0516c741LLU,
                                      0x351d94668d7c35eLLU,  0x73160e24fc8768daLLU,
                                      0x562ca11220c31698LLU, 0x3dba336235c48913LLU };
  std::vector<ullint> xrsb_result(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsb_result[i] = xrs256ppb.revealBitString();
    xrs256ppb.uniformRandomNumber();
  }
  std::vector<int> xrsb_diff(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsb_diff[i] = xrsb_result[i] - orbit[i];
  }
  check(xrsb_diff, RelationalOperator::EQUAL, std::vector<ullint>(orbit.size(), 0LLU),
        "The Xoshiro256++ random number generator did not produce the expected bit strings when "
        "initialized with a specific state.");
  Xoshiro256ppGenerator xrs256ppc({ 0x7743a154e17a5e9bLLU, 0x7823a1cd9453899bLLU,
                                    0x976589eefbb1c7f5LLU, 0x702cf168260fa29eLLU });
  xrs256ppc.jump();
  const std::vector<ullint> orbit_ii = { 0x39c4396d8759c874LLU, 0x4b948d9de69752ecLLU,
                                         0x871591604b03d9a6LLU, 0x444d6d471322d17bLLU,
                                         0xb0a9eb9383bf0803LLU, 0x481f6c796c1d0ecaLLU,
                                         0xb89a346b480341bfLLU, 0x1494bad1d1b19126LLU,
                                         0xa2f5ca0a0ab0805LLU,  0x75a4de1da308cc8fLLU };
  std::vector<ullint> xrsc_result(orbit_ii.size());
  for (size_t i = 0; i < orbit_ii.size(); i++) {
    xrsc_result[i] = xrs256ppc.revealBitString();
    xrs256ppc.uniformRandomNumber();
  }
  std::vector<int> xrsc_diff(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsc_diff[i] = xrsc_result[i] - orbit_ii[i];
  }
  check(xrsc_diff, RelationalOperator::EQUAL, std::vector<ullint>(orbit_ii.size(), 0LLU),
        "The Xoshiro256++ random number generator did not produce the expected bit strings after "
        "traversing a jump.");
  
  // Verify rounding results
  section(3);
  const size_t szt_a = 159283;
  const size_t szt_a_round = roundUp<size_t>(szt_a, 32);
  check(szt_a_round, RelationalOperator::EQUAL, Approx(159296).margin(1.0e-6), "Rounding upwards "
        "to the nearest increment of 32 failed.");
  const int szt_b_round = roundUp<size_t>(szt_a, 128);
  check(szt_b_round, RelationalOperator::EQUAL, Approx(159360).margin(1.0e-6), "Rounding upwards "
        "to the nearest increment of 128 failed.");
  const ulint prime_composite = 2 * 2 * 5 * 3 * 7 * 19;
  const std::vector<ulint> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
  const std::vector<ulint> factors_i = primeFactors(prime_composite, primes, 7);
  const std::vector<ulint> factors_ii = primeFactors(prime_composite, primes, 9);
  const std::vector<ulint> answer_i = { 2, 1, 1, 1, 0, 0, 0 };
  const std::vector<ulint> answer_ii = { 2, 1, 1, 1, 0, 0, 0, 1, 0 };
  check(factors_i, RelationalOperator::EQUAL, Approx(answer_i).margin(1.0e-8), "Prime "
        "factorization with numbers up to 7 fails.");

  // Check matrix math from the lightweight, onboard library
  section(4);
  const int rows_a = 7;
  const int cols_a = 5;
  const int rows_b = 5;
  const int cols_b = 7;
  Hybrid<double> tr_mat_a(rows_a * cols_a, "random matrix");
  Hybrid<double> tr_mat_b(rows_b * cols_b, "random matrix");
  Hybrid<double> tr_mat_c(rows_a * cols_b, "posdef matrix");
  Hybrid<double> tr_mat_d(cols_a * rows_b, "posdef matrix");
  double* tr_a_ptr = tr_mat_a.data();
  double* tr_b_ptr = tr_mat_b.data();
  for (int i = 0; i < rows_a * cols_a; i++) {
    tr_a_ptr[i] = prng_c.gaussianRandomNumber();
  }
  for (int i = 0; i < rows_b * cols_b; i++) {
    tr_b_ptr[i] = prng_c.gaussianRandomNumber();
  }

  // Check that the GEMM matrix multiplication snapshot file exists
  const std::string matrices_snp = oe.getStormmSourcePath() + osSeparator() + "test" +
                                   osSeparator() + "Math" + osSeparator() + "matrices.m";
  snp_found = (getDrivePathType(matrices_snp) == DrivePathType::FILE) ? TestPriority::CRITICAL :
                                                                        TestPriority::ABORT;
  if (snp_found == TestPriority::ABORT && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("The snapshot file " + matrices_snp + "was not found.  Make sure that the $STORMM_SOURCE "
           "environment variable is set to the root soruce directory, where src/ and test/ "
           "subdirectories can be found.  A number of subsequent tests will be skipped.",
           "test_math");
  }

  // Try a direct matrix-matrix multiplication, producing a non-symmetric [7 by 7] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_b.data(), rows_b, cols_b,
                 tr_mat_c.data());
  snapshot(matrices_snp, polyNumericVector(tr_mat_c.readHost()), "tr_ab", 1.0e-5, "Matrix "
           "multiply result for [7 by 5] x [5 by 7] matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE, snp_found);

  // Overwrite the result with a positive definite [7 by 7] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_a.data(), rows_a, cols_a, tr_mat_c.data(),
                 1.0, 1.0, 0.0, TransposeState::AS_IS, TransposeState::TRANSPOSE);
  snapshot(matrices_snp, polyNumericVector(tr_mat_c.readHost()), "tr_aat", 1.0e-5, "Matrix "
           "multiply result for [7 by 5] x [7 by 5](T) matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Obtain a new, positive definite [5 by 5] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_a.data(), rows_a, cols_a, tr_mat_d.data(),
                 1.0, 1.0, 0.0, TransposeState::TRANSPOSE, TransposeState::AS_IS);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "tr_ata", 1.0e-5, "Matrix "
           "multiply result for [7 by 5](T) x [7 by 5] matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);
  
  // Try inverting a positive definite matrix
  Hybrid<double> tr_mat_e(cols_a * rows_b, "inverse matrix");
  invertSquareMatrix(tr_mat_d.data(), tr_mat_e.data(), cols_a);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "inv_ata", 1.0e-5, "Matrix "
           "inversion result for a [5 by 5] matrix is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Overwrite the positive definite matrix with a non-symmetric [5 by 5] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_b.data(), rows_b, cols_b, tr_mat_d.data(),
                 1.0, 1.0, 0.0, TransposeState::TRANSPOSE, TransposeState::TRANSPOSE);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "tr_atbt", 1.0e-5, "Matrix "
           "multiply result for [7 by 5](T) x [5 by 7](T) matrices is incorrect.",
           oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           snp_found);

  // Create a positive-definite matrix, then compute its eigenvalues and eigenvectors (this is
  // better accomplished by a routine like BLAS dsyevd, as it is more than just real and symmetric,
  // but the slower jacobi routine is all STORMM has got without real BLAS compiled).
  const size_t rank = 8;
  Hybrid<double> base_matrix(rank * rank);
  Hybrid<double> posdef_matrix(rank * rank);
  Hybrid<double> eigenvectors(rank * rank);
  Hybrid<double> eigenvalues(rank);
  double* dbase = base_matrix.data();
  double* dposdef = posdef_matrix.data();
  for (size_t i = 0; i < rank * rank; i++) {
    dbase[i] = prng_c.gaussianRandomNumber();
  }
  matrixMultiply(dbase, rank, rank, dbase, rank, rank, dposdef, 1.0, 1.0, 0.0,
                 TransposeState::TRANSPOSE, TransposeState::AS_IS);
  jacobiEigensolver(&posdef_matrix, &eigenvectors, &eigenvalues, rank);
  snapshot(matrices_snp, polyNumericVector(eigenvectors.readHost()), "eigvec", 1.0e-5,
           "Eigenvectors for a rank-8 positive definite matrix are incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);
  snapshot(matrices_snp, polyNumericVector(eigenvalues.readHost()), "eigval", 1.0e-5,
           "Eigenvalues for a rank-8 positive definite matrix are incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Try a much bigger eigenvalue problem and check its results
  const size_t big_rank = 95;
  HpcMatrix<double> mtrx_base(big_rank, big_rank, MatrixFillMode::RANDN, &prng_c);
  HpcMatrix<double> mtrx_posdef(big_rank, big_rank);
  HpcMatrix<double> mtrx_eigvec(big_rank, big_rank);
  eigenvalues.resize(big_rank);
  matrixMultiply(mtrx_base, mtrx_base, &mtrx_posdef, 1.0, 1.0, 0.0, TransposeState::TRANSPOSE,
                 TransposeState::AS_IS);
  dposdef = mtrx_posdef.memptr();
  double bad_value = mtrx_posdef(big_rank / 2, 0) * 1.098;
  std::swap(bad_value, dposdef[big_rank / 2]);
  CHECK_THROWS(jacobiEigensolver(&mtrx_posdef, &mtrx_eigvec, &eigenvalues, ExceptionResponse::DIE),
               "The jacobiEigensolver() routine attempted to work on a non-symmetric matrix.");
  std::swap(bad_value, dposdef[big_rank / 2]);
  HpcMatrix<double> copy_posdef = mtrx_posdef;
  jacobiEigensolver(&mtrx_posdef, &mtrx_eigvec, &eigenvalues, ExceptionResponse::DIE);
  std::vector<double> eigvsums(big_rank, 0.0);
  Hybrid<double> mtrx_eigtest(big_rank, "test_eig");
  for (size_t i = 0; i < big_rank; i++) {
    Hybrid<double> icol_ptr = mtrx_eigvec.colptr(i);
    matrixVectorMultiply(copy_posdef, icol_ptr, &mtrx_eigtest, 1.0, 1.0, 0.0,
                         TransposeState::AS_IS);
    eigvsums[i] = sum<double>(mtrx_eigtest) - sum<double>(icol_ptr) * eigenvalues.readHost(i);
  }
  check(eigvsums, RelationalOperator::EQUAL,
        Approx(std::vector<double>(big_rank, 0.0)).margin(1.0e-8), "Eigenvalues and eigenvectors "
        "produced by jacobiEigensolver() are incorrect.");

  // Try computing box transformation matrices
  const double lx = 64.1;
  const double ly = 60.3;
  const double lz = 48.9;
  const double alpha =  98.7 * pi / 180.0;
  const double beta  = 103.2 * pi / 180.0;
  const double gamma =  85.4 * pi / 180.0;
  HpcMatrix<double> umat(3, 3);
  HpcMatrix<double> invu(3, 3);
  HpcMatrix<double> xfrm_prod(3, 3);
  computeBoxTransform(lx, ly, lz, alpha, beta, gamma, umat.memptr(), invu.memptr());
  matrixMultiply(umat, invu, &xfrm_prod);
  HpcMatrix<double> ident(3, 3, MatrixFillMode::EYE);
  std::vector<double> linear_umat(9, 0.0);
  std::vector<double> linear_invu(9, 0.0);
  std::vector<double> linear_xfrm_prod(9, 0.0);
  std::vector<double> linear_ident(9, 0.0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      linear_umat[(j * 3) + i] = umat(i, j);
      linear_invu[(j * 3) + i] = invu(i, j);
      linear_xfrm_prod[(j * 3) + i] = xfrm_prod(i, j);
      linear_ident[(j * 3) + i] = ident(i, j);
    }
  }

  // Data entered in column-major format may look like a transpose at first.  Check
  // the box transformation matrices against results from an external implementation.
  const std::vector<double> umat_answer = {
    1.5600624025e-02,  0.0000000000e+00,  0.0000000000e+00,
   -1.2551964058e-03,  1.6637338818e-02,  0.0000000000e+00,
    3.5203270719e-03,  2.3009524689e-03,  2.1204798362e-02 };
  const std::vector<double> invu_answer = {
    6.4100000000e+01,  0.0000000000e+00,  0.0000000000e+00,
    4.8359951370e+00,  6.0105766371e+01,  0.0000000000e+00,
   -1.1166357548e+01, -6.5221328286e+00,  4.7159137423e+01 };
  check(linear_umat, RelationalOperator::EQUAL, Approx(umat_answer).margin(1.0e-8),
        "The box transformation matrix is incorrect.");
  check(linear_invu, RelationalOperator::EQUAL, Approx(invu_answer).margin(1.0e-8),
        "The inverse transformation matrix is incorrect.");
  check(linear_xfrm_prod, RelationalOperator::EQUAL, linear_ident, "The product of box "
        "transformation and inverse transformation matrices is not the identity matrix.");

  // Further checks on the mean, standard deviation, variance, correlation, dot product,
  // magnitude, and projection operations
  section(1);
  check(mean(invu_answer), RelationalOperator::EQUAL, Approx(17.6124898394).margin(tiny),
        "The mean value of a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::VARIANCE), RelationalOperator::EQUAL,
        Approx(816.0346461040).margin(tiny), "The variance of a double-precision vector is not "
        "correct.");
  check(variance(invu_answer, VarianceMethod::STANDARD_DEVIATION), RelationalOperator::EQUAL,
        Approx(30.2991580224).margin(tiny), "The standard deviation of a double-precision vector "
        "is not correct.");
  check(variance(invu_answer, VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION),
        RelationalOperator::EQUAL, Approx(28.5663201359).margin(tiny), "The root mean squared "
        "deviation of a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::COEFFICIENT_OF_VARIATION),
        RelationalOperator::EQUAL, Approx(1.7203222428).margin(tiny), "The coefficient of "
        "variation for a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::NORMALIZED_RMSD),
        RelationalOperator::EQUAL, Approx(1.6219353650).margin(tiny), "The normalized root mean "
        "squared deviation for a double-precision vector is not correct.");
  check(magnitude(umat_answer), RelationalOperator::EQUAL, Approx(3.1449747049e-02).margin(tiny),
        "The magnitude of a double-precision vector is not computed correctly.");
  check(dot(invu_answer, umat_answer), RelationalOperator::EQUAL,
        Approx(2.9396135279).margin(tiny), "The dot product of two double-precision vectors is "
        "not computed correctly.");
  std::vector<double> invu_perturb(invu_answer);
  std::vector<double> invu_partial(invu_answer);
  for (size_t i = 6; i < invu_partial.size(); i++) {
    invu_partial[i] = 0.0;
  }
  for (size_t i = 0; i < invu_perturb.size(); i++) {
    invu_perturb[i] += xrs256pp.uniformRandomNumber() - 0.5;
    invu_partial[i] += xrs256pp.uniformRandomNumber() - 0.5;
  }
  std::vector<double> invu_project(invu_answer.size());
  project(invu_answer, invu_partial, &invu_project);
  const std::vector<double> invu_project_answer = { 63.5808430141, -0.2998731521,  0.4369816989,
                                                     5.2640142928, 60.3620556749,  0.0677375859,
                                                    -0.2593021730, -0.1475558457, -0.3946802908 };
  check(invu_project, RelationalOperator::EQUAL, Approx(invu_project_answer).margin(tiny),
        "The projection of one vector onto another does not meet expectations.");
  const std::vector<double> corr_x1 = {  0.5, -1.0, -5.4,  6.7,  7.2,  3.8, -4.1, 9.3 };
  const std::vector<double> corr_x2 = {  0.7, -0.8, -4.9,  9.2,  5.3,  4.0, -5.9, 7.9 };
  check(pearson(corr_x1, corr_x2), RelationalOperator::EQUAL, Approx(0.9651506029).margin(tiny),
        "The correlation coefficient computed for two vectors is incorrect.");

  // Check single-precision random number generation
  section(2);
  Xoroshiro128pGenerator test_128;
  Xoroshiro128pGenerator test_128f;
  Xoshiro256ppGenerator test_256;
  Xoshiro256ppGenerator test_256f;
  const size_t npts = 64;
  std::vector<double> t128_uni(npts), t128f_uni(npts), t128_gss(npts), t128f_gss(npts);
  std::vector<double> t256_uni(npts), t256f_uni(npts), t256_gss(npts), t256f_gss(npts);
  for (size_t i = 0; i < npts; i++) {
    t128_uni[i]  = test_128.uniformRandomNumber();
    t128f_uni[i] = test_128f.spUniformRandomNumber();
    t128_gss[i]  = test_128.gaussianRandomNumber();
    t128f_gss[i] = test_128f.spGaussianRandomNumber();
    t256_uni[i]  = test_256.uniformRandomNumber();
    t256f_uni[i] = test_256f.spUniformRandomNumber();
    t256_gss[i]  = test_256.gaussianRandomNumber();
    t256f_gss[i] = test_256f.spGaussianRandomNumber();
  }
  check(t128f_uni, RelationalOperator::EQUAL, Approx(t128_uni).margin(3.0e-7), "Random numbers "
        "pulled from a uniform distribution with the Xoroshiro128+ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t128f_gss, RelationalOperator::EQUAL, Approx(t128_gss).margin(6.0e-6), "Random numbers "
        "pulled from a normal distribution with the Xoroshiro128+ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t256f_uni, RelationalOperator::EQUAL, Approx(t256_uni).margin(3.0e-7), "Random numbers "
        "pulled from a uniform distribution with the Xoshiro256++ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t256f_gss, RelationalOperator::EQUAL, Approx(t256_gss).margin(6.0e-6), "Random numbers "
        "pulled from a normal distribution with the Xoshiro256++ generator do not have the same "
        "values when produced in single- versus double-precision.");

  // Check the behavior of the Xoroshiro128+ series object
  Xoroshiro128pGenerator leader_xrs128p(9183084, 75);
  const ullint2 init_state = leader_xrs128p.revealState();
  Xoroshiro128pGenerator test_xrs128p(105892, 95);
  test_xrs128p.setState(init_state);
  std::vector<double> orig_uni_output(npts), rstr_uni_output(npts);
  for (int i = 0; i < npts; i++) {
    orig_uni_output[i] = leader_xrs128p.uniformRandomNumber();
    rstr_uni_output[i] = test_xrs128p.uniformRandomNumber();
  }
  check(rstr_uni_output, RelationalOperator::EQUAL, orig_uni_output, "Resetting the state of a "
        "Xoroshiro128+ generator failed to restart the sequence as expected.");
  const int ngen = 1024;
  const int bank_depth = 8;
  RandomNumberMill<double> my_128p_series(init_state, ngen, bank_depth);
  CHECK_THROWS(my_128p_series.getBankValue(ngen + 3, 7), "Invalid random number bank access was "
               "permitted in a Xoroshiro128+ generator series object.");
  int mseries_state = 0;
  const int stride_length = my_128p_series.getRefreshStride();
  check(stride_length, RelationalOperator::EQUAL, ngen / bank_depth, "The length of a refresh "
        "stride in the Xoroshiro128+ generator series should be such that running through a "
        "number of refresh cycles equal to the depth of the object's banks should cover all "
        "generators.");
  const int sample_a = 0;
  const int sample_b = ngen - 53;
  const int sample_c = (ngen / 2) - 37;
  const int sample_d = ngen / 4;
  const int maxcyc = 4;
  std::vector<double> a_path(bank_depth * maxcyc), b_path(bank_depth * maxcyc);
  std::vector<double> c_path(bank_depth * maxcyc), d_path(bank_depth * maxcyc);
  for (int cyc = 0; cyc < maxcyc; cyc++) {
    for (int i = 0; i < bank_depth; i++) {
      a_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_a, i);
      b_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_b, i);
      c_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_c, i);
      d_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_d, i);
      my_128p_series.gaussianRandomNumbers(i * stride_length, (i + 1) * stride_length);
    }
  }
  const std::vector<double> a_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_a,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> b_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_b,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> c_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_c,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> d_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_d,
                                                              ngen, bank_depth, maxcyc);
  check(a_path, RelationalOperator::EQUAL, a_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_a) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(b_path, RelationalOperator::EQUAL, b_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_b) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(c_path, RelationalOperator::EQUAL, c_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_c) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(d_path, RelationalOperator::EQUAL, d_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_d) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
