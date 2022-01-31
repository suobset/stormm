#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/benchmark.h"
#include "../../src/Random/random.h"

using omni::data_types::int4;
using omni::diskutil::osSeparator;
using omni::diskutil::PrintSituation;
using omni::errors::rtWarn;
using omni::parse::polyNumericVector;
using omni::random::Ran2Generator;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);

  // Section 1
  section("Structs and methods");
  
  // Section 2
  section("Snapshotting");

  // Section 3
  section("Timings");

  // Perform basic checks of the Approx object
  StopWatch section_timer("Trial Timer");
  section_timer.addCategory(unitTestSectionName(1));
  section(1);
  Approx gray_number(9.6, ComparisonType::ABSOLUTE, 0.11);
  check(gray_number.test(9.7), "Approx object fails to perform a real-to-real scalar comparison.");
  Approx gray_vector(std::vector<double>{9.8, 8.1, 8.0, 4.5, 7.3}, ComparisonType::ABSOLUTE,
                     0.011);
  CHECK_THROWS(gray_vector.test(9.7), "Approx object performs a scalar-to-vector "
               "comparison without throwing an exception.");
  check(gray_vector.test(std::vector<double>{9.79, 8.11, 7.99, 4.49, 7.31}), "Approx object fails "
        "to pass a correct vector-to-vector comparison.");
  check(gray_vector.test(std::vector<double>{9.79, 8.08, 7.99, 4.49, 7.31}) == false,
        "Approx object fails to reject an unacceptable real vector-to-vector comparison.");
  CHECK_THROWS(gray_vector.test(std::vector<double>{9.8, 8.1, 8.0}), "Approx object performs a "
               "comparison on vectors of different lengths without throwing an exception.");
  check(gray_number, RelationalOperator::GREATER_THAN, 9.4, "Approximate comparison fails to "
        "process a scalar greater-than inequality.");
  check((gray_number > 9.8) == false, "Approximate comparison fails to properly declare a "
        "scalar greater-than inequality to be false.");
  check(4.5, RelationalOperator::LESS_THAN, gray_number, "Approximate comparison fails to "
        "process a scalar less-than inequality.");
  check((9.8 < gray_number) == false, "Approximate comparison fails to properly declare a "
        "scalar less-than inequality to be false.");
  check(4, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(6, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(4, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(6, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(8), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(std::vector<double>{9.79, 8.08, 7.99, 4.49, 7.31}, RelationalOperator::LE,
        gray_vector, "Approximate less-than-or-equal comparisons should give the benefit of the "
        "doubt, and let tolerances work in favor of any inequality comparisons.  This was not the "
        "case for a comparison of two vectors.");
  check((std::vector<double>{11.79, 8.08, 7.99, 4.49, 7.31} < gray_vector) == false,
        "Approximate comparison fails to properly identify a less-than comparison of two vectors "
        "as false.  Any element violating the inequality should render the entire comparison "
        "false.");
  check(std::vector<double>{11.79, 80.08, 71.99, 40.49, 17.31}, RelationalOperator::GREATER_THAN,
        gray_vector, "Approximate comparison fails to process a greater-than comparison of two "
        "vectors.");
  check(6, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(4, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(6, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(8), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(908153, RelationalOperator::GE, 908153, "Inferred approximate comparison is too harsh in "
        "comparing two medium-sized integer values by greater-than-or-equal inequality.");
  check(908153, RelationalOperator::LE, 908153, "Inferred approximate comparison is too harsh in "
        "comparing two medium-sized integer values by less-than-or-equal inequality.");
  check(-90153, RelationalOperator::LT, 90153, "Inferred approximate comparison misses an obvious "
        "comparison of a significant number with its own negative value.");
  const std::vector<double> short_vector = {11.79, 80.08, 71.99, 40.49};
  CHECK_THROWS(bool test_a = (short_vector > gray_vector), "Approximate comparison processes a "
               "greater-than comparison of two vectors despite them being of different lengths.");
  CHECK_THROWS(bool test_b = (short_vector < gray_vector), "Approximate comparison processes a "
               "less-than comparison of two vectors despite them being of different lengths.");
  CHECK_THROWS(bool test_c = (short_vector >= gray_vector), "Approximate comparison processes a "
               "greater-than-or-equal comparison of two vectors despite them being of different "
               "lengths.");
  CHECK_THROWS(bool test_d = (short_vector <= gray_vector), "Approximate comparison processes a "
               "less-than-or-equal comparison of two vectors despite them being of different "
               "lengths.");
  section_timer.assignTime(1);

  // Examine file snapshotting by first recording data, then re-reading it.
  section_timer.addCategory(unitTestSectionName(2));
  section(2);
  const int n_pts = 100;
  Ran2Generator prng(71277);
  std::vector<double> x(n_pts, 0.0);
  std::vector<double> y(n_pts, 0.0);
  std::vector<double> z(n_pts, 0.0);
  for (int i = 0; i < n_pts; i++) {
    x[i] = prng.gaussianRandomNumber();
    y[i] = prng.gaussianRandomNumber();
    z[i] = prng.gaussianRandomNumber();
  }
  TestPriority snp_tests = (oe.getTemporaryDirectoryAccess()) ? TestPriority::CRITICAL :
                                                                TestPriority::ABORT;
  if (snp_tests == TestPriority::ABORT) {
    rtWarn("The temporary directory " + oe.getTemporaryDirectoryPath() + " is not writeable and "
           "will therefore not support some subsequent tests.  Make sure that the $OMNI_TMPDIR "
           "environment variable is set to a directory where you have write permissions.",
           "test_unit_test");
  }
  const std::string snp_file = oe.getTemporaryDirectoryPath() + osSeparator() + "xyz_randoms.m";
  snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::OPEN_NEW, snp_tests);
  snapshot(snp_file, polyNumericVector(y), "y_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(z), "z_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4, "Failed to read x_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(y), "y_randoms", 1.0e-4, "Failed to read y_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(z), "z_randoms", 1.0e-4, "Failed to read z_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  CHECK_THROWS_SOFT(snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4,
                             "This should not work.", SnapshotOperation::SNAPSHOT, 1.0e-8,
                             NumberFormat::STANDARD_REAL, PrintSituation::OPEN_NEW),
                    "Snapshotting was able to open a new file " + snp_file +
                    " despite it already existing.", snp_tests);
  if (snp_tests == TestPriority::CRITICAL) {
    oe.logFileCreated(snp_file);
  }
  section_timer.assignTime(2);

  // Test the benchmarking 
  section_timer.addCategory(unitTestSectionName(3) + ", Part A");
  section(3);
  double t = 0.0;
  for (int i = 0; i < 7500 * n_pts; i++) {
    t += prng.gaussianRandomNumber();
  }
  section_timer.assignTime(3);
  section_timer.addCategory(unitTestSectionName(3) + ", Part B");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 7500 * n_pts; j++) {
      t += prng.gaussianRandomNumber();
    }
    section_timer.assignTime(4);
  }
  section_timer.addCategory(unitTestSectionName(3) + ", Part C");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 15000 * n_pts; j++) {
      t += prng.gaussianRandomNumber();
    }
    section_timer.assignTime(5);
  }
  check(section_timer.getCategorySamples(4), RelationalOperator::EQUAL, 4, "An incorrect number "
        "of samples is recorded for " + section_timer.getCategoryName(4) + ".");
  check(section_timer.getCategorySamples(3), RelationalOperator::EQUAL, 1, "An incorrect number "
        "of samples is recorded for " + section_timer.getCategoryName(3) + ".");
  check(section_timer.getCategoryMinimumTime(3), RelationalOperator::EQUAL,
        section_timer.getCategoryMaximumTime(3), "The minimum and maximum intervals recorded for "
        "a timing category with only one entry should be identical.");
  check(section_timer.getCategoryAverageInterval(5) /
        section_timer.getCategoryAverageInterval(4), RelationalOperator::EQUAL,
        Approx(2.0).margin(0.1), "Timings should roughly double for producing twice the quantity "
        "of random numbers.", TestPriority::NON_CRITICAL);

  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
