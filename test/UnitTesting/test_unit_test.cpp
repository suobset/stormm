#include "copyright.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/Random/random.h"

using stormm::constants::warp_size_int;
#ifndef STORMM_USE_HPC
using stormm::data_types::int4;
#endif
using stormm::errors::rtWarn;
using stormm::parse::polyNumericVector;
using stormm::random::Ran2Generator;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using namespace stormm::diskutil;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Structs and methods");
  
  // Section 2
  section("Snapshotting");

  // Section 3
  section("Test system manager");
  
  // Section 4
  section("Timings");

  // Perform basic checks of the Approx object
  StopWatch section_timer("Trial Timer");
  section_timer.addCategory(unitTestSectionName(1));
  section(1);
  Approx gray_number(9.6, ComparisonType::ABSOLUTE, 0.11);
  check(gray_number.test(9.7), "Approx object fails to perform a real-to-real scalar comparison.");
  Approx gray_vector(std::vector<double>{9.8, 8.1, 8.0, 4.5, 7.3}, ComparisonType::ABSOLUTE,
                     0.011);
  check(gray_vector.test(9.7) == false, "A scalar-to-vector approximate comparison was judged "
        "successful.");
  check(gray_vector.test(std::vector<double>{9.79, 8.11, 7.99, 4.49, 7.31}), "Approx object fails "
        "to pass a correct vector-to-vector comparison.");
  check(gray_vector.test(std::vector<double>{9.79, 8.08, 7.99, 4.49, 7.31}) == false,
        "Approx object fails to reject an unacceptable real vector-to-vector comparison.");
  check(gray_vector.test(std::vector<double>{9.8, 8.1, 8.0}) == false, "Approx object performs a "
        "comparison on vectors of different lengths and returns success.");
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
  check((short_vector > gray_vector) == false, "Approximate comparison processes a greater-than "
        "comparison of two vectors with different lengths and returns success.");
  check((short_vector < gray_vector) == false, "Approximate comparison processes a less-than "
        "comparison of two vectors with two different lengths and returns success.");
  check((short_vector >= gray_vector) == false, "Approximate comparison processes a "
        "greater-than-or-equal comparison of two vectors with different lengths, returning "
        "success.");
  check((short_vector <= gray_vector) == false, "Approximate comparison processes a "
        "less-than-or-equal comparison of two vectors, returning success despite them being of "
        "different lengths.");
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
           "will therefore not support some subsequent tests.  Make sure that the $STORMM_TMPDIR "
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

  // Test the test system manager
  section(3);
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string top_name_ext = "top";
  const std::string crd_name_ext = "inpcrd";
  const std::vector<std::string> sysnames = { "symmetry_L1", "stereo_L1", "med_1", "med_4" };
  TestSystemManager tsm(base_top_name, top_name_ext, sysnames, base_crd_name, crd_name_ext,
                        sysnames);
  std::vector<std::string> all_top_names, all_crd_names;
  all_top_names.reserve(sysnames.size());
  all_crd_names.reserve(sysnames.size());
  for (size_t i = 0; i < sysnames.size(); i++) {
    all_top_names.push_back(base_top_name + osc + sysnames[i] + "." + top_name_ext);
    all_crd_names.push_back(base_crd_name + osc + sysnames[i] + "." + crd_name_ext);
  }
  bool systems_exist = true;
  for (size_t i = 0; i < sysnames.size(); i++) {
    systems_exist = (systems_exist &&
                     getDrivePathType(all_top_names[i]) == DrivePathType::FILE &&
                     getDrivePathType(all_crd_names[i]) == DrivePathType::FILE);
  }
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  check(tsm.getTestingStatus() == do_tests, "A hand-coded search of the files composing the "
        "TestSystemManager object does not produce the same testing status.");
  check(tsm.getTopologyFile(1), RelationalOperator::EQUAL, all_top_names[1],
        "The TestSystemManager does not report the expected name for one of its topology files.",
        tsm.getTestingStatus());
  check(tsm.getCoordinateFile(1), RelationalOperator::EQUAL, all_crd_names[1],
        "The TestSystemManager does not report the expected name for one of its topology files.",
        tsm.getTestingStatus());
  check(tsm.getSystemCount(), RelationalOperator::EQUAL, 4, "The number of systems in the "
        "TestSystemManager is incorrect.", tsm.getTestingStatus());
  const CoordinateFrame cf_one = tsm.exportCoordinateFrame(1);
  check(cf_one.getAtomCount(), RelationalOperator::EQUAL, 78, "The TestSystemManager does not "
        "export CoordinateFrame objects correctly.", tsm.getTestingStatus());
  const PhaseSpace ps_two = tsm.exportPhaseSpace(2);
  check(ps_two.getAtomCount(), RelationalOperator::EQUAL, 44, "The TestSystemManager does not "
        "export PhaseSpace objects correctly.", tsm.getTestingStatus());
  const PhaseSpaceSynthesis all_coords = tsm.exportPhaseSpaceSynthesis({ 0, 1, 2, 3 }, 0.0,
                                                                       1, 24);
  const CoordinateFrame cf_two = all_coords.exportCoordinates(2);
  const std::vector<double> cf_two_xyz = cf_two.getInterlacedCoordinates();
  const std::vector<double> ps_two_xyz = ps_two.getInterlacedCoordinates();
  const double roundoff_err = meanUnsignedError(cf_two_xyz, ps_two_xyz);
  check(roundoff_err, RelationalOperator::LESS_THAN, 5.0e-8, "The TestSystemManager's synthesis "
        "products contain unexpected roundoff error, or perhaps the wrong order of systems.",
        tsm.getTestingStatus());
  const double perturbation_sigma = 0.1;
  const int igseed = 915087;
  const CoordinateSeries<float> cs_three = tsm.exportCoordinateSeries<float>(3, 4,
                                                                             perturbation_sigma,
                                                                             igseed);
  Xoshiro256ppGenerator main_xrs(igseed);
  CoordinateSeries<double> cs_three_main(tsm.exportCoordinateFrame(3), 4);
  CoordinateSeriesWriter<double> csr = cs_three_main.data();
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < csr.natom; j++) {
      const size_t atom_ij = (i * roundUp(csr.natom, warp_size_int)) + j;
      csr.xcrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
      csr.ycrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
      csr.zcrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
    }
  }
  const CoordinateFrame cs_outcome      = cs_three.exportFrame(1);
  const CoordinateFrame cs_outcome_main = cs_three_main.exportFrame(1);
  const std::vector<double> frm_one_xyz      = cs_outcome.getInterlacedCoordinates();
  const std::vector<double> frm_one_xyz_main = cs_outcome_main.getInterlacedCoordinates();
  const double noise_roundoff = meanUnsignedError(frm_one_xyz, frm_one_xyz_main);
  check(noise_roundoff, RelationalOperator::LESS_THAN, 1.0e-7, "The TestSystemManager does not "
        "introduce the expected noise into a CoordinateSeries when requested.",
        tsm.getTestingStatus());  
  
  // Test the benchmarking 
  section_timer.addCategory(unitTestSectionName(3) + ", Part A");
  section(4);
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

  return countGlobalTestFailures();
}
