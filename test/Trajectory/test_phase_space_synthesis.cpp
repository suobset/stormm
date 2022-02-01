#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::random::Ran2Generator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using namespace omni::cuda;
using namespace omni::trajectory;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test PhaseSpaceSynthesis layout");

  // Create some vectors of random numbers, then upload them and test what happens when perturbing
  // atomic coordinates by these numbers.
  Ran2Generator my_prng(oe.getRandomSeed());
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_crd_name = base_crd_name + osc + "tip3p.inpcrd";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_crd_name = base_crd_name + osc + "tip4p.inpcrd";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const bool files_exist = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip3p_ps.buildFromFile(tip3p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    tip4p_ps.buildFromFile(tip4p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    trpcage_ps.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtWarn("The topology and coordinate files for the TIP3P and TIP4P water boxes as well as the "
           "Trp-cage miniprotein in water must be available in ${OMNI_SOURCE}/test/ "
           "subdirectories Topology and Trajectory, respectively.  Check the $OMNI_SOURCE "
           "environment variable to make sure that it is set properly.  A number of tests will "
           "be skipped.", "test_precision_model");
  }

  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
