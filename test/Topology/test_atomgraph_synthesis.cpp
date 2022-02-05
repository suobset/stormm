#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_synthesis.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::random::Ran2Generator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using namespace omni::topology;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test AtomGraphSynthesis layout");

  // Create some vectors of random numbers, then upload them and test what happens when perturbing
  // atomic coordinates by these numbers.
  const char osc = osSeparator();
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_nbfix_top_name = base_top_name + osc + "trpcage_in_water_nbfix.top";
  const std::string ubiquitin_top_name = base_top_name + osc + "ubiquitin.top";
  const bool files_exist = (getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_nbfix_top_name) == DrivePathType::FILE &&
                            getDrivePathType(ubiquitin_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag, nbfix_ag, ubiquitin_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    nbfix_ag.buildFromPrmtop(trpcage_nbfix_top_name);
    ubiquitin_ag.buildFromPrmtop(ubiquitin_top_name);
  }
  else {
    rtWarn("The topology files for the TIP3P and TIP4P water boxes as well as two versions of the "
           "solvated Trp-cage miniprotein and ubiquitin must be available in ${OMNI_SOURCE}/test/ "
           "subdirectories Topology and Trajectory, respectively.  Check the $OMNI_SOURCE "
           "environment variable to make sure that it is set properly.  A number of tests will "
           "be skipped.", "test_atomgraph_synthesis");
  }

  // Create the synthesis
  const std::vector<AtomGraph*> all_tops = { &tip3p_ag, &tip4p_ag, &trpcage_ag, &nbfix_ag,
                                             &ubiquitin_ag };
  const std::vector<int> system_ids = { 0, 1, 2, 3, 4, 3, 3, 4, 2, 1, 1, 3 };
  AtomGraphSynthesis synth(all_tops, system_ids);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
