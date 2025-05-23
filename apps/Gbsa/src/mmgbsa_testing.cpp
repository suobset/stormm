#include "../../../src/copyright.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/FileManagement/file_util.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Parsing/parsing_enumerators.h"
#include "../../../src/Parsing/textfile.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/summary_file.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/unit_test_enumerators.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "mmgbsa_problem_set.h"
#include "mmgbsa_testing.h"

namespace mmgbsa {

using stormm::testing::TestEnvironment;
using stormm::testing::countGlobalTestFailures;
using stormm::diskutil::getBaseName;
using stormm::diskutil::osSeparator;
using stormm::errors::rtErr;
using stormm::namelist::FilesControls;
using stormm::parse::TextFile;
using stormm::parse::TextOrigin;
using stormm::review::stormmWatermark;
using stormm::synthesis::SystemCache;
using stormm::testing::check;
using stormm::testing::RelationalOperator;
using stormm::testing::TestVerbosity;

//-------------------------------------------------------------------------------------------------
int runUnitTests() {
  const char* mock_argv[] = { "nothing" };
  TestEnvironment oe(1, mock_argv);

  // Checks on the class for staging the calculations, based on mock input
  const char osc = osSeparator();
  const std::string input_str("&files\n    -sys { -p " + oe.getStormmSourcePath() + osc + "test" +
                              osc + "Topology" + osc + "symmetry_C1.top\n"
                              "        -c " + oe.getStormmSourcePath() + osc + "test" + osc + 
                              "Trajectory" + osc + "symmetry_C1.inpcrd -label ligand }\n"
                              "    -sys { -p " + oe.getStormmSourcePath() + osc + "test" + osc +
                              "Topology" + osc + "symmetry_C4.top\n"
                              "        -c " + oe.getStormmSourcePath() + osc + "test" + osc +
                              "Trajectory" + osc + "symmetry_C4.inpcrd -label ligand }\n"
                              "    -sys { -p " + oe.getStormmSourcePath() + osc + "test" + osc +
                              "Topology" + osc + "symmetry_C3.top\n"
                              "        -c " + oe.getStormmSourcePath() + osc + "test" + osc +
                              "Trajectory" + osc + "symmetry_C3.inpcrd -label receptor }\n"
                              "&end\n\n"
                              "&solvent\n  igb = 7,\n&end\n");
  TextFile tf(input_str, TextOrigin::RAM);
  int start_line = 0;
  FilesControls ficon(tf, &start_line);
  SystemCache sysche(ficon);
  try {
    MMGBSAProblemSet mockset(sysche.getSelfPointer());
    check(mockset.getProblemCount(), RelationalOperator::EQUAL, 5, "The number of problems "
          "reported by a MMGBSAProblemSet does not meet expectations.");
    const int rec_cache_idx = mockset.getReceptorTopologyCacheIndex();
    check(sysche.getSystemTopology(rec_cache_idx).getAtomCount(), RelationalOperator::EQUAL, 25,
          "The number of atoms found in the receptor of a mock system does not meet "
          "expectations.  This may indicate a confusion as to which system is the receptor and "
          "which is the ligand.");
    check(mockset.getProblemTopologyPointer(4)->getAtomCount(), RelationalOperator::EQUAL,
          mockset.getProblemTopologyPointer(0)->getAtomCount() +
          mockset.getProblemTopologyPointer(2)->getAtomCount(), "The number of atoms in the "
          "complex of a mock receptor and ligand system, consisting of test systems " +
          getBaseName(mockset.getProblemTopologyPointer(0)->getFileName()) + "(R) and " +
          getBaseName(mockset.getProblemTopologyPointer(0)->getFileName()) + "(L) does not meet "
          "expectations.");
  }
  catch (std::runtime_error) {
    rtErr(std::string("Errors occurred when trying to assemble the files for the unit testing.  "
                      "Check the STORMM source path environment variable ${STORMM_SOURCE} to "
                      "ensure that the test directories are in the expected places, e.g. "
                      "${STORMM_SOURCE}") +
          osc + "test" + osc + " contains Trajectory and Topology subdirectories, with "
          "symmetry_C*.inpcrd and symmetry_C*.top files in them.", "mmgbsa", "runUnitTests");
  }

  // Report the testing results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}

//-------------------------------------------------------------------------------------------------
int runRegressionTests() {
  return countGlobalTestFailures();
}

} // namespace mmgbsa
