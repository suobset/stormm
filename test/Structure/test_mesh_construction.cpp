#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::review::stormmSplash;
using namespace stormm::structure;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Establish various mesh objects");

  // Create various meshes.  Test constructor overloads and traps.
  section(1);
  BackgroundMesh<double> bgm_a;
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
