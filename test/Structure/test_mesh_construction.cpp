#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::review::stormmSplash;
using namespace stormm::diskutil;
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
  char osc = osSeparator();
  const std::string crd_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string top_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  TestSystemManager tsm(top_base_name, "top", { "trpcage" },
                        crd_base_name, "inpcrd", { "trpcage" });
  const CoordinateFrame trpi_cf = tsm.exportCoordinateFrame(0);
  BackgroundMesh<double> bgm_a;
  bgm_a.setTopologyPointer(tsm.getTopologyPointer(0));
  bgm_a.setCoordinatePointer(&trpi_cf);
  bgm_a.setMeshParameters(10.0, 1.0);
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
