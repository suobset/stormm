#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::data_types::ullint;
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
  BackgroundMesh<ullint> bgm_a(GridDetail::OCCLUSION, NonbondedPotential::CLASH);
  bgm_a.setTopologyPointer(tsm.getTopologyPointer(0));
  bgm_a.setCoordinatePointer(&trpi_cf);
  bgm_a.setMeshParameters(10.0, 1.0);
  bgm_a.setProbeRadius(1.4);
  bgm_a.computeField();
  BackgroundMeshWriter<double, ullint> bgmw_a = bgm_a.dpData();
  const std::vector<double> trpi_mesh_dims = { static_cast<double>(bgmw_a.dims.na),
                                               static_cast<double>(bgmw_a.dims.nb),
                                               static_cast<double>(bgmw_a.dims.nc),
                                               bgmw_a.dims.invu[0], bgmw_a.dims.invu[1],
                                               bgmw_a.dims.invu[2], bgmw_a.dims.invu[3],
                                               bgmw_a.dims.invu[4], bgmw_a.dims.invu[5],
                                               bgmw_a.dims.invu[6], bgmw_a.dims.invu[7],
                                               bgmw_a.dims.invu[8] };
  const std::vector<double> trpi_mesh_dims_ans = { 45.0, 43.0, 37.0,  1.0,  0.0,  0.0,
                                                    0.0,  1.0,  0.0,  0.0,  0.0,  1.0 };
  check(trpi_mesh_dims, RelationalOperator::EQUAL, trpi_mesh_dims_ans, "The mesh dimensions "
        "computed for Trp-cage do not meet expectations.", tsm.getTestingStatus());
  BackgroundMesh<double> bgm_b(GridDetail::NONBONDED_FIELD, NonbondedPotential::ELECTROSTATIC,
                               tsm.getTopologyPointer(0), &trpi_cf, 10.0, 1.0);
  BackgroundMeshWriter<double, double> bgmw_b = bgm_b.dpData();
  
  // CHECK
  for (int i = 0; i < 45; i++) {
    printf("%12.4lf\n", bgmw_b.coeffs[64 * ((((23 * 43) + 22) * 45) + i)]);
  }
  // END CHECK

  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
