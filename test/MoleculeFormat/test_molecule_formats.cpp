#include "copyright.h"
#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MoleculeFormat/mdl_mol_obj.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::structure;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  const char osc = osSeparator();

  // Section 1
  section("Test the MDL MOL entry");

  // Open an MDL MOL file with a single entry in V2000 format
  const std::string base_name = oe.getStormmSourcePath() + osc + "test" + osc + "MoleculeFormat";
  const std::string mdl_name      = base_name + osc + "sulfonamide.mol";
  const std::string sdf_name      = base_name + osc + "sulfonamide_rots.sdf";
  const std::string chemaxon_name = base_name + osc + "sdf_chemaxon.sdf";
  const bool files_ready = (getDrivePathType(mdl_name) == DrivePathType::FILE &&
                            getDrivePathType(sdf_name) == DrivePathType::FILE &&
                            getDrivePathType(chemaxon_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_ready) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (files_ready == false) {
    rtWarn("Files for MDL MOL format molecules were not found.  Check the STORMM source path, "
           "currently " + oe.getStormmSourcePath() + ", to ensure that it contains a valid "
           "subdirectory test/MoleculeFormat.  Subsequent tests will be skipped.");
  }
  TextFile mdl_a(mdl_name);
  MdlMolObj sulf_a(mdl_a);
  check(sulf_a.getAtomCount(), RelationalOperator::EQUAL, 16, "The atom count from an MDL MOL "
        "file is incorrect.", do_tests);
  check(sulf_a.getBondCount(), RelationalOperator::EQUAL, 15, "The bond count from an MDL MOL "
        "file is incorrect.", do_tests);
  const std::vector<double> sulf_a_ycrds = { -1.0845,  0.1506,  0.2311,  1.6123, -0.3787, -0.8405,
                                             -0.0173, -1.8467, -0.8392, -1.4937,  0.6911, -1.5531,
                                             -1.3843,  0.5307, -0.6815,  0.7010 };
  check(sulf_a.getCoordinates(CartesianDimension::Y), RelationalOperator::EQUAL,
        Approx(sulf_a_ycrds).margin(verytiny), "Coordinates for an MDL MOL entry were not taken "
        "up as expected.  The Y coordinate values do not match.", do_tests);
  const std::vector<int> sulf_a_znumbers = { 6, 7, 16, 8, 8, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  check(sulf_a.getAtomicNumbers(), RelationalOperator::EQUAL, sulf_a_znumbers, "Atomic numbers "
        "inferred for an MDL MOL entry do notmeet expectations.", do_tests);
  std::vector<MdlMolObj> sulf_rtm = readStructureDataFile(sdf_name);
  const std::vector<double> sulf_rtm_c_xcrds = {  2.2071,  1.4262, -0.0626, -0.5187,  0.1066,
                                                 -1.1940, -1.9484,  1.8501,  3.2455,  2.1797,
                                                  1.4358, -1.8772, -0.5739, -2.5472, -2.6168,
                                                 -1.2508 };
  check(sulf_rtm[2].getCoordinates(CartesianDimension::X), RelationalOperator::EQUAL,
        Approx(sulf_rtm_c_xcrds).margin(verytiny), "Coordinates for an MDL MOL entry within an "
        "SD file were not taken up as expected.  The X coordinate values do not match.", do_tests);

  // Read an SD file with V2000 MDL MOL entries containing some properties
  std::vector<MdlMolObj> chemaxon_mols = readStructureDataFile(chemaxon_name);
  check(chemaxon_mols[1].getPropertiesCount(), RelationalOperator::EQUAL, 1, "An incorrect number "
        "of properties were found in the second V2000 MDL MOL entry of " + chemaxon_name + ".",
        do_tests);
  const std::vector<int> fc_result = { chemaxon_mols[1].getFormalCharge(22),
                                       chemaxon_mols[1].getFormalCharge(22) };
  check(fc_result, RelationalOperator::EQUAL, std::vector<int>(2, 1), "Formal charges were not "
        "properly interpreted from an MDL MOL V2000 format property.", do_tests);

  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
