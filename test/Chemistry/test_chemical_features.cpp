#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"

using omni::constants::ExceptionResponse;
using omni::data_types::int2;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::topology::AtomGraph;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::PhaseSpace;

using namespace omni::chemistry;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Lay out the testing environment
  TestEnvironment oe(argc, argv);

  // Test the existence of topology and coordinate files
  const char osc = osSeparator();
  const std::string base_chem_name = oe.getOmniSourcePath() + osc + "test" + osc + "Chemistry";
  const std::string base_nml_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Namelists";
  const std::string base_nmlp_name = base_nml_name + osc + "topol";
  const std::string base_nmlc_name = base_nml_name + osc + "coord";
  const std::string base_top_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";

  // Ligands with an eight-membered ring that will avoid anti-aromaticity by puckering
  const std::string mol1_top_name  = base_chem_name + osc + "lig1_c8h8.top";
  const std::string mol1_crd_name  = base_chem_name + osc + "lig1_c8h8.inpcrd";
  const std::string mol2_top_name  = base_chem_name + osc + "lig2_c8h8.top";
  const std::string mol2_crd_name  = base_chem_name + osc + "lig2_c8h8.inpcrd";

  // Additional ligands, with and without water
  const std::string mol3_top_name    = base_chem_name + osc + "lig3_1fsj.top";
  const std::string mol3_crd_name    = base_chem_name + osc + "lig3_1fsj.inpcrd";
  const std::string drug_vs_top_name = base_top_name + osc + "drug_example_vs.top";
  const std::string drug_vs_crd_name = base_crd_name + osc + "drug_example_vs.inpcrd";

  // Small peptides and proteins
  const std::string ala_top_name     = base_nmlp_name + osc + "ala.top";
  const std::string ala_crd_name     = base_nmlc_name + osc + "ala.inpcrd";
  const std::string gly_top_name     = base_nmlp_name + osc + "gly_gly.top";
  const std::string gly_crd_name     = base_nmlc_name + osc + "gly_gly.inpcrd";
  const std::string phe_top_name     = base_nmlp_name + osc + "phe.top";
  const std::string phe_crd_name     = base_nmlc_name + osc + "phe.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string ubiquit_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string ubiquit_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";

  // Water boxes
  const std::string tip5p_top_name = base_top_name + osc + "tip5p.top";
  const std::string tip5p_crd_name = base_crd_name + osc + "tip5p.rst";
  
  // Check the existence of all files
  const std::vector<std::string> top_files = { mol1_top_name, mol2_top_name, mol3_top_name,
                                               drug_vs_top_name, ala_top_name, gly_top_name,
                                               phe_top_name, trpcage_top_name, ubiquit_top_name };
  const std::vector<std::string> crd_files = { mol1_crd_name, mol2_crd_name, mol3_crd_name,
                                               drug_vs_crd_name, ala_crd_name, gly_crd_name,
                                               phe_crd_name, trpcage_crd_name, ubiquit_crd_name };
  bool files_exist = true;
  const size_t nsys = top_files.size();
  for (size_t i = 0; i < nsys; i++) {
    files_exist = (files_exist && getDrivePathType(top_files[i]) == DrivePathType::FILE);
    files_exist = (files_exist && getDrivePathType(crd_files[i]) == DrivePathType::FILE);
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (files_exist == false) {
    rtWarn("Some required files were not found.  Check the ${OMNI_SOURCE} environment variable, "
           "currently set to " + oe.getOmniSourcePath() + ", for validity.  Subsequent tests will "
           "be skipped.", "test_chemical_features");
  }
  std::vector<AtomGraph> sys_ag(nsys);
  std::vector<PhaseSpace> sys_ps(nsys);
  if (files_exist) {
    for (size_t i = 0; i < nsys; i++) {
      sys_ag[i].buildFromPrmtop(top_files[i], ExceptionResponse::SILENT);
      sys_ps[i].buildFromFile(crd_files[i], CoordinateFileKind::UNKNOWN);
    }
  }
  std::vector<int> first_mol_size(nsys);
  std::vector<ChemicalFeatures> sys_chem;
  sys_chem.reserve(nsys);
  const MapRotatableGroups mapg_yes = MapRotatableGroups::YES;
  const MapRotatableGroups mapg_no  = MapRotatableGroups::NO;
  if (files_exist) {
    for (size_t i = 0; i < nsys; i++) {
      const int2 first_mol_lims = sys_ag[i].getMoleculeLimits(0);
      first_mol_size[i] = first_mol_lims.y - first_mol_lims.x;
      sys_chem.emplace_back(&sys_ag[i], CoordinateFrameReader(sys_ps[i]),
                            (first_mol_size[i] < 120) ? mapg_yes : mapg_no);
    }
  }

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
