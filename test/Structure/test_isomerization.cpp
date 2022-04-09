#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Structure/isomerization.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/unit_test.h"

// CHECK
#include "../../src/FileManagement/file_util.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
using omni::trajectory::TrajectoryKind;
using omni::trajectory::CoordinateFileKind;
using omni::diskutil::PrintSituation;
// END CHECK

using omni::chemistry::ChemicalFeatures;
using omni::chemistry::MapRotatableGroups;
using omni::chemistry::RotatorGroup;
using omni::constants::tiny;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::math::computeBoxTransform;
using omni::energy::evaluateBondTerms;
using omni::energy::evaluateAngleTerms;
using omni::energy::ScoreCard;
using omni::topology::AtomGraph;
using omni::topology::ChemicalDetailsKit;
using omni::topology::UnitCellType;
using omni::topology::ValenceKit;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFrameReader;
using omni::trajectory::PhaseSpace;
using namespace omni::structure;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Check rotatable bonds throughout a structure.  If the topology supplied pertains to a system
// with many molecules, only the first will be analyzed.  Options are provided to take snapshots
// of detailed particle positions.
//-------------------------------------------------------------------------------------------------
void checkRotationalSampling(const AtomGraph &ag, const PhaseSpace &ps,
                             const ChemicalFeatures &chemfe, const TestEnvironment &oe,
                             const TestPriority do_tests,
                             const std::string &snp_var_name = std::string("")) {
  ScoreCard sc(1);
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const CoordinateFrameReader cfr(ps);
  const double orig_bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
  const double orig_angl_e = evaluateAngleTerms(vk, cfr, &sc, 0);
  section(1);
  
  // Rotate about various bonds.  This will generate all sorts of clashes, but bond and angle
  // energies should be unaffected.
  const std::vector<RotatorGroup> rt_grp = chemfe.getRotatableBondGroups();
  const int nrt = rt_grp.size();
  PhaseSpace ps_copy(ps);
  PhaseSpaceWriter psw = ps_copy.data();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double> rot_crd(3 * (cdk.mol_limits[1] - cdk.mol_limits[0]) * nrt);
  std::vector<double> repos_dev(nrt), ubond_dev(nrt), uangl_dev(nrt);
  int rcpos = 0;
  for (int i = 0; i < nrt; i++) {
    rotateAboutBond(&ps_copy, rt_grp[i].root_atom, rt_grp[i].pivot_atom,
                    rt_grp[i].rotatable_atoms, 2.0 * omni::symbols::pi / 3.0);

    // Record the positions of atoms in the first molecule
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      rot_crd[rcpos] = psw.xcrd[cdk.mol_contents[j]];
      rcpos++;
      rot_crd[rcpos] = psw.ycrd[cdk.mol_contents[j]];
      rcpos++;
      rot_crd[rcpos] = psw.zcrd[cdk.mol_contents[j]];
      rcpos++;
    }

    // Record the bond and angle energies
    const double new_bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
    const double new_angl_e = evaluateAngleTerms(vk, cfr, &sc, 0);
    ubond_dev[i] = fabs(new_bond_e - orig_bond_e);
    uangl_dev[i] = fabs(new_angl_e - orig_angl_e);

    // Reverse the rotation
    rotateAboutBond(&ps_copy, rt_grp[i].root_atom, rt_grp[i].pivot_atom,
                    rt_grp[i].rotatable_atoms, -2.0 * omni::symbols::pi / 3.0);

    // Check that the molecule was returned to its original state
    double rmsd = 0.0;
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      const size_t jatom = cdk.mol_contents[j];
      const double dx = psw.xcrd[jatom] - psr.xcrd[jatom];
      const double dy = psw.ycrd[jatom] - psr.ycrd[jatom];
      const double dz = psw.zcrd[jatom] - psr.zcrd[jatom];
      rmsd += (dx * dx) + (dy * dy) + (dz * dz);
    }
    rmsd = sqrt(rmsd / static_cast<double>(cdk.mol_limits[1] - cdk.mol_limits[0]));
    repos_dev[i] = rmsd;
  }
  const char osc = osSeparator();
  const std::string base_iso_path = oe.getOmniSourcePath() + osc + "test" + osc + "Structure";
  const std::string rcrd_snapshot = base_iso_path + osc + "rotated_coords.m";
  const bool snps_exist = (getDrivePathType(rcrd_snapshot) == DrivePathType::FILE);
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_var_name.size() > 0) {
    snapshot(rcrd_snapshot, polyNumericVector(rot_crd), snp_var_name, 1.0e-6, "Coordinate "
             "obtained from rotation of selected bonds do not meet expectations.",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE,
             do_snps);
  }
  const Approx target(std::vector<double>(nrt, 0.0), ComparisonType::ABSOLUTE, 1.0e-6); 
  check(repos_dev, RelationalOperator::EQUAL, target, "Reversing the rotational operations of "
        "selected bonds does not return the molecule to its original state.", do_tests);
  check(ubond_dev, RelationalOperator::EQUAL, target, "Bond energies are changed by rotation "
        "about a bond.", do_tests);
  check(uangl_dev, RelationalOperator::EQUAL, target, "Angle energies are changed by rotation "
        "about a bond.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Rotation about a bond");

  // Section 2
  section("Selected chiral inversions");
  
  // Get a handful of realistic systems
  const char osc = osSeparator();
  const std::string base_top_path = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string drug_top_path = base_top_path + osc + "drug_example.top";
  const std::string drug_crd_path = base_crd_path + osc + "drug_example.inpcrd";
  const std::string lig1_top_path = base_top_path + osc + "stereo_L1.top";
  const std::string lig1_crd_path = base_crd_path + osc + "stereo_L1.inpcrd";
  const std::string lig2_top_path = base_top_path + osc + "symmetry_L1.top";
  const std::string lig2_crd_path = base_crd_path + osc + "symmetry_L1.inpcrd";
  const std::string trpc_top_path = base_top_path + osc + "trpcage.top";
  const std::string trpc_crd_path = base_crd_path + osc + "trpcage.inpcrd";
  const bool files_exist = (getDrivePathType(drug_top_path) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(lig1_top_path) == DrivePathType::FILE &&
                            getDrivePathType(lig1_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(lig2_top_path) == DrivePathType::FILE &&
                            getDrivePathType(lig2_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_top_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_crd_path) == DrivePathType::FILE);
  if (files_exist == false) {
    rtWarn("Files for various drug molecules and a miniprotein were not found.  Check the "
           "$OMNI_SOURCE environment variable to ensure that " + drug_top_path + " and " +
           drug_crd_path + " become valid paths.  Some tests will be skipped",
           "test_local_arrangement");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph drug_ag = (files_exist) ? AtomGraph(drug_top_path) : AtomGraph();
  PhaseSpace drug_ps = (files_exist) ? PhaseSpace(drug_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures drug_feat = (files_exist) ? ChemicalFeatures(&drug_ag, drug_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph trpc_ag = (files_exist) ? AtomGraph(trpc_top_path) : AtomGraph();
  PhaseSpace trpc_ps = (files_exist) ? PhaseSpace(trpc_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures trpc_feat = (files_exist) ? ChemicalFeatures(&trpc_ag, trpc_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();

  // CHECK
  AtomGraph lig1_ag = (files_exist) ? AtomGraph(lig1_top_path) : AtomGraph();
  PhaseSpace lig1_ps = (files_exist) ? PhaseSpace(lig1_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures lig1_feat = (files_exist) ? ChemicalFeatures(&lig1_ag, lig1_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph lig2_ag = (files_exist) ? AtomGraph(lig2_top_path) : AtomGraph();
  PhaseSpace lig2_ps = (files_exist) ? PhaseSpace(lig2_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures lig2_feat = (files_exist) ? ChemicalFeatures(&lig2_ag, lig2_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  printf("Ligand 1 has %2d rotatable bonds and %d chiral centers.\n",
         lig1_feat.getRotatableBondCount(), lig1_feat.getChiralCenterCount());
  printf("Ligand 2 has %2d rotatable bonds and %d chiral centers.\n",
         lig2_feat.getRotatableBondCount(), lig2_feat.getChiralCenterCount());
  // END CHECK
  
  // Rotate bonds within each system
  checkRotationalSampling(drug_ag, drug_ps, drug_feat, oe, do_tests, "drug_rot_iso");
  checkRotationalSampling(trpc_ag, trpc_ps, trpc_feat, oe, do_tests);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
