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
#include "../../src/Structure/global_manipulation.h"
#include "../../src/Structure/isomerization.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/rmsd.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::MapRotatableGroups;
using stormm::chemistry::RotatorGroup;
using stormm::constants::tiny;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::math::computeBoxTransform;
using stormm::energy::evaluateBondTerms;
using stormm::energy::evaluateAngleTerms;
using stormm::energy::ScoreCard;
using stormm::structure::rotateCoordinates;
using stormm::structure::rmsd;
using stormm::structure::RmsdMethod;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::UnitCellType;
using stormm::topology::ValenceKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::PhaseSpace;
using namespace stormm::structure;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Check rotatable bonds throughout a structure.  If the topology supplied pertains to a system
// with many molecules, only the first will be analyzed.  Options are provided to take snapshots
// of detailed particle positions.
//
// Arguments:
//   ag:            System topology
//   ps:            Coordinates of the system
//   chemfe:        Chemical features of the system
//   oe:            System environment variables (for snapshotting behavior)
//   snp_var_name:  The variable name under which to store or access snapshot data
//   expectation:   The manner in which to approach the snapshot file, if a variable name has been
//                  supplied to direct the flow of data into it (when updating snapshot files)
//-------------------------------------------------------------------------------------------------
void checkRotationalSampling(const AtomGraph &ag, const PhaseSpace &ps,
                             const ChemicalFeatures &chemfe, const TestEnvironment &oe,
                             const TestPriority do_tests,
                             const std::string &snp_var_name = std::string(""),
                             const PrintSituation expectation = PrintSituation::APPEND) {
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
  PhaseSpace rotation_copy(ps);
  PhaseSpaceWriter psw = rotation_copy.data();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double> rot_crd(3 * (cdk.mol_limits[1] - cdk.mol_limits[0]) * nrt);
  std::vector<double> repos_dev(nrt), ubond_dev(nrt), uangl_dev(nrt);
  int rcpos = 0;
  for (int i = 0; i < nrt; i++) {
    rotateAboutBond(&rotation_copy, rt_grp[i].root_atom, rt_grp[i].pivot_atom,
                    rt_grp[i].rotatable_atoms, 2.0 * stormm::symbols::pi / 3.0);

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
    const double new_bond_e = evaluateBondTerms(vk, CoordinateFrameReader(rotation_copy), &sc, 0);
    const double new_angl_e = evaluateAngleTerms(vk, CoordinateFrameReader(rotation_copy), &sc, 0);
    ubond_dev[i] = fabs(new_bond_e - orig_bond_e);
    uangl_dev[i] = fabs(new_angl_e - orig_angl_e);

    // Reverse the rotation
    rotateAboutBond(&rotation_copy, rt_grp[i].root_atom, rt_grp[i].pivot_atom,
                    rt_grp[i].rotatable_atoms, -2.0 * stormm::symbols::pi / 3.0);

    // Check that the molecule was returned to its original state
    repos_dev[i] = rmsd(ps, rotation_copy, ag, RmsdMethod::ALIGN_GEOM, cdk.mol_limits[0],
                        cdk.mol_limits[1]);
  }
  const char osc = osSeparator();
  const std::string base_iso_path = oe.getStormmSourcePath() + osc + "test" + osc + "Structure";
  const std::string rcrd_snapshot = base_iso_path + osc + "rotated_coords.m";
  const bool snps_exist = (getDrivePathType(rcrd_snapshot) == DrivePathType::FILE);
  if (snps_exist == false && oe.takeSnapshot() == SnapshotOperation::COMPARE) {
    rtWarn("The snapshot file " + rcrd_snapshot + " was not found.  Check the ${STORMM_SOURCE} "
           "environment variable, currently set to " + oe.getStormmSourcePath() + ", for validity.  "
           "Subsequent tests will be skipped.", "test_isomerization");
  }
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_var_name.size() > 0) {
    snapshot(rcrd_snapshot, polyNumericVector(rot_crd), snp_var_name, 1.0e-6, "Coordinates "
             "obtained from rotation of selected bonds do not meet expectations.",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, expectation, do_snps);
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
// Check chiral centers throughout a structure.  If the topology supplied pertains to a system
// with many molecules, only the first will be analyzed.  Options are provided to take snapshots
// of detailed particle positions resulting from inverting the various chiral centers.
//
// Arguments:
//   ag:            System topology
//   ps:            Coordinates of the system
//   chemfe:        Chemical features of the system
//   oe:            System environment variables (for snapshotting behavior)
//   snp_var_name:  The variable name under which to store or access snapshot data
//   expectation:   The manner in which to approach the snapshot file, if a variable name has been
//                  supplied to direct the flow of data into it (when updating snapshot files)
//-------------------------------------------------------------------------------------------------
void checkChiralSampling(const AtomGraph &ag, const PhaseSpace &ps,
                         const ChemicalFeatures &chemfe, const TestEnvironment &oe,
                         const TestPriority do_tests,
                         const std::string &snp_var_name = std::string(""),
                         const PrintSituation expectation = PrintSituation::APPEND) {
  ScoreCard sc(1);
  const std::vector<RotatorGroup> inv_grp = chemfe.getChiralInversionGroups();
  const std::vector<ChiralInversionProtocol> protocols = chemfe.getChiralInversionMethods();
  const std::vector<int> centers = chemfe.getChiralCenters();
  const int nchiral = protocols.size();
  if (nchiral == 0) {
    return;
  }
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const CoordinateFrameReader cfr(ps);
  const double orig_bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
  PhaseSpace inversion_copy(ps);
  PhaseSpaceWriter psw = inversion_copy.data();
  const PhaseSpaceReader psr = ps.data();
  int invcpos = 0;
  std::vector<double> inverted_crd(3 * nchiral * (cdk.mol_limits[1] - cdk.mol_limits[0]));
  std::vector<double> repos_dev(nchiral), ubond_dev(nchiral);
  for (int i = 0; i < nchiral; i++) {
    flipChiralCenter(&inversion_copy, i, centers, protocols, inv_grp);

    // Record the positions of atoms in the first molecule
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      inverted_crd[invcpos] = psw.xcrd[cdk.mol_contents[j]];
      invcpos++;
      inverted_crd[invcpos] = psw.ycrd[cdk.mol_contents[j]];
      invcpos++;
      inverted_crd[invcpos] = psw.zcrd[cdk.mol_contents[j]];
      invcpos++;
    }
    
    // Record the bond and angle energies
    const double new_bond_e = evaluateBondTerms(vk, CoordinateFrameReader(inversion_copy), &sc, 0);
    ubond_dev[i] = fabs(new_bond_e - orig_bond_e);
    
    // Reverse the inversion to check that the molecule recovers its initial state
    flipChiralCenter(&inversion_copy, i, centers, protocols, inv_grp);
    repos_dev[i] = rmsd(ps, inversion_copy, ag, RmsdMethod::ALIGN_GEOM, cdk.mol_limits[0],
                        cdk.mol_limits[1]);
  }
  const char osc = osSeparator();
  const std::string base_iso_path = oe.getStormmSourcePath() + osc + "test" + osc + "Structure";
  const std::string invcrd_snapshot = base_iso_path + osc + "inverted_coords.m";
  const bool snps_exist = (getDrivePathType(invcrd_snapshot) == DrivePathType::FILE);
  if (snps_exist == false && oe.takeSnapshot() == SnapshotOperation::COMPARE) {
    rtWarn("The snapshot file " + invcrd_snapshot + " was not found.  Check the ${STORMM_SOURCE} "
           "environment variable, currently set to " + oe.getStormmSourcePath() + ", for validity.  "
           "Subsequent tests will be skipped.", "test_isomerization");
  }
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_var_name.size() > 0) {
    snapshot(invcrd_snapshot, polyNumericVector(inverted_crd), snp_var_name, 1.0e-6, "Coordinates "
             "obtained from inversion of selected chiral centers do not meet expectations.",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, expectation, do_snps);
  }
  const Approx target(std::vector<double>(nchiral, 0.0), ComparisonType::ABSOLUTE, 1.0e-6);
  check(repos_dev, RelationalOperator::EQUAL, target, "Reversing the inversion operations of "
        "selected chiral centers does not return the molecule to its original state.", do_tests);
  check(ubond_dev, RelationalOperator::EQUAL, target, "Bond energies are changed by inversion "
        "of a chiral center.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Rotation about a bond");

  // Section 2
  section("Selected chiral inversions");

  // Section 3
  section("Test RMSD calculations");
  
  // Get a handful of realistic systems
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
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
           "$STORMM_SOURCE environment variable to ensure that " + drug_top_path + " and " +
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
  
  // Rotate bonds within each system
  checkRotationalSampling(drug_ag, drug_ps, drug_feat, oe, do_tests, "drug_rot_iso",
                          PrintSituation::OVERWRITE);
  checkRotationalSampling(trpc_ag, trpc_ps, trpc_feat, oe, do_tests);
  checkRotationalSampling(lig1_ag, lig1_ps, lig1_feat, oe, do_tests, "lig1_rot_iso");
  checkRotationalSampling(lig2_ag, lig2_ps, lig2_feat, oe, do_tests, "lig2_rot_iso");

  // Isomerize chiral centers
  checkChiralSampling(drug_ag, drug_ps, drug_feat, oe, do_tests, "drug_chir_iso",
                      PrintSituation::OVERWRITE);
  checkChiralSampling(trpc_ag, trpc_ps, trpc_feat, oe, do_tests);
  checkChiralSampling(lig1_ag, lig1_ps, lig1_feat, oe, do_tests, "lig1_chir_iso");
  checkChiralSampling(lig2_ag, lig2_ps, lig2_feat, oe, do_tests, "lig2_chir_iso");

  // Test RMSD computations on simple structures
  section(3);
  CoordinateFrame starfish_a(7);
  CoordinateFrameWriter strfa = starfish_a.data();
  strfa.xcrd[1] =  1.0;
  strfa.ycrd[2] =  1.0;
  strfa.xcrd[3] = -1.0;
  strfa.ycrd[4] = -1.0;
  strfa.zcrd[5] =  1.0;
  strfa.zcrd[6] = -1.0;
  CoordinateFrame starfish_b(starfish_a);
  CoordinateFrameWriter strfb = starfish_b.data();
  rotateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 0.0, 0.0, stormm::symbols::pi / 4.0, 0,
                    strfb.natom);
  double rms_no_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd,
                                             strfb.ycrd, strfb.zcrd, nullptr,
                                             RmsdMethod::NO_ALIGN_GEOM, 0, strfa.natom);
  double rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd,
                                          strfb.ycrd, strfb.zcrd, nullptr, RmsdMethod::ALIGN_GEOM,
                                          0, strfa.natom);
  check(rms_no_align, RelationalOperator::EQUAL, Approx(0.578562967).margin(1.0e-8),
        "Positional (non-aligned) RMSD computed for coordinates pre-shifted to their respective "
        "centers of mass does not produce the expected result.");
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates pre-shifted to their respective centers of mass "
        "does not produce the expected result.");
  translateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 5.0, 4.8, 9.7, 0, strfb.natom);
  rms_no_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                      strfb.zcrd, nullptr, RmsdMethod::NO_ALIGN_GEOM, 0,
                                      strfa.natom);
  rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                   strfb.zcrd, nullptr, RmsdMethod::ALIGN_GEOM, 0, strfa.natom);
  check(rms_no_align, RelationalOperator::EQUAL, Approx(11.935859211).margin(1.0e-8),
        "Positional (non-aligned) RMSD computed for coordinates differing in their respective "
        "centers of mass does not produce the expected result.");
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates differing in their respective centers of mass "
        "does not produce the expected result.");
  rotateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 0.1, -0.3, -0.25, 0, strfb.natom);
  rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                   strfb.zcrd, nullptr, RmsdMethod::ALIGN_GEOM, 0, strfa.natom);
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates differing in their respective centers of mass, "
        "one rotated a second time for more frustration, does not produce the expected result.");
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
