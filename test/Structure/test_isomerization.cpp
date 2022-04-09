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
  const std::string trpc_top_path = base_top_path + osc + "trpcage.top";
  const std::string trpc_crd_path = base_crd_path + osc + "trpcage.inpcrd";
  const std::string base_iso_path = oe.getOmniSourcePath() + osc + "test" + osc + "Structure";
  const bool files_exist = (getDrivePathType(drug_top_path) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_top_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_crd_path) == DrivePathType::FILE);
  if (files_exist == false) {
    rtWarn("Files for a drug molecule, in water and inside a periodic box, were not found.  Check "
           "the $OMNI_SOURCE environment variable to ensure that " + drug_top_path + " and " +
           drug_crd_path + " become valid paths.  Some tests will be skipped",
           "test_local_arrangement");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph drug_ag = (files_exist) ? AtomGraph(drug_top_path) : AtomGraph();
  PhaseSpace drug_ps = (files_exist) ? PhaseSpace(drug_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  CoordinateFrameReader drug_cfr = (files_exist) ? CoordinateFrameReader(drug_ps) :
                                                   CoordinateFrameReader(CoordinateFrame().data());
  const ChemicalFeatures drug_feat = (files_exist) ? ChemicalFeatures(&drug_ag, drug_cfr,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph trpc_ag = (files_exist) ? AtomGraph(trpc_top_path) : AtomGraph();
  PhaseSpace trpc_ps = (files_exist) ? PhaseSpace(trpc_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  CoordinateFrameReader trpc_cfr = (files_exist) ? CoordinateFrameReader(trpc_ps) :
                                                   CoordinateFrameReader(CoordinateFrame().data());
  const ChemicalFeatures trpc_feat = (files_exist) ? ChemicalFeatures(&trpc_ag, trpc_cfr,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();

  // Compute the bond and angle energies.  These quantities should not change when manipulating
  // rotatable bonds.
  ScoreCard sc(1);
  const ValenceKit<double> vk = drug_ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = drug_ag.getChemicalDetailsKit();
  const double drug_bond_e = evaluateBondTerms(vk, drug_cfr, &sc, 0);
  const double drug_angl_e = evaluateAngleTerms(vk, drug_cfr, &sc, 0);
  section(1);
  
  // Rotate about various bonds.  This will generate all sorts of clashes, but bond and angle
  // energies should be unaffected.
  const std::vector<RotatorGroup> drug_rt = drug_feat.getRotatableBondGroups();
  const int ndrug_rt = drug_rt.size();
  const PhaseSpace drug_ps_copy(drug_ps);
  PhaseSpaceWriter dpsw = drug_ps.data();
  const PhaseSpaceReader dpsr = drug_ps_copy.data();
  std::vector<double> drug_rot_crd(3 * (cdk.mol_limits[1] - cdk.mol_limits[0]) * ndrug_rt);
  std::vector<double> drug_repos_dev(ndrug_rt);
  int rcpos = 0;
  for (int i = 0; i < ndrug_rt; i++) {
    rotateAboutBond(&drug_ps, drug_rt[i].root_atom, drug_rt[i].pivot_atom,
                    drug_rt[i].rotatable_atoms, 2.0 * omni::symbols::pi / 3.0);

    // Record the positions of atoms in the first molecule
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      drug_rot_crd[rcpos] = dpsw.xcrd[cdk.mol_contents[j]];
      rcpos++;
      drug_rot_crd[rcpos] = dpsw.ycrd[cdk.mol_contents[j]];
      rcpos++;
      drug_rot_crd[rcpos] = dpsw.zcrd[cdk.mol_contents[j]];
      rcpos++;
    }
    rotateAboutBond(&drug_ps, drug_rt[i].root_atom, drug_rt[i].pivot_atom,
                    drug_rt[i].rotatable_atoms, -2.0 * omni::symbols::pi / 3.0);

    // Check that the molecule was returned to its original state
    double rmsd = 0.0;
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      const size_t jatom = cdk.mol_contents[j];
      const double dx = dpsw.xcrd[jatom] - dpsr.xcrd[jatom];
      const double dy = dpsw.ycrd[jatom] - dpsr.ycrd[jatom];
      const double dz = dpsw.zcrd[jatom] - dpsr.zcrd[jatom];
      rmsd += (dx * dx) + (dy * dy) + (dz * dz);
    }
    rmsd = sqrt(rmsd / static_cast<double>(cdk.mol_limits[1] - cdk.mol_limits[0]));
    drug_repos_dev[i] = rmsd;
  }
  const std::string rcrd_snapshot = base_iso_path + osc + "rotated_coords.m";
  const bool snps_exist = (getDrivePathType(rcrd_snapshot) == DrivePathType::FILE);
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  snapshot(rcrd_snapshot, polyNumericVector(drug_rot_crd), std::string("rot_iso"), 1.0e-6,
           "Coordinate obtained from rotation of selected bonds do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE,
           do_snps);
  check(drug_repos_dev, RelationalOperator::EQUAL, std::vector<double>(ndrug_rt, 0.0), "Reversing "
        "the rotational operations of selected bonds does not return the molecule to its original "
        "state.", do_tests);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
