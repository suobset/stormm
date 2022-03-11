#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::double2;
using omni::constants::ExceptionResponse;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::parse::NumberFormat;
using omni::parse::polyNumericVector;
using omni::topology::AtomGraph;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::PhaseSpace;
using omni::trajectory::TrajectoryKind;
using namespace omni::energy;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Valence energy evaluation");

  // Section 2
  section("Valence force evaluation");

  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getOmniSourcePath() + osc + "test" + osc + "Potential";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const bool systems_exist = (getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein, the DHFR globular protein (with CHARMM potential "
           "details), and alanine dipeptide (with ff19SB) were not found.  These files should be "
           "found in the ${OMNI_SOURCE}/test/Topology and ${OMNI_SOURCE}/test/Trajectory "
           "directories.  Check the $OMNI_SOURCE environment variable.  A number of tests will be "
           "skipped.", "test_valence_evaluation");
  }
  AtomGraph trpcage_ag, dhfr_ag, alad_ag;
  PhaseSpace trpcage_ps, dhfr_ps, alad_ps;
  if (systems_exist) {
    trpcage_ag.buildFromPrmtop(trpcage_top_name, ExceptionResponse::SILENT);
    trpcage_ps.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  ScoreCard all_systems_sc(3);
  const int trpcage_idx = 0;
  const int dhfr_idx = 1;
  const int alad_idx = 2;
  const TrajectoryKind tkind = TrajectoryKind::FORCES;
  trpcage_ps.initializeForces();
  const double trpcage_bond_e = evaluateBondTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_bond_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_angl_e = evaluateAngleTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                   EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_angl_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double2 trpcage_dihe_e = evaluateDihedralTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                       EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_dihe_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_ubrd_e = evaluateUreyBradleyTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                         EvaluateForce::YES, trpcage_idx);
  const double trpcage_cimp_e = evaluateCharmmImproperTerms(trpcage_ag, &trpcage_ps,
                                                            &all_systems_sc, EvaluateForce::YES,
                                                            trpcage_idx);
  const double trpcage_cmap_e = evaluateCmapTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  dhfr_ps.initializeForces();
  const double dhfr_bond_e = evaluateBondTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                               EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_bond_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_angl_e = evaluateAngleTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_angl_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double2 dhfr_dihe_e = evaluateDihedralTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                    EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_dihe_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_ubrd_e = evaluateUreyBradleyTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                      EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_ubrd_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_cimp_e = evaluateCharmmImproperTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                         EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_cimp_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_cmap_e = evaluateCmapTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                               EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_cmap_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_bond_e = evaluateBondTerms(alad_ag, &alad_ps, &all_systems_sc,
                                               EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_bond_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_angl_e = evaluateAngleTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_angl_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double2 alad_dihe_e = evaluateDihedralTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                    EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_dihe_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_ubrd_e = evaluateUreyBradleyTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                      EvaluateForce::YES, alad_idx);
  const double alad_cimp_e = evaluateCharmmImproperTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                         EvaluateForce::YES, alad_idx);
  const double alad_cmap_e = evaluateCmapTerms(alad_ag, &alad_ps, &all_systems_sc,
                                               EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_cmap_frc = alad_ps.getInterlacedCoordinates(tkind);

  // Collect pertinent energy results
  const std::vector<double> trpcage_acc = all_systems_sc.reportInstantaneousStates(trpcage_idx);
  const std::vector<double> dhfr_acc    = all_systems_sc.reportInstantaneousStates(dhfr_idx);
  const std::vector<double> alad_acc    = all_systems_sc.reportInstantaneousStates(alad_idx);

  const std::vector<double> trpcage_valence_energy_answer =
    {   12.38440829,   69.48736586, -142.41607189,    0.77975874,    0.00000000,    0.00000000,
         0.00000000 };
  const std::vector<double> dhfr_valence_energy_answer =
    {  147.47163346,  439.91217087,  754.04779195,    0.00000000,   31.87732134,   18.58590078,
       -85.10501247 };
  const std::vector<double> alad_valence_energy_answer =
    {    0.38340254,    0.54479669,    2.42761153,    0.00708592,    0.00000000,    0.00000000,
        -0.38861139 };
  const std::vector<double> trpcage_valence_energy_result =
    { trpcage_bond_e, trpcage_angl_e, trpcage_dihe_e.x, trpcage_dihe_e.y, trpcage_ubrd_e,
      trpcage_cimp_e, trpcage_cmap_e };
  const std::vector<double> dhfr_valence_energy_result =
    { dhfr_bond_e, dhfr_angl_e, dhfr_dihe_e.x, dhfr_dihe_e.y, dhfr_ubrd_e, dhfr_cimp_e,
      dhfr_cmap_e };
  const std::vector<double> alad_valence_energy_result =
    { alad_bond_e, alad_angl_e, alad_dihe_e.x, alad_dihe_e.y, alad_ubrd_e, alad_cimp_e,
      alad_cmap_e };
  section(1);
  check(trpcage_valence_energy_result, RelationalOperator::EQUAL, trpcage_valence_energy_answer,
        "Trp-cage valence energies do not meet expectations.", do_tests);
  check(dhfr_valence_energy_result, RelationalOperator::EQUAL, dhfr_valence_energy_answer,
        "DHFR valence energies do not meet expectations.", do_tests);
  check(alad_valence_energy_result, RelationalOperator::EQUAL, alad_valence_energy_answer,
        "Alanine dipeptide valence energies do not meet expectations.", do_tests);

  // Re-compute valence energies with a CoordinateFrame abstract (no force computations)
  ScoreCard secondary_sc(3);
  const CoordinateFrameReader trpcage_cfr(trpcage_ps);
  const CoordinateFrameReader dhfr_cfr(dhfr_ps);
  const CoordinateFrameReader alad_cfr(alad_ps);
  const double2 trpcage_dihe_energy_ii = evaluateDihedralTerms(trpcage_ag, trpcage_cfr,
                                                               &secondary_sc, trpcage_idx);
  const std::vector<double> trpcage_valence_energy_ii = {
    evaluateBondTerms(trpcage_ag, trpcage_cfr, &secondary_sc, trpcage_idx),
    evaluateAngleTerms(trpcage_ag, trpcage_cfr, &secondary_sc, trpcage_idx),
    trpcage_dihe_energy_ii.x, trpcage_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(trpcage_ag, trpcage_cfr, &secondary_sc, trpcage_idx),
    evaluateCharmmImproperTerms(trpcage_ag, trpcage_cfr, &secondary_sc, trpcage_idx),
    evaluateCmapTerms(trpcage_ag, trpcage_cfr, &secondary_sc, trpcage_idx)
  };
  const double2 dhfr_dihe_energy_ii = evaluateDihedralTerms(dhfr_ag, dhfr_cfr, &secondary_sc,
                                                            dhfr_idx);
  const std::vector<double> dhfr_valence_energy_ii = {
    evaluateBondTerms(dhfr_ag, dhfr_cfr, &secondary_sc, dhfr_idx),
    evaluateAngleTerms(dhfr_ag, dhfr_cfr, &secondary_sc, dhfr_idx),
    dhfr_dihe_energy_ii.x, dhfr_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(dhfr_ag, dhfr_cfr, &secondary_sc, dhfr_idx),
    evaluateCharmmImproperTerms(dhfr_ag, dhfr_cfr, &secondary_sc, dhfr_idx),
    evaluateCmapTerms(dhfr_ag, dhfr_cfr, &secondary_sc, dhfr_idx)
  };
  const double2 alad_dihe_energy_ii = evaluateDihedralTerms(alad_ag, alad_cfr, &secondary_sc,
                                                            alad_idx);
  const std::vector<double> alad_valence_energy_ii = {
    evaluateBondTerms(alad_ag, alad_cfr, &secondary_sc, alad_idx),
    evaluateAngleTerms(alad_ag, alad_cfr, &secondary_sc, alad_idx),
    alad_dihe_energy_ii.x, alad_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(alad_ag, alad_cfr, &secondary_sc, alad_idx),
    evaluateCharmmImproperTerms(alad_ag, alad_cfr, &secondary_sc, alad_idx),
    evaluateCmapTerms(alad_ag, alad_cfr, &secondary_sc, alad_idx)
  };
  check(trpcage_valence_energy_ii, RelationalOperator::EQUAL, trpcage_valence_energy_answer,
        "Trp-cage valence energies do not meet expectations when computed with a CoordinateFrame "
        "abstract.", do_tests);
  check(dhfr_valence_energy_ii, RelationalOperator::EQUAL, dhfr_valence_energy_answer,
        "DHFR valence energies do not meet expectations when computed with a CoordinateFrame "
        "abstract.", do_tests);
  check(alad_valence_energy_ii, RelationalOperator::EQUAL, alad_valence_energy_answer,
        "Alanine dipeptide valence energies do not meet expectations when computed with a "
        "CoordinateFrame abstract.", do_tests);
  
  // Check for the existence of snapshot files
  const std::string trpcage_snapshot(base_ptl_name + osc + "trpcage_details.m");
  const std::string dhfr_snapshot(base_ptl_name + osc + "dhfr_details.m");
  const std::string alad_snapshot(base_ptl_name + osc + "ala_dipeptide_details.m");
  const bool snps_exist = (getDrivePathType(trpcage_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(alad_snapshot) == DrivePathType::FILE);
  const TestPriority snap_check = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot files " + alad_snapshot + " were not found.  These files contain reference "
           "forces for checking the valence energy derivative calculations.  Check that the "
           "${OMNI_SOURCE} environment variable is set properly so that these snapshot files may "
           "be found in ${OMNI_SOURCE}/test/Potential/.  Subsequent tests will be skipped until "
           "these reference files are available.", "test_valence_evaluation");
  }
  section(2);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_bond_frc), "trpcage_bond",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic bond stretching interactions in the "
           "Trp-cage system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_angl_frc), "trpcage_angl",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic angle bending interactions in the "
           "Trp-cage system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::APPEND, snap_check);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_dihe_frc), "trpcage_dihe",
           NumberFormat::SCIENTIFIC, "Forces due to cosine-based dihedral (proper and improper) "
           "interactions in the Trp-cage system do not meet expectations.", oe.takeSnapshot(),
           1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_bond_frc), "dhfr_bond", NumberFormat::SCIENTIFIC,
           "Forces due to harmonic bond stretching interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::OVERWRITE,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_angl_frc), "dhfr_angl", NumberFormat::SCIENTIFIC,
           "Forces due to harmonic angle bending interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_dihe_frc), "dhfr_dihe", NumberFormat::SCIENTIFIC,
           "Forces due to cosine-based dihedral interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_ubrd_frc), "dhfr_ubrd", NumberFormat::SCIENTIFIC,
           "Forces due to Urey-Bradley interactions in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_cimp_frc), "dhfr_cimp", NumberFormat::SCIENTIFIC,
           "Forces due to CHARMM improper terms in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_cmap_frc), "dhfr_cmap", NumberFormat::SCIENTIFIC,
           "Forces due to CMAP terms in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_bond_frc), "ala_dipeptide_bond_forces",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic bond stretching interactions in the "
           "Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_angl_frc), "ala_dipeptide_angl_forces",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic angle bending interactions in the "
           "Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_dihe_frc), "ala_dipeptide_dihe_forces",
           NumberFormat::SCIENTIFIC, "Forces due to cosine-based dihedral (proper and improper) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_cmap_frc), "ala_dipeptide_cmap_forces",
           NumberFormat::SCIENTIFIC, "Forces due to ff19SB CMAP terms in the Ala dipeptide system "
           "do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);

  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;
}
