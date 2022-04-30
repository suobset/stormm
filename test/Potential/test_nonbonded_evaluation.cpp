#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Potential/nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::double2;
using omni::int2;
using omni::int3;
using omni::constants::ExceptionResponse;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::energy::StateVariable;
using omni::energy::StaticExclusionMask;
using omni::errors::rtWarn;
using omni::parse::char4ToString;
using omni::parse::NumberFormat;
using omni::parse::polyNumericVector;
using omni::topology::AtomGraph;
using omni::topology::NonbondedKit;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::PhaseSpace;
using omni::trajectory::TrajectoryKind;
using namespace omni::energy;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Check that the sending atom exclusion masks match the receiving atom masks for each tile of
// a StaticExclusionMask object.
//
// Arguments:
//   se:  Non-bonded all-to-all exclusion mask for the system
//-------------------------------------------------------------------------------------------------
int2 checkReflexiveMarks(const StaticExclusionMask &se) {
  int2 result = { 0, 0 };
  const int natom = se.getAtomCount();
  const int nsptile = (natom + 255) / 256;
  for (int sti = 0; sti < nsptile; sti++) {
    const int ni_tile = (sti < nsptile - 1) ? 16 : (natom - (sti * 256)) / 16;
    for (int stj = 0; stj < nsptile; stj++) {
      const int nj_tile = (stj < nsptile - 1) ? 16 : (natom - (stj * 256)) / 16;
      for (int ti = 0; ti < ni_tile; ti++) {
        for (int tj = 0; tj < nj_tile; tj++) {

          // Get the mask for this tile
          std::vector<uint> cmask = se.getTileExclusions(sti, stj, ti, tj);

          // Loop over the first 16 atoms, assemble the corresponding unsigned ints for the
          // second 16, and compare.
          for (int i = 0; i < 16; i++) {
            uint rec_mask = 0x0;
            for (int j = 0; j < 16; j++) {
              rec_mask |= (((cmask[16 + j] >> i) & 0x1) << i);
            }
            if (cmask[i] != rec_mask) {
              for (int j = 0; j < 16; j++) {
                const bool send_excl = ((cmask[i] >> j) & 0x1);
                const bool recv_excl = ((cmask[16 + j] >> i) & 0x1);
                result.x += (send_excl && (recv_excl == false));
                result.y += ((send_excl == false) && recv_excl);
              }
            }
          }
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check all exclusions from the original topology against those marked in a StaticExclusionMask
// object.
//
// Arguments:
//   se:  Non-bonded all-to-all exclusion mask for the system
//-------------------------------------------------------------------------------------------------
int3 checkMarkedExclusions(const StaticExclusionMask &se) {
  const int natom = se.getTopologyPointer()->getAtomCount();
  const NonbondedKit<double> nbk = se.getTopologyPointer()->getDoublePrecisionNonbondedKit();
  int3 n_errors = { 0, 0, 0 };
  for (int i = 0; i < natom; i++) {
    for (int j = 0; j <= i; j++) {
      const bool is_excl = se.testExclusion(i, j);
      if (j == i && is_excl) {
        n_errors.x += 1;
        continue;
      }
      bool topol_excl = false;
      for (int k = nbk.nb11_bounds[i]; k < nbk.nb11_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb11x[k]);
      }
      for (int k = nbk.nb12_bounds[i]; k < nbk.nb12_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb12x[k]);
      }
      for (int k = nbk.nb13_bounds[i]; k < nbk.nb13_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb13x[k]);
      }
      for (int k = nbk.nb14_bounds[i]; k < nbk.nb14_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb14x[k]);
      }
      n_errors.y += (topol_excl && (is_excl == false));
      n_errors.z += ((topol_excl == false) && is_excl);
    }
  }
  return n_errors;
}

//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  const int input_timings  = timer.addCategory("Input file parsing");
  const int sem_timings    = timer.addCategory("Static exclusion mask");
  const int attn14_timings = timer.addCategory("Compute 1:4 electrostatics");
  const int nb_timings     = timer.addCategory("Compute non-bonded interactions");
  
  // Section 1
  section("Tests of the StaticExclusionMask object");

  // Section 2
  section("Non-bonded 1:4 energy evaluation");

  // Section 3
  section("Non-bonded energy evaluation");

  // Section 4
  section("Non-bonded force evaluation");

  // Section 5
  section("Fixed-precision energy accumulator performance");

  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getOmniSourcePath() + osc + "test" + osc + "Potential";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string trpw_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpw_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string trpp_top_name = base_top_name + osc + "trpcage_no_z.top";
  const bool systems_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein (with and without water), the DHFR globular "
           "protein (with CHARMM potential details), and alanine dipeptide (with ff19SB) were not "
           "found.  These files should be found in the ${OMNI_SOURCE}/test/Topology and "
           "${OMNI_SOURCE}/test/Trajectory directories.  Check the $OMNI_SOURCE environment "
           "variable.  A number of tests will be skipped.", "test_nonbonded_evaluation");
  }

  // Read topologies and coordinates
  AtomGraph trpi_ag, dhfr_ag, alad_ag, trpw_ag, trpp_ag;
  PhaseSpace trpi_ps, dhfr_ps, alad_ps, trpw_ps;
  if (systems_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpw_ag.buildFromPrmtop(trpw_top_name, ExceptionResponse::SILENT);
    trpw_ps.buildFromFile(trpw_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpp_ag.buildFromPrmtop(trpw_top_name, ExceptionResponse::SILENT);
  }
  timer.assignTime(input_timings);
  
  // Prepare exclusion masks for each system (the constructor will work even if the
  // system has zero atoms)
  const StaticExclusionMask trpi_semask(&trpi_ag);
  const StaticExclusionMask dhfr_semask(&dhfr_ag);
  const StaticExclusionMask alad_semask(&alad_ag);
  const StaticExclusionMask trpw_semask(&trpw_ag);
  timer.assignTime(sem_timings);

  // Check exclusions for three systems against the original topologies
  section(1);
  const std::vector<const StaticExclusionMask*> my_masks = { &trpi_semask, &alad_semask,
                                                             &trpw_semask };
  const std::vector<std::string> my_systems = { "Trp-cage (unsolvated)", "alanine dipeptide",
                                                "Trp-cage (solvated)" };
  for (size_t i = 0; i < my_masks.size(); i++) {
    const int3 n_errors = checkMarkedExclusions(*(my_masks[i]));
    check(n_errors.x == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions along the diagonal in the " + my_systems[i] +
          " system.  These should be handled implicitly by any non-bonded loop.", do_tests);
    check(n_errors.y == 0, "The StaticExclusionMask object failed to record " +
          std::to_string(n_errors.y) + " exclusions marked in the topology for the " +
          my_systems[i] + " system.", do_tests);
    check(n_errors.z == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.z) + " exclusions not in the original topology for the " +
          my_systems[i] + " system.", do_tests);
  }

  // Check the reflexive nature of the exclusion mask tiles
  for (size_t i = 0; i < my_masks.size(); i++) {
    const int2 n_errors = checkReflexiveMarks(*(my_masks[i]));
    check(n_errors.x == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions in the sending atoms but not in the receiving "
          "atoms for the " + my_systems[i] + " system.", do_tests);
    check(n_errors.y == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions in the receiving atoms but not in the sending "
          "atoms for the " + my_systems[i] + " system.", do_tests);
  }

  // Accumulate energies and forces
  ScoreCard all_systems_sc(5), secondary_sc(5);
  const int trpi_idx = 0;
  const int dhfr_idx = 1;
  const int alad_idx = 2;
  const int trpw_idx = 3;
  const int trpp_idx = 4;
  TrajectoryKind tkind = TrajectoryKind::FORCES;

  // Check for the existence of snapshot files
  const std::string trpi_snapshot(base_ptl_name + osc + "trpcage_nb_details.m");
  const std::string dhfr_snapshot(base_ptl_name + osc + "dhfr_nb_details.m");
  const std::string alad_snapshot(base_ptl_name + osc + "ala_dipeptide_nb_details.m");
  const bool snps_exist = (getDrivePathType(trpi_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(alad_snapshot) == DrivePathType::FILE);
  const TestPriority snap_check = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot files " + trpi_snapshot + ", " + dhfr_snapshot + ", and " + alad_snapshot +
           " were not found.  These files contain reference forces for checking the valence "
           "energy derivative calculations.  Check that the ${OMNI_SOURCE} environment variable "
           "is set properly so that these snapshot files may be found in "
           "${OMNI_SOURCE}/test/Potential/.  Subsequent tests will be skipped until these "
           "reference files are available.", "test_nonbonded_evaluation");
  }

  // Compute 1:4 electrostatic and Lennard-Jones forces on the isolated Trp-cage system
  section(2);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpi_14_e = evaluateAttenuated14Terms(trpi_ag, &trpi_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpi_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpi_14_elec_frc = trpi_ps.getInterlacedCoordinates(tkind);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(trpi_ag, &trpi_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, trpi_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpi_14_vdw_frc = trpi_ps.getInterlacedCoordinates(tkind);
  check(trpi_14_e.x, RelationalOperator::EQUAL, Approx(1458.0998129).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpi_14_e.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6), "van-der Waals "
        "1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);
  snapshot(trpi_snapshot, polyNumericVector(trpi_14_elec_frc), "trpcage_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the Trp-cage (isolated boundary conditions) system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::OVERWRITE, snap_check);
  snapshot(trpi_snapshot, polyNumericVector(trpi_14_vdw_frc), "trpcage_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the Trp-cage (isolated boundary conditions) system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);
  timer.assignTime(0);
  const double2 trpi_14_e_ii = evaluateAttenuated14Terms(trpi_ag, CoordinateFrameReader(trpi_ps),
                                                         &secondary_sc, trpi_idx);
  timer.assignTime(attn14_timings);
  check(trpi_14_e_ii.x, RelationalOperator::EQUAL, Approx(1458.0998129).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly when evaluating energy only with a coordinate frame object.", do_tests);
  check(trpi_14_e_ii.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6),
        "van-der Waals 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly when evaluating energy only with a coordinate frame object.", do_tests);
  
  // Compute 1:4 electrostatic and Lennard-Jones forces on the DHFR system
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  const double2 dhfr_14_e = evaluateAttenuated14Terms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      dhfr_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> dhfr_14_elec_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(dhfr_ag, &dhfr_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, dhfr_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> dhfr_14_vdw_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  check(dhfr_14_e.x, RelationalOperator::EQUAL, Approx(6507.3376751).margin(1.0e-6),
        "Electrostatic 1:4 energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  check(dhfr_14_e.y, RelationalOperator::EQUAL, Approx(367.0925927).margin(1.0e-6),
        "van-der Waals 1:4 energy for DHFR (CHARMM force field, with special 1:4 LJ coefficients) "
        "was not computed correctly.", do_tests);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_14_elec_frc), "dhfr_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_14_vdw_frc), "dhfr_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12,
           PrintSituation::APPEND, snap_check);

  // Compute 1:4 electrostatic and Lennard-Jones forces on the alanine dipeptide system
  alad_ps.initializeForces();
  timer.assignTime(0);
  const double2 alad_14_e = evaluateAttenuated14Terms(alad_ag, &alad_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      alad_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> alad_14_elec_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(alad_ag, &alad_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, alad_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> alad_14_vdw_frc = alad_ps.getInterlacedCoordinates(tkind);
  check(alad_14_e.x, RelationalOperator::EQUAL, Approx(46.8072425).margin(1.0e-6), "Electrostatic "
        "1:4 energy for alanine dipeptide (ff19SB force field) was not computed correctly.",
        do_tests);
  check(alad_14_e.y, RelationalOperator::EQUAL, Approx(3.1589264).margin(1.0e-6), "van-der Waals "
        "1:4 energy for alanine dipeptide (ff19SB force field) was not computed correctly.",
        do_tests);
  snapshot(alad_snapshot, polyNumericVector(alad_14_elec_frc), "ala_dipeptide_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::OVERWRITE, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_14_vdw_frc), "ala_dipeptide_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);

  // Compute 1:4 electrostatic and Lennard-Jones forces on the solvated Trp-cage system
  trpw_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpw_14_e = evaluateAttenuated14Terms(trpw_ag, &trpw_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpw_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpw_14_elec_frc = trpw_ps.getInterlacedCoordinates(tkind);
  trpw_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(trpw_ag, &trpw_ps, &secondary_sc, EvaluateForce::YES,
                            EvaluateForce::NO, trpw_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpw_14_vdw_frc = trpw_ps.getInterlacedCoordinates(tkind);
  check(trpw_14_e.x, RelationalOperator::EQUAL, Approx(1458.0996855).margin(1.0e-6),
        "Electrostatic 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly.", do_tests);
  check(trpw_14_e.y, RelationalOperator::EQUAL, Approx(62.6481797).margin(1.0e-6), "van-der Waals "
        "1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);
  timer.assignTime(0);
  const double2 trpw_14_e_ii = evaluateAttenuated14Terms(trpw_ag, CoordinateFrameReader(trpw_ps),
                                                         &all_systems_sc, trpw_idx);
  timer.assignTime(attn14_timings);
  check(trpw_14_e_ii.x, RelationalOperator::EQUAL, Approx(1458.0996855).margin(1.0e-6),
        "Electrostatic 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly when evaluating energy only with a CoordinateFrame abstract.",
        do_tests);
  check(trpw_14_e_ii.y, RelationalOperator::EQUAL, Approx(62.6481797).margin(1.0e-6),
        "van-der Waals 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly when evaluating energy only with a CoordinateFrame abstract.",
        do_tests);

  // Compute 1:4 electrostatic energies on the Trp-cage system with no explicit scaling factors in
  // the topology
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpp_14_e = evaluateAttenuated14Terms(trpp_ag, &trpi_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpp_idx);
  timer.assignTime(attn14_timings);
  trpi_ps.initializeForces();
  check(trpp_14_e.x, RelationalOperator::EQUAL, Approx(1458.0996843).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpp_14_e.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6), "van-der Waals "
        "1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);

  // Compute non-bonded energy for the isolated Trp-cage system
  section(3);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpi_nonb_e = evaluateNonbondedEnergy(trpi_ag, trpi_semask, &trpi_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, trpi_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> trpi_elec_frc = trpi_ps.getInterlacedCoordinates(tkind);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(trpi_ag, trpi_semask, &trpi_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, trpi_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> trpi_vdw_frc = trpi_ps.getInterlacedCoordinates(tkind);
  check(trpi_nonb_e.x, RelationalOperator::EQUAL, Approx(-1816.1478854).margin(1.0e-6),
        "Electrostatic energy for isolated Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpi_nonb_e.y, RelationalOperator::EQUAL, Approx(-119.3500215).margin(1.0e-6),
        "van-der Waals energy for isolated Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);

  // Compute non-bonded energy for the DHFR system
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  const double2 dhfr_nonb_e = evaluateNonbondedEnergy(dhfr_ag, dhfr_semask, &dhfr_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, dhfr_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> dhfr_elec_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(dhfr_ag, dhfr_semask, &dhfr_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, dhfr_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> dhfr_vdw_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  check(dhfr_nonb_e.x, RelationalOperator::EQUAL, Approx(-10036.4655324).margin(1.0e-6),
        "Electrostatic energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  check(dhfr_nonb_e.y, RelationalOperator::EQUAL, Approx(-1009.1558078).margin(1.0e-6),
        "van-der Waals energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  timer.assignTime(0);
  const double2 dhfr_nonb_e_ii = evaluateNonbondedEnergy(dhfr_ag, dhfr_semask,
                                                         CoordinateFrameReader(dhfr_ps),
                                                         &secondary_sc, dhfr_idx);
  timer.assignTime(nb_timings);
  check(dhfr_nonb_e_ii.x, RelationalOperator::EQUAL, Approx(-10036.4655324).margin(1.0e-6),
        "Electrostatic energy for DHFR (CHARMM force field) was not computed correctly when "
        "evaluating energy only with a CoordinateFrame abstract.", do_tests);
  check(dhfr_nonb_e_ii.y, RelationalOperator::EQUAL, Approx(-1009.1558078).margin(1.0e-6),
        "van-der Waals energy for DHFR (CHARMM force field) was not computed correctly when "
        "evaluating energy only with a CoordinateFrame abstract.", do_tests);

  // Compute non-bonded energy for the alanine dipeptide system
  alad_ps.initializeForces();
  timer.assignTime(0);
  const double2 alad_nonb_e = evaluateNonbondedEnergy(alad_ag, alad_semask, &alad_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, alad_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> alad_elec_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(alad_ag, alad_semask, &alad_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, alad_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> alad_vdw_frc = alad_ps.getInterlacedCoordinates(tkind);
  check(alad_nonb_e.x, RelationalOperator::EQUAL, Approx(-78.8463310).margin(1.0e-6),
        "Electrostatic energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly.", do_tests);
  check(alad_nonb_e.y, RelationalOperator::EQUAL, Approx(-1.2452480).margin(1.0e-6),
        "van-der Waals energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly.", do_tests);
  timer.assignTime(0);
  const double2 alad_nonb_e_ii = evaluateNonbondedEnergy(alad_ag, alad_semask,
                                                         CoordinateFrameReader(alad_ps),
                                                         &secondary_sc, alad_idx);
  timer.assignTime(nb_timings);
  check(alad_nonb_e_ii.x, RelationalOperator::EQUAL, Approx(-78.8463310).margin(1.0e-6),
        "Electrostatic energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly when evaluating energy only with a CoordinateFrame abstract.", do_tests);
  check(alad_nonb_e_ii.y, RelationalOperator::EQUAL, Approx(-1.2452480).margin(1.0e-6),
        "van-der Waals energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly when evaluating energy only with a CoordinateFrame abstract.", do_tests);  

  // Check forces for the above three systems (skip the Trp-cage system with no Z numbers or
  // explicit scaling factors (old form of Amber prmtop), as well as the Trp-cage system in water
  // (it's really a periodic box, useful for computing a meaningful 1:4 non-bonded energy sum).
  section(4);
  snapshot(trpi_snapshot, polyNumericVector(trpi_elec_frc), "trpcage_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the Trp-cage (isolated boundary conditions) system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(trpi_snapshot, polyNumericVector(trpi_vdw_frc), "trpcage_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the Trp-cage (isolated boundary conditions) system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_elec_frc), "dhfr_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_vdw_frc), "dhfr_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_elec_frc), "ala_dipeptide_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_vdw_frc), "ala_dipeptide_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);

  // Check the relevant, accumulated energies against the double-precision standard
  section(5);
  const std::vector<double> trpi_acc = all_systems_sc.reportInstantaneousStates(trpi_idx);
  const std::vector<double> dhfr_acc = all_systems_sc.reportInstantaneousStates(dhfr_idx);
  const std::vector<double> alad_acc = all_systems_sc.reportInstantaneousStates(alad_idx);
  const std::vector<double> trpi_nbe_answer = {   1458.0998129,     62.6481811,
                                                 -1816.1478854,   -119.3500215 };
  const std::vector<double> dhfr_nbe_answer = {   6507.3376751,    367.0925927,
                                                -10036.4655324,  -1009.1558078 };
  const std::vector<double> alad_nbe_answer = {     46.8072425,      3.1589264,
                                                   -78.8463310,     -1.2452480 };
  const int qq14_idx = static_cast<size_t>(StateVariable::ELECTROSTATIC_ONE_FOUR);
  const int lj14_idx = static_cast<size_t>(StateVariable::VDW_ONE_FOUR);
  const int qqnb_idx = static_cast<size_t>(StateVariable::ELECTROSTATIC);
  const int ljnb_idx = static_cast<size_t>(StateVariable::VDW);
  const std::vector<double> trpi_acc_result = { trpi_acc[qq14_idx], trpi_acc[lj14_idx],
                                                trpi_acc[qqnb_idx], trpi_acc[ljnb_idx] };
  const std::vector<double> dhfr_acc_result = { dhfr_acc[qq14_idx], dhfr_acc[lj14_idx],
                                                dhfr_acc[qqnb_idx], dhfr_acc[ljnb_idx] };
  const std::vector<double> alad_acc_result = { alad_acc[qq14_idx], alad_acc[lj14_idx],
                                                alad_acc[qqnb_idx], alad_acc[ljnb_idx] };
  check(trpi_acc_result, RelationalOperator::EQUAL, Approx(trpi_nbe_answer).margin(2.0e-5),
        "Trp-cage (isolated boundary conditions) non-bonded energy accumulators do not reproduce "
        "the double-precision targets to within tolerances.", do_tests);
  check(dhfr_acc_result, RelationalOperator::EQUAL, Approx(dhfr_nbe_answer).margin(1.0e-3),
        "DHFR non-bonded energy accumulators do not reproduce the double-precision targets to "
        "within tolerances.", do_tests);
  check(alad_acc_result, RelationalOperator::EQUAL, Approx(alad_nbe_answer).margin(1.0e-5),
        "Alanine dipeptide non-bonded energy accumulators do not reproduce the double-precision "
        "targets to within tolerances.", do_tests);
  timer.assignTime(0);

  // Print results
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
