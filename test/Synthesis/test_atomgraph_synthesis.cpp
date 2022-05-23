#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/eval_synthesis.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::constants::verytiny;
using omni::data_types::double2;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::energy::evalSyValenceEnergy;
using omni::energy::EvaluateForce;
using omni::energy::ScoreCard;
using omni::energy::StateVariable;
using omni::errors::rtWarn;
using omni::random::Ran2Generator;
using omni::restraints::RestraintApparatus;
using namespace omni::synthesis;
using namespace omni::topology;
using namespace omni::trajectory;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Compare the forces from a single system and one member of a multi-system representation.
// Return the maximum absolute difference in Cartesian force components on any atom.
//
// Arguments:
//   ps:     A single system's coordinates, velocities, and forces
//   psy:    The multi-system object
//   sysid:  System ID within the multi-system object
//-------------------------------------------------------------------------------------------------
double getForceDeviation(const PhaseSpace &ps, const PhaseSpaceSynthesis *psy, int sysid) {
  const CoordinateFrame sys_frc = psy->exportCoordinates(sysid, TrajectoryKind::FORCES);
  const std::vector<double> frc_a = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> frc_b = sys_frc.getInterlacedCoordinates();
  return maxAbsoluteDifference(frc_a, frc_b);
}

//-------------------------------------------------------------------------------------------------
// Evaluate basic force field terms and CHARMM force field extensions evaluated with a PhaseSpace-
// and AtomGraph-Synthesis representation of many systems.  Check energies and forces reported
// against a simpler calculation of the corresponding quantity.
//
// Arguments:
//   poly_ag:   The synthesis of topologies
//   poly_ps:   The synthesis of coordinate / velocity / force representations
//   do_tests:  Indication that tests are possible
//-------------------------------------------------------------------------------------------------
void checkSynthesis(const AtomGraphSynthesis &poly_ag, PhaseSpaceSynthesis *poly_ps,
                    const TestPriority do_tests) {
  
  // Get the valence abstract and prepare for energy calculations
  SyValenceKit<double> syvk = poly_ag.getDoublePrecisionValenceKit();
  ScoreCard sc(poly_ps->getSystemCount(), 1, 32);
  poly_ps->initializeForces();
  evalSyValenceEnergy<double>(syvk, poly_ps->data(), &sc, EvaluateForce::YES, VwuTask::BOND,
                              VwuGoal::ACCUMULATE_FORCES, 0);
  const int nsys = poly_ps->getSystemCount();
  Approx error_limits(std::vector<double>(nsys, 0.0), ComparisonType::ABSOLUTE, verytiny);
  std::vector<double> bond_nrg, bond_nrg_answer, bond_frc_deviations;
  ScoreCard tmp_sc(1, 1, 32);    
  for (int i = 0; i < nsys; i++) {
    bond_nrg.push_back(sc.reportInstantaneousStates(StateVariable::BOND, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    bond_nrg_answer.push_back(evaluateBondTerms(poly_ag.getTopologyPointer(i), &psi, &tmp_sc,
                                                EvaluateForce::YES));
    bond_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.0e-5));
  check(bond_nrg, RelationalOperator::EQUAL, Approx(bond_nrg_answer).margin(5.5e-7), "Bond "
        "energies computed using the synthesis methods are inconsistent with those computed using "
        "a simpler approach.", do_tests);
  check(bond_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to bond "
        "interactions are inconsistent with those computed using a simpler approach.", do_tests);
  poly_ps->initializeForces();
  evalSyValenceEnergy<double>(syvk, poly_ps->data(), &sc, EvaluateForce::YES, VwuTask::ANGL,
                              VwuGoal::ACCUMULATE_FORCES, 0);
  std::vector<double> angl_nrg, angl_nrg_answer, angl_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    angl_nrg.push_back(sc.reportInstantaneousStates(StateVariable::ANGLE, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    angl_nrg_answer.push_back(evaluateAngleTerms(poly_ag.getTopologyPointer(i), &psi, &tmp_sc,
                                                 EvaluateForce::YES));
    angl_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 5.0e-6));
  check(angl_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to harmonic "
        "angle interactions are inconsistent with those computed using a simpler approach.",
        do_tests);
  check(angl_nrg, RelationalOperator::EQUAL, Approx(angl_nrg_answer).margin(3.1e-7),
        "Harmonic angle energies computed using the synthesis methods are inconsistent with those "
        "computed using a simpler approach.", do_tests);
  poly_ps->initializeForces();
  evalSyValenceEnergy<double>(syvk, poly_ps->data(), &sc, EvaluateForce::YES, VwuTask::DIHE,
                              VwuGoal::ACCUMULATE_FORCES, 0);
  std::vector<double> dihe_nrg, impr_nrg, dihe_nrg_answer, impr_nrg_answer, dihe_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    dihe_nrg.push_back(sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, i));
    impr_nrg.push_back(sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    const double2 du = evaluateDihedralTerms(poly_ag.getTopologyPointer(i), &psi, &tmp_sc,
                                             EvaluateForce::YES);
    dihe_nrg_answer.push_back(du.x);
    impr_nrg_answer.push_back(du.y);
    dihe_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.5e-6));
  check(dihe_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "cosine-based dihedral interactions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(dihe_nrg, RelationalOperator::EQUAL, Approx(dihe_nrg_answer).margin(1.0e-7),
        "Cosine-based dihedral energies computed using the synthesis methods are inconsistent "
        "with those computed using a simpler approach.", do_tests);
  check(impr_nrg, RelationalOperator::EQUAL, Approx(impr_nrg_answer).margin(2.5e-9),
        "Cosine-based improper dihedral energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);
  poly_ps->initializeForces();
  evalSyValenceEnergy<double>(syvk, poly_ps->data(), &sc, EvaluateForce::YES, VwuTask::INFR14,
                              VwuGoal::ACCUMULATE_FORCES, 0);
  std::vector<double> qq14_nrg, lj14_nrg, qq14_nrg_answer, lj14_nrg_answer, attn14_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    qq14_nrg.push_back(sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC_ONE_FOUR, i));
    lj14_nrg.push_back(sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    const double2 du = evaluateAttenuated14Terms(poly_ag.getTopologyPointer(i), &psi, &tmp_sc,
                                                 EvaluateForce::YES, EvaluateForce::YES);
    qq14_nrg_answer.push_back(du.x);
    lj14_nrg_answer.push_back(du.y);
    attn14_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 4.5e-6));
  check(attn14_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "attenuated 1:4 interactions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(qq14_nrg, RelationalOperator::EQUAL, Approx(qq14_nrg_answer).margin(1.5e-7),
        "Attenuated 1:4 electrostatic energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);
  check(lj14_nrg, RelationalOperator::EQUAL, Approx(lj14_nrg_answer).margin(2.5e-7),
        "Attenuated 1:4 van-der Waals energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);
  poly_ps->initializeForces();
  evalSyValenceEnergy<double>(syvk, poly_ps->data(), &sc, EvaluateForce::YES, VwuTask::CMAP,
                              VwuGoal::ACCUMULATE_FORCES, 0);
  std::vector<double> cmap_nrg, cmap_nrg_answer, cmap_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    cmap_nrg.push_back(sc.reportInstantaneousStates(StateVariable::CMAP, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    cmap_nrg_answer.push_back(evaluateCmapTerms(poly_ag.getTopologyPointer(i), &psi, &tmp_sc,
                                                EvaluateForce::YES));
    cmap_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.0e-6));
  check(cmap_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to CMAP "
        "energy surface contributions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(cmap_nrg, RelationalOperator::EQUAL, Approx(cmap_nrg_answer).margin(3.1e-7),
        "CMAP energies computed using the synthesis methods are inconsistent with those "
        "computed using a simpler approach.", do_tests);
}

//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);
  StopWatch timer;

  // Section 1
  section("Test AtomGraphSynthesis layout");

  // Section 2
  section("Test molecular mechanics potential calculations");

  // Section 3
  section("Traps for bad input");
  
  // Create some vectors of random numbers, then upload them and test what happens when perturbing
  // atomic coordinates by these numbers.
  const char osc = osSeparator();
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_nbfix_top_name = base_top_name + osc + "trpcage_in_water_nbfix.top";
  const std::string ubiquitin_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string drug_top_name = base_top_name + osc + "drug_example_vs.top";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const bool files_exist = (getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_nbfix_top_name) == DrivePathType::FILE &&
                            getDrivePathType(ubiquitin_top_name) == DrivePathType::FILE &&
                            getDrivePathType(drug_top_name) == DrivePathType::FILE &&
                            getDrivePathType(brbz_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag, trpcage2_ag, trpcage3_ag, nbfix_ag, ubiquitin_ag;
  AtomGraph drug_ag, brbz_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    trpcage2_ag.buildFromPrmtop(trpcage_top_name);
    trpcage3_ag.buildFromPrmtop(trpcage_top_name);
    nbfix_ag.buildFromPrmtop(trpcage_nbfix_top_name);
    ubiquitin_ag.buildFromPrmtop(ubiquitin_top_name);
    drug_ag.buildFromPrmtop(drug_top_name);
    brbz_ag.buildFromPrmtop(brbz_top_name);
  }
  else {
    rtWarn("The topology files for the TIP3P and TIP4P water boxes as well as two versions of the "
           "solvated Trp-cage miniprotein, ubiquitin, and two drug molecules must be available in "
           "${OMNI_SOURCE}/test/ subdirectories Topology and Trajectory, respectively.  Check the "
           "$OMNI_SOURCE environment variable to make sure that it is set properly.  A number of "
           "tests will be skipped.", "test_atomgraph_synthesis");
  }

  // Set one of the Trp-cage systems to have a different topology source name.  This will trick
  // the synthesis generator into thinking that this is a unique system, whereas the third copy
  // is not, leaving six (not seven) total topologies in the synthesis.
  trpcage2_ag.setSource(trpcage_top_name + "_2");

  // Create the synthesis
  const std::vector<AtomGraph*> all_tops = { &tip3p_ag, &tip4p_ag, &trpcage_ag, &trpcage2_ag,
                                             &trpcage3_ag, &nbfix_ag, &ubiquitin_ag, &drug_ag,
                                             &brbz_ag};
  const std::vector<int> system_ids = { 0, 1, 2, 3, 4, 3, 3, 5, 2, 1, 1, 3, 6, 7, 8 };
  AtomGraphSynthesis poly_ag(all_tops, system_ids, ExceptionResponse::SILENT, &timer);
  
  // Check various descriptors
  section(1);
  check(poly_ag.getAtomCount(), RelationalOperator::EQUAL, 50570, "The topology synthesis does "
        "not contain the expected number of atoms.", do_tests);
  check(poly_ag.getVirtualSiteCount(), RelationalOperator::EQUAL, 1264, "The topology synthesis "
        "does not contain the expected number of virtual sites.", do_tests);
  std::vector<int> valence_term_counts;
  if (files_exist) {
    valence_term_counts = { poly_ag.getBondTermCount(), poly_ag.getAngleTermCount(),
                            poly_ag.getDihedralTermCount(), poly_ag.getUreyBradleyTermCount(),
                            poly_ag.getCharmmImproperTermCount(), poly_ag.getCmapTermCount() };
  }
  else {
    valence_term_counts = { 0, 0, 0, 0, 0, 0 };
  }
  const std::vector<int> valence_term_answer = { 50608, 6823, 16765, 0, 0, 71 };
  check(valence_term_counts, RelationalOperator::EQUAL, valence_term_answer, "The topology "
        "synthesis contains incorrect numbers of some valence terms.", do_tests);

  // Get the coordinates for all structures
  section(2);
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps, trpcage2_ps, trpcage3_ps, nbfix_ps, ubiquitin_ps;
  PhaseSpace drug_ps, brbz_ps;
  const std::string base_crd_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string tip3p_crd_name     = base_crd_name + osc + "tip3p.inpcrd";
  const std::string tip4p_crd_name     = base_crd_name + osc + "tip4p.inpcrd";
  const std::string trpcage_crd_name   = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string ubiquitin_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";
  const std::string drug_crd_name      = base_crd_name + osc + "drug_example_vs.inpcrd";
  const std::string brbz_crd_name      = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  const bool coords_exist = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(tip4p_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(ubiquitin_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(drug_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(brbz_crd_name) == DrivePathType::FILE);
  if (coords_exist) {
    tip3p_ps.buildFromFile(tip3p_crd_name);
    tip4p_ps.buildFromFile(tip4p_crd_name);
    trpcage_ps.buildFromFile(trpcage_crd_name);
    ubiquitin_ps.buildFromFile(ubiquitin_crd_name);
    drug_ps.buildFromFile(drug_crd_name);
    brbz_ps.buildFromFile(brbz_crd_name);
  }
  else {
    rtWarn("Coordinates for the periodic systems needed to accompany the first AtomGraphSynthesis "
           "were not found.  Check the installation and the ${OMNI_SOURCE} environment variable.  "
           "Subsequent tests will be skipped.\n");
  }
  const TestPriority do_per_eval = (coords_exist && files_exist) ? TestPriority::CRITICAL :
                                                                   TestPriority::ABORT;
  
  const std::string base_pept_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Namelist" +
                                         osc + "topol";
  const std::string base_pept_crd_name = oe.getOmniSourcePath() + osc + "test" + osc + "Namelist" +
                                         osc + "coord";
  const std::string brbz_vs_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const std::string brbz_vs_crd_name = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  std::vector<AtomGraph*> ag_list;
  std::vector<PhaseSpace> ps_list;
  ag_list.reserve(poly_ag.getSystemCount());
  ps_list.reserve(poly_ag.getSystemCount());
  for (int i = 0; i < poly_ag.getSystemCount(); i++) {
    AtomGraph *ag_ptr = poly_ag.getTopologyPointer(i);
    ag_list.push_back(ag_ptr);
    if (ag_ptr == &tip3p_ag) {
      ps_list.push_back(tip3p_ps);
    }
    else if (ag_ptr == &tip4p_ag) {
      ps_list.push_back(tip4p_ps);
    }
    else if (ag_ptr == &trpcage_ag || ag_ptr == &trpcage2_ag || ag_ptr == &trpcage3_ag ||
             ag_ptr == &nbfix_ag) {
      ps_list.push_back(trpcage_ps);
    }
    else if (ag_ptr == &ubiquitin_ag) {
      ps_list.push_back(ubiquitin_ps);
    }
    else if (ag_ptr == &drug_ag) {
      ps_list.push_back(drug_ps);
    }
    else if (ag_ptr == &brbz_ag) {
      ps_list.push_back(brbz_ps);
    }
  }
  PhaseSpaceSynthesis poly_ps(ps_list, ag_list);
  checkSynthesis(poly_ag, &poly_ps, do_tests);

  // Prepare some more systems
  const std::string tiso_top_name = base_top_name + osc + "trpcage.top";
  const std::string tiso_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string lig1_top_name = base_top_name + osc + "stereo_L1_vs.top";
  const std::string lig1_crd_name = base_crd_name + osc + "stereo_L1_vs.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const bool new_exist = (getDrivePathType(tiso_top_name) == DrivePathType::FILE &&
                          getDrivePathType(lig1_top_name) == DrivePathType::FILE &&
                          getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                          getDrivePathType(tiso_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(lig1_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(dhfr_crd_name) == DrivePathType::FILE);
  AtomGraph tiso_ag, lig1_ag, dhfr_ag;
  PhaseSpace tiso_ps, lig1_ps, dhfr_ps;
  if (new_exist) {
    tiso_ag.buildFromPrmtop(tiso_top_name, ExceptionResponse::SILENT);
    lig1_ag.buildFromPrmtop(lig1_top_name, ExceptionResponse::SILENT);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    tiso_ps.buildFromFile(tiso_crd_name);
    lig1_ps.buildFromFile(lig1_crd_name);
    dhfr_ps.buildFromFile(dhfr_crd_name);
  }
  else {
    rtWarn("Files corresponding to various systems in isolated boundary conditions were not "
           "found.  Check the ${OMNI_SOURCE} environment variable.  The necessary directories "
           "are the same as for other files needed by this test program.  Subsequent tests will "
           "be skipped.", "test_atomgraph_synthesis");
  }
  std::vector<AtomGraph*> agn_list = { &tiso_ag, &lig1_ag, &dhfr_ag };
  std::vector<PhaseSpace> psn_list = { tiso_ps, lig1_ps, dhfr_ps };
  AtomGraphSynthesis poly_agn(agn_list, { 0, 1, 2 }, ExceptionResponse::SILENT, &timer);
  PhaseSpaceSynthesis poly_psn(psn_list, agn_list);
  checkSynthesis(poly_agn, &poly_psn, do_tests);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
