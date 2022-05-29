#include "../../src/FileManagement/file_listing.h"
#include "../../src/ForceField/forcefield_element.h"
#include "../../src/ForceField/forcefield_enumerators.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::errors;
using namespace omni::mm;
using namespace omni::modeling;
using namespace omni::namelist;
using namespace omni::parse;
using namespace omni::topology;
using namespace omni::trajectory;
using namespace omni::testing;

int main(const int argc, const char* argv[]) {

  // Initialize the test environment
  TestEnvironment oe(argc, argv);

  // Section 1: test bond modifications
  section("Test bond parameter modifications");

  // Section 2: test angle modifications
  section("Test angle parameter modifications");

  // Section 3: test dihedral modifications
  section("Test dihedral parameter modifications");

  // Load a series of topologies, then modify them a bit at a time
  const char osc = osSeparator();
  const std::string topology_home   = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string trajectory_home = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string chemistry_home  = oe.getOmniSourcePath() + osc + "test" + osc + "Chemistry";
  const std::string  mol1_top_name  = topology_home + osc + "stereo_L1.top";
  const std::string  mol2_top_name  = topology_home + osc + "stereo_L1_vs.top";
  const std::string  mol3_top_name  = topology_home + osc + "symmetry_L1.top";
  const std::string  mol4_top_name  = topology_home + osc + "symmetry_L1_vs.top";
  const std::string  mol5_top_name  = topology_home + osc + "bromobenzene.top";
  const std::string  mol6_top_name  = topology_home + osc + "bromobenzene_vs.top";
  const std::string  mol7_top_name  = chemistry_home + osc + "lig1_c8h8.top";
  const std::string  mol8_top_name  = chemistry_home + osc + "lig2_c8h8.top";
  const std::string  mol9_top_name  = chemistry_home + osc + "lig3_1fsj.top";
  const std::string mol10_top_name  = chemistry_home + osc + "morphine_like.top";
  const std::string  mol1_crd_name  = trajectory_home + osc + "stereo_L1.inpcrd";
  const std::string  mol2_crd_name  = trajectory_home + osc + "stereo_L1_vs.inpcrd";
  const std::string  mol3_crd_name  = trajectory_home + osc + "symmetry_L1.inpcrd";
  const std::string  mol4_crd_name  = trajectory_home + osc + "symmetry_L1_vs.inpcrd";
  const std::string  mol5_crd_name  = trajectory_home + osc + "bromobenzene.inpcrd";
  const std::string  mol6_crd_name  = trajectory_home + osc + "bromobenzene_vs.inpcrd";
  const std::string  mol7_crd_name  = chemistry_home + osc + "lig1_c8h8.inpcrd";
  const std::string  mol8_crd_name  = chemistry_home + osc + "lig2_c8h8.inpcrd";
  const std::string  mol9_crd_name  = chemistry_home + osc + "lig3_1fsj.inpcrd";
  const std::string mol10_crd_name  = chemistry_home + osc + "morphine_like.inpcrd";
  const std::vector<std::string> all_topologies  = {
    mol1_top_name,  mol2_top_name,  mol3_top_name, mol4_top_name,  mol5_top_name,
    mol6_top_name,  mol7_top_name,  mol8_top_name, mol9_top_name, mol10_top_name };
  const std::vector<std::string> all_coordinates = {
    mol1_crd_name,  mol2_crd_name,  mol3_crd_name,  mol4_crd_name,  mol5_crd_name,
    mol6_crd_name,  mol7_crd_name,  mol8_crd_name,  mol9_crd_name, mol10_crd_name };
  bool files_exist = true;
  const size_t system_count = all_topologies.size();
  for (size_t i = 0; i < system_count; i++) {
    files_exist = (files_exist && (getDrivePathType(all_topologies[i]) == DrivePathType::FILE));
    files_exist = (files_exist && (getDrivePathType(all_coordinates[i]) == DrivePathType::FILE));
  }
  if (files_exist == false) {
    rtWarn("Files for various small molecule structures in " + topology_home + ", " +
           trajectory_home + ", and " + chemistry_home + " were not found.  Check the "
           "${OMNI_SOURCE} environment variable, currently set to " + oe.getOmniSourcePath() +
           ", for validity.  Tests will be skipped.", "test_topology_mods");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<AtomGraph> ag_list;
  std::vector<PhaseSpace> ps_list;
  std::vector<StaticExclusionMask> se_list;
  if (files_exist) {
    ag_list.reserve(system_count);
    ps_list.reserve(system_count);
    for (size_t i = 0; i < system_count; i++) {
      ag_list.emplace_back(all_topologies[i], ExceptionResponse::SILENT);
      ps_list.emplace_back(all_coordinates[i]);
    }
    se_list.reserve(system_count);
    for (size_t i = 0; i < system_count; i++) {    
      se_list.emplace_back(&ag_list[i]);
    }
  }
  else {
    ag_list.resize(system_count);
    ps_list.resize(system_count);
    se_list.resize(system_count);
  }

  // Evaluate energies with the canonical topologies
  ScoreCard sc_orig(system_count);
  MinimizeControls mincon;
  mincon.setSteepestDescentCycles(25);
  mincon.setTotalCycles(50);
  for (size_t i = 0; i < system_count; i++) {
    RestraintApparatus ra(&ag_list[i]);
    ScoreCard min_sc = minimize(&ps_list[i], &ag_list[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_list[i], &sc_orig, ag_list[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_e0 = sc_orig.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_e0 = sc_orig.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);
  ForceFieldElement ca_ha_bond(ParameterKind::BOND, stringToChar4("ca"), stringToChar4("ha"));
  ForceFieldElement bl_bm_bond(ParameterKind::BOND, stringToChar4("bl"), stringToChar4("bm"));
  ca_ha_bond.setStiffness(100.0);
  ca_ha_bond.setEquilibrium(1.5);
  bl_bm_bond.setStiffness(120.0);
  bl_bm_bond.setEquilibrium(1.7);
  std::vector<AtomGraph> ag_mods;
  std::vector<PhaseSpace> ps_mods;
  for (size_t i = 0; i < system_count; i++) {
    ag_mods.emplace_back(ag_list[i]);
    ps_mods.emplace_back(ps_list[i]);
  }
  ScoreCard sc_mods(system_count);
  for (size_t i = 0; i < system_count; i++) {
    ca_ha_bond.apply(&ag_mods[i]);
    bl_bm_bond.apply(&ag_mods[i]);
    RestraintApparatus ra(&ag_mods[i]);
    ScoreCard min_sc = minimize(&ps_mods[i], &ag_mods[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_mods[i], &sc_mods, ag_mods[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_em = sc_mods.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_em= sc_mods.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_em =
    sc_mods.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_em =
    sc_mods.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
