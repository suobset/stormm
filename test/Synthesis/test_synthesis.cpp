#include <algorithm>
#include <string>
#include <vector>
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Chemistry/atommask.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/mixed_types.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Math/sorting.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Potential/eval_valence_workunit.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::chemistry::AtomMask;
using omni::chemistry::ChemicalFeatures;
using omni::chemistry::MapRotatableGroups;
using omni::constants::ExceptionResponse;
using omni::constants::verytiny;
using omni::data_types::int2;
using omni::data_types::uint2;
using omni::data_types::uint3;
using omni::data_types::double2;
using omni::data_types::double3;
using omni::data_types::double4;
using omni::data_types::ValueWithCounter;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::openOutputFile;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::math::readBitFromMask;
using omni::math::reduceUniqueValues;
using omni::math::UniqueValueHandling;
using omni::namelist::FilesControls;
using omni::parse::TextFile;
using omni::random::Xoroshiro128pGenerator;
using omni::restraints::applyHydrogenBondPreventors;
using omni::restraints::applyPositionalRestraints;
using omni::restraints::BoundedRestraint;
using omni::restraints::RestraintApparatus;
using omni::restraints::RestraintKit;
using omni::structure::distance;
using omni::structure::angle;
using omni::structure::dihedral_angle;
using namespace omni::energy;
using namespace omni::numerics;
using namespace omni::synthesis;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

//-------------------------------------------------------------------------------------------------
// Check the coverage of a simple task form a ValenceWorkUnit.
//
// Arguments:
//   accd:              Accumulator directives: coverage is only counted if the accumulator bit is
//                      set to 1, which means energy is accumulated by the work unit
//   taskid:            List of simple task topological index numbers
//   coverage:          Coverage array (accumulated and returned)
//   range_problem_in:  Current state of the array range indexing problem detector
//-------------------------------------------------------------------------------------------------
bool checkSimpleTaskCoverage(const std::vector<uint> &accd, const std::vector<int> &taskid,
                             std::vector<int> *coverage, const bool range_problem_in) {
  bool range_problem = range_problem_in;
  const int n_items = coverage->size();
  int* cov_ptr = coverage->data();
  for (size_t j = 0; j < taskid.size(); j++) {
    if (readBitFromMask(accd, j) == 0) {
      continue;
    }
    if (taskid[j] >= 0 && taskid[j] < n_items) {
      cov_ptr[taskid[j]] += 1;
    }
    else {
      range_problem = true;
    }
  }
  return range_problem;
}

//-------------------------------------------------------------------------------------------------
// Check that the naive distance between two particles is the distance computed upon re-imaging.
// Return TRUE is this is so, FALSE otherwise.
//
// Arguments:
//   i:     The first particle
//   j:     The second particle
//   cfr:   Coordinates of the particles, plus box information
//-------------------------------------------------------------------------------------------------
bool checkNaiveDistance(const int i, const int j, const CoordinateFrameReader &cfr) {
  const double reim_dist  = distance(i, j, cfr);
  const double naive_dist = distance<double, double>(i, j, cfr.xcrd, cfr.ycrd, cfr.zcrd, nullptr,
                                                     nullptr, UnitCellType::NONE);
  return (fabs(reim_dist - naive_dist) < omni::constants::tiny);
}

//-------------------------------------------------------------------------------------------------
// Run a series of tests using valence work units.
//
// Arguments:
//   top_name:  Name of the topology to use
//   crd_name:  Name of the coordinate file to use
//   oe:        Operating environment information (for error reporting)
//   my_prng:   Random number generator (modified by use inside this function)
//-------------------------------------------------------------------------------------------------
void runValenceWorkUnitTests(const std::string &top_name, const std::string &crd_name,
                             const TestEnvironment &oe, Xoroshiro128pGenerator *my_prng) {
  const bool files_exist = (getDrivePathType(top_name) == DrivePathType::FILE &&
                            getDrivePathType(crd_name) == DrivePathType::FILE);
  if (files_exist == false) {
    rtWarn("The topology and input coordinates for the DHFR system appear to be missing.  Check "
           "the ${OMNI_SOURCE} variable (currently " + oe.getOmniSourcePath() + ") to make sure "
           "that " + top_name + " and " + crd_name + " valid paths.", "test_synthesis");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph ag  = (files_exist) ? AtomGraph(top_name, ExceptionResponse::SILENT) :
                                  AtomGraph();
  PhaseSpace ps = (files_exist) ? PhaseSpace(crd_name) : PhaseSpace();
  const CoordinateFrameReader cfr(ps);

  // Create a set of restraints on the structure, starting from backbone positional restraints
  // and adding things that randomly keep various two-, three-, and four-point meansurements near
  // their observed values.
  std::vector<BoundedRestraint> mol_rstr = applyPositionalRestraints(&ag, cfr, "@N,CA,C,O", 1.4);
  for (size_t i = 0; i < mol_rstr.size(); i++) {
    const double3 trgt = mol_rstr[i].getTargetSite();
    mol_rstr[i].setTargetSite({ trgt.x + 0.5 - my_prng->uniformRandomNumber(),
                                 trgt.y + 0.5 - my_prng->uniformRandomNumber(),
                                 trgt.z + 0.5 - my_prng->uniformRandomNumber() });
  }
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const int solute_extent = ag.getLastSoluteAtom();
  const int min_iinc = std::max(solute_extent / 8, 4);
  for (int i = min_iinc; i < solute_extent; i += min_iinc) {
    for (int j = 0; j < i; j += min_iinc) {
      const double orig_distance = distance(i, j, ps) + 0.5 - my_prng->uniformRandomNumber();
      if (checkNaiveDistance(i, j, cfr) == false) {
        continue;
      }
      mol_rstr.emplace_back(i, j, &ag, 1.0, 1.4, 0.0, orig_distance - 0.25, orig_distance + 0.25,
                            orig_distance + 100.0);
      if (checkNaiveDistance(i, i - 2, cfr) && checkNaiveDistance(j, i - 2, cfr)) {
        const double orig_angle = angle(i, j, i - 2, ps) +
                                  (0.1 * (0.5 - my_prng->uniformRandomNumber()));
        if (orig_angle > 0.1 && orig_angle < 3.0) {
          mol_rstr.emplace_back(i, j, i - 2, &ag, 2.3, 0.9, 0.0, orig_angle - 0.05,
                                orig_angle + 0.02, omni::symbols::pi - omni::constants::tiny);
        }
      }
      if (checkNaiveDistance(i, i - 1, cfr) && checkNaiveDistance(j, i - 1, cfr) &&
          checkNaiveDistance(j, j + 1, cfr)) {
        const double orig_dihedral = dihedral_angle(i, i - 1, j, j + 1, ps) +
                                     (0.8 * (0.5 - my_prng->uniformRandomNumber()));
        mol_rstr.emplace_back(i, i - 1, j, j + 1, &ag, 2.3, 0.9, orig_dihedral - 1.0,
                              orig_dihedral - 0.03, orig_dihedral + 0.03, orig_dihedral + 1.0);
      }
    }
  }

  // Create the restraint apparatus, then the valence work units.
  RestraintApparatus ra(mol_rstr);
  ValenceDelegator vdel(&ag, &ra);
  const std::vector<ValenceWorkUnit> all_vwu = buildValenceWorkUnits(&vdel);
  const int n_vwu = all_vwu.size();

  // Check the coverage in valence work units with an independent tally of each atom and term
  std::vector<int> bond_coverage(vk.nbond, 0);
  std::vector<int> angl_coverage(vk.nangl, 0);
  std::vector<int> dihe_coverage(vk.ndihe, 0);
  std::vector<int> ubrd_coverage(vk.nubrd, 0);
  std::vector<int> cimp_coverage(vk.ncimp, 0);
  std::vector<int> cmap_coverage(vk.ncmap, 0);
  std::vector<int> infr_coverage(vk.ninfr14, 0);
  bool bond_range_problem = false;
  bool angl_range_problem = false;
  bool dihe_range_problem = false;
  bool ubrd_range_problem = false;
  bool cimp_range_problem = false;
  bool cmap_range_problem = false;
  bool infr_range_problem = false;
  for (int i = 0; i < n_vwu; i++) {
    const std::vector<uint2> cbnd_insr = all_vwu[i].getCompositeBondInstructions();
    const std::vector<uint> cbnd_accd  = all_vwu[i].getAccumulationFlags(VwuTask::CBND);
    const std::vector<int> cbnd_taskid = all_vwu[i].getCompositeBondTaskList();
    for (size_t j = 0; j < cbnd_insr.size(); j++) {
      if (readBitFromMask(cbnd_accd, j) == 0) {
        continue;
      }
      if ((cbnd_insr[j].x >> 20) & 0x1) {
        if (cbnd_taskid[j] >= 0 && cbnd_taskid[j] < vk.nubrd) {
          ubrd_coverage[cbnd_taskid[j]] += 1;
        }
        else {
          ubrd_range_problem = true;
        }
      }
      else {
        if (cbnd_taskid[j] >= 0 && cbnd_taskid[j] < vk.nbond) {
          bond_coverage[cbnd_taskid[j]] += 1;
        }
        else {
          bond_range_problem = true;
        }
      }
    }
    angl_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::ANGL),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::ANGL),
                                                 &angl_coverage, angl_range_problem);
    const std::vector<uint3> cdhe_insr = all_vwu[i].getCompositeDihedralInstructions();
    const std::vector<uint> cdhe_accd  = all_vwu[i].getAccumulationFlags(VwuTask::CDHE);
    const std::vector<int2> cdhe_taskid = all_vwu[i].getCompositeDihedralTaskList();
    for (size_t j = 0; j < cdhe_insr.size(); j++) {
      if (readBitFromMask(cdhe_accd, j) == 0) {
        continue;
      }
      if ((cdhe_insr[j].x >> 30) & 0x1) {
        if (cdhe_taskid[j].x >= 0 && cdhe_taskid[j].x < vk.ncimp) {
          cimp_coverage[cdhe_taskid[j].x] += 1;
        }
        else {
          cimp_range_problem = true;
        }
      }
      else {
        if (cdhe_taskid[j].x >= 0 && cdhe_taskid[j].x < vk.ndihe) {
          dihe_coverage[cdhe_taskid[j].x] += 1;
        }
        else {
          dihe_range_problem = true;          
        }
        if (cdhe_taskid[j].y >= 0) {
          if (cdhe_taskid[j].y < vk.ndihe) {
            dihe_coverage[cdhe_taskid[j].y] += 1;
          }
          else {
            dihe_range_problem = true;          
          }
        }
      }
    }
    cmap_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::CMAP),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::CMAP),
                                                 &cmap_coverage, cmap_range_problem);
    infr_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::INFR14),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::INFR14),
                                                 &infr_coverage, infr_range_problem);
  }
  const std::vector<int> bond_coverage_answer(vk.nbond, 1);
  const std::vector<int> angl_coverage_answer(vk.nangl, 1);
  const std::vector<int> dihe_coverage_answer(vk.ndihe, 1);
  const std::vector<int> ubrd_coverage_answer(vk.nubrd, 1);
  const std::vector<int> cimp_coverage_answer(vk.ncimp, 1);
  const std::vector<int> cmap_coverage_answer(vk.ncmap, 1);
  const std::vector<int> infr_coverage_answer(vk.ninfr14, 1);
  if (vk.nbond > 0) {
    check(bond_range_problem == false, "Composite bond instructions reference a bad topology bond "
          "index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(bond_coverage, RelationalOperator::EQUAL, bond_coverage_answer, "Bond accumulation is "
          "incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.nangl > 0) {
    check(angl_range_problem == false, "Composite angle instructions reference a bad topology "
          "index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(angl_coverage, RelationalOperator::EQUAL, angl_coverage_answer, "Angle accumulation is "
          "incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.ndihe > 0) {
    check(dihe_range_problem == false, "Composite dihedral instructions reference a bad topology "
          "dihedral index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(dihe_coverage, RelationalOperator::EQUAL, dihe_coverage_answer, "Dihedral accumulation "
          "is incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.nubrd > 0) {
    check(ubrd_range_problem == false, "Composite bond instructions reference a bad topology "
          "Urey-Bradley index in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
    check(ubrd_coverage, RelationalOperator::EQUAL, ubrd_coverage_answer, "Urey-Bradley "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ncimp > 0) {
    check(cimp_range_problem == false, "Composite dihedral instructions reference a bad CHARMM "
          "improper topology index in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
    check(cimp_coverage, RelationalOperator::EQUAL, cimp_coverage_answer, "CHARMM improper "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ncmap > 0) {
    check(cmap_range_problem == false, "CMAP instructions reference a bad topology index in "
          "valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(cmap_coverage, RelationalOperator::EQUAL, cmap_coverage_answer, "CMAP surface term "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ninfr14 > 0) {
    check(infr_range_problem == false, "Inferred 1:4 instructions reference a bad topology index "
          "in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(infr_coverage, RelationalOperator::EQUAL, infr_coverage_answer, "Inferred 1:4 "
          "non-bonded interaction accumulation is incorrect in valence work units for topology " +
          ag.getFileName() + ".", do_tests);
  }

  // Check each component of the energies the work units evaluate
  ScoreCard sc(2);
  PhaseSpace ps_vwu(ps);
  evaluateBondTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::BOND);
  const std::vector<double> bond_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> bond_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateAngleTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::ANGL);
  const std::vector<double> angl_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> angl_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateDihedralTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::DIHE);
  const std::vector<double> dihe_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> dihe_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateUreyBradleyTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::UBRD);
  const std::vector<double> ubrd_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> ubrd_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateCharmmImproperTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::CIMP);
  const std::vector<double> cimp_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> cimp_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateCmapTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::CMAP);
  const std::vector<double> cmap_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> cmap_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateAttenuated14Terms(&ag, &ps, &sc, EvaluateForce::YES, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::INFR14);
  const std::vector<double> attn_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> attn_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  
  // Restraint energies need to be plucked from the results immediately, as the reference
  // calculations only evaluate restraints all at once whereas the ValenceWorkUnit evaluation
  // implements them one type at a time.
  evaluateRestraints(&ra, &ps, &sc, EvaluateForce::YES, 0);
  const double rest_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RPOSN);
  const double rposn_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RBOND);
  const double rbond_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RANGL);
  const double rangl_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RDIHE);
  const double rdihe_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  const std::vector<double> rstr_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> rstr_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  
  // Read off other energy terms
  const std::vector<double> bond_e = sc.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_e = sc.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_e = sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_e =
    sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);
  const std::vector<double> ubrd_e = sc.reportInstantaneousStates(StateVariable::UREY_BRADLEY);
  const std::vector<double> cimp_e = sc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER);
  const std::vector<double> cmap_e = sc.reportInstantaneousStates(StateVariable::CMAP);
  const std::vector<double> qq14_e =
    sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC_ONE_FOUR);
  const std::vector<double> lj14_e =
    sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR);

  // Check that energy terms computed with ValenceWorkUnits match the reference calculations
  check(bond_e[0], RelationalOperator::EQUAL, Approx(bond_e[1]).margin(1.0e-6), "Bond energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(angl_e[0], RelationalOperator::EQUAL, Approx(angl_e[1]).margin(1.0e-6), "Angle energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(dihe_e[0], RelationalOperator::EQUAL, Approx(dihe_e[1]).margin(1.0e-6), "Cosine-based "
        "proper dihedral energies computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(impr_e[0], RelationalOperator::EQUAL, Approx(impr_e[1]).margin(1.0e-6), "Cosine-based "
        "improper dihedral energies computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);
  check(ubrd_e[0], RelationalOperator::EQUAL, Approx(ubrd_e[1]).margin(1.0e-6), "Urey-Bradley "
        "energies computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(cimp_e[0], RelationalOperator::EQUAL, Approx(cimp_e[1]).margin(1.0e-6), "CHARMM improper "
        "energies computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(cmap_e[0], RelationalOperator::EQUAL, Approx(cmap_e[1]).margin(1.0e-6), "CMAP energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(qq14_e[0], RelationalOperator::EQUAL, Approx(qq14_e[1]).margin(1.0e-6), "Electrostatic "
        "1:4 interactions computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(lj14_e[0], RelationalOperator::EQUAL, Approx(lj14_e[1]).margin(1.0e-6), "Electrostatic "
        "1:4 interactions computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(rest_e, RelationalOperator::EQUAL, Approx(rposn_e + rbond_e + rangl_e + rdihe_e, 1.0e-6),
        "Restraint penalty energies computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);

  // Check that forces computed with ValenceWorkUnits match the reference calculations
  check(bond_frc, RelationalOperator::EQUAL, Approx(bond_frc_ref).margin(1.0e-6), "Forces due to "
        "bonds computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(angl_frc, RelationalOperator::EQUAL, Approx(angl_frc_ref).margin(1.0e-6), "Forces due to "
        "angles computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(dihe_frc, RelationalOperator::EQUAL, Approx(dihe_frc_ref).margin(1.0e-6), "Forces due to "
        "proper and improper (cosine-based) dihedral terms computed with ValenceWorkUnits do not "
        "agree with the reference calculations (topology " + top_name + ").", do_tests);
  check(attn_frc, RelationalOperator::EQUAL, Approx(attn_frc_ref).margin(1.0e-6), "Forces due to "
        "attenuated 1:4 interactions computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);
  check(ubrd_frc, RelationalOperator::EQUAL, Approx(ubrd_frc_ref).margin(1.0e-6), "Forces due to "
        "Urey-Bradley terms computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(cimp_frc, RelationalOperator::EQUAL, Approx(cimp_frc_ref).margin(1.0e-6), "Forces due to "
        "CHARMM improper dihedrals computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(cmap_frc, RelationalOperator::EQUAL, Approx(cmap_frc_ref).margin(1.0e-6), "Forces due to "
        "CMAP potential surfaces computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(rstr_frc, RelationalOperator::EQUAL, Approx(rstr_frc_ref).margin(1.0e-6), "Forces due to "
        "restraint potentials computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Accumulate contributing atoms base on a series of force terms.
//
// Arguments:
//   item_count:   The number of energy terms of some type, the trusted length of {i,j,...}_atoms
//   result:       Array of all atoms that influence any one atom (accumulated and returned)
//   i_atoms:      The first list of atoms associated with whatever energy term
//   j_atoms:      The second list of atoms associated with whatever energy term
//   k_atoms:      The third list of atoms associated with whatever energy term
//   l_atoms:      The fourth list of atoms associated with whatever energy term
//   m_atoms:      The fifth list of atoms associated with whatever energy term
//-------------------------------------------------------------------------------------------------
void accumulateContributingAtoms(const int item_count, std::vector<std::vector<int>> *result,
                                 const int* i_atoms, const int* j_atoms = nullptr,
                                 const int* k_atoms = nullptr, const int* l_atoms = nullptr,
                                 const int* m_atoms = nullptr) {
  std::vector<int>* result_data = result->data();
  for (int pos = 0; pos < item_count; pos++) {
    std::vector<int> atom_list;
    atom_list.push_back(i_atoms[pos]);
    int order = 1;
    if (j_atoms != nullptr) {
      atom_list.push_back(j_atoms[pos]);
      order++;
    }
    if (k_atoms != nullptr) {
      atom_list.push_back(k_atoms[pos]);
      order++;
    }
    if (l_atoms != nullptr) {
      atom_list.push_back(l_atoms[pos]);
      order++;
    }
    if (m_atoms != nullptr) {
      atom_list.push_back(m_atoms[pos]);
      order++;
    }
    for (int i = 0; i < order - 1; i++) {
      for (int j = i + 1; j < order; j++) {
        result_data[atom_list[i]].push_back(atom_list[j]);
        result_data[atom_list[j]].push_back(atom_list[i]);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Determine a list of all atoms that contribute to forces on any particular atom through valence
// terms, restraints, or virtual sites.  This is done by looping over all topology terms and
// restraints while extending lists pertaining to each atom of the topology.
//
// Arguments:
//   ag:      The system topology
//   ra:      The system restraint apparatus, a sort of supplemental energy surface
//-------------------------------------------------------------------------------------------------
std::vector<std::vector<int>> getAtomForceContributors(const AtomGraph &ag,
                                                       const RestraintApparatus &ra) {
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const RestraintKit<double, double2, double4> rar = ra.getDoublePrecisionAbstract();

  // Initialize the result.  Every atom list includes the atom itself.
  std::vector<std::vector<int>> result(vk.natom, std::vector<int>());
  for (int i = 0; i < vk.natom; i++) {
    result[i].push_back(i);
  }
  accumulateContributingAtoms(vk.nbond, &result, vk.bond_i_atoms, vk.bond_j_atoms);
  accumulateContributingAtoms(vk.nangl, &result, vk.angl_i_atoms, vk.angl_j_atoms,
                              vk.angl_k_atoms);
  accumulateContributingAtoms(vk.ndihe, &result, vk.dihe_i_atoms, vk.dihe_j_atoms, vk.dihe_k_atoms,
                              vk.dihe_l_atoms);
  accumulateContributingAtoms(vk.nubrd, &result, vk.ubrd_i_atoms, vk.ubrd_k_atoms);
  accumulateContributingAtoms(vk.ncimp, &result, vk.cimp_i_atoms, vk.cimp_j_atoms, vk.cimp_k_atoms,
                              vk.cimp_l_atoms);
  accumulateContributingAtoms(vk.ncmap, &result, vk.cmap_i_atoms, vk.cmap_j_atoms, vk.cmap_k_atoms,
                              vk.cmap_l_atoms, vk.cmap_m_atoms);
  accumulateContributingAtoms(vk.ninfr14, &result, vk.infr14_i_atoms, vk.infr14_l_atoms);
  accumulateContributingAtoms(rar.nposn, &result, rar.rposn_atoms);
  accumulateContributingAtoms(rar.nbond, &result, rar.rbond_i_atoms, rar.rbond_j_atoms);
  accumulateContributingAtoms(rar.nangl, &result, rar.rangl_i_atoms, rar.rangl_j_atoms,
                              rar.rangl_k_atoms);
  accumulateContributingAtoms(rar.ndihe, &result, rar.rdihe_i_atoms, rar.rdihe_j_atoms,
                              rar.rdihe_k_atoms, rar.rdihe_l_atoms);

  // All atoms of each virtual site frame should be included with every atom of the frame
  for (int pos = 0; pos < vsk.nsite; pos++) {
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos]);
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos],
                                  &vsk.frame3_idx[pos]);
      break;
    case VirtualSiteKind::FIXED_4:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos],
                                  &vsk.frame3_idx[pos], &vsk.frame4_idx[pos]);
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Finally, all atoms that contribute to forces on any virtual site should be included with
  // every frame atom.
  for (int pos = 0; pos < vsk.nsite; pos++) {
    const int vatom = vsk.vs_atoms[pos];
    result[vsk.frame1_idx[pos]].insert(result[vsk.frame1_idx[pos]].end(),
                                       result[vatom].begin(), result[vatom].end());
    result[vsk.frame2_idx[pos]].insert(result[vsk.frame2_idx[pos]].end(),
                                       result[vatom].begin(), result[vatom].end());
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      result[vsk.frame3_idx[pos]].insert(result[vsk.frame3_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      break;
    case VirtualSiteKind::FIXED_4:
      result[vsk.frame3_idx[pos]].insert(result[vsk.frame3_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      result[vsk.frame4_idx[pos]].insert(result[vsk.frame4_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Reduce each vector to unique values and return the result
  for (int i = 0; i < vk.natom; i++) {
    reduceUniqueValues(&result[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);

  // Section 1
  section("SystemCache construction");

  // Section 2
  section("PhaseSpaceSynthesis layout");

  // Section 3
  section("Valence work unit construction");
  
  section(1);
  const char osc = osSeparator();
  const TestPriority test_sysc = (oe.getTemporaryDirectoryAccess()) ? TestPriority::CRITICAL :
                                                                      TestPriority::ABORT;
  const std::string mdin_name = oe.getTemporaryDirectoryPath() + osc + "cachetest.in";
  if (oe.getTemporaryDirectoryAccess()) {
    std::ofstream foutp = openOutputFile(mdin_name, PrintSituation::OVERWRITE, "Prepare a brief "
                                         "input file for testing the SystemCache machinery.");
    std::string buffer("&files\n  -p ");
    buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "topol" + osc +
              ".*.top\n  -c ";
    buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "coord" + osc +
              ".*.inpcrd\n&end\n";
    foutp.write(buffer.c_str(), buffer.size());
    foutp.close();
    oe.logFileCreated(mdin_name);
  }
  const bool mdin_exists = (getDrivePathType(mdin_name) == DrivePathType::FILE);
  const TextFile tf = (mdin_exists) ? TextFile(mdin_name) : TextFile();
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  SystemCache sysc(fcon, ExceptionResponse::SILENT);
  const int nsys = sysc.getSystemCount();
  check(nsys, RelationalOperator::EQUAL, 16, "The number of systems detected with regular "
        "expression searching did not meet expectations.  This may indicate a problem with the "
        "${OMNI_SOURCE} environment variable that did not show up in test program setup.",
        test_sysc);
  const std::vector<AtomGraph*> ag_ptr = sysc.getSystemTopologyPointer();
  const std::vector<PhaseSpace*> ps_ptr = sysc.getCoordinatePointer();
  std::vector<int> ag_atoms(nsys);
  std::vector<int> ps_atoms(nsys);
  for (int i = 0; i < nsys; i++) {
    ag_atoms[i] = ag_ptr[i]->getAtomCount();
    ps_atoms[i] = ps_ptr[i]->getAtomCount();
  }
  check(ag_atoms, RelationalOperator::EQUAL, ps_atoms, "The numbers of atoms found with pointers "
        "matched to topologies and coordinate sets of the same SystemCache object disagree.",
        test_sysc);
  std::vector<int> ag_ref_atoms(nsys);
  std::vector<int> ps_ref_atoms(nsys);
  for (int i = 0; i < nsys; i++) {
    const AtomGraph  &ag_ref = sysc.getSystemTopologyReference(i);
    const PhaseSpace &ps_ref = sysc.getCoordinateReference(i);
    ag_ref_atoms[i] = ag_ref.getAtomCount();
    ps_ref_atoms[i] = ps_ref.getAtomCount();
  }
  check(ag_ref_atoms, RelationalOperator::EQUAL, ps_ref_atoms, "The numbers of atoms found with "
        "references to matching topologies and coordinate sets of the same SystemCache object "
        "disagree.", test_sysc);
  check(ag_ref_atoms, RelationalOperator::EQUAL, ps_atoms, "The numbers of atoms found with "
        "pointers and references matched to topologies and coordinate sets of the same "
        "SystemCache object disagree.", test_sysc);
  
  // Create some topologies and coordinate sets.
  Xoroshiro128pGenerator my_prng(oe.getRandomSeed());
  const std::string base_crd_name = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_crd_name = base_crd_name + osc + "tip3p.inpcrd";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_crd_name = base_crd_name + osc + "tip4p.inpcrd";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const bool files_exist = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps, tip3p_ps_2, tip4p_ps_2, trpcage_ps_2;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip3p_ps.buildFromFile(tip3p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    tip4p_ps.buildFromFile(tip4p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    trpcage_ps.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtWarn("The topology and coordinate files for the TIP3P and TIP4P water boxes as well as the "
           "Trp-cage miniprotein in water must be available in ${OMNI_SOURCE}/test/ "
           "subdirectories Topology and Trajectory, respectively.  Check the $OMNI_SOURCE "
           "environment variable to make sure that it is set properly.  A number of tests will "
           "be skipped.", "test_phase_space_synthesis");
  }
  
  // Make a PhaseSpaceSynthesis the meticulous way
  section(2);
  std::vector<PhaseSpace>  psv = { tip3p_ps, tip4p_ps, trpcage_ps, tip3p_ps, tip4p_ps,
                                   tip3p_ps, tip3p_ps, tip4p_ps, trpcage_ps };
  const std::vector<AtomGraph*> agv = { &tip3p_ag, &tip4p_ag, &trpcage_ag, &tip3p_ag, &tip4p_ag,
                                        &tip3p_ag, &tip3p_ag, &tip4p_ag, &trpcage_ag };
  for (size_t i = 3; i < psv.size(); i++) {
    PhaseSpaceWriter pswi = psv[i].data();
    for (int j = 0; j < pswi.natom; j++) {
      pswi.xcrd[j] += 0.02 * my_prng.gaussianRandomNumber();
      pswi.ycrd[j] += 0.02 * my_prng.gaussianRandomNumber();
      pswi.zcrd[j] += 0.02 * my_prng.gaussianRandomNumber();
    }
  }
  PhaseSpaceSynthesis psynth(psv, agv);   

  // Try extracting a system from it
  PhaseSpace tip3p_ps_copy(tip3p_ps.getAtomCount(), tip3p_ps.getUnitCellType());
  psynth.extractPhaseSpace(&tip3p_ps_copy, 3);
  PhaseSpaceWriter tip3p_orig_writer = tip3p_ps.data();
  PhaseSpaceWriter tip3p_muta_writer = psv[3].data();
  PhaseSpaceWriter tip3p_copy_writer = tip3p_ps_copy.data();
  std::vector<double> y_orig(tip3p_orig_writer.natom);
  std::vector<double> y_muta(tip3p_orig_writer.natom);
  std::vector<double> y_copy(tip3p_orig_writer.natom);
  for (int i = 0; i < tip3p_orig_writer.natom; i++) {
    y_orig[i] = tip3p_orig_writer.ycrd[i];
    y_muta[i] = tip3p_muta_writer.ycrd[i];
    y_copy[i] = tip3p_copy_writer.ycrd[i];
  }
  check(y_muta, RelationalOperator::EQUAL, y_copy, "The PhaseSpaceSynthesis object was not able "
        "to return the correct image of one of its systems.", do_tests);
  const std::vector<double> scaling_answer = {
    default_globalpos_scale_lf, default_localpos_scale_lf, default_velocity_scale_lf,
    default_force_scale_lf, default_inverse_globalpos_scale_lf, default_inverse_localpos_scale_lf,
    default_inverse_velocity_scale_lf, default_inverse_force_scale_lf };
  const std::vector<int> scale_bits_answer = {
    default_globalpos_scale_bits, default_localpos_scale_bits, default_velocity_scale_bits,
    default_force_scale_bits };
  PsSynthesisWriter psynth_w = psynth.data();
  const std::vector<double> scaling_result = { psynth_w.gpos_scale, psynth_w.lpos_scale,
                                               psynth_w.vel_scale, psynth_w.frc_scale,
                                               psynth_w.inv_gpos_scale, psynth_w.inv_lpos_scale,
                                               psynth_w.inv_vel_scale, psynth_w.inv_frc_scale };
  const std::vector<float> scaling_result_f = { psynth_w.gpos_scale_f, psynth_w.lpos_scale_f,
                                                psynth_w.vel_scale_f, psynth_w.frc_scale_f,
                                                psynth_w.inv_gpos_scale_f,
                                                psynth_w.inv_lpos_scale_f,
                                                psynth_w.inv_vel_scale_f,
                                                psynth_w.inv_frc_scale_f };
  const std::vector<int> scale_bits_result = { psynth_w.gpos_bits, psynth_w.lpos_bits,
                                               psynth_w.vel_bits, psynth_w.frc_bits };
  check(scaling_result, RelationalOperator::EQUAL,
        Approx(scaling_answer, ComparisonType::RELATIVE, omni::constants::verytiny),
        "Double-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.");
  check(scaling_result_f, RelationalOperator::EQUAL,
        Approx(scaling_answer, ComparisonType::RELATIVE, omni::constants::verytiny),
        "Single-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.");
  check(scale_bits_result, RelationalOperator::EQUAL, scale_bits_answer, "Fixed-precision bit "
        "counts found in the PhaseSpaceSynthesis object's writer do not meet expectations.");
  check(psynth_w.time_step, RelationalOperator::EQUAL, 1.0, "The time step was not correctly "
        "copied into a PhaseSpaceSynthesis writeable abstract.");
  PhaseSpaceSynthesis psynth2(psv, agv, std::vector<Thermostat>(1), std::vector<Barostat>(1),
                              2.5, 24, 25, 40, 28);
  PsSynthesisWriter psynth_w2 = psynth2.data();
  const std::vector<double> scaling_answer2 = { pow(2.0, 24), pow(2.0, 25), pow(2.0, 40),
                                                pow(2.0, 28), 1.0 / pow(2.0, 24),
                                                1.0 / pow(2.0, 25), 1.0 / pow(2.0, 40),
                                                1.0 / pow(2.0, 28) };
  const std::vector<float> scaling_answer2_f(scaling_answer2.begin(), scaling_answer2.end());
  const std::vector<int> scale_bits_answer2 = { 24, 25, 40, 28 };
  const std::vector<double> scaling_result2 = { psynth_w2.gpos_scale, psynth_w2.lpos_scale,
                                                psynth_w2.vel_scale, psynth_w2.frc_scale,
                                                psynth_w2.inv_gpos_scale, psynth_w2.inv_lpos_scale,
                                                psynth_w2.inv_vel_scale, psynth_w2.inv_frc_scale };
  const std::vector<float> scaling_result2_f = { psynth_w2.gpos_scale_f, psynth_w2.lpos_scale_f,
                                                 psynth_w2.vel_scale_f, psynth_w2.frc_scale_f,
                                                 psynth_w2.inv_gpos_scale_f,
                                                 psynth_w2.inv_lpos_scale_f,
                                                 psynth_w2.inv_vel_scale_f,
                                                 psynth_w2.inv_frc_scale_f };
  const std::vector<int> scale_bits_result2 = { psynth_w2.gpos_bits, psynth_w2.lpos_bits,
                                                psynth_w2.vel_bits, psynth_w2.frc_bits };
  check(scaling_result2, RelationalOperator::EQUAL,
        Approx(scaling_answer2, ComparisonType::RELATIVE, omni::constants::verytiny),
        "Double-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.");
  check(scaling_result2_f, RelationalOperator::EQUAL,
        Approx(scaling_answer2, ComparisonType::RELATIVE, omni::constants::verytiny),
        "Single-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.");
  check(scale_bits_result2, RelationalOperator::EQUAL, scale_bits_answer2, "Fixed-precision bit "
        "counts found in the PhaseSpaceSynthesis object's writer do not meet expectations.");
  check(psynth_w2.time_step, RelationalOperator::EQUAL, 2.5, "The time step was not correctly "
        "copied into a PhaseSpaceSynthesis writeable abstract.");
  psynth2.extractPhaseSpace(&tip3p_ps_copy, 3);
  for (int i = 0; i < tip3p_orig_writer.natom; i++) {
    y_orig[i] = tip3p_orig_writer.ycrd[i];
    y_muta[i] = tip3p_muta_writer.ycrd[i];
    y_copy[i] = tip3p_copy_writer.ycrd[i];
  }
  check(y_muta, RelationalOperator::EQUAL, Approx(y_copy).margin(1.0e-6),
        "The PhaseSpaceSynthesis object returns an incorrect image of one of its systems, even "
        "after compensating for a lower global position resolution", do_tests);
  check(y_muta, RelationalOperator::NOT_EQUAL, Approx(y_copy).margin(1.0e-8),
        "The PhaseSpaceSynthesis object returns an image of one of its systems with higher "
        "fidelity to the original than expected.", do_tests);

  // Prepare valence work units for the array of topologies
  section(3);
  std::vector<RestraintApparatus> ra_vec;
  for (int i = 0; i < sysc.getTopologyCount(); i++) {
    const int example_system_idx = sysc.getCoordinateExample(i);
    const AtomGraph *ag_i = sysc.getTopologyPointer(i);
    const PhaseSpace &ps_i = sysc.getCoordinateReference(example_system_idx);
    const CoordinateFrameReader cfr_i(ps_i);
    const ChemicalFeatures chemfe_i(ag_i, ps_i, MapRotatableGroups::YES);
    const AtomMask bkbn_i(":* & @CA,N,C,O", ag_i, &chemfe_i, cfr_i);
    ra_vec.emplace_back(applyHydrogenBondPreventors(ag_i, chemfe_i, 64.0, 3.1));
    ra_vec[i].addRestraints(applyPositionalRestraints(ag_i, cfr_i, bkbn_i, 16.0));
  }
  bool force_partner_counts_match = true;
  bool force_partners_match = true;
  for (int i = 0; i < sysc.getTopologyCount(); i++) {
    ValenceDelegator vdel(sysc.getTopologyPointer(i), &ra_vec[i]);
    const std::vector<ValenceWorkUnit> vwu_i = buildValenceWorkUnits(&vdel);
    const AtomGraph &ag_i = sysc.getTopologyReference(i);
    const std::vector<std::vector<int>> fcontrib = getAtomForceContributors(ag_i, ra_vec[i]);
    const ValenceKit<double> vk = ag_i.getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> vsk = ag_i.getDoublePrecisionVirtualSiteKit();
    for (int j = 0; j < vk.natom; j++) {
      const std::vector<int> vwu_deps = vdel.findForcePartners(j);
      force_partner_counts_match = (force_partner_counts_match &&
                                    vwu_deps.size() == fcontrib[j].size());
      const std::vector<ValueWithCounter<int>> missing_dependencies =
        findUnmatchedValues(vwu_deps, fcontrib[j], UniqueValueHandling::CONFIRM_ALL_COPIES);
      force_partners_match = (missing_dependencies.size() == 0LLU);
    }
  }
  check(force_partner_counts_match, "The counts of force-relevant partners in ValenceWorkUnit "
        "objects made for a series of amino acid dipeptides do not agree with numbers computed "
        "through an alternative method.", do_tests);
  check(force_partners_match, "Lists of force-relevant partners in ValenceWorkUnit objects made "
        "for a series of amino acid dipeptides do not agree with those assembled through an "
        "alternative method.", do_tests);

  // Run diagnotics of the valence work units on a simple system, one with only a single work unit
  runValenceWorkUnitTests(sysc.getSystemTopologyReference(0).getFileName(),
                          sysc.getCoordinateReference(0).getFileName(), oe, &my_prng);

  // Read a larger topology that will be forced to split its contents among several work units
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  runValenceWorkUnitTests(dhfr_top_name, dhfr_crd_name, oe, &my_prng);

  // Read a solvated topology that will be forced to split its contents among work units, and
  // include water molecules with SETTLE constraint groups.
  const std::string ubiq_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";
  const std::string ubiq_top_name = base_top_name + osc + "ubiquitin.top";
  runValenceWorkUnitTests(ubiq_top_name, ubiq_crd_name, oe, &my_prng);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
