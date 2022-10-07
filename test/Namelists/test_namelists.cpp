#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/ForceField/forcefield_enumerators.h"
#include "../../src/Namelists/input.h"
#include "../../src/Namelists/nml_ffmorph.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Namelists/nml_random.h"
#include "../../src/Namelists/nml_solvent.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::chemistry::ChemicalFeatures;
using stormm::constants::ExceptionResponse;
using stormm::constants::tiny;
using stormm::diskutil::osSeparator;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::errors::rtWarn;
using stormm::modeling::ForceFieldElement;
using stormm::modeling::ParameterKind;
using stormm::parse::separateText;
using stormm::parse::strcmpCased;
using stormm::parse::TextOrigin;
using stormm::topology::AtomGraph;
using stormm::topology::AtomicRadiusSet;
using stormm::topology::ImplicitSolventModel;
using stormm::trajectory::CoordinateFrame;
using namespace stormm::namelist;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Test what should be bad namelist input and make sure that it throws.
//
// Arguments:
//   nml_name:       Name of the namelist to test (i.e. &files)
//   content:        Content of the erroneous namelist.  This function add the trailing &end card.
//   error_message:  Error message to display if the content does NOT throw an exception.  This
//                   function adds "in the (namelist name) namelist" to the end.
//-------------------------------------------------------------------------------------------------
void testBadNamelist(const std::string &nml_name, const std::string &content,
                     const std::string &error_message) {
  const std::string bad_idea = std::string("&") + nml_name + " " + content + " &end";
  const TextFile bad_input(bad_idea, TextOrigin::RAM);
  int start_line = 0;
  bool found_nml;
  
  // Redirect stdout and record the number of test failures thus far.  This is to suppress Alert
  // messages that occur leading up to an ultimate namelist failure, which is supposed to happen
  // and thus should not call for any attention.
  FILE fp_old = *stdout;
  *stdout = *fopen("/dev/null", "w");
  const int initial_failures = gbl_test_results.getFailureCount();
  const std::string updated_error = error_message + " in the &" + nml_name + " namelist.";
  
  // Run the test, which may send messages to stdout
  if (strcmpCased(nml_name, "minimize")) {
    CHECK_THROWS(MinimizeControls t_mincon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "random")) {
    CHECK_THROWS(RandomControls t_rngcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "solvent")) {
    CHECK_THROWS(SolventControls t_watcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "ffmorph")) {
    CHECK_THROWS(FFMorphControls t_ffmcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else {
    rtErr("The namelist &" + nml_name + " does not pair with any known case.", "test_namelists");
  }

  // Reset stdout and reprint the error message if there was a new failure.
  *stdout = fp_old;
  if (gbl_test_results.getFailureCount() > initial_failures) {
    const std::string parsed_msg = terminalFormat(error_message + " in the &" + nml_name +
                                                  " namelist.", "", "", 14, 0, 14);
    printf("Check FAILED: %s\n", parsed_msg.c_str());
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test the &files namelist");

  // Section 2
  section("Test the &minimize namelist");

  // Section 3
  section("Test the &random namelist");

  // Section 4
  section("Test the &solvent namelist");

  // Section 5
  section("Test the &ffmorph namelist");

  // Section 6
  section("Test the &restraint namelist");
  
  // The files namelist is perhaps the most complex due to its interchangeable defaults, and
  // will be critical to the operation of any STORMM app
  section(1);
  const char osc = osSeparator();
  const std::string input_base = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists";
  const std::string main_file = input_base + osc + "testrun.in";
  const bool input_exists = (getDrivePathType(main_file) == DrivePathType::FILE);
  const TestPriority do_tests = (input_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const TextFile tf = (input_exists) ? TextFile(main_file) : TextFile();
  int start_line = 0;
  FilesControls filcon(tf, &start_line);
  const int nfreetop = filcon.getFreeTopologyCount();
  const int nfreecrd = filcon.getFreeCoordinatesCount();
  check(nfreetop, RelationalOperator::EQUAL, 6, "The &files namelist did not record the expected "
        "number of free topologies.", do_tests);
  check(nfreecrd, RelationalOperator::EQUAL, 6, "The &files namelist did not record the expected "
        "number of free input coordinate sets.", do_tests);
  check(filcon.getSystemDefinitionCount(), RelationalOperator::EQUAL, 2, "The &files namelist did "
        "not record the expected number of defined systems.", do_tests);
  std::vector<AtomGraph> ags;
  ags.reserve(nfreetop);
  std::vector<int> top_atom_counts(nfreetop);
  for (int i = 0; i < nfreetop; i++) {
    ags.push_back(AtomGraph(filcon.getFreeTopologyName(i), ExceptionResponse::SILENT));
    top_atom_counts[i] = ags[i].getAtomCount();
  }
  check(sum<int>(top_atom_counts), RelationalOperator::EQUAL, 213, "The free topologies named in "
        "the &files namelist do not contain the expected, combined number of atoms.", do_tests);
  std::vector<CoordinateFrame> cfs;
  ags.reserve(nfreecrd);
  std::vector<int> crd_atom_counts(nfreecrd);
  for (int i = 0; i < nfreecrd; i++) {
    cfs.push_back(CoordinateFrame(filcon.getFreeCoordinateName(i), CoordinateFileKind::UNKNOWN,
                                  0));
    crd_atom_counts[i] = cfs[i].getAtomCount();
  }
  check(sum<int>(crd_atom_counts), RelationalOperator::EQUAL, 213, "The free topologies named in "
        "the &files namelist do not contain the expected, combined number of atoms.", do_tests);

  // The minimization namelist contains integer and real-valued numbers
  std::string bad;
  section(2);
  start_line = 0;
  bool found_nml;
  MinimizeControls mincon(tf, &start_line, &found_nml);
  check(mincon.getTotalCycles(), RelationalOperator::EQUAL, 1000, "The &minimize namelist did "
        "not convey the correct total number of cycles.", do_tests);
  check(mincon.getSteepestDescentCycles(), RelationalOperator::EQUAL, 100, "The &minimize "
        "namelist did not convey the correct number of steepest descent cycles.", do_tests);
  check(mincon.getElectrostaticCutoff(), RelationalOperator::EQUAL, Approx(9.0).margin(tiny),
        "The &minimize namelist did not convey the correct electrostatic cutoff.", do_tests);
  check(mincon.getLennardJonesCutoff(), RelationalOperator::EQUAL, Approx(9.0).margin(tiny),
        "The &minimize namelist did not convey the correct van der-Waals cutoff.", do_tests);
  check(mincon.getInitialStep(), RelationalOperator::EQUAL, Approx(0.02).margin(tiny),
        "The &minimize namelist did not convey the correct initial step size.", do_tests);
  check(mincon.getConvergenceTarget(), RelationalOperator::EQUAL, Approx(0.00008).margin(tiny),
        "The &minimize namelist did not convey the correct convergence criterion.", do_tests);
  testBadNamelist("minimize", "maxcyc = -1", "A negative cycle count was accepted");
  testBadNamelist("minimize", "ncyc = 50, maxcyc = 49", "The number of steepest descent cycles "
                  "was allowed to exceed the number of all cycles");
  testBadNamelist("minimize", "ncyc = 5, maxcyc = 9, dx0 = 0.0", "An initial step size of zero "
                  "was accepted");
  testBadNamelist("minimize", "drms = 0.00000000001", "An exceedingly small convergence step "
                  "was permitted");
  
  // The random number control namelist provides a number of useful options
  section(3);
  start_line = 0;
  RandomControls rngcon(tf, &start_line, &found_nml);
  check(rngcon.getRandomSeed(), RelationalOperator::EQUAL, 67108863, "The &random namelist did "
        "not convey the correct random seed.", do_tests);
  check(rngcon.getStreamCount(), RelationalOperator::EQUAL, 128, "The &random namelist did "
        "not convey the correct stream count.", do_tests);
  check(rngcon.getProductionStride(), RelationalOperator::EQUAL, 32, "The &random namelist did "
        "not convey the correct production stride.", do_tests);
  check(rngcon.getWarmupCycleCount(), RelationalOperator::EQUAL, 63, "The &random namelist did "
        "not convey the correct warmup cycle count.", do_tests);
  testBadNamelist("random", "igseed=0", "A random seed of zero was accepted");
  testBadNamelist("random", "igstreams = -5", "A negative number of random streams was accepted");
  testBadNamelist("random", "igstreams = 0", "A count of zero random streams was accepted");
  testBadNamelist("random", "igstride = -2", "A negative random production stride was accepted");

  // The solvent namelist provides a means of specifying the solvent model
  section(4);
  start_line = 0;
  SolventControls watcon(tf, &start_line, &found_nml);
  check(watcon.getBornRadiiCutoff(), RelationalOperator::EQUAL, Approx(25.0).margin(tiny),
        "The &solvent namelist reports the wrong Born radius cutoff (rgbmax).", do_tests);
  check(watcon.getInternalDielectric(), RelationalOperator::EQUAL, Approx(1.5).margin(tiny),
        "The &solvent namelist reports the wrong internal dieletric (intdiel).", do_tests);
  check(watcon.getExternalDielectric(), RelationalOperator::EQUAL, Approx(77.9).margin(tiny),
        "The &solvent namelist reports the wrong external dieletric (extdiel).", do_tests);
  check(watcon.getPBRadiiSet() == AtomicRadiusSet::MBONDI2, "The &solvent namelist reports the "
        "wrong PB radii set (pbradii)", do_tests);
  check(watcon.getImplicitSolventModel() == ImplicitSolventModel::OBC_GB_II, "The &solvent "
        "namelist reports the wrong implicit solvent model (igb)", do_tests);
  testBadNamelist("solvent", "pbradii = MBONDI4", "An invalid radius set was accepted");
  testBadNamelist("solvent", "pbradii = MBONDI2, pbradii = MBONDI", "Duplicate specification of "
                  "the radius set keyword was allowed");
  testBadNamelist("solvent", "rgbmax = -16.0", "A negative Born radius cutoff was allowed");
  testBadNamelist("solvent", "intdiel = 0.001", "A ridiculously small internal dielectric was "
                  "allowed");
  testBadNamelist("solvent", "extdiel = 0.005", "A ridiculously small external dielectric was "
                  "allowed");
  
  // The force field morphing namelist allows parameter modifications at run time
  section(5);
  start_line = 0;
  FFMorphControls ffmcon(tf, &start_line, &found_nml);
  check(ffmcon.getEditCount(ParameterKind::BOND), RelationalOperator::EQUAL, 3, "The &ffmorph "
        "namelist does not convey the expected number of bond parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::ANGLE), RelationalOperator::EQUAL, 2, "The &ffmorph "
        "namelist does not convey the expected number of angle parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::DIHEDRAL), RelationalOperator::EQUAL, 3, "The &ffmorph "
        "namelist does not convey the expected number of dihedral parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::UREY_BRADLEY), RelationalOperator::EQUAL, 1,
        "The &ffmorph namelist does not convey the expected number of Urey-bradley interaction "
        "edits.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::BOND, 0).testStiffnessModification(), "The stiffness "
        "parameter of a bond parameter is not properly marked for modification.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::BOND, 1).testEquilibriumModification() == false,
        "The equilibrium parameter of a bond parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::ANGLE, 1).testStiffnessModification() == false,
        "The stiffness constant of an angle parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 1).testAmplitudeModification(),
        "The amplitude of a dihedral parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 1).getAmplitude(), RelationalOperator::EQUAL,
        Approx(2.07).margin(tiny), "The amplitude of a dihedral parameter is not properly marked "
        "for modification.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 0).getPhaseAngle(), RelationalOperator::EQUAL,
        Approx(0.24 * stormm::symbols::pi / 180.0).margin(tiny), "The amplitude of a dihedral "
        "parameter is not properly marked for modification.", do_tests);
  testBadNamelist("ffmorph", "bond { -ti CX -tj CT }", "An incomplete bond keyword entry was "
                  "accepted");
  testBadNamelist("ffmorph", "bond { -ti CX -tj CT -blah }", "A bond keyword entry with a "
                  "spurious subkey was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tj CT -k 57.0 }", "An angle keyword entry with "
                  "too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tk CT -theta0 107.8 }", "An angle keyword entry "
                  "with too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -tj CT -tk CT -theta0 107.8 }", "An angle keyword "
                  "entry with too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tj CT -tk CT }", "An angle keyword "
                  "entry with no modifications was accepted");
  testBadNamelist("ffmorph", "dihedral { -ti CG -tk CV -tl CN -n 2 -amp 0.97 }", "A dihedral "
                  "keyword entry with too few atom types was accepted");
  testBadNamelist("ffmorph", "dihedral { -ti CG -tk CV -tj OS -tl CN -tl H1 -n 2 -amp 0.97 }",
                  "A dihedral keyword entry with repeated atom types was accepted");

  // The restraint namelist can build individual restraints as well as ensembles of them
  section(6);
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  TestSystemManager tsm(base_top_name, "top",
                        { "stereo_L1_vs", "symmetry_L1_vs", "drug_example_vs_iso",
                          "bromobenzene_vs_iso", "med_1", "med_2", "med_3", "med_4", "med_5" },
                        base_crd_name, "inpcrd",
                        { "stereo_L1_vs", "symmetry_L1_vs", "drug_example_vs_iso",
                          "bromobenzene_vs_iso", "med_1", "med_2", "med_3", "med_4", "med_5" });
  
  start_line = 0;
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
