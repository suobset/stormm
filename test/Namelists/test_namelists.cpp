#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/input.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Namelists/nml_random.h"
#include "../../src/Namelists/nml_solvent.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::constants::tiny;
using omni::diskutil::osSeparator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtWarn;
using omni::parse::separateText;
using omni::topology::AtomGraph;
using omni::trajectory::CoordinateFrame;
using namespace omni::namelist;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

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
  
  // The files namelist is perhaps the most complex due to its interchangeable defaults, and
  // will be critical to the operation of any OMNI app
  section(1);
  const char osc = osSeparator();
  const std::string input_base = oe.getOmniSourcePath() + osc + "test" + osc + "Namelists";
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

  // The random number control namelist provides a number of useful options
  start_line = 0;
  RandomControls rngcon(tf, &start_line);
  check(rngcon.getRandomSeed(), RelationalOperator::EQUAL, 67108863, "The &random namelist did "
        "not convey the correct random seed.");
  check(rngcon.getStreamCount(), RelationalOperator::EQUAL, 128, "The &random namelist did "
        "not convey the correct stream count.");
  check(rngcon.getProductionStride(), RelationalOperator::EQUAL, 32, "The &random namelist did "
        "not convey the correct production stride.");
  check(rngcon.getWarmupCycleCount(), RelationalOperator::EQUAL, 63, "The &random namelist did "
        "not convey the correct warmup cycle count.");
  start_line = 0;
  SolventControls watcon(tf, &start_line);

  
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
