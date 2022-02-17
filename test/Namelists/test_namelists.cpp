#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/input.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Namelists/nml_random.h"
#include "../../src/Namelists/nml_solvent.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::diskutil::osSeparator;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtWarn;
using omni::parse::separateText;
using omni::topology::AtomGraph;
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
  ags.reserve(filcon.getFreeTopologyCount());
  std::vector<int> atom_counts(nfreetop);
  for (int i = 0; i < filcon.getFreeTopologyCount(); i++) {
    ags.push_back(AtomGraph(filcon.getFreeTopologyName(i), ExceptionResponse::SILENT));
    atom_counts[i] = ags[i].getAtomCount();
  }
  check(sum<int>(atom_counts), RelationalOperator::EQUAL, 213, "The free topologies named in "
        "the &files namelist do not contain the expected, combined number of atoms.");

  // The minimization namelist contains integer and real-valued numbers
  section(2);
  start_line = 0;
  MinimizeControls mincon(tf, &start_line);
  start_line = 0;
  RandomControls rngcon(tf, &start_line);
  start_line = 0;
  SolventControls watcon(tf, &start_line);

  
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
