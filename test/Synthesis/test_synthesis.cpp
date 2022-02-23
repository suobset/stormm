#include <string>
#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::openOutputFile;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::namelist::FilesControls;
using omni::parse::TextFile;
using omni::random::Ran2Generator;
using omni::synthesis::SystemCache;
using namespace omni::card;
using namespace omni::trajectory;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test SystemCache construction");

  // Section 2
  section("Test PhaseSpaceSynthesis layout");

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
  const std::vector<AtomGraph*> ag_ptr = sysc.getTopologyPointer();
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
    const AtomGraph  &ag_ref = sysc.getTopologyReference(i);
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
  
  // Create some vectors of random numbers, then upload them and test what happens when perturbing
  // atomic coordinates by these numbers.
  Ran2Generator my_prng(oe.getRandomSeed());
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
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps;
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

  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
