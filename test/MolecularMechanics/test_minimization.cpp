#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/behavior.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Restraints/bounded_restraint.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace omni::chemistry;
using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::errors;
using namespace omni::mm;
using namespace omni::namelist;
using namespace omni::restraints;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  StopWatch timer;

  // Section 1
  section("Minimize a dipeptide");

  // Read files
  const char osc = osSeparator();
  const std::string base_top_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene.top";
  const std::string brbz_crd_name = base_crd_name + osc + "bromobenzene.inpcrd";
  const std::string bbvs_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const std::string bbvs_crd_name = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  const std::vector<std::string> all_top = { alad_top_name, brbz_top_name, bbvs_top_name };
  const std::vector<std::string> all_crd = { alad_crd_name, brbz_crd_name, bbvs_crd_name };
  const int system_count = all_top.size();
  bool files_exist = true;
  for (int i = 0; i < system_count; i++) {
    files_exist = (getDrivePathType(all_top[i]) == DrivePathType::FILE &&
                   getDrivePathType(all_crd[i]) == DrivePathType::FILE && files_exist);
  }
  if (files_exist == false) {
    rtWarn("Topology and coordinate files must be available for tests to proceed.  Check the "
           "${OMNI_SOURCE} environment variable, currently set to " + oe.getOmniSourcePath() +
           ", to verify that it is valid.  Subsequent tests will be skipped.",
           "test_minimization");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<AtomGraph> all_ag;
  std::vector<PhaseSpace> all_ps;
  std::vector<ChemicalFeatures> all_chemfe;
  std::vector<RestraintApparatus> all_ra;
  std::vector<StaticExclusionMask> all_se;
  all_ag.reserve(system_count);
  all_ps.reserve(system_count);
  all_chemfe.reserve(system_count);
  all_ra.reserve(system_count);
  all_se.reserve(system_count);
  if (files_exist) {
    for (int i = 0; i < system_count; i++) {
      all_ag.emplace_back(all_top[i], ExceptionResponse::SILENT);
      all_ps.emplace_back(all_crd[i]);
      all_chemfe.emplace_back(&all_ag[i], all_ps[i], MapRotatableGroups::YES, 300.0, &timer);
      const std::vector<BoundedRestraint> brs = applyHydrogenBondPreventors(&all_ag[i],
                                                                            all_chemfe[i], 16.0); 
      all_ra.emplace_back(brs, &all_ag[i]);
      all_se.emplace_back(&all_ag[i]);
    }
  }
  else {
    all_ag.resize(system_count);
    all_ps.resize(system_count);
    all_chemfe.resize(system_count);
    all_ra.resize(system_count);
    all_se.resize(system_count);
  }
  
  // Try the dipeptide--this systems contains CMAPs in addition to basic Amber force field terms
  section(1);
  MinimizeControls mincon;
  minimize(&all_ps[0], all_ag[0], all_ra[0], all_se[0], mincon);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
