#include <string>
#include <vector>
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::ExceptionResponse;
using omni::constants::verytiny;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::openOutputFile;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::namelist::FilesControls;
using omni::parse::TextFile;
using omni::random::Xoroshiro128pGenerator;
using omni::synthesis::SystemCache;
using namespace omni::numerics;
using namespace omni::synthesis;
using namespace omni::trajectory;
using namespace omni::testing;

int main(int argc, char* argv[]) {

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
  for (int i = 0; i < sysc.getTopologyCount(); i++) {
    
  }
       
  // Summary evaluation
  printTestSummary(oe.getVerbosity());

  return 0;
}
