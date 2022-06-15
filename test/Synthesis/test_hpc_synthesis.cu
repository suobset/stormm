// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/rounding.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_valence_potential.cuh"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"

using namespace omni::card;
using namespace omni::constants;
using namespace omni::errors;
using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::math;
using namespace omni::mm;
using namespace omni::parse;
using namespace omni::synthesis;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  
  // Section 1
  section("Coordinate compilation and staging");
  
  // Section 2
  section("Topology compilation and staging");

  // Get the GPU specs.  Set of parameters for the work units and launch grids.
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  int nblocks = gpu.getSMPCount();
  int nthreads = gpu.getMaxThreadsPerBlock();
  int max_vwu_atoms = maximum_valence_work_unit_atoms;
  if (gpu.getArchMajor() == 6 && gpu.getArchMinor() == 1) {
    nblocks *= 2;
    nthreads /= 2;
    max_vwu_atoms /= 2;
  }

  // Configure the relevant kernels for this executable.
  valenceKernelSetup();
  
  // Collect coordinates and topologies
  const char osc = osSeparator();
  std::string buffer("&files\n  -p ");
  buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "topol" + osc +
            ".*.top\n  -c ";
  buffer += oe.getOmniSourcePath() + osc + "test" + osc + "Namelists" + osc + "coord" + osc +
            ".*.inpcrd\n&end\n";
  const TextFile tf(buffer, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  const SystemCache sysc(fcon, ExceptionResponse::SILENT, MapRotatableGroups::NO, &timer);

  // Form the syntheses of topologies and coordinates
  section(1);
  const int nsys = sysc.getSystemCount();
  const TestPriority do_tests = (nsys > 0) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (nsys == 0) {
    rtWarn("No system topology and coordinate pairs were detected.  Subsequent tests will be "
           "skipped.", "test_hpc_synthesis");
  }
  std::vector<int> topology_indices(nsys, 0);
  for (int i = 0; i < nsys; i++) {
    topology_indices[i] = i;
  }
  AtomGraphSynthesis poly_ag(sysc.getSystemTopologyPointerCC(), topology_indices,
                             ExceptionResponse::WARN, max_vwu_atoms, &timer);
  PhaseSpaceSynthesis poly_ps(sysc);
  check(poly_ag.getSystemCount(), RelationalOperator::EQUAL, poly_ps.getSystemCount(),
        "PhaseSpaceSynthesis and AtomGraphSynthesis objects formed from the same SystemCache have "
        "different numbers of systems inside of them.", do_tests);

  // Upload the compiled systems and check the results
  PsSynthesisWriter psw = poly_ps.data();
  poly_ag.upload();
  poly_ps.upload();
  std::vector<double> gpu_charges = poly_ag.getPartialCharges<double>(HybridTargetLevel::DEVICE);
  int padded_atom_count = 0;
  for (int i = 0; i < nsys; i++) {
    if (i == nsys - 1) {
      padded_atom_count += poly_ag.getSystemTopologyPointer(i)->getAtomCount();
    }
    else {
      padded_atom_count += roundUp(poly_ag.getSystemTopologyPointer(i)->getAtomCount(),
                                   warp_size_int);
    }
  }
  std::vector<double> rbt_charges(padded_atom_count, 0.0);
  int atom_offset = 0;
  for (int i = 0; i < nsys; i++) {
    const AtomGraph *iag_ptr = poly_ag.getSystemTopologyPointer(i);
    std::vector<double> ichg = iag_ptr->getPartialCharge<double>();
    const int natom = iag_ptr->getAtomCount();
    for (int j = 0; j < natom; j++) {
      rbt_charges[atom_offset + j] = ichg[j];
    }
    const int padded_natom = roundUp(natom, warp_size_int);
    if (i < nsys - 1) {
      for (int j = natom; j < padded_natom; j++) {
        gpu_charges[atom_offset + j] = 0.0;
      }
    }
    atom_offset += padded_natom;
  }
  check(gpu_charges, RelationalOperator::EQUAL, rbt_charges, "Charges pulled from the GPU in an "
        "AtomGraphSynthesis object do not meet expectations.", do_tests);

  // Allocate resources for various kernels
  CacheResource tb_space(nblocks, max_vwu_atoms);
  MolecularMechanicsControls mmctrl;
  mmctrl.primeWorkUnitCounters(gpu, poly_ag);
  ScoreCard sc(nsys, 1, 32);

  // Launch the valence evaluation kernel for small systems with only bonds, angles, dihedrals,
  // and 1:4 attenuated interactions.
  launchValenceSp(poly_ag, &mmctrl, &poly_ps, &sc, &tb_space, EvaluateForce::YES,
                  EvaluateEnergy::NO, VwuGoal::ACCUMULATE, ForceAccumulationMethod::SPLIT, gpu);
  std::vector<double> frc_mues(nsys);
  const std::vector<double> frc_mue_tolerance(nsys, 3.5e-5);
  std::vector<double> frc_max_errors(nsys);
  const std::vector<double> frc_max_error_tolerance(nsys, 2.0e-4);
  for (int i = 0; i < nsys; i++) {
    PhaseSpace devc_result = poly_ps.exportSystem(i, HybridTargetLevel::DEVICE);
    PhaseSpace host_result = poly_ps.exportSystem(i, HybridTargetLevel::HOST);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    evalValeMM(&host_result, &isc, poly_ag.getSystemTopologyPointer(i), EvaluateForce::YES, 0);
    const TrajectoryKind frcid = TrajectoryKind::FORCES;
    const std::vector<double> devc_frc = devc_result.getInterlacedCoordinates(frcid);
    const std::vector<double> host_frc = host_result.getInterlacedCoordinates(frcid);
    frc_mues[i] = meanUnsignedError(devc_frc, host_frc);
    frc_max_errors[i] = maxAbsoluteDifference(devc_frc, host_frc);
  }
  check(frc_mues, RelationalOperator::LESS_THAN, frc_mue_tolerance, "Forces obtained by the "
        "valence interaction kernel, operating on systems without external restraints, exceed the "
        "tolerance for mean unsigned errors in their vector components.", do_tests);
  check(frc_max_errors, RelationalOperator::LESS_THAN, frc_max_error_tolerance, "Forces obtained "
        "by the valence interaction kernel, operating on systems without external restraints, "
        "exceed the maximum allowed errors for forces acting on any one particle.", do_tests);
  poly_ps.initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl.incrementStep();
  launchValenceSp(poly_ag, &mmctrl, &poly_ps, &sc, &tb_space, EvaluateForce::YES,
                  EvaluateEnergy::NO, VwuGoal::ACCUMULATE, ForceAccumulationMethod::WHOLE, gpu);
  for (int i = 0; i < nsys; i++) {
    PhaseSpace devc_result = poly_ps.exportSystem(i, HybridTargetLevel::DEVICE);
    PhaseSpace host_result = poly_ps.exportSystem(i, HybridTargetLevel::HOST);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    evalValeMM(&host_result, &isc, poly_ag.getSystemTopologyPointer(i), EvaluateForce::YES, 0);
    const TrajectoryKind frcid = TrajectoryKind::FORCES;
    const std::vector<double> devc_frc = devc_result.getInterlacedCoordinates(frcid);
    const std::vector<double> host_frc = host_result.getInterlacedCoordinates(frcid);
    frc_mues[i] = meanUnsignedError(devc_frc, host_frc);
    frc_max_errors[i] = maxAbsoluteDifference(devc_frc, host_frc);
  }
  check(frc_mues, RelationalOperator::LESS_THAN, frc_mue_tolerance, "Forces obtained by the "
        "valence interaction kernel, operating on systems without external restraints, exceed the "
        "tolerance for mean unsigned errors in their vector components when accumulated in int64.",
        do_tests);
  check(frc_max_errors, RelationalOperator::LESS_THAN, frc_max_error_tolerance, "Forces obtained "
        "by the valence interaction kernel, operating on systems without external restraints, "
        "exceed the maximum allowed errors for forces acting on any one particle with "
        "fixed-precision accumulation in int64.", do_tests);
  
  // Create a set of larger systems, now involving CMAPs and other CHARMM force field terms
  const std::string topology_base = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpi_top_name = topology_base + osc + "trpcage.top";
  const std::string dhfr_top_name = topology_base + osc + "dhfr_cmap.top";
  const std::string alad_top_name = topology_base + osc + "ala_dipeptide.top";
  const std::string coordinate_base = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpi_crd_name = coordinate_base + osc + "trpcage.inpcrd";
  const std::string dhfr_crd_name = coordinate_base + osc + "dhfr_cmap.inpcrd";
  const std::string alad_crd_name = coordinate_base + osc + "ala_dipeptide.inpcrd";  
  const bool files_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                            getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                            getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(alad_crd_name) == DrivePathType::FILE);
  AtomGraph trpi_ag, dhfr_ag, alad_ag;
  PhaseSpace trpi_ps, dhfr_ps, alad_ps;
  if (files_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name);
    dhfr_ps.buildFromFile(dhfr_crd_name);
    alad_ps.buildFromFile(alad_crd_name);
  }
  else {
    rtWarn("Files for several systems in implicit solvent were not found.  Check the "
           "${OMNI_SOURCE} environment variable for validity.  Subsequent tests will be skipped.");
  }
#if 0
  const std::vector<AtomGraph*> bigger_tops = { &trpi_ag, &trpi_ag, &dhfr_ag, &dhfr_ag, &alad_ag,
                                                &trpi_ag, &dhfr_ag, &dhfr_ag, &alad_ag, &trpi_ag };
  const std::vector<PhaseSpace> bigger_crds = { trpi_ps, trpi_ps, dhfr_ps, dhfr_ps, alad_ps,
                                                trpi_ps, dhfr_ps, dhfr_ps, alad_ps, trpi_ps };
#endif
  for (int len = 4; len < 336; len += 4) {
    const std::vector<AtomGraph*> bigger_tops(len, &trpi_ag);
    const std::vector<PhaseSpace> bigger_crds(len, trpi_ps);
    std::vector<int> trpi_indices(len);
    for (int i = 0; i < len; i++) {
      trpi_indices[i] = i;
    }
    PhaseSpaceSynthesis big_poly_ps(bigger_crds, bigger_tops);
    AtomGraphSynthesis big_poly_ag(bigger_tops, trpi_indices, ExceptionResponse::SILENT,
                                   max_vwu_atoms, &timer);
    big_poly_ag.upload();
    big_poly_ps.upload();
    timer.assignTime(0);
    const int i_timings = timer.addCategory("GPU VWU evaluation " + std::to_string(len));
    const int cpu_timings = timer.addCategory("CPU VWU evaluation " + std::to_string(len));
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < 1000; j++) {
        mmctrl.incrementStep();
        launchValenceSp(big_poly_ag, &mmctrl, &big_poly_ps, &sc, &tb_space, EvaluateForce::YES,
                        EvaluateEnergy::YES, VwuGoal::ACCUMULATE, ForceAccumulationMethod::SPLIT,
                        gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(i_timings);
    }
  }
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
