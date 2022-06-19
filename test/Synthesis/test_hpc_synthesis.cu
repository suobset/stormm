// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Constants/fixed_precision.h"
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
#include "../../src/Restraints/restraint_apparatus.h"
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
#include "assemble_restraints.h"

using namespace omni::card;
using namespace omni::constants;
using namespace omni::errors;
using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::math;
using namespace omni::mm;
using namespace omni::numerics;
using namespace omni::parse;
using namespace omni::restraints;
using namespace omni::synthesis;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

//-------------------------------------------------------------------------------------------------
// Compute forces due to valence interactions acting on a series of systems using topology and
// coordinate compilations.  Check the results against the accumulations for individual systems.
// This function assumes that the topology and coordinate compilations have already been uploaded
// to the device.
//
// Arguments:
//   poly_ps:        Coordinates for many systems
//   mmctrl:         Molecular mechanics progress counters
//   tb_space:       Thread block resources, pre-allocated on the GPU
//   poly_ag:        Topologies for many systems
//   facc_method:    Force accumulation method
//   prec:           Precision level at which to perform the calculations (may not be compatible
//                   with all force accumulation methods)
//   gpu:            Details of the GPU to use
//   mue_tol:        Tolerance for mean unsigned error in forces 
//   max_error_tol:  Tolerance for maximum unsigned error in forces 
//-------------------------------------------------------------------------------------------------
void checkCompilationForces(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                            CacheResource *tb_space, const AtomGraphSynthesis &poly_ag,
                            const ForceAccumulationMethod facc_method, const PrecisionLevel prec,
                            const GpuDetails &gpu, const double mue_tol,
                            const double max_error_tol, const TestPriority do_tests) {
  const int nsys = poly_ps->getSystemCount();
  std::vector<double> frc_mues(nsys);
  const std::vector<double> frc_mue_tolerance(nsys, mue_tol);
  std::vector<double> frc_max_errors(nsys);
  const std::vector<double> frc_max_error_tolerance(nsys, max_error_tol);
  ScoreCard sc(nsys, 1, 32);
  poly_ps->initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl->incrementStep();
  switch (prec) {
  case PrecisionLevel::SINGLE:
  case PrecisionLevel::SINGLE_PLUS:
    launchValenceSp(poly_ag, mmctrl, poly_ps, &sc, tb_space, EvaluateForce::YES,
                    EvaluateEnergy::NO, VwuGoal::ACCUMULATE, facc_method, gpu);
    break;
  case PrecisionLevel::DOUBLE:
    launchValenceDp(poly_ag, mmctrl, poly_ps, &sc, tb_space, EvaluateForce::YES,
                    EvaluateEnergy::NO, VwuGoal::ACCUMULATE, gpu);    
    break;
  }
  for (int i = 0; i < nsys; i++) {
    PhaseSpace host_result = poly_ps->exportSystem(i, HybridTargetLevel::HOST);
    PhaseSpace devc_result = poly_ps->exportSystem(i, HybridTargetLevel::DEVICE);
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
        "tolerance for mean unsigned errors in their vector components.  Force accumulation "
        "method: " + getForceAccumulationMethodName(facc_method) + ".  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(frc_max_errors, RelationalOperator::LESS_THAN, frc_max_error_tolerance, "Forces obtained "
        "by the valence interaction kernel, operating on systems without external restraints, "
        "exceed the maximum allowed errors for forces acting on any one particle.  Force "
        "accumulation method: " + getForceAccumulationMethodName(facc_method) + ".  Precision "
        "level in the calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Compute energies of a series of systems, using topology and coordinate compilations on the GPU,
// due to valence interactions.  Check the results against the accumulations for individual
// systems.  This function assumes that the topology and coordinate compilations have already been
// uploaded to the device.
//
// Arguments:
//-------------------------------------------------------------------------------------------------
void checkCompilationEnergies(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                              CacheResource *tb_space, const AtomGraphSynthesis &poly_ag,
                              const PrecisionLevel prec, const GpuDetails &gpu,
                              const double bond_tol, const double angl_tol, const double dihe_tol,
                              const double impr_tol, const double ubrd_tol, const double cimp_tol,
                              const double cmap_tol, const double lj14_tol, const double qq14_tol,
                              const double rstr_tol, const TestPriority do_tests) {
  const int nsys = poly_ps->getSystemCount();
  ScoreCard sc(nsys, 1, 32);
  poly_ps->initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl->incrementStep();
  switch (prec) {
  case PrecisionLevel::SINGLE:
  case PrecisionLevel::SINGLE_PLUS:
    launchValenceSp(poly_ag, mmctrl, poly_ps, &sc, tb_space, EvaluateForce::NO,
                    EvaluateEnergy::YES, VwuGoal::ACCUMULATE, ForceAccumulationMethod::SPLIT, gpu);
    break;
  case PrecisionLevel::DOUBLE:
    launchValenceDp(poly_ag, mmctrl, poly_ps, &sc, tb_space, EvaluateForce::NO,
                    EvaluateEnergy::YES, VwuGoal::ACCUMULATE, gpu);
    break;
  }
  sc.download();
  std::vector<double> cpu_bond(nsys), gpu_bond(nsys), cpu_angl(nsys), gpu_angl(nsys);
  std::vector<double> cpu_dihe(nsys), gpu_dihe(nsys), cpu_impr(nsys), gpu_impr(nsys);
  std::vector<double> cpu_ubrd(nsys), gpu_ubrd(nsys), cpu_cimp(nsys), gpu_cimp(nsys);
  std::vector<double> cpu_cmap(nsys), gpu_cmap(nsys), cpu_qq14(nsys), gpu_qq14(nsys);
  std::vector<double> cpu_lj14(nsys), gpu_lj14(nsys);
  for (int i = 0; i < nsys; i++) {
    PhaseSpace devc_result = poly_ps->exportSystem(i, HybridTargetLevel::DEVICE);
    PhaseSpace host_result = poly_ps->exportSystem(i, HybridTargetLevel::HOST);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    evalValeMM(&host_result, &isc, poly_ag.getSystemTopologyPointer(i), EvaluateForce::NO, 0);
    gpu_bond[i] =  sc.reportInstantaneousStates(StateVariable::BOND, i);
    cpu_bond[i] = isc.reportInstantaneousStates(StateVariable::BOND, 0);
    gpu_angl[i] =  sc.reportInstantaneousStates(StateVariable::ANGLE, i);
    cpu_angl[i] = isc.reportInstantaneousStates(StateVariable::ANGLE, 0);
    gpu_dihe[i] =  sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, i);
    cpu_dihe[i] = isc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, 0);
    gpu_impr[i] =  sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, i);
    cpu_impr[i] = isc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, 0);
    gpu_ubrd[i] =  sc.reportInstantaneousStates(StateVariable::UREY_BRADLEY, i);
    cpu_ubrd[i] = isc.reportInstantaneousStates(StateVariable::UREY_BRADLEY, 0);
    gpu_cimp[i] =  sc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER, i);
    cpu_cimp[i] = isc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER, 0);
    gpu_cmap[i] =  sc.reportInstantaneousStates(StateVariable::CMAP, i);
    cpu_cmap[i] = isc.reportInstantaneousStates(StateVariable::CMAP, 0);
    gpu_qq14[i] =  sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC_ONE_FOUR, i);
    cpu_qq14[i] = isc.reportInstantaneousStates(StateVariable::ELECTROSTATIC_ONE_FOUR, 0);
    gpu_lj14[i] =  sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, i);
    cpu_lj14[i] = isc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, 0);
  }
  check(gpu_bond, RelationalOperator::EQUAL, Approx(cpu_bond).margin(bond_tol), "Bond energies "
        "computed on the CPU and GPU do not agree.  Precision level in the calculation: " +
        getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_angl, RelationalOperator::EQUAL, Approx(cpu_angl).margin(angl_tol), "Angle energies "
        "computed on the CPU and GPU do not agree.  Precision level in the calculation: " +
        getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_dihe, RelationalOperator::EQUAL, Approx(cpu_dihe).margin(dihe_tol), "Proper "
        "dihedral energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_impr, RelationalOperator::EQUAL, Approx(cpu_impr).margin(impr_tol), "Improper "
        "dihedral energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_ubrd, RelationalOperator::EQUAL, Approx(cpu_ubrd).margin(ubrd_tol), "Urey-Bradley "
        "energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_cimp, RelationalOperator::EQUAL, Approx(cpu_cimp).margin(cimp_tol), "CHARMM "
        "improper dihedral energies computed on the CPU and GPU do not agree.  Precision level in "
        "the calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_cmap, RelationalOperator::EQUAL, Approx(cpu_cmap).margin(cmap_tol), "CMAP "
        "energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_qq14, RelationalOperator::EQUAL, Approx(cpu_qq14).margin(qq14_tol), "Electrostatic "
        "1:4 energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
  check(gpu_lj14, RelationalOperator::EQUAL, Approx(cpu_lj14).margin(lj14_tol), "Lennard-Jones "
        "1:4 energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getPrecisionLevelName(prec) + ".", do_tests);
}

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
  checkCompilationForces(&poly_ps, &mmctrl, &tb_space, poly_ag, ForceAccumulationMethod::WHOLE,
                         PrecisionLevel::DOUBLE, gpu, 3.5e-6, 2.0e-5, do_tests);
  checkCompilationForces(&poly_ps, &mmctrl, &tb_space, poly_ag, ForceAccumulationMethod::SPLIT,
                         PrecisionLevel::SINGLE, gpu, 3.5e-5, 2.0e-4, do_tests);
  checkCompilationForces(&poly_ps, &mmctrl, &tb_space, poly_ag, ForceAccumulationMethod::WHOLE,
                         PrecisionLevel::SINGLE, gpu, 3.5e-5, 2.0e-4, do_tests);
  checkCompilationEnergies(&poly_ps, &mmctrl, &tb_space, poly_ag, PrecisionLevel::DOUBLE,
                           gpu, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                           1.0e-6, 1.0e-6, do_tests);
  checkCompilationEnergies(&poly_ps, &mmctrl, &tb_space, poly_ag, PrecisionLevel::SINGLE,
                           gpu, 1.5e-5, 1.5e-5, 5.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 6.0e-6,
                           2.2e-5, 1.0e-6, do_tests);
  
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

  // Read some larger topologies, with CHARMM CMAP and other force field terms
  const std::vector<AtomGraph*> bigger_tops = { &trpi_ag, &dhfr_ag, &alad_ag };
  const std::vector<PhaseSpace> bigger_crds = { trpi_ps, dhfr_ps, alad_ps };
  PhaseSpaceSynthesis big_poly_ps(bigger_crds, bigger_tops);
  const std::vector<int> big_top_indices = { 0, 1, 2 };
  AtomGraphSynthesis big_poly_ag(bigger_tops, big_top_indices, ExceptionResponse::SILENT,
                                 max_vwu_atoms, &timer);
  big_poly_ag.upload();
  big_poly_ps.upload();
  timer.assignTime(0);
  checkCompilationForces(&big_poly_ps, &mmctrl, &tb_space, big_poly_ag,
                         ForceAccumulationMethod::WHOLE, PrecisionLevel::DOUBLE, gpu, 3.5e-6,
                         2.0e-5, do_tests);
  checkCompilationForces(&big_poly_ps, &mmctrl, &tb_space, big_poly_ag,
                         ForceAccumulationMethod::SPLIT, PrecisionLevel::SINGLE, gpu, 7.5e-5,
                         3.0e-3, do_tests);
  checkCompilationForces(&big_poly_ps, &mmctrl, &tb_space, big_poly_ag,
                         ForceAccumulationMethod::WHOLE, PrecisionLevel::SINGLE, gpu, 7.5e-5,
                         3.0e-3, do_tests);
  checkCompilationEnergies(&big_poly_ps, &mmctrl, &tb_space, big_poly_ag, PrecisionLevel::DOUBLE,
                           gpu, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 6.0e-6, 1.0e-6,
                           1.0e-6, 1.0e-6, do_tests);
  checkCompilationEnergies(&big_poly_ps, &mmctrl, &tb_space, big_poly_ag, PrecisionLevel::SINGLE,
                           gpu, 1.5e-4, 2.2e-5, 9.0e-5, 1.5e-5, 6.0e-5, 3.0e-5, 6.0e-6, 7.5e-5,
                           2.2e-4, 1.0e-6, do_tests);  

  // Read some topologies with virtual sites.  First, test the forces that appear to act on the
  // virtual sites.  Add restraints to these ligands.
  const std::string brbz_top_name = topology_base + osc + "bromobenzene_iso.top";
  const std::string lig1_top_name = topology_base + osc + "stereo_L1.top";
  const std::string lig2_top_name = topology_base + osc + "symmetry_L1.top";
  const std::string brbz_crd_name = coordinate_base + osc + "bromobenzene_iso.inpcrd";
  const std::string lig1_crd_name = coordinate_base + osc + "stereo_L1.inpcrd";
  const std::string lig2_crd_name = coordinate_base + osc + "symmetry_L1.inpcrd";
  const bool ligands_exist = (getDrivePathType(brbz_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig1_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig2_top_name) == DrivePathType::FILE &&
                              getDrivePathType(brbz_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig1_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig2_crd_name) == DrivePathType::FILE);
  AtomGraph brbz_ag, lig1_ag, lig2_ag;
  PhaseSpace brbz_ps, lig1_ps, lig2_ps;
  if (ligands_exist) {
    brbz_ag.buildFromPrmtop(brbz_top_name);
    lig1_ag.buildFromPrmtop(lig1_top_name);
    lig2_ag.buildFromPrmtop(lig2_top_name);
    brbz_ps.buildFromFile(brbz_crd_name);
    lig1_ps.buildFromFile(lig1_crd_name);
    lig2_ps.buildFromFile(lig2_crd_name);
  }
  RestraintApparatus brbz_ra = assembleRestraints(&brbz_ag, brbz_ps);
  RestraintApparatus lig1_ra = assembleRestraints(&lig1_ag, lig1_ps);
  RestraintApparatus lig2_ra = assembleRestraints(&lig2_ag, lig2_ps);
  const std::vector<AtomGraph*> ligand_ag_list = { &brbz_ag, &lig1_ag, &lig2_ag };
  const std::vector<PhaseSpace> ligand_ps_list = {  brbz_ps,  lig1_ps,  lig2_ps };
  const std::vector<RestraintApparatus*> ligand_ra_list = { &brbz_ra, &lig1_ra, &lig2_ra };
  const std::vector<int> ligand_tiling = { 0, 1, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 1, 0 };
  PhaseSpaceSynthesis ligand_poly_ps(ligand_ps_list, ligand_ag_list, ligand_tiling);
  AtomGraphSynthesis ligand_poly_ag(ligand_ag_list, //ligand_ra_list, ligand_tiling,
                                    ligand_tiling,
                                    ExceptionResponse::WARN, max_vwu_atoms, &timer);
  ligand_poly_ag.upload();
  ligand_poly_ps.upload();
  timer.assignTime(0);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &tb_space, ligand_poly_ag,
                         ForceAccumulationMethod::WHOLE, PrecisionLevel::DOUBLE, gpu, 3.5e-6,
                         2.0e-5, do_tests);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &tb_space, ligand_poly_ag,
                         ForceAccumulationMethod::SPLIT, PrecisionLevel::SINGLE, gpu, 7.5e-5,
                         3.0e-3, do_tests);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &tb_space, ligand_poly_ag,
                         ForceAccumulationMethod::WHOLE, PrecisionLevel::SINGLE, gpu, 7.5e-5,
                         3.0e-3, do_tests);
  checkCompilationEnergies(&ligand_poly_ps, &mmctrl, &tb_space, ligand_poly_ag,
                           PrecisionLevel::DOUBLE, gpu, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                           1.0e-6, 6.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, do_tests);
  checkCompilationEnergies(&ligand_poly_ps, &mmctrl, &tb_space, ligand_poly_ag,
                           PrecisionLevel::SINGLE, gpu, 1.5e-4, 2.2e-5, 9.0e-5, 1.5e-5, 6.0e-5,
                           3.0e-5, 6.0e-6, 7.5e-5, 2.2e-4, 1.0e-6, do_tests);
  
#if 0
  for (int len = 4; len < 36; len += 4) {
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
#endif

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
