// -*-c++-*-
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/namelist_emulator.h"
#include "../../src/Numerics/numeric_enumerators.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/hub_and_spoke.h"
#include "../../src/Structure/settle.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/integration.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::symbols;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

#include "../../src/Math/rounding.cui"

//-------------------------------------------------------------------------------------------------
// Copy the initial velocities and forces into the coordinate synthesis, preparing the system for
// another force calculation.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
initSynthesisForceVelocity(PsSynthesisWriter poly_psw, const llint* frc_init,
                           const int* frc_init_ovrf, const llint* vel_init,
                           const int* vel_init_ovrf) {
  const int last_sys = poly_psw.system_count - 1;
  const int padded_natom = poly_psw.atom_starts[last_sys] +
                           devcRoundUp(poly_psw.atom_counts[last_sys], warp_size_int);
  const int adv = blockDim.x * gridDim.x;
  for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < padded_natom; i += adv) {
    poly_psw.xfrc[i] = frc_init[i                     ];
    poly_psw.yfrc[i] = frc_init[i +      padded_natom ];
    poly_psw.zfrc[i] = frc_init[i + (2 * padded_natom)];
    if (poly_psw.frc_bits > force_scale_nonoverflow_bits) {
      poly_psw.xfrc_ovrf[i] = frc_init_ovrf[i                     ];
      poly_psw.yfrc_ovrf[i] = frc_init_ovrf[i +      padded_natom ];
      poly_psw.zfrc_ovrf[i] = frc_init_ovrf[i + (2 * padded_natom)];
    }
    poly_psw.xvel[i] = vel_init[i                     ];
    poly_psw.yvel[i] = vel_init[i +      padded_natom ];
    poly_psw.zvel[i] = vel_init[i + (2 * padded_natom)];
    if (poly_psw.vel_bits > velocity_scale_nonoverflow_bits) {
      poly_psw.xvel_ovrf[i] = vel_init_ovrf[i                     ];
      poly_psw.yvel_ovrf[i] = vel_init_ovrf[i +      padded_natom ];
      poly_psw.zvel_ovrf[i] = vel_init_ovrf[i + (2 * padded_natom)];
    }
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Baseline variables
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif

  // Lay out time categories for program profiling
  const int time_load  = timer.addCategory("Load input files");
  const int time_build = timer.addCategory("Class object setup");
  const int time_test  = timer.addCategory("CPU testing");
  const int time_init  = timer.addCategory("Initialization");
  const int time_valn  = timer.addCategory("Valence kernel with initialization");

  // Take in command line inputs
  CommandLineParser clip("valence_interactions", "A benchmarking program for measuring the rate "
                         "at which various GPU kernels can compute all valence interactions and "
                         "then, if requested, move particles.", { "-timings" });
  clip.addStandardAmberInputs("-c", "-p", "-ig_seed", "-temp0");
  clip.addStandardBenchmarkingInputs({ "-iter", "-trials", "-replicas", "-precision",
                                       "-skip_cpu_check", "-eval_nrg", "-omit_frc", "-fp_bits",
                                       "-blur" });

  // Custom inputs for this benchmarking program
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-move", NamelistType::BOOLEAN);
  t_nml->addHelp("-move", "Request that atom positions be updated based on the computed forces");
  t_nml->addKeyword("-rigid_geom", NamelistType::BOOLEAN);
  t_nml->addHelp("-rigid_geom", "Request that geometric constraints (rigid bonds to hydrogen "
                 "atoms) be enforced");
  t_nml->addKeyword("-acc", NamelistType::STRING, std::string("split"));
  t_nml->addHelp("-acc", "Specify how to accumulate forces.  Fixed-precision is always in effect, "
                 "and the \"split\" fixed-precision accumulation is obligatory for "
                 "\"double\"-precision computations.  For most architectures in single-precision "
                 "mode, \"split\" accumulation (two int32_t accumulators) is faster than "
                 "\"whole\" accumulation (one int64_t accumulator).");
  
  // Initialize the testing environment such that it cooperates with this program's own
  // CommandLineParser to read user input.
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Read command line instructions
  clip.parseUserInput(argc, argv);
  const std::string inpcrd_file = t_nml->getStringValue("-c");
  const std::string topology_file = t_nml->getStringValue("-p");
  const bool skip_cpu_check = t_nml->getBoolValue("-skip_cpu_check");
  const bool enforce_rigid_geom = t_nml->getBoolValue("-rigid_geom");
  const int n_trials = t_nml->getIntValue("-trials");
  const int n_repeats = t_nml->getIntValue("-iter");
  const int n_replicas = t_nml->getIntValue("-replicas");
  const int fp_bits = t_nml->getIntValue("-fp_bits");
  const int ig_seed = t_nml->getIntValue("-ig_seed");
  const double blur = t_nml->getRealValue("-blur");
  const double temperature = t_nml->getRealValue("-temp0");
  const PrecisionModel prec = translatePrecisionModel(t_nml->getStringValue("-precision"));
  AccumulationMethod acc_meth;
  switch (translateAccumulationMethod(t_nml->getStringValue("-acc"))) {
  case AccumulationMethod::AUTOMATIC:
  case AccumulationMethod::SPLIT:
    acc_meth = AccumulationMethod::SPLIT;
    break;
  case AccumulationMethod::WHOLE:
    acc_meth = AccumulationMethod::WHOLE;
    break;
  }
  EvaluateEnergy eval_nrg = t_nml->getBoolValue("-eval_nrg") ? EvaluateEnergy::YES:
                                                               EvaluateEnergy::NO;
  EvaluateForce eval_frc;
  if (t_nml->getBoolValue("-omit_frc")) {

    // Force the evaluation of energy if forces are to be omitted.
    eval_nrg = EvaluateEnergy::YES;
    eval_frc = EvaluateForce::NO;
  }
  else {
    eval_frc = EvaluateForce::YES;
  }
  const VwuGoal purpose = (t_nml->getBoolValue("-move")) ? VwuGoal::MOVE_PARTICLES :
                                                           VwuGoal::ACCUMULATE;

  // A Hybrid object was created to engage the GPU.  Absorb any bootup time into "miscellaneous."
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
  }
#ifdef STORMM_USE_HPC
  DynamicsControls dyncon;
  if (enforce_rigid_geom) {
    dyncon.setGeometricConstraints(ApplyConstraints::YES);
  }
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), dyncon, ExceptionResponse::DIE);
  timer.assignTime(time_load);
  MolecularMechanicsControls mmctrl;
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(n_replicas, 0),
                                                              blur, ig_seed, 36, 44, fp_bits);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(n_replicas, 0),
                                                            ExceptionResponse::DIE, gpu);
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    {
      const StaticExclusionMaskSynthesis poly_se(poly_ps.getSystemTopologyPointer(),
                                                 std::vector<int>(poly_ps.getSystemCount(), 0));
      poly_ag.loadNonbondedWorkUnits(poly_se);
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
  ScoreCard sc(poly_ps.getSystemCount(), 1, 32);
  const CoreKlManager launcher(gpu, poly_ag);
  const int2 valence_lp = launcher.getValenceKernelDims(prec, eval_frc, eval_nrg, acc_meth,
                                                        purpose, ClashResponse::NONE);
  CacheResource valence_tb_space(valence_lp.x, maximum_valence_work_unit_atoms);
  mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                               purpose, prec, prec, poly_ag);
  Thermostat tstat(poly_ag, ThermostatKind::NONE, temperature);
  if (enforce_rigid_geom) {
    tstat.setGeometryConstraints(ApplyConstraints::YES);
  }

  // Create some vectors of random forces and velocities to copy into the synthesis in order to
  // "initialize" the non-bonded forces acting on all particles and their velocities.  This takes
  // advantage of the contiguity of force and velocity arrays in the PhaseSpaceSynthesis.
  const size_t nvalue = 3 * poly_ag.getPaddedAtomCount();
  Hybrid<llint> vel_init(nvalue, "vel_init_x");
  Hybrid<int> vel_init_ovrf(nvalue, "vel_init_x_ovrf");
  Hybrid<llint> frc_init(nvalue, "frc_init_x");
  Hybrid<int> frc_init_ovrf(nvalue, "frc_init_x_ovrf");
  Xoshiro256ppGenerator xrs(ig_seed);
  addRandomNoise(&xrs, &frc_init, &frc_init_ovrf, 10.0, pow(2.0, fp_bits));
  llint* vel_ptr = vel_init.data();
  int* vel_ovrf_ptr = vel_init_ovrf.data();
  const SyAtomUpdateKit<double,
                        double2, double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit();
  PsSynthesisWriter host_psw = poly_ps.data();
  const double ebeta = sqrt(boltzmann_constant_gafs * temperature) * host_psw.vel_scale_f;
  for (size_t i = 0; i < nvalue; i++) {
    const double mss = poly_auk.inv_masses[i];

    // The fixed-precision velocity scaling is folded into ebeta
    const double dv = (mss > small) ? ebeta * sqrt(mss) * xrs.gaussianRandomNumber() : 0.0;
    const int95_t idv = hostDoubleToInt95(dv);
    vel_ptr[i] = idv.x;
    vel_ovrf_ptr[i] = idv.y;
  }

  // Idealize the geometry of water molecules in each structure, after the blurring.  Load the
  // initial forces (intended to mimic non-bonded interactions) and velocities into each
  // individual system as appropriate.
  std::vector<PhaseSpace> ps_vec;
  const std::vector<AtomGraph*> ag_vec = poly_ps.getSystemTopologyPointer();
  ps_vec.reserve(poly_ps.getSystemCount());
  const double inv_frc_scale = pow(2.0, -poly_ps.getForceAccumulationBits());
  const double inv_vel_scale = pow(2.0, -poly_ps.getVelocityBits());
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    ps_vec.emplace_back(poly_ps.exportSystem(i), HybridFormat::HOST_ONLY);
    idealizeSettleReference(&ps_vec.back(), ag_vec[i], PrecisionModel::DOUBLE);
    PhaseSpaceWriter psw = ps_vec.back().data();
    const int atom_offset = poly_ps.getAtomOffset(i);
    const llint* xfrc_prim_ptr = frc_init.data();
    const llint* yfrc_prim_ptr = &xfrc_prim_ptr[poly_ag.getPaddedAtomCount()];
    const llint* zfrc_prim_ptr = &yfrc_prim_ptr[poly_ag.getPaddedAtomCount()];
    const int* xfrc_ovrf_ptr = frc_init_ovrf.data();
    const int* yfrc_ovrf_ptr = &xfrc_ovrf_ptr[poly_ag.getPaddedAtomCount()];
    const int* zfrc_ovrf_ptr = &yfrc_ovrf_ptr[poly_ag.getPaddedAtomCount()];
    const llint* xvel_prim_ptr = vel_init.data();
    const llint* yvel_prim_ptr = &xvel_prim_ptr[poly_ag.getPaddedAtomCount()];
    const llint* zvel_prim_ptr = &yvel_prim_ptr[poly_ag.getPaddedAtomCount()];
    const int* xvel_ovrf_ptr = vel_init_ovrf.data();
    const int* yvel_ovrf_ptr = &xvel_ovrf_ptr[poly_ag.getPaddedAtomCount()];
    const int* zvel_ovrf_ptr = &yvel_ovrf_ptr[poly_ag.getPaddedAtomCount()];
    for (int j = 0; j < psw.natom; j++) {
      const size_t jo = j + atom_offset;
      psw.xfrc[j] = hostInt95ToDouble(xfrc_prim_ptr[jo], xfrc_ovrf_ptr[jo]) * inv_frc_scale;
      psw.yfrc[j] = hostInt95ToDouble(yfrc_prim_ptr[jo], yfrc_ovrf_ptr[jo]) * inv_frc_scale;
      psw.zfrc[j] = hostInt95ToDouble(zfrc_prim_ptr[jo], zfrc_ovrf_ptr[jo]) * inv_frc_scale;
      psw.xvel[j] = hostInt95ToDouble(xvel_prim_ptr[jo], xvel_ovrf_ptr[jo]) * inv_vel_scale;
      psw.yvel[j] = hostInt95ToDouble(yvel_prim_ptr[jo], yvel_ovrf_ptr[jo]) * inv_vel_scale;
      psw.zvel[j] = hostInt95ToDouble(zvel_prim_ptr[jo], zvel_ovrf_ptr[jo]) * inv_vel_scale;
    }
    poly_ps.importSystem(ps_vec.back(), i);
  }
  
  // Upload the critical components
  poly_ag.upload();
  poly_ps.upload();
  frc_init.upload();
  frc_init_ovrf.upload();
  vel_init.upload();
  vel_init_ovrf.upload();
  timer.assignTime(time_build);
  
  // Benchmark the force and velocity initialization process
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter devc_psw = poly_ps.data(devc_tier);
  const llint* devc_frc_init = frc_init.data(devc_tier);
  const llint* devc_vel_init = vel_init.data(devc_tier);
  const int* devc_frc_init_ovrf = frc_init_ovrf.data(devc_tier);
  const int* devc_vel_init_ovrf = vel_init_ovrf.data(devc_tier);
  timer.assignTime(0);
  for (int i = 0; i < n_trials; i++) {
    for (int j = 0; j < n_repeats; j++) {
      initSynthesisForceVelocity<<<gpu.getSMPCount(),
                                   large_block_size>>>(devc_psw, devc_frc_init, devc_frc_init_ovrf,
                                                       devc_vel_init, devc_vel_init_ovrf);
    }
    cudaDeviceSynchronize();  
    timer.assignTime(time_init);
  }
  
  // Check the topology synthesis for the boundary conditions on each system.  If periodic boundary
  // conditions are in effect, create a 

  // Run the valence kernel
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyRestraintKit<double,
                           double2,
                           double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit(devc_tier);
      const SyAtomUpdateKit<double,
                            double2,
                            double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      MMControlKit<double> ctrl = mmctrl.dpData(devc_tier);
      ThermostatWriter<double> tstw = tstat.dpData(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<double> gmem_rval = valence_tb_space.dpData(devc_tier);
      timer.assignTime(0);
      for (int i = 0; i < n_trials; i++) {
        for (int j = 0; j < n_repeats; j++) {
          initSynthesisForceVelocity<<<gpu.getSMPCount(),
                                       large_block_size>>>(devc_psw, devc_frc_init,
                                                           devc_frc_init_ovrf, devc_vel_init,
                                                           devc_vel_init_ovrf);
          ctrl.step += 1;
          launchValence(poly_vk, poly_rk, &ctrl, &devc_psw, poly_auk, &tstw, &scw, &gmem_rval,
                        eval_frc, eval_nrg, purpose, valence_lp);
        }
        cudaDeviceSynchronize();
        timer.assignTime(time_valn);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyRestraintKit<float,
                           float2,
                           float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit(devc_tier);
      const SyAtomUpdateKit<float,
                            float2,
                            float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      MMControlKit<float> ctrl = mmctrl.spData(devc_tier);
      ThermostatWriter<float> tstw = tstat.spData(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<float> gmem_rval = valence_tb_space.spData(devc_tier);
      timer.assignTime(0);
      for (int i = 0; i < n_trials; i++) {
        for (int j = 0; j < n_repeats; j++) {
          initSynthesisForceVelocity<<<gpu.getSMPCount(),
                                       large_block_size>>>(devc_psw, devc_frc_init,
                                                           devc_frc_init_ovrf, devc_vel_init,
                                                           devc_vel_init_ovrf);
          ctrl.step += 1;
          launchValence(poly_vk, poly_rk, &ctrl, &devc_psw, poly_auk, &tstw, &scw, &gmem_rval,
                        eval_frc, eval_nrg, purpose, acc_meth, valence_lp);
        }
        cudaDeviceSynchronize();
        timer.assignTime(time_valn);
      }
    }
    break;
  }
  
  // Run a separate check on the validity of the forces
  if (skip_cpu_check == false) {
    double frc_tol, vel_tol, pos_tol;
    switch (prec) {
    case PrecisionModel::DOUBLE:
      frc_tol = 1.0e-6;
      vel_tol = 1.0e-6;
      break;
    case PrecisionModel::SINGLE:
      frc_tol = 1.8e-3;
      vel_tol = 1.0e-6;
      break;
    }
    const TrajectoryKind tfrc = TrajectoryKind::FORCES;
    poly_ps.download();
    switch (purpose) {
    case VwuGoal::MOVE_PARTICLES:

      // Check that the positions and velocities of particles are consistent--the forces will be
      // unavailable, but evaluating them becomes necessary, and the energy is produced as a
      // consequence of any CPU-based calculation.
      for (int i = 0; i < poly_ps.getSystemCount(); i++) {
        evalValeMM(&ps_vec[i], &sc, ag_vec[i], eval_frc, i);
        PhaseSpaceWriter psw_cpu = ps_vec[i].data();
        const VirtualSiteKit<double> vsk = ag_vec[i]->getDoublePrecisionVirtualSiteKit();
        const ConstraintKit<double> cnk = ag_vec[i]->getDoublePrecisionConstraintKit();
        const ChemicalDetailsKit cdk = ag_vec[i]->getChemicalDetailsKit();
        transmitVirtualSiteForces(&psw_cpu, vsk);
        Thermostat tstat_i(psw_cpu.natom, ThermostatKind::NONE, temperature);
        if (enforce_rigid_geom) {
          tstat_i.setGeometryConstraints(ApplyConstraints::YES);
        }
        const ThermostatReader tstr = tstat_i.dpData();
        velocityVerletVelocityUpdate(&psw_cpu, cdk, tstr);
        if (enforce_rigid_geom) {
          rattleVelocities(&psw_cpu, cnk, tstr.dt, tstr.rattle_tol, 100, RattleMethod::CENTER_SUM);
          settleVelocities(&psw_cpu, cnk);
        }
        velocityVerletCoordinateUpdate(&psw_cpu, cdk, tstr);
        if (enforce_rigid_geom) {
          shakePositions(&psw_cpu, cnk, tstr.dt, tstr.rattle_tol, 100, RattleMethod::CENTER_SUM);
          settlePositions<double, double3>(&psw_cpu, cnk, tstr.dt);
        }
        const PhaseSpace psi_gpu = poly_ps.exportSystem(i);
        const PhaseSpaceReader psr_gpu = psi_gpu.data();
        const CoordinateCycle bcyc = CoordinateCycle::BLACK;
        const TrajectoryKind tvel = TrajectoryKind::VELOCITIES;
        const std::vector<double> cpu_vel = ps_vec[i].getInterlacedCoordinates(bcyc, tvel);
        const std::vector<double> gpu_vel = psi_gpu.getInterlacedCoordinates(bcyc, tvel);
        check(gpu_vel, RelationalOperator::EQUAL, Approx(cpu_vel).margin(vel_tol), "Velocities "
              "computed by CPU and GPU-based routines do not agree in replica " +
              std::to_string(i) + ".");
        const TrajectoryKind tpos = TrajectoryKind::POSITIONS;
        const std::vector<double> cpu_loc = ps_vec[i].getInterlacedCoordinates(bcyc, tpos);
        const std::vector<double> gpu_loc = psi_gpu.getInterlacedCoordinates(bcyc, tpos);
        check(gpu_loc, RelationalOperator::EQUAL, Approx(cpu_loc).margin(vel_tol), "Velocities "
              "computed by CPU and GPU-based routines do not agree in replica " +
              std::to_string(i) + ".");
      }
      break;
    case VwuGoal::ACCUMULATE:

      // Check the forces, positions, and velocities of particles.
      switch (eval_frc) {
      case EvaluateForce::YES:
        for (int i = 0; i < poly_ps.getSystemCount(); i++) {
          evalValeMM(&ps_vec[i], &sc, ag_vec[i], eval_frc, i);
          const PhaseSpace psi_gpu = poly_ps.exportSystem(i);
          const std::vector<double> cpu_frc = ps_vec[i].getInterlacedCoordinates(tfrc);
          const std::vector<double> gpu_frc = psi_gpu.getInterlacedCoordinates(tfrc);
          check(gpu_frc, RelationalOperator::EQUAL, Approx(cpu_frc).margin(frc_tol), "Forces "
                "computed by CPU and GPU-based routines do not agree in replica " +
                std::to_string(i) + ".");
        }
        break;
      case EvaluateForce::NO:
        break;
      }
      break;
    }

  }
#else // STORMM_USE_HPC
  rtWarn("This benchmarking program requires GPU support to run.", "valence_interactions");
#endif // STORMM_USE_HPC
    
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return 0;
}

