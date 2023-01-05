// -*-c++-*-
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/kernel_manager.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_intake.h"
#include "../../src/Trajectory/coordinate_intake.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::random;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Initialize the test environment
  const TestEnvironment oe(argc, argv);
  StopWatch timer;

  // Section 1: test the thermostat mechanism
  section("Thermostat diagnostics");

  // Section 2: test the propagation of atoms against CPU results
  section("Compare GPU and CPU dynamics");

  // Section 3: test consistency of system propagation
  section("Self-consistency of GPU dynamics");
  
  // Read topology and starting coordinate files
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  bool files_ok;
  AtomGraph trpi_ag = loadTopology(trpi_top_name, &files_ok);
  const ImplicitSolventModel born_model = ImplicitSolventModel::HCT_GB;
  trpi_ag.setImplicitSolventModel(born_model);
  
  // Prep the GPU
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);

  // Try a rather weird means of initializing a vector of PhaseSpace objects, to test the
  // first-classness of the PhaseSpace object itself.
  std::vector<PhaseSpace> trpi_ps_vec(1, loadPhaseSpace(trpi_crd_name, &files_ok));
  const std::vector<AtomGraph*> trpi_ag_vec(1, &trpi_ag);
  const std::vector<int> trpi_tiling(1, 0);
  AtomGraphSynthesis trpi_poly_ag(trpi_ag_vec, trpi_tiling, ExceptionResponse::WARN, gpu, &timer);
  StaticExclusionMaskSynthesis trpi_poly_se(trpi_poly_ag.getTopologyPointers(),
                                            trpi_poly_ag.getTopologyIndices());
  trpi_poly_ag.loadNonbondedWorkUnits(trpi_poly_se, InitializationTask::GB_LANGEVIN_DYNAMICS, 15,
                                      gpu);
  PhaseSpaceSynthesis trpi_poly_ps(trpi_ps_vec, trpi_ag_vec, trpi_tiling);
  trpi_poly_ag.upload();
  trpi_poly_se.upload();
  trpi_poly_ps.upload();
  
  // Test the thermostat construction
  Thermostat trpi_heat_bath(ThermostatKind::LANGEVIN, 305.0);
  trpi_heat_bath.setAtomCount(trpi_poly_ag.getPaddedAtomCount());
  trpi_heat_bath.setRandomCacheDepth(trpi_poly_ag.getRandomCacheDepth());
  trpi_heat_bath.initializeRandomStates(PrecisionModel::SINGLE, 9815734, 25, gpu);
  std::vector<double> gpu_result, cpu_result;
  std::vector<ullint4> gpu_gstate, cpu_gstate;
  Xoshiro256ppGenerator xrs(9815734, 25);
  for (int i = 0; i < trpi_heat_bath.getAtomCount(); i+= 1024) {
    xrs.setState(trpi_heat_bath.getGeneratorState(i, HybridTargetLevel::HOST));
    for (int j = 0; j < trpi_heat_bath.getRandomCacheDepth() * 3; j++) {
      gpu_result.push_back(trpi_heat_bath.getCachedRandomResult(PrecisionModel::SINGLE, i, j,
                                                                HybridTargetLevel::DEVICE));
      cpu_result.push_back(xrs.spGaussianRandomNumber());
    }
    gpu_gstate.push_back(trpi_heat_bath.getGeneratorState(i, HybridTargetLevel::DEVICE));
    cpu_gstate.push_back(xrs.revealState());
  }
  check(gpu_result, RelationalOperator::EQUAL, Approx(cpu_result).margin(1.0e-6),
        "Random numbers generated by the CPU and GPU do not agree for a Langevin thermostat.");

  // Additional data structures needed for dynamics
  const KernelManager trpi_launcher(gpu, trpi_poly_ag);
  MolecularMechanicsControls trpi_mmctrl;
  trpi_mmctrl.primeWorkUnitCounters(trpi_launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                    PrecisionModel::SINGLE, trpi_poly_ag);
  const int2 vale_lp = trpi_launcher.getValenceKernelDims(PrecisionModel::SINGLE,
                                                          EvaluateForce::YES, EvaluateEnergy::NO,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::MOVE_PARTICLES,
                                                          ClashResponse::NONE);
  const int2 nonb_lp =
    trpi_launcher.getNonbondedKernelDims(PrecisionModel::SINGLE,
                                         trpi_poly_ag.getNonbondedWorkType(), EvaluateForce::YES,
                                         EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                                         born_model, ClashResponse::NONE);
  CacheResource trpi_vale_tb_space(vale_lp.x, maximum_valence_work_unit_atoms);
  CacheResource trpi_nonb_tb_space(nonb_lp.x, small_block_max_atoms);
  ImplicitSolventWorkspace trpi_isw(trpi_poly_ag.getSystemAtomOffsets(),
                                    trpi_poly_ag.getSystemAtomCounts(), PrecisionModel::SINGLE);
  ScoreCard trpi_sc(trpi_tiling.size(), 8, 32);
  const int ts_timings = timer.addCategory("Lazy time step");
  timer.assignTime(0);
  for (int i = 0; i < 10000; i++) {
    launchNonbonded(PrecisionModel::SINGLE, trpi_poly_ag, trpi_poly_se, &trpi_mmctrl,
                    &trpi_poly_ps, &trpi_heat_bath, &trpi_sc, &trpi_nonb_tb_space, &trpi_isw,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                    trpi_launcher);
    launchValence(PrecisionModel::SINGLE, trpi_poly_ag, &trpi_mmctrl, &trpi_poly_ps,
                  &trpi_heat_bath, &trpi_sc, &trpi_vale_tb_space, EvaluateForce::YES,
                  EvaluateEnergy::NO, VwuGoal::MOVE_PARTICLES, AccumulationMethod::SPLIT,
                  trpi_launcher);
    trpi_poly_ps.initializeForces(gpu, HybridTargetLevel::DEVICE);
    trpi_isw.initialize(HybridTargetLevel::DEVICE, CoordinateCycle::PRESENT, gpu);
    trpi_mmctrl.incrementStep();
    trpi_heat_bath.incrementStep();
  }
  cudaDeviceSynchronize();
  timer.assignTime(ts_timings);
  
  // Pre-construct the necessary abstracts and other information, then re-run the dynamics
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const SyValenceKit<float> trpi_poly_vk = trpi_poly_ag.getSinglePrecisionValenceKit(devc);
  const SeMaskSynthesisReader trpi_poly_ser = trpi_poly_se.data(devc);
  const SyNonbondedKit<float,
                       float2> trpi_poly_nbk = trpi_poly_ag.getSinglePrecisionNonbondedKit(devc);
  const SyRestraintKit<float, float2, float4> trpi_poly_rk =
    trpi_poly_ag.getSinglePrecisionRestraintKit(devc);
  const SyAtomUpdateKit<float, float2, float4> trpi_poly_auk =
    trpi_poly_ag.getSinglePrecisionAtomUpdateKit(devc);
  MMControlKit<float> trpi_ctrl = trpi_mmctrl.spData(devc);
  PsSynthesisWriter trpi_poly_psw = trpi_poly_ps.data(devc);
  ScoreCardWriter trpi_scw = trpi_sc.data(devc);
  CacheResourceKit<float> trpi_vale_r = trpi_vale_tb_space.spData(devc);
  CacheResourceKit<float> trpi_nonb_r = trpi_nonb_tb_space.spData(devc);
  ThermostatWriter<float> trpi_tstw = trpi_heat_bath.spData(devc);
  ISWorkspaceKit trpi_iswk = trpi_isw.spData(devc);
  const NbwuKind nbwu_type = trpi_poly_ag.getNonbondedWorkType();
  const int2 gbr_lp = trpi_launcher.getBornRadiiKernelDims(PrecisionModel::SINGLE, nbwu_type,
                                                           AccumulationMethod::SPLIT, born_model);
  const int2 gbd_lp = trpi_launcher.getBornDerivativeKernelDims(PrecisionModel::SINGLE, nbwu_type,
                                                                AccumulationMethod::SPLIT,
                                                                born_model);
  const int xts_timings = timer.addCategory("Fast time step");
  timer.assignTime(0);
  for (int i = 0; i < 10000; i++) {
    launchNonbonded(trpi_poly_ag.getNonbondedWorkType(), trpi_poly_nbk, trpi_poly_ser,
                    &trpi_ctrl, &trpi_poly_psw, &trpi_tstw, &trpi_scw, &trpi_nonb_r, &trpi_iswk,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT, gbr_lp,
                    nonb_lp, gbd_lp);
    launchValence(trpi_poly_vk, trpi_poly_rk, &trpi_ctrl, &trpi_poly_psw, trpi_poly_auk,
                  &trpi_tstw, &trpi_scw, &trpi_vale_r, EvaluateForce::YES, EvaluateEnergy::NO,
                  VwuGoal::MOVE_PARTICLES, AccumulationMethod::SPLIT, vale_lp);
    trpi_ctrl.step += 1;
  }
  cudaDeviceSynchronize();
  timer.assignTime(xts_timings);
  
  // Display timings and test results
  if (oe.getDisplayTimingsOrder()) {
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
