// -*-c++-*-
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/kernel_manager.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/scorecard.h"
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
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"

using namespace stormm::card;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::mm;
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
  
  // Read topology and starting coordinate files
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpi_crd_name = base_crd_name + osc + "myoglobin.inpcrd";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpi_top_name = base_top_name + osc + "myoglobin.top";
  bool files_ok;
  AtomGraph trpi_ag = loadTopology(trpi_top_name, ExceptionResponse::WARN, TopologyKind::AMBER,
                                   &files_ok);
  const ImplicitSolventModel born_model = ImplicitSolventModel::HCT_GB;
  trpi_ag.setImplicitSolventModel(born_model);

  // Prep the GPU
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);

  // Try a rather weird means of initializing a vector of PhaseSpace objects, to test the
  // first-classness of the PhaseSpace object itself.
  std::vector<PhaseSpace> trpi_ps_vec(1, loadPhaseSpace(trpi_crd_name, ExceptionResponse::WARN,
                                                        CoordinateFileKind::AMBER_ASCII_RST,
                                                        &files_ok));
  const std::vector<AtomGraph*> trpi_ag_vec(1, &trpi_ag);
  const std::vector<int> trpi_tiling(100, 0);
  AtomGraphSynthesis trpi_poly_ag(trpi_ag_vec, trpi_tiling, ExceptionResponse::WARN, gpu, &timer);
  StaticExclusionMaskSynthesis trpi_poly_se(trpi_poly_ag.getTopologyPointers(),
                                            trpi_poly_ag.getTopologyIndices());
  trpi_poly_ag.loadNonbondedWorkUnits(trpi_poly_se);
  PhaseSpaceSynthesis trpi_poly_ps(trpi_ps_vec, trpi_ag_vec, trpi_tiling);
  trpi_poly_ag.upload();
  trpi_poly_se.upload();
  trpi_poly_ps.upload();

  // Additional Data structures needed for dynamics
  const KernelManager trpi_launcher(gpu, trpi_poly_ag);
  MolecularMechanicsControls trpi_mmctrl;
  trpi_mmctrl.primeWorkUnitCounters(trpi_launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                    PrecisionModel::SINGLE, trpi_poly_ag);
  const int2 vale_lp = trpi_launcher.getValenceKernelDims(PrecisionModel::SINGLE,
                                                          EvaluateForce::YES, EvaluateEnergy::NO,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::MOVE_PARTICLES);
  const int2 nonb_lp =
    trpi_launcher.getNonbondedKernelDims(PrecisionModel::SINGLE,
                                         trpi_poly_ag.getNonbondedWorkType(), EvaluateForce::YES,
                                         EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                                         born_model);
  CacheResource trpi_vale_tb_space(vale_lp.x, maximum_valence_work_unit_atoms);
  CacheResource trpi_nonb_tb_space(nonb_lp.x, small_block_max_atoms);
  ImplicitSolventWorkspace trpi_isw(trpi_poly_ag.getSystemAtomOffsets(),
                                    trpi_poly_ag.getSystemAtomCounts(), PrecisionModel::SINGLE);
  ScoreCard trpi_sc(trpi_tiling.size(), 8, 32);
  const int ts_timings = timer.addCategory("Lazy time step");
  timer.assignTime(0);
  for (int i = 0; i < 1000; i++) {
    launchNonbonded(PrecisionModel::SINGLE, trpi_poly_ag, trpi_poly_se, &trpi_mmctrl,
                    &trpi_poly_ps, &trpi_sc, &trpi_nonb_tb_space, &trpi_isw, EvaluateForce::YES,
                    EvaluateEnergy::NO, AccumulationMethod::SPLIT, trpi_launcher);
    launchValence(PrecisionModel::SINGLE, trpi_poly_ag, &trpi_mmctrl, &trpi_poly_ps, &trpi_sc,
                  &trpi_vale_tb_space, EvaluateForce::YES, EvaluateEnergy::NO,
                  VwuGoal::MOVE_PARTICLES, AccumulationMethod::SPLIT, trpi_launcher);
    trpi_poly_ps.initializeForces(gpu, HybridTargetLevel::DEVICE);
    trpi_isw.initialize(HybridTargetLevel::DEVICE, CoordinateCycle::PRESENT, gpu);
    trpi_mmctrl.incrementStep();
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
  MMControlKit<float> trpi_ctrl = trpi_mmctrl.spData(devc);
  PsSynthesisWriter trpi_poly_psw = trpi_poly_ps.data(devc);
  ScoreCardWriter trpi_scw = trpi_sc.data(devc);
  CacheResourceKit<float> trpi_vale_r = trpi_vale_tb_space.spData(devc);
  CacheResourceKit<float> trpi_nonb_r = trpi_nonb_tb_space.spData(devc);
  ISWorkspaceKit trpi_iswk = trpi_isw.spData(devc);
  const NbwuKind nbwu_type = trpi_poly_ag.getNonbondedWorkType();
  const int2 gbr_lp = trpi_launcher.getBornRadiiKernelDims(PrecisionModel::SINGLE, nbwu_type,
                                                           AccumulationMethod::SPLIT, born_model);
  const int2 gbd_lp = trpi_launcher.getBornDerivativeKernelDims(PrecisionModel::SINGLE, nbwu_type,
                                                                AccumulationMethod::SPLIT,
                                                                born_model);
  const int xts_timings = timer.addCategory("Fast time step");
  timer.assignTime(0);
  for (int i = 0; i < 1000; i++) {
    launchNonbonded(trpi_poly_ag.getNonbondedWorkType(), trpi_poly_nbk, trpi_poly_ser,
                    &trpi_ctrl, &trpi_poly_psw, &trpi_scw, &trpi_nonb_r, &trpi_iswk,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT, gbr_lp,
                    nonb_lp, gbd_lp);
    launchValence(trpi_poly_vk, trpi_poly_rk, &trpi_ctrl, &trpi_poly_psw, &trpi_scw, &trpi_vale_r,
                  EvaluateForce::YES, EvaluateEnergy::NO, VwuGoal::MOVE_PARTICLES,
                  AccumulationMethod::SPLIT, vale_lp);
    trpi_poly_ps.initializeForces(gpu, HybridTargetLevel::DEVICE);
    trpi_isw.initialize(HybridTargetLevel::DEVICE, CoordinateCycle::PRESENT, gpu);
    trpi_ctrl.step += 1;
  }
  cudaDeviceSynchronize();
  timer.assignTime(xts_timings);

  // CHECK
  printf("There are %4d non-bonded work units.\n", trpi_poly_ag.getNonbondedWorkUnitCount());
  printf("There are %4d valence    work units.\n", trpi_poly_ag.getValenceWorkUnitCount());
  // END CHECK
  
  // Display time
  if (oe.getDisplayTimingsOrder()) {
    timer.printResults();
  }
  
  return 0;
}
