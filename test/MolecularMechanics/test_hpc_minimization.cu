// -*-c++-*-
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/kernel_manager.h"
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Constants/behavior.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/reduction_abstracts.h"
#include "../../src/Math/reduction_bridge.h"
#include "../../src/Math/reduction_enumerators.h"
#include "../../src/Math/hpc_reduction.h"
#include "../../src/MolecularMechanics/hpc_minimization.h"
#include "../../src/MolecularMechanics/line_minimization.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace omni::card;
using namespace omni::constants;
using namespace omni::diskutil;
using namespace omni::energy;
using namespace omni::errors;
using namespace omni::math;
using namespace omni::mm;
using namespace omni::synthesis;
using namespace omni::testing;
using namespace omni::topology;
using namespace omni::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  
  // Section 1
  section("Minimize a collection of drug molecules and dipeptides");

  // Read small molecules and compile them into a synthesis
  const char osc = osSeparator();
  const std::string base_top_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene_iso.top";
  const std::string brbz_crd_name = base_crd_name + osc + "bromobenzene_iso.inpcrd";
  const std::string lig1_top_name = base_top_name + osc + "stereo_L1.top";
  const std::string lig1_crd_name = base_crd_name + osc + "stereo_L1.inpcrd";
  const std::string lig2_top_name = base_top_name + osc + "symmetry_L1.top";
  const std::string lig2_crd_name = base_crd_name + osc + "symmetry_L1.inpcrd";
  const std::vector<std::string> all_top = { alad_top_name, brbz_top_name, lig1_top_name,
                                             lig2_top_name };
  const std::vector<std::string> all_crd = { alad_crd_name, brbz_crd_name, lig1_crd_name,
                                             lig2_crd_name };
  const int small_mol_count = all_top.size();
  bool files_exist = true;
  for (int i = 0; i < small_mol_count; i++) {
    files_exist = (getDrivePathType(all_top[i]) == DrivePathType::FILE &&
                   getDrivePathType(all_crd[i]) == DrivePathType::FILE && files_exist);
  }
  std::vector<AtomGraph> small_mol_ag;
  std::vector<AtomGraph*> small_mol_ag_ptr;
  std::vector<PhaseSpace> small_mol_ps;
  if (files_exist) {
    small_mol_ag.reserve(small_mol_count);
    small_mol_ps.reserve(small_mol_count);
    small_mol_ag_ptr.resize(small_mol_count);
    for (int i = 0; i < small_mol_count; i++) {
      small_mol_ag.emplace_back(all_top[i], ExceptionResponse::SILENT);
      small_mol_ps.emplace_back(all_crd[i]);
      small_mol_ag_ptr[i] = &small_mol_ag[i];
    }
  }
  else {
    small_mol_ag.resize(small_mol_count);
    small_mol_ag.resize(small_mol_count);
    rtWarn("Topology and coordinate files for a number of small molecules and dipeptides were not "
           "found.  Check the ${OMNI_SOURCE} environment variable, currently set to " +
           oe.getOmniSourcePath() + ", for validity.  Subsequent tests will be skipped.",
           "test_hpc_minimization");
  }
  std::vector<int> small_mol_id = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 3, 1, 2, 2, 1, 3, 0 };
  //std::vector<int> small_mol_id = { 2, 3, 2, 3, 0, 1, 2, 3, 2, 3, 3, 2, 2, 3, 3, 2 };
  small_mol_id.resize(8192);
  for (int i = 0; i < 512; i++) {
    for (int j = 0; j < 16; j++) {
      small_mol_id[(16 * i) + j] = small_mol_id[j];
    }
  }
  AtomGraphSynthesis small_poly_ag(small_mol_ag_ptr, small_mol_id, ExceptionResponse::WARN, gpu,
                                   &timer);
  StaticExclusionMaskSynthesis small_poly_se(small_poly_ag.getTopologyPointers(), small_mol_id);
  small_poly_ag.loadNonbondedWorkUnits(small_poly_se);
  PhaseSpaceSynthesis small_poly_ps(small_mol_ps, small_mol_ag_ptr, small_mol_id);
  
  // Create the minimization instructions
  MinimizeControls mincon;
  mincon.setTotalCycles(500);
  
  // Create a molecular mechanics control object based on the minimization operations
  MolecularMechanicsControls mmctrl(mincon);

  // Track energies in the systems
  ScoreCard sc(small_mol_id.size(), mincon.getTotalCycles(), 32);

  // Obtain kernel launch parameters for the workload
  KernelManager launcher(gpu, small_poly_ag);
  
  // Lay out GPU cache resources
#if 0
  const int2 vale_fe_lp = launcher.getValenceKernelDims(PrecisionModel::SINGLE, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(PrecisionModel::SINGLE, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 nonb_lp = launcher.getNonbondedKernelDims(PrecisionModel::SINGLE,
                                                       NbwuKind::TILE_GROUPS,
                                                       EvaluateForce::YES, EvaluateEnergy::YES,
                                                       ForceAccumulationMethod::SPLIT);
  const int2 redu_lp = launcher.getReductionKernelDims(PrecisionModel::SINGLE,
                                                       ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
#endif
  const int2 vale_fe_lp = { 588, 128 };
  const int2 vale_xe_lp = { 672, 128 };
  const int2 nonb_lp    = { 420, 256 };
  const int2 redu_lp    = { 336, 256 };
  
  // CHECK
  printf("Launch valence (FE) kernel on %4d blocks, %4d threads\n", vale_fe_lp.x, vale_fe_lp.y);
  printf("Launch valence (xE) kernel on %4d blocks, %4d threads\n", vale_xe_lp.x, vale_xe_lp.y);
  printf("Launch nonbonded    kernel on %4d blocks, %4d threads\n", nonb_lp.x, nonb_lp.y);
  printf("Launch reduction    kernel on %4d blocks, %4d threads\n", redu_lp.x, redu_lp.y);
  // END CHECK
  
  CacheResource valence_tb_reserve(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_reserve(nonb_lp.x, small_block_max_atoms);

  // Upload the synthesis and prime the pumps
  small_poly_ag.upload();
  small_poly_ps.upload();
  small_poly_se.upload();
  mmctrl.primeWorkUnitCounters(launcher, PrecisionModel::SINGLE, small_poly_ag);
  
  // Obtain the appropriate abstracts
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyValenceKit<float> small_poly_vk = small_poly_ag.getSinglePrecisionValenceKit(tier);
  const SyNonbondedKit<float> small_poly_nbk = small_poly_ag.getSinglePrecisionNonbondedKit(tier);
  const SeMaskSynthesisReader small_poly_ser = small_poly_se.data(tier);
  const SyRestraintKit<float, float2, float4> small_poly_rk =
    small_poly_ag.getSinglePrecisionRestraintKit(tier);
  const NbwuKind nb_work_type = small_poly_ag.getNonbondedWorkType();
  MMControlKit<float> ctrl = mmctrl.spData(tier);
  PsSynthesisWriter small_poly_psw = small_poly_ps.data(tier);
  ScoreCardWriter scw = sc.data(tier);
  CacheResourceKit vale_tbk = valence_tb_reserve.spData(tier);
  CacheResourceKit nonb_tbk = nonbond_tb_reserve.spData(tier);
  ReductionKit small_poly_redk(small_poly_ag, tier);
  ReductionBridge small_poly_rbg(small_poly_ag.getReductionWorkUnitCount());
  ConjGradSubstrate cgsbs(&small_poly_ps, &small_poly_rbg, tier);
  LineMinimization line_record(small_poly_ag.getSystemCount());
  line_record.primeMoveLengths(mmctrl.getInitialMinimizationStep());
  LinMinWriter lmw = line_record.data();
  
  // Run minimizations
  const int min_timings = timer.addCategory("Minimization of small molecules");
  timer.assignTime(0);
  for (int i = 0; i < mincon.getTotalCycles(); i++) {

    // First stage of the cycle: compute forces and obtain the conjugate gradient move.
    small_poly_ps.initializeForces(gpu, HybridTargetLevel::DEVICE);
    sc.initialize(HybridTargetLevel::DEVICE, gpu);
    launchNonbonded(nb_work_type, small_poly_nbk, small_poly_ser, &ctrl, &small_poly_psw,
                    &scw, &nonb_tbk, EvaluateForce::YES, EvaluateEnergy::YES,
                    ForceAccumulationMethod::SPLIT, nonb_lp);
    launchValence(small_poly_vk, small_poly_rk, &ctrl, &small_poly_psw,
                  &scw, &vale_tbk, EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                  ForceAccumulationMethod::SPLIT, vale_fe_lp);
    if (i == 0) {
      small_poly_ps.primeConjugateGradientCalculation(gpu, tier);
    }
    launchConjugateGradient(small_poly_redk, &cgsbs, &ctrl, redu_lp);

    // Second stage of the cycle: advance once along the line and recompute the energy.
    launchLineAdvance(PrecisionModel::SINGLE, &small_poly_psw, small_poly_redk, scw,
                      &lmw, 0, redu_lp);
    ctrl.step += 1;
    sc.initialize(HybridTargetLevel::DEVICE, gpu);
    launchNonbonded(nb_work_type, small_poly_nbk, small_poly_ser, &ctrl, &small_poly_psw,
                    &scw, &nonb_tbk, EvaluateForce::NO, EvaluateEnergy::YES,
                    ForceAccumulationMethod::SPLIT, nonb_lp);
    launchValence(small_poly_vk, small_poly_rk, &ctrl, &small_poly_psw,
                  &scw, &vale_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                  ForceAccumulationMethod::SPLIT, vale_xe_lp);

    // Third stage of the cycle: advance once more along the line and recompute the energy.
    launchLineAdvance(PrecisionModel::SINGLE, &small_poly_psw, small_poly_redk, scw,
                      &lmw, 1, redu_lp);
    ctrl.step += 1;
    sc.initialize(HybridTargetLevel::DEVICE, gpu);
    launchNonbonded(nb_work_type, small_poly_nbk, small_poly_ser, &ctrl, &small_poly_psw,
                    &scw, &nonb_tbk, EvaluateForce::NO, EvaluateEnergy::YES,
                    ForceAccumulationMethod::SPLIT, nonb_lp);
    launchValence(small_poly_vk, small_poly_rk, &ctrl, &small_poly_psw,
                  &scw, &vale_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                  ForceAccumulationMethod::SPLIT, vale_xe_lp);

    // Final stage of the cycle: advance a final time along the line, recompute the energy, fit
    // a cubic polynomial to guess the best overall advancement, and place the system there.
    launchLineAdvance(PrecisionModel::SINGLE, &small_poly_psw, small_poly_redk, scw,
                      &lmw, 2, redu_lp);
    ctrl.step += 1;
    sc.initialize(HybridTargetLevel::DEVICE, gpu);
    launchNonbonded(nb_work_type, small_poly_nbk, small_poly_ser, &ctrl, &small_poly_psw,
                    &scw, &nonb_tbk, EvaluateForce::NO, EvaluateEnergy::YES,
                    ForceAccumulationMethod::SPLIT, nonb_lp);
    launchValence(small_poly_vk, small_poly_rk, &ctrl, &small_poly_psw,
                  &scw, &vale_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                  ForceAccumulationMethod::SPLIT, vale_xe_lp);
    launchLineAdvance(PrecisionModel::SINGLE, &small_poly_psw, small_poly_redk, scw,
                      &lmw, 3, redu_lp);
    ctrl.step += 1;

    // CHECK
#if 0
    sc.initialize(HybridTargetLevel::DEVICE, gpu);
    launchNonbonded(nb_work_type, small_poly_nbk, small_poly_ser, &ctrl, &small_poly_psw,
                    &scw, &nonb_tbk, EvaluateForce::NO, EvaluateEnergy::YES,
                    ForceAccumulationMethod::SPLIT, nonb_lp);
    launchValence(small_poly_vk, small_poly_rk, &ctrl, &small_poly_psw,
                  &scw, &vale_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                  ForceAccumulationMethod::SPLIT, vale_xe_lp);
    for (int j = 15; j < 30; j += 173) {
      PhaseSpace chkj_ps = small_poly_ps.exportSystem(j, HybridTargetLevel::DEVICE);
      const std::vector<double> gpu_frc = chkj_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
      chkj_ps.initializeForces();
      ScoreCard tmp_sc(1, 1, 32);
      StaticExclusionMask chkj_se(small_poly_ag.getSystemTopologyPointer(j));
      evalNonbValeMM(&chkj_ps, &tmp_sc, small_poly_ag.getSystemTopologyPointer(j), chkj_se,
                     EvaluateForce::YES, 0);
      const std::vector<double> cpu_frc = chkj_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
      const double total_e = tmp_sc.reportTotalEnergy(0);
      printf("  %12.6lf\n", total_e);
    }
#endif
    // END CHECK
  }
  cudaDeviceSynchronize();
  timer.assignTime(min_timings);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());

  return 0;
}
