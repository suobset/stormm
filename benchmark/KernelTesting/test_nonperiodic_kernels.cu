// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include <string>
#include <vector>
#include "../../src/Accelerator/hpc_config.cuh"
#include "../../src/Accelerator/kernel_manager.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_valence_potential.cuh"
#include "../../src/Potential/hpc_nonbonded_potential.cuh"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"


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
// Get a SystemCache object containing all topologies and coordinates in a pair of directories.
//
// Arguments:
//   topol_path:  A series of strings that will be joined into the topology directory name
//   coord_path:  A series of strings that will be joined into the coordinate directory name
//   oe:          Contains critical shell variables such as the $OMNI source path where the
//                named directories are expected to reside
//-------------------------------------------------------------------------------------------------
SystemCache directorySweep(const std::vector<std::string> &topol_path,
                           const std::vector<std::string> &coord_path, const TestEnvironment &oe) {
  
  // Collect coordinates and topologies
  const char osc = osSeparator();
  std::string buffer("&files\n  -p ");
  buffer += oe.getOmniSourcePath() + osc + "benchmark";
  for (size_t i = 0; i < topol_path.size(); i++) {
    buffer += osc + topol_path[i];
  }
  buffer += "\n  -c ";
  buffer += oe.getOmniSourcePath() + osc + "benchmark";
  for (size_t i = 0; i < topol_path.size(); i++) {
    buffer += osc + coord_path[i];
  }
  buffer += "\n&end\n";
  const TextFile tf(buffer, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  return SystemCache(fcon, ExceptionResponse::SILENT, MapRotatableGroups::NO);
}

//-------------------------------------------------------------------------------------------------
// Replicate a single topology and coordinate system many times, then run kernels to obtain
// timings.
//
// Arguments:
//   ag:       The topology to replicate
//   ps:       The coordinates to replicate
//   nrep:     The number of replicas to make
//   mmctrl:   Step counter and progress counters for all work units
//   gpu:      Details of the GPU to use in calculations
//   timer:    Object to record the timings
//   desc:     Description of the system to run (this will factor into the timings section names)
//-------------------------------------------------------------------------------------------------
void replicaProcessing(AtomGraph *ag, const PhaseSpace &ps, const int nrep,
                       MolecularMechanicsControls *mmctrl, const GpuDetails &gpu,
                       StopWatch *timer, const PrecisionModel prec, const EvaluateForce eval_frc,
                       const EvaluateEnergy eval_nrg, const ForceAccumulationMethod acc_meth,
                       const VwuGoal purpose) {
  std::vector<AtomGraph*> ag_vec(1, ag);
  std::vector<PhaseSpace> ps_vec(1, ps);
  std::vector<int> ag_idx(nrep, 0);
  AtomGraphSynthesis poly_ag(ag_vec, ag_idx, ExceptionResponse::SILENT,
                             maximum_valence_work_unit_atoms);
  PhaseSpaceSynthesis poly_ps(ps_vec, ag_vec, ag_idx);
  StaticExclusionMaskSynthesis poly_se(ag_vec, ag_idx);
  SeMaskSynthesisReader poly_ser = poly_se.data();
  poly_ag.loadNonbondedWorkUnits(poly_se);
  KernelManager launcher(gpu, poly_ag);
  ScoreCard sc(nrep, 1, 32);
  const int2 valence_lp = launcher.getValenceKernelDims(prec, eval_frc, eval_nrg, acc_meth,
                                                        purpose);
  const int2 nonbond_lp = launcher.getNonbondedKernelDims(prec, NbwuKind::TILE_GROUPS, eval_frc,
                                                          eval_nrg, acc_meth);
  CacheResource valence_tb_space(valence_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_space(nonbond_lp.x, small_block_max_atoms);
  mmctrl->primeWorkUnitCounters(launcher, prec, poly_ag);
  
  // Upload the critical components
  poly_ag.upload();
  poly_se.upload();
  poly_ps.upload();

  // Some common variables for either branch
  const std::string valk_name = valenceKernelKey(prec, eval_frc, eval_nrg, acc_meth, purpose);
  const std::string nnbk_name = nonbondedKernelKey(prec, NbwuKind::TILE_GROUPS, eval_frc,
                                                   eval_nrg, acc_meth);
  const int sys_val_time = timer->addCategory(getBaseName(ag->getFileName()) + " on " +
                                              valk_name + " (" + std::to_string(nrep) + ")");
  const int sys_nb_time = timer->addCategory(getBaseName(ag->getFileName()) + " on " +
                                             nnbk_name + " (" + std::to_string(nrep) + ")");
  
  // Obtain abstracts outside the inner loop, in case this is a significant contributor to the
  // run time.  Forces will only be initialized once, and thereafter calculated repeatedly to test
  // only the run time of the one kernel.
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<double> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<double,
                           double2,
                           double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      PsSynthesisWriter poly_psw = poly_ps.data(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<double> gmem_rval = valence_tb_space.dpData(devc_tier);
      CacheResourceKit<double> gmem_rnnb = nonbond_tb_space.dpData(devc_tier);
      poly_ps.initializeForces(gpu, devc_tier);
      timer->assignTime(0);

      // Test the valence kernel
      for (int i = 0; i < 1000; i++) {
        ctrl.step += 1;
        launchValenceDp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_rval, eval_frc, eval_nrg,
                        purpose, launcher);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_val_time);

      // Test the non-bonded kernel
      poly_ps.initializeForces(gpu, devc_tier);
      timer->assignTime(0);
      for (int i = 0; i < 1000; i++) {
        ctrl.step += 1;
        launchNonbondedTileGroupsDp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_rnnb,
                                    eval_frc, eval_nrg, launcher);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_nb_time);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<float> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<float,
                           float2,
                           float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      PsSynthesisWriter poly_psw = poly_ps.data(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<float> gmem_rval = valence_tb_space.spData(devc_tier);
      CacheResourceKit<float> gmem_rnnb = nonbond_tb_space.spData(devc_tier);
      poly_ps.initializeForces(gpu, devc_tier);      
      timer->assignTime(0);

      // Test the valence kernel
      for (int i = 0; i < 1000; i++) {
        ctrl.step += 1;
        launchValenceSp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_rval, eval_frc, eval_nrg,
                        purpose, acc_meth, launcher);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_val_time);

      // Test the non-bonded kernel
      poly_ps.initializeForces(gpu, devc_tier);      
      timer->assignTime(0);
      for (int i = 0; i < 1000; i++) {
        ctrl.step += 1;
        launchNonbondedTileGroupsSp(poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_rnnb,
                                    eval_frc, eval_nrg, acc_meth, launcher);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_nb_time);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);

  // Configure the relevant kernels for this executable.
  valenceKernelSetup();

  // Read dipeptides
  const std::vector<std::string> dipeptide_topols = { "Topologies", "Dipeptides",
                                                      ".*_ff1.*SB.top" };
  const std::vector<std::string> dipeptide_coords = { "Coordinates", "Dipeptides",
                                                      ".*_ff1.*SB.inpcrd"};
  SystemCache dipeptide_sc = directorySweep(dipeptide_topols, dipeptide_coords, oe);
  MolecularMechanicsControls mmctrl;
  
  // Loop over the dipeptides one at a time, make syntheses of each of them individually, and
  // test kernels.
  const int ndipeptides = dipeptide_sc.getSystemCount();
  const std::vector<int> batch_multiplier = { 1, 3, 5, 10, 20, 40 };
  for (int i = 0; i < ndipeptides; i++) {
    AtomGraph *iag_ptr = dipeptide_sc.getSystemTopologyPointer(i);
    for (int j = 0; j < 6; j++) {
      const int ncopy = gpu.getSMPCount() * batch_multiplier[j];
      replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                        ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                        ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                        ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                        ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE);
      
      // Only do double-precision calculations for low replica numbers--this can be strenuous on
      // many architectures, particularly in the non-bonded kernel.
      if (ncopy < 10) {
        replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                          &timer, PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                          ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
        replicaProcessing(iag_ptr, dipeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                          &timer, PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                          ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      }
    }
  }

  // Read dipeptides
  const std::vector<std::string> tripeptide_topols = { "Topologies", "Tripeptides",
                                                       ".*_ff1.*SB.top" };
  const std::vector<std::string> tripeptide_coords = { "Coordinates", "Tripeptides",
                                                       ".*_ff1.*SB.inpcrd"};
  SystemCache tripeptide_sc = directorySweep(tripeptide_topols, tripeptide_coords, oe);
  
  // Loop over the tripeptides one at a time, make syntheses of each of them individually, and
  // test kernels.
  const int ntripeptides = tripeptide_sc.getSystemCount();
  for (int i = 0; i < ntripeptides; i++) {
    AtomGraph *iag_ptr = tripeptide_sc.getSystemTopologyPointer(i);
    for (int j = 0; j < 6; j++) {
      const int ncopy = gpu.getSMPCount() * batch_multiplier[j];
      replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                        ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                        ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                        ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE);
      replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                        &timer, PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                        ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE);
      if (ncopy < 10) {
        replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                          &timer, PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                          ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
        replicaProcessing(iag_ptr, tripeptide_sc.getCoordinateReference(i), ncopy, &mmctrl, gpu,
                          &timer, PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                          ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE);
      }
    }
  }

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
}
