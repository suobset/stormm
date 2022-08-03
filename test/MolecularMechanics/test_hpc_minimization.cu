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
// Check that all forces on equivalent systems are equal, to within a reasonable tolerance.  It is
// assumed that systems described by the same topology are identical--no restraints, or at least
// no unique restraints, are present.
//
// Arguments:
//   poly_ps:  Compilation of coordinates and forces
//   poly_ag:  Compilation of topologies for all systems
//   err_msg:  Error to display if forces between equivalent systems do not match up
//   c_limit:  Tolerance for coordinates to disagree before being deemed inconsistent
//   f_limit:  Tolerance for forces to disagree before being deemed inconsistent
//   step_no:  Number of the minimization step
//-------------------------------------------------------------------------------------------------
bool checkConsistency(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                      const std::string &err_msg, const double c_limit, const double f_limit,
                      const int step_no) {
  int n_coord_mismatch = 0;
  int n_force_mismatch = 0;
  int n_syscrd_mismatch = 0;
  int n_sysfrc_mismatch = 0;
  PsSynthesisReader poly_psr = poly_ps.data();
  const std::vector<int> ag_indices = poly_ag.getTopologyIndices();
  std::vector<bool> covered(poly_psr.system_count, false);
  std::vector<bool> syscrd_ok(poly_psr.system_count, true);
  std::vector<bool> sysfrc_ok(poly_psr.system_count, true);
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (covered[i]) {
      continue;
    }
    covered[i] = true;
    const int iag_no = ag_indices[i];
    for (int j = i + 1; j < poly_psr.system_count; j++) {
      if (covered[j] == false && ag_indices[j] == iag_no) {
        covered[j] = true;
        int itrack = poly_psr.atom_starts[i];
        const int jtrack_lim = poly_psr.atom_starts[j] + poly_psr.atom_counts[j];
        for (int jtrack = poly_psr.atom_starts[j]; jtrack < jtrack_lim; jtrack++) {
          if (fabs(poly_psr.xcrd[itrack] - poly_psr.xcrd[jtrack]) > f_limit ||
              fabs(poly_psr.ycrd[itrack] - poly_psr.ycrd[jtrack]) > f_limit ||
              fabs(poly_psr.zcrd[itrack] - poly_psr.zcrd[jtrack]) > f_limit) {
            n_coord_mismatch++;
            if (syscrd_ok[j]) {
              n_syscrd_mismatch++;
              syscrd_ok[j] = false;
            }
          }
          if (fabs(poly_psr.xfrc[itrack] - poly_psr.xfrc[jtrack]) > f_limit ||
              fabs(poly_psr.yfrc[itrack] - poly_psr.yfrc[jtrack]) > f_limit ||
              fabs(poly_psr.zfrc[itrack] - poly_psr.zfrc[jtrack]) > f_limit) {
            n_force_mismatch++;
            if (sysfrc_ok[j]) {
              n_sysfrc_mismatch++;
              sysfrc_ok[j] = false;
            }
          }
          itrack++;
        }
      }
    }
  }
  if (n_coord_mismatch > 0) {
    rtWarn("A total of " + std::to_string(n_coord_mismatch) + " atoms' coordinates were "
           "inconsistent among the first instance of a system and subsequent replicas with the "
           "same topology.  In all, " + std::to_string(n_syscrd_mismatch) + " systems displayed "
           "errors.  Checked at: " + err_msg + ", step " + std::to_string(step_no) + ".");
  }
  if (n_force_mismatch > 0) {
    rtWarn("A total of " + std::to_string(n_force_mismatch) + " atoms' forces were inconsistent "
           "among the first instance of a system and subsequent replicas with the same topology.  "
           "In all, " + std::to_string(n_sysfrc_mismatch) + " systems displayed errors.  Checked "
           "at: " + err_msg + ", step " + std::to_string(step_no) + ".");
  }
  return (n_coord_mismatch == 0 && n_force_mismatch == 0);
}

//-------------------------------------------------------------------------------------------------
// Set forces and coordinates for all systems to be exactly the same.
//
// Arguments:
//   poly_ps:  Compilation of coordinates and forces
//   poly_ag:  Compilation of topologies for all systems
//-------------------------------------------------------------------------------------------------
void mandateEquality(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag) {
  PsSynthesisWriter poly_psr = poly_ps->data();
  const std::vector<int> ag_indices = poly_ag.getTopologyIndices();
  std::vector<bool> covered(poly_psr.system_count, false);
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (covered[i]) {
      continue;
    }
    covered[i] = true;
    const int iag_no = ag_indices[i];
    for (int j = i + 1; j < poly_psr.system_count; j++) {
      if (ag_indices[j] == iag_no) {
        covered[j] = true;
        int itrack = poly_psr.atom_starts[i];
        const int jtrack_lim = poly_psr.atom_starts[j] + poly_psr.atom_counts[j];
        for (int jtrack = poly_psr.atom_starts[j]; jtrack < jtrack_lim; jtrack++) {
          poly_psr.xfrc[jtrack] = poly_psr.xfrc[itrack];
          poly_psr.yfrc[jtrack] = poly_psr.yfrc[itrack];
          poly_psr.zfrc[jtrack] = poly_psr.zfrc[itrack];
          itrack++;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Perform energy minimization with meticulous checks to ensure that the process is consistent and
// reproducible.
//-------------------------------------------------------------------------------------------------
void metaMinimization(const std::vector<AtomGraph*> &ag_ptr_vec,
                      const std::vector<PhaseSpace> &ps_vec, const std::vector<int> &mol_id_vec,
                      const GpuDetails &gpu, const TestPriority do_tests,
                      const PrecisionModel prec = PrecisionModel::DOUBLE, const int gpos_bits = 40,
                      const int frc_bits = 40, const int maxcyc = 500,
                      const bool enforce_same_track = true, const bool check_mm = true,
                      StopWatch *timer = nullptr) {
  AtomGraphSynthesis poly_ag(ag_ptr_vec, mol_id_vec, ExceptionResponse::WARN, gpu, timer);
  StaticExclusionMaskSynthesis poly_se(poly_ag.getTopologyPointers(), mol_id_vec);
  poly_ag.loadNonbondedWorkUnits(poly_se);
  PhaseSpaceSynthesis poly_ps(ps_vec, ag_ptr_vec, mol_id_vec, gpos_bits, 24, 40, frc_bits);
  
  // Create the minimization instructions
  MinimizeControls mincon;
  mincon.setTotalCycles(maxcyc);
  
  // Create separate molecular mechanics control objects based on the minimization operations
  // for each of the ways that the valence kernel gets subdivided. 
  MolecularMechanicsControls mmctrl_fe(mincon);
  MolecularMechanicsControls mmctrl_xe(mincon);

  // Track energies in the systems
  ScoreCard sc(mol_id_vec.size(), mincon.getTotalCycles(), 32);

  // Obtain kernel launch parameters for the workload
  KernelManager launcher(gpu, poly_ag);
  
  // Lay out GPU cache resources
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        ForceAccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE);
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec,
                                                       NbwuKind::TILE_GROUPS,
                                                       EvaluateForce::YES, EvaluateEnergy::YES,
                                                       ForceAccumulationMethod::SPLIT);
  const int2 redu_lp = launcher.getReductionKernelDims(prec,
                                                       ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  CacheResource valence_fe_tb_reserve(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource valence_xe_tb_reserve(vale_xe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_reserve(nonb_lp.x, small_block_max_atoms);

  // CHECK
  printf("Launch valence (FE) kernel on %4d blocks, %4d threads per block.\n", vale_fe_lp.x,
         vale_fe_lp.y);
  printf("Launch valence ( E) kernel on %4d blocks, %4d threads per block.\n", vale_xe_lp.x,
         vale_xe_lp.y);
  printf("Launch nonbonded    kernel on %4d blocks, %4d threads per block.\n", nonb_lp.x,
         nonb_lp.y);
  printf("Launch reduction    kernel on %4d blocks, %4d threads per block.\n", redu_lp.x,
         redu_lp.y);
  // END CHECK
  
  // Upload the synthesis and prime the pumps
  poly_ag.upload();
  poly_ps.upload();
  poly_se.upload();
  mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES, prec,
                                  poly_ag);
  mmctrl_xe.primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES, prec, poly_ag);
  
  // Obtain the appropriate abstracts
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const SyValenceKit<double> d_poly_vk = poly_ag.getDoublePrecisionValenceKit(devc);
  const SyValenceKit<float>  f_poly_vk = poly_ag.getSinglePrecisionValenceKit(devc);
  const SyNonbondedKit<double> d_poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc);
  const SyNonbondedKit<float>  f_poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  const SeMaskSynthesisReader poly_ser = poly_se.data(devc);
  const SyRestraintKit<double,
                       double2, double4> d_poly_rk = poly_ag.getDoublePrecisionRestraintKit(devc);
  const SyRestraintKit<float,
                       float2, float4> f_poly_rk = poly_ag.getSinglePrecisionRestraintKit(devc);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  MMControlKit<double> d_ctrl_fe = mmctrl_fe.dpData(devc);
  MMControlKit<double> d_ctrl_xe = mmctrl_xe.dpData(devc);
  MMControlKit<float>  f_ctrl_fe = mmctrl_fe.spData(devc);
  MMControlKit<float>  f_ctrl_xe = mmctrl_xe.spData(devc);
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  ScoreCardWriter scw = sc.data(devc);
  CacheResourceKit<double> d_vale_fe_tbk = valence_fe_tb_reserve.dpData(devc);
  CacheResourceKit<double> d_vale_xe_tbk = valence_xe_tb_reserve.dpData(devc);
  CacheResourceKit<double> d_nonb_tbk = nonbond_tb_reserve.dpData(devc);
  CacheResourceKit<float> f_vale_fe_tbk = valence_fe_tb_reserve.spData(devc);
  CacheResourceKit<float> f_vale_xe_tbk = valence_xe_tb_reserve.spData(devc);
  CacheResourceKit<float> f_nonb_tbk = nonbond_tb_reserve.spData(devc);
  ReductionKit poly_redk(poly_ag, devc);
  ReductionBridge poly_rbg(poly_ag.getReductionWorkUnitCount());
  ConjGradSubstrate cgsbs(&poly_ps, &poly_rbg, devc);
  LineMinimization line_record(poly_ag.getSystemCount());
  line_record.primeMoveLengths(mmctrl_fe.getInitialMinimizationStep());
  LinMinWriter lmw = line_record.data(devc);

  // Run minimizations
  const int meta_timings = (timer == nullptr) ?
                           0 : timer->addCategory("Minimization of small molecules");
  if (timer != nullptr) {
    timer->assignTime(0);
  }
  std::vector<double> cpu_total_e(mincon.getTotalCycles(), 0.0);
  std::vector<double> gpu_total_e(mincon.getTotalCycles(), 0.0);
  std::vector<double> force_mue(mincon.getTotalCycles(), 0.0);
  int consistency_failures = 0;
  for (int i = 0; i < mincon.getTotalCycles(); i++) {
    
    // First stage of the cycle: compute forces and obtain the conjugate gradient move.
    poly_ps.initializeForces(gpu, devc);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_fe, &poly_psw, &scw, &d_nonb_tbk,
                      EvaluateForce::YES, EvaluateEnergy::YES, nonb_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_fe, &poly_psw, &scw, &d_vale_fe_tbk,
                    EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_fe_lp);
      d_ctrl_fe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_fe, &poly_psw, &scw, &f_nonb_tbk,
                      EvaluateForce::YES, EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                      nonb_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_fe, &poly_psw, &scw, &f_vale_fe_tbk,
                    EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    ForceAccumulationMethod::SPLIT, vale_fe_lp);
      f_ctrl_fe.step += 1;
      break;
    }

    // Check the forces computed for a couple of systems.  This is somewhat redundant, but serves
    // as a sanity check in case other aspects of the energy minimization show problems.
    const int jlim = (check_mm) ? (3 * i) + 1 : 3 * i;
    for (int j = 3 * i; j < jlim; j++) {
      const int jmod = j % poly_ag.getSystemCount();
      PhaseSpace chkj_ps = poly_ps.exportSystem(jmod, devc);
      const std::vector<double> gpu_frc = chkj_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
      chkj_ps.initializeForces();
      ScoreCard tmp_sc(1, 1, 32);
      StaticExclusionMask chkj_se(poly_ag.getSystemTopologyPointer(jmod));
      evalNonbValeMM(&chkj_ps, &tmp_sc, poly_ag.getSystemTopologyPointer(jmod), chkj_se,
                     EvaluateForce::YES, 0);
      const std::vector<double> cpu_frc = chkj_ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
      cpu_total_e[i] = tmp_sc.reportTotalEnergy(0);
      gpu_total_e[i] = sc.reportTotalEnergy(jmod, devc);
      force_mue[i] = meanUnsignedError(cpu_frc, gpu_frc);
    }
    
    // Download and check the forces for each system to verify consistency.  If the forces are
    // consistent enough, set them to be exactly consistent, and do the same with the coordinates,
    // to avoid miniscule roundoff errors that could otherwis creep in over hundreds of steps.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Force computation", 5.0e-7, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }

    // Perform the conjugate gradient transformation
    if (i == 0) {
      poly_ps.primeConjugateGradientCalculation(gpu, devc);
    }
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchConjugateGradient(poly_redk, &cgsbs, &d_ctrl_fe, redu_lp);
      break;
    case PrecisionModel::SINGLE:
      launchConjugateGradient(poly_redk, &cgsbs, &f_ctrl_fe, redu_lp);
      break;
    }
    
    // Download and check the conjugate gradient transformation.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Conjugate gradient transformation",
                           5.0e-7, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Second stage of the cycle: advance once along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 0, redu_lp);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &scw, &d_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, &scw, &d_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &scw, &f_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                      nonb_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, &scw, &f_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    ForceAccumulationMethod::SPLIT, vale_xe_lp);
      f_ctrl_xe.step += 1;
      break;
    }
    
    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance I", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Third stage of the cycle: advance once more along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 1, redu_lp);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &scw, &d_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, &scw, &d_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &scw, &f_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                      nonb_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, &scw, &f_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    ForceAccumulationMethod::SPLIT, vale_xe_lp);
      f_ctrl_xe.step += 1;
      break;
    }
    
    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance II", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Final stage of the cycle: advance a final time along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 2, redu_lp);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &scw, &d_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, &scw, &d_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &scw, &f_nonb_tbk,
                      EvaluateForce::NO, EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                      nonb_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, &scw, &f_vale_xe_tbk,
                    EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    ForceAccumulationMethod::SPLIT, vale_xe_lp);
      f_ctrl_xe.step += 1;
      break;
    }

    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance III", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      } 
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Fit a cubic polynomial to guess the best overall advancement, and place the system there.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 3, redu_lp);

    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance IV", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
  }
  cudaDeviceSynchronize();
  if (timer != nullptr) {
    timer->assignTime(meta_timings);
  }
  if (enforce_same_track) {
    check(consistency_failures == 0, "Consistency failures occurred during the minimizations.  "
          "The results in identical systems diverged after approximately " +
          std::to_string((consistency_failures + 5) / 6) + " steps.", do_tests);
  }
  if (check_mm) {
    const int ave_atom_count = poly_ag.getAtomCount() / poly_ag.getSystemCount();
    const double etol = (prec == PrecisionModel::DOUBLE) ? 5.0e-7 : 5.0e-4;
    check(cpu_total_e, RelationalOperator::EQUAL, Approx(gpu_total_e).margin(etol), "Energies of "
          "relaxed structures did not agree.  Average structure size: " +
          std::to_string(ave_atom_count) + ".", do_tests);
    check(force_mue, RelationalOperator::EQUAL,
          Approx(std::vector<double>(mincon.getTotalCycles(), 0.0)).margin(1.0e-6), "Snapshots of "
          "forces taken during energy minimization on the GPU do not agree with their CPU-derived "
          "counterparts.", do_tests);
  }
  
  // Verify that the energies for all systems meet the expected values
  const std::vector<double> final_e = sc.reportTotalEnergies(devc);
  std::vector<double> sampled_e;
  const int nsys = poly_ps.getSystemCount();
  for (int i = 0; i < nsys; i++) {
    
  }
       
  // CHECK
  printf("final_e = [\n");
  int j = 0;
  for (size_t i = 0; i < mol_id_vec.size(); i++) {
    printf("  %12.8lf", final_e[i]);
    j++;
    if (j == 16) {
      printf("\n");
      j = 0;
    }
  }
  if (j > 0) {
    printf("\n");
  }
  printf("];\n");
  // END CHECK
}

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
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<int> small_mol_id = { 0, 1, 2, 3, 0, 1, 2, 3, 0, 3, 1, 2, 2, 1, 3, 0 };
  small_mol_id.resize(8192);
  for (int i = 0; i < small_mol_id.size() / 16; i++) {
    for (int j = 0; j < 16; j++) {
      small_mol_id[(16 * i) + j] = small_mol_id[j];
    }
  }
#if 0
  metaMinimization(small_mol_ag_ptr, small_mol_ps, small_mol_id, gpu, do_tests,
                   PrecisionModel::DOUBLE);
#endif
  metaMinimization(small_mol_ag_ptr, small_mol_ps, small_mol_id, gpu, do_tests,
                   PrecisionModel::SINGLE, 28, 24, 500, false, false, &timer);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  
  return 0;
}
