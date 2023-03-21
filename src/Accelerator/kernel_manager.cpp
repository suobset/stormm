// -*-c++-*-
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#ifdef STORMM_USE_HPC
#include "Math/hpc_reduction.h"
#include "Math/reduction_workunit.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Structure/hpc_rmsd.h"
#include "Structure/hpc_virtual_site_handling.h"
#endif
#include "kernel_manager.h"

namespace stormm {
namespace card {

using constants::ExceptionResponse;
#ifdef STORMM_USE_HPC
using energy::queryValenceKernelRequirements;
using energy::queryNonbondedKernelRequirements;
using energy::queryBornRadiiKernelRequirements;
using energy::queryBornDerivativeKernelRequirements;
using stmath::queryReductionKernelRequirements;
using stmath::optReductionKernelSubdivision;
using structure::queryRMSDKernelRequirements;
using structure::queryVirtualSiteKernelRequirements;
using synthesis::optValenceKernelSubdivision;
using synthesis::optVirtualSiteKernelSubdivision;
#endif
using stmath::findBin;
using parse::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat() :
  block_size_limit{1}, shared_usage{0}, block_dimension{1}, grid_dimension{1}, register_usage{0},
  kernel_name{std::string("")}
{}
  
//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                           const int register_usage_in, const int shared_usage_in,
                           const int block_subdivision, const GpuDetails &gpu,
                           const std::string &kernel_name_in) :
    block_size_limit{lb_max_threads_per_block},
    shared_usage{shared_usage_in},
    block_dimension{(block_size_limit / block_subdivision / warp_size_int) * warp_size_int},
    grid_dimension{block_subdivision * lb_min_blocks_per_smp * gpu.getSMPCount()},
    register_usage{register_usage_in},
    kernel_name{kernel_name_in}
{
  // Refine the register usage (this is unreliable with current cudart function calls, and should
  // not be trusted even after this step).
  const std::vector<int> register_break_points = {  0, 40, 48, 56, 64, 72, 80, 128, 256 };
  const std::vector<int> register_warp_counts  = { 48, 40, 36, 32, 28, 24, 16,   8 };
  const int register_bin = findBin(register_break_points, register_usage, ExceptionResponse::WARN);
  register_usage = register_warp_counts[register_bin];
}

//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                           const int register_usage_in, const int shared_usage_in,
                           const GpuDetails &gpu, const std::string &kernel_name_in) :
    KernelFormat(lb_max_threads_per_block, lb_min_blocks_per_smp, register_usage_in,
                 shared_usage_in, 1, gpu, kernel_name_in)
{}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const cudaFuncAttributes &attr, const int lb_min_blocks_per_smp,
                           const int block_subdivision, const GpuDetails &gpu,
                           const std::string &kernel_name_in) :
    KernelFormat(attr.maxThreadsPerBlock, lb_min_blocks_per_smp, attr.numRegs,
                 attr.sharedSizeBytes, block_subdivision, gpu, kernel_name_in)
{}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
int2 KernelFormat::getLaunchParameters() const {
  return { grid_dimension, block_dimension };
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getRegisterUsage() const {
  return register_usage;
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getBlockSizeLimit() const {
  return block_size_limit;
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getSharedMemoryRequirement() const {
  return shared_usage;
}

//-------------------------------------------------------------------------------------------------
const std::string& KernelFormat::getKernelName() const {
  return kernel_name;
}

//-------------------------------------------------------------------------------------------------
KernelManager::KernelManager(const GpuDetails &gpu_in, const AtomGraphSynthesis &poly_ag) :
    gpu{gpu_in},
    valence_block_multiplier{valenceBlockMultiplier()},
    nonbond_block_multiplier_dp{nonbondedBlockMultiplier(gpu_in, poly_ag.getUnitCellType(),
                                                         PrecisionModel::DOUBLE,
                                                         poly_ag.getImplicitSolventModel())},
    nonbond_block_multiplier_sp{nonbondedBlockMultiplier(gpu_in, poly_ag.getUnitCellType(),
                                                         PrecisionModel::SINGLE,
                                                         poly_ag.getImplicitSolventModel())},
    gbradii_block_multiplier_dp{gbRadiiBlockMultiplier(gpu_in, PrecisionModel::DOUBLE)},
    gbradii_block_multiplier_sp{gbRadiiBlockMultiplier(gpu_in, PrecisionModel::SINGLE)},
    gbderiv_block_multiplier_dp{gbDerivativeBlockMultiplier(gpu_in, PrecisionModel::DOUBLE)},
    gbderiv_block_multiplier_sp{gbDerivativeBlockMultiplier(gpu_in, PrecisionModel::SINGLE)},
    reduction_block_multiplier{reductionBlockMultiplier()},
    virtual_site_block_multiplier_dp{virtualSiteBlockMultiplier(PrecisionModel::DOUBLE)},
    virtual_site_block_multiplier_sp{virtualSiteBlockMultiplier(PrecisionModel::SINGLE)},
    k_dictionary{}
{
#ifdef STORMM_USE_HPC
  // Valence kernel entries
  const int valence_d_div = optValenceKernelSubdivision(poly_ag.getSystemAtomCounts(),
                                                        PrecisionModel::DOUBLE,
                                                        EvaluateForce::YES, gpu_in);
  const int valence_fxe_div = optValenceKernelSubdivision(poly_ag.getSystemAtomCounts(),
                                                          PrecisionModel::SINGLE,
                                                          EvaluateForce::NO, gpu_in);
  const int valence_ffx_div = optValenceKernelSubdivision(poly_ag.getSystemAtomCounts(),
                                                          PrecisionModel::SINGLE,
                                                          EvaluateForce::YES, gpu_in);
  const int valence_ffe_div = optValenceKernelSubdivision(poly_ag.getSystemAtomCounts(),
                                                          PrecisionModel::SINGLE,
                                                          EvaluateForce::YES, gpu_in);
  const std::vector<ClashResponse> clash_policy = { ClashResponse::NONE, ClashResponse::FORGIVE };
  const std::vector<std::string> clash_ext = { "", "NonClash" };
  for (int i = 0; i < 2; i++) {
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::NO, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_d_div, "kdsValenceEnergyAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_d_div, "kdsValenceAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_d_div, "kdsValenceForceAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_d_div, "kdsValenceEnergyAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_d_div, "kdsValenceForceEnergyAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::NO, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_fxe_div, "kfsValenceEnergyAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_ffx_div, "kfsValenceAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_ffx_div, "kfValenceAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_ffx_div, "kfsValenceForceAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                         AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_ffx_div, "kfValenceForceAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_ffe_div, "kfValenceEnergyAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES, clash_policy[i],
                         valence_ffe_div, "kfsValenceEnergyAtomUpdate" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_ffe_div, "kfValenceForceEnergyAccumulation" + clash_ext[i]);
    catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                         AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE, clash_policy[i],
                         valence_ffe_div, "kfsValenceForceEnergyAccumulation" + clash_ext[i]);
  
    // Non-bonded kernel entries
    const std::vector<ImplicitSolventModel> is_models = { ImplicitSolventModel::NONE,
                                                          ImplicitSolventModel::HCT_GB,
                                                          ImplicitSolventModel::NECK_GB };
    const std::vector<std::string> is_model_names = { "Vacuum", "GB", "GBNeck" };
    for (int j = 0; j < 3; j++) {
      switch (poly_ag.getUnitCellType()) {
      case UnitCellType::NONE:
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgd" + is_model_names[j] + "Energy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgds" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgds" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgf" + is_model_names[j] + "Energy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::NO, AccumulationMethod::WHOLE, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "Force" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::SPLIT, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);
        catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                               EvaluateEnergy::YES, AccumulationMethod::WHOLE, is_models[j],
                               clash_policy[i],
                               "ktgfs" + is_model_names[j] + "ForceEnergy" + clash_ext[i]);

        // Generalized Born radii and radial derivative kernel entries
        if (is_models[j] != ImplicitSolventModel::NONE && clash_policy[i] == ClashResponse::NONE) {
          catalogBornRadiiKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::SPLIT, is_models[j],
                                 "ktgdsCalculate" + is_model_names[j] + "Radii");
          catalogBornRadiiKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::SPLIT, is_models[j],
                                 "ktgfsCalculate" + is_model_names[j] + "Radii");
          catalogBornRadiiKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                 AccumulationMethod::WHOLE, is_models[j],
                                 "ktgfCalculate" + is_model_names[j] + "Radii");
          catalogBornDerivativeKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::SPLIT, is_models[j],
                                      "ktgdsCalculate" + is_model_names[j] + "Derivatives");
          catalogBornDerivativeKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::SPLIT, is_models[j],
                                      "ktgfsCalculate" + is_model_names[j] + "Derivatives");
          catalogBornDerivativeKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS,
                                      AccumulationMethod::WHOLE, is_models[j],
                                      "ktgfCalculate" + is_model_names[j] + "Derivatives");
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
  }
  
  // Reduction kernel entries
  const int reduction_div = optReductionKernelSubdivision(poly_ag.getSystemAtomCounts(), gpu_in);
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, reduction_div, "kdgtConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, reduction_div, "kdscConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, reduction_div, "kdrsConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, reduction_div, "kdrdConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, reduction_div, "kfgtConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, reduction_div, "kfscConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, reduction_div, "kfrsConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, reduction_div, "kfrdConjGrad");

  // Virtual site kernel entries
  const int vsite_div = optVirtualSiteKernelSubdivision(poly_ag.getValenceWorkUnitAbstracts(),
                                                        poly_ag.getValenceWorkUnitCount());
  catalogVirtualSiteKernel(PrecisionModel::DOUBLE, VirtualSiteActivity::PLACEMENT, vsite_div,
                           "kdPlaceVirtualSites");
  catalogVirtualSiteKernel(PrecisionModel::SINGLE, VirtualSiteActivity::PLACEMENT, vsite_div,
                           "kfPlaceVirtualSites");
  catalogVirtualSiteKernel(PrecisionModel::DOUBLE, VirtualSiteActivity::TRANSMIT_FORCES, vsite_div,
                           "kdTransmitVSiteForces");
  catalogVirtualSiteKernel(PrecisionModel::SINGLE, VirtualSiteActivity::TRANSMIT_FORCES, vsite_div,
                           "kfTransmitVSiteForces");

  // RMSD kernel entries
  catalogRMSDKernel(PrecisionModel::DOUBLE, RMSDTask::REFERENCE, "kdComputeRMSDToReference");
  catalogRMSDKernel(PrecisionModel::SINGLE, RMSDTask::REFERENCE, "kfComputeRMSDToReference");
  catalogRMSDKernel(PrecisionModel::DOUBLE, RMSDTask::MATRIX, "kdComputeRMSDMatrix");
  catalogRMSDKernel(PrecisionModel::SINGLE, RMSDTask::MATRIX, "kfComputeRMSDMatrix");
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogValenceKernel(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const VwuGoal purpose,
                                         const ClashResponse collision_handling,
                                         const int subdivision, const std::string &kernel_name) {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose,
                                             collision_handling);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogValenceKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryValenceKernelRequirements(prec, eval_force, eval_nrg,
                                                                 acc_meth, purpose,
                                                                 collision_handling);
  k_dictionary[k_key] = KernelFormat(attr, valence_block_multiplier, subdivision, gpu,
                                     kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogNonbondedKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const ClashResponse collision_handling,
                                           const std::string &kernel_name) {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth, igb,
                                               collision_handling);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogNonbondedKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryNonbondedKernelRequirements(prec, kind, eval_force,
                                                                   eval_nrg, acc_meth, igb,
                                                                   collision_handling);
  int nonbond_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    nonbond_block_multiplier = nonbond_block_multiplier_dp;
     break;
  case PrecisionModel::SINGLE:
    nonbond_block_multiplier = nonbond_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, nonbond_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogBornRadiiKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const std::string &kernel_name) {
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth, igb);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Born radii kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogBornRadiiKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryBornRadiiKernelRequirements(prec, kind, acc_meth, igb);
  int gbradii_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gbradii_block_multiplier = gbradii_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    gbradii_block_multiplier = gbradii_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, gbradii_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogBornDerivativeKernel(const PrecisionModel prec, const NbwuKind kind,
                                                const AccumulationMethod acc_meth,
                                                const ImplicitSolventModel igb,
                                                const std::string &kernel_name) {
  const std::string k_key = bornDerivativeKernelKey(prec, kind, acc_meth, igb);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Born radii derivative kernel identifier " + k_key + " already exists in the kernel "
          "map.", "KernelManager", "catalogBornRadiiKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryBornDerivativeKernelRequirements(prec, kind, acc_meth, igb);
  int gbderiv_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gbderiv_block_multiplier = gbderiv_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    gbderiv_block_multiplier = gbderiv_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, gbderiv_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogReductionKernel(const PrecisionModel prec, const ReductionGoal purpose,
                                           const ReductionStage process, const int subdivision,
                                           const std::string &kernel_name) {
  const std::string k_key = reductionKernelKey(prec, purpose, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Reduction kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogReductionKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryReductionKernelRequirements(prec, purpose, process);
  k_dictionary[k_key] = KernelFormat(attr, reduction_block_multiplier, subdivision, gpu,
                                     kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogVirtualSiteKernel(const PrecisionModel prec,
                                             const VirtualSiteActivity purpose,
                                             const int subdivision,
                                             const std::string &kernel_name) {
  const std::string k_key = virtualSiteKernelKey(prec, purpose);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Virtual site handling kernel identifier " + k_key + " already exists in the kernel "
          "map.", "KernelManager", "catalogVirtualSiteKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryVirtualSiteKernelRequirements(prec, purpose);
  int virtual_site_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    virtual_site_block_multiplier = virtual_site_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    virtual_site_block_multiplier = virtual_site_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, virtual_site_block_multiplier, subdivision, gpu,
                                     kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogRMSDKernel(const PrecisionModel prec, const RMSDTask order,
                                      const std::string &kernel_name) {
  const std::string k_key = rmsdKernelKey(prec, order);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("RMSD calculation kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogRMSDKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryRMSDKernelRequirements(prec, order);
  int rmsd_block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    rmsd_block_multiplier = rmsd_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    rmsd_block_multiplier = rmsd_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, rmsd_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}
                                             
//-------------------------------------------------------------------------------------------------
int2 KernelManager::getValenceKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const VwuGoal purpose,
                                         const ClashResponse collision_handling) const {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose,
                                             collision_handling);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getValenceKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getNonbondedKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb,
                                           const ClashResponse collision_handling) const {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth, igb,
                                               collision_handling);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getNonbondedKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getBornRadiiKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth,
                                           const ImplicitSolventModel igb) const {
  if (igb == ImplicitSolventModel::NONE) {
    return { 0, 0 };
  }
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth, igb);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Born radii computation kernel identifier " + k_key + " was not found in the kernel "
          "map.", "KernelManager", "getBornRadiiKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getBornDerivativeKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                                const AccumulationMethod acc_meth,
                                                const ImplicitSolventModel igb) const {
  if (igb == ImplicitSolventModel::NONE) {
    return { 0, 0 };
  }
  const std::string k_key = bornDerivativeKernelKey(prec, kind, acc_meth, igb);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Born radii derivative computation kernel identifier " + k_key + " was not found in "
          "the kernel map.", "KernelManager", "getBornRadiiKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getReductionKernelDims(const PrecisionModel prec, const ReductionGoal purpose,
                                           const ReductionStage process) const {
  const std::string k_key = reductionKernelKey(prec, purpose, process);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Reduction kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getReductionKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getVirtualSiteKernelDims(const PrecisionModel prec,
                                             const VirtualSiteActivity purpose) const {
  const std::string k_key = virtualSiteKernelKey(prec, purpose);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Virtual site handling kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getVirtualSiteKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getRMSDKernelDims(const PrecisionModel prec, const RMSDTask order) const {
  const std::string k_key = rmsdKernelKey(prec, order);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("RMSD calculation kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getRMSDKerenlDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
const GpuDetails& KernelManager::getGpu() const {
  return gpu;
}

//-------------------------------------------------------------------------------------------------
void KernelManager::printLaunchParameters(const std::string &k_key) const {
  if (k_key.size() == 0) {
    return;
  }
  if (strcmpCased(k_key, "all", CaseSensitivity::NO)) {
    for (auto it = k_dictionary.begin(); it != k_dictionary.end(); it++) {
      printLaunchParameters(it->first);
    }
    return;
  }
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("No kernel with identifier " + k_key + " is known.", "KernelManager",
          "printLaunchParameters");
  }
  const int2 lp = k_dictionary.at(k_key).getLaunchParameters();
  const int mtpb = k_dictionary.at(k_key).getBlockSizeLimit();
  const int nreg = k_dictionary.at(k_key).getRegisterUsage();
  printf("  %12.12s :: %4d blocks, %4d threads, %3d registers (%4d SMP, %4d max threads)\n",
         k_key.c_str(), lp.x, lp.y, nreg, gpu.getSMPCount(), mtpb);
}

//-------------------------------------------------------------------------------------------------
int valenceBlockMultiplier() {
#ifdef STORMM_USE_HPC
  return 2;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int nonbondedBlockMultiplier(const GpuDetails &gpu, const UnitCellType unit_cell,
                             const PrecisionModel prec, const ImplicitSolventModel igb) {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  switch (unit_cell) {
  case UnitCellType::NONE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
      }
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      return 2;
    case PrecisionModel::SINGLE:
      return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 2 : 3;
    }
  }
#  else
  // Other vendors are not known to make GPUs that have special requirements
  switch (unit_cell) {
  case UnitCellType::NONE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return 5;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        return 3;
      case PrecisionModel::SINGLE:
        return 5;
      }
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    return 3;
  }
#  endif
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
int gbRadiiBlockMultiplier(const GpuDetails &gpu, const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  }
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
int gbDerivativeBlockMultiplier(const GpuDetails &gpu, const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  }
  __builtin_unreachable();
#else
  return 1;
#endif  
}

//-------------------------------------------------------------------------------------------------
int reductionBlockMultiplier() {
#ifdef STORMM_USE_HPC
  return 4;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int virtualSiteBlockMultiplier(const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 3;
  case PrecisionModel::SINGLE:
    return 4;
  }
  __builtin_unreachable();
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int rmsdBlockMultiplier(const PrecisionModel prec) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return 4;
  case PrecisionModel::SINGLE:
    return 5;
  }
  __builtin_unreachable();
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
std::string valenceKernelKey(const PrecisionModel prec, const EvaluateForce eval_force,
                             const EvaluateEnergy eval_nrg, const AccumulationMethod acc_meth,
                             const VwuGoal purpose, const ClashResponse collision_handling) {
  std::string k_key("vale_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (eval_force) {
  case EvaluateForce::YES:
    k_key += "f";
    break;
  case EvaluateForce::NO:
    k_key += "x";
    break;
  }
  switch (eval_nrg) {
  case EvaluateEnergy::YES:
    k_key += "e";
    break;
  case EvaluateEnergy::NO:
    k_key += "x";
    break;
  }
  if (eval_force == EvaluateForce::YES) {
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case AccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case AccumulationMethod::AUTOMATIC:
      break;
    }
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
      k_key += "a";
      break;
    case VwuGoal::MOVE_PARTICLES:
      k_key += "m";
      break;
    }
  }
  switch (collision_handling) {
  case ClashResponse::NONE:
    k_key += "_cl";
    break;
  case ClashResponse::FORGIVE:
    k_key += "_nc";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string nonbondedKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const EvaluateForce eval_force, const EvaluateEnergy eval_nrg,
                               const AccumulationMethod acc_meth, const ImplicitSolventModel igb,
                               const ClashResponse collision_handling) {
  std::string k_key("nonb_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    k_key += "tg";
    break;
  case NbwuKind::SUPERTILES:
    k_key += "st";
    break;
  case NbwuKind::HONEYCOMB:
    k_key += "hc";
    break;
  case NbwuKind::UNKNOWN:
    break;
  }
  switch (igb) {
  case ImplicitSolventModel::NONE:
    k_key += "_vac_";
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    k_key += "_gbs_";
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    k_key += "_gbn_";
    break;
  }
  switch (eval_force) {
  case EvaluateForce::YES:
    k_key += "f";
    break;
  case EvaluateForce::NO:
    k_key += "x";
    break;
  }
  switch (eval_nrg) {
  case EvaluateEnergy::YES:
    k_key += "e";
    break;
  case EvaluateEnergy::NO:
    k_key += "x";
    break;
  }
  if (eval_force == EvaluateForce::YES) {
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case AccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case AccumulationMethod::AUTOMATIC:
      break;
    }
  }
  switch (collision_handling) {
  case ClashResponse::NONE:
    k_key += "_cl";
    break;
  case ClashResponse::FORGIVE:
    k_key += "_nc";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string appendBornKernelKey(const PrecisionModel prec, const NbwuKind kind,
                                const AccumulationMethod acc_meth,
                                const ImplicitSolventModel igb) {
  std::string app_key("");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    app_key += "d";
    break;
  case PrecisionModel::SINGLE:
    app_key += "f";
    break;
  }
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    app_key += "tg";
    break;
  case NbwuKind::SUPERTILES:
    app_key += "st";
    break;
  case NbwuKind::HONEYCOMB:
    app_key += "hc";
    break;
  case NbwuKind::UNKNOWN:
    break;
  }
  switch (acc_meth) {
  case AccumulationMethod::SPLIT:
    app_key += "s";
    break;
  case AccumulationMethod::WHOLE:
    app_key += "w";
    break;
  case AccumulationMethod::AUTOMATIC:
    break;
  }
  switch (igb) {
  case ImplicitSolventModel::NONE:
    app_key += "_vac";
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    app_key += "_gbs";
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    app_key += "_gbn";
    break;
  }
  return app_key;
}

//-------------------------------------------------------------------------------------------------
std::string bornRadiiKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const AccumulationMethod acc_meth, const ImplicitSolventModel igb) {
  std::string k_key("gbrd_");
  return k_key + appendBornKernelKey(prec, kind, acc_meth, igb);
}

//-------------------------------------------------------------------------------------------------
std::string bornDerivativeKernelKey(const PrecisionModel prec, const NbwuKind kind,
                                    const AccumulationMethod acc_meth,
                                    const ImplicitSolventModel igb) {
  std::string k_key("gbdv_");
  return k_key + appendBornKernelKey(prec, kind, acc_meth, igb);
}

//-------------------------------------------------------------------------------------------------
std::string reductionKernelKey(const PrecisionModel prec, const ReductionGoal purpose,
                               const ReductionStage process) {
  std::string k_key("redc_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (purpose) {
  case ReductionGoal::NORMALIZE:
  case ReductionGoal::CENTER_ON_ZERO:
    break;
  case ReductionGoal::CONJUGATE_GRADIENT:
    k_key += "_cg";
    break;
  }
  switch (process) {
  case ReductionStage::GATHER:
    k_key += "_gt";
    break;
  case ReductionStage::SCATTER:
    k_key += "_sc";
    break;
  case ReductionStage::RESCALE:
    k_key += "_rs";
    break;
  case ReductionStage::ALL_REDUCE:
    k_key += "_rd";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string virtualSiteKernelKey(const PrecisionModel prec, const VirtualSiteActivity purpose) {
  std::string k_key("vste_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (purpose) {
  case VirtualSiteActivity::PLACEMENT:
    k_key += "_pl";
    break;
  case VirtualSiteActivity::TRANSMIT_FORCES:
    k_key += "_xm";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string rmsdKernelKey(const PrecisionModel prec, const RMSDTask order) {
  std::string k_key("rmsd_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (order) {
  case RMSDTask::REFERENCE:
    k_key += "_r";
    break;
  case RMSDTask::MATRIX:
    k_key += "_m";
    break;
  }
  return k_key;
}

} // namespace card
} // namespace stormm
