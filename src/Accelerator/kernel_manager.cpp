// -*-c++-*-
#include "Constants/behavior.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#ifdef OMNI_USE_HPC
#include "Potential/hpc_nonbonded_potential.cuh"
#include "Potential/hpc_valence_potential.cuh"
#endif
#include "kernel_manager.h"

namespace omni {
namespace card {

using constants::ExceptionResponse;
#ifdef OMNI_USE_HPC
using energy::queryValenceKernelRequirements;
#endif
using math::findBin;
using parse::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat() :
    block_dimension{0}, grid_dimension{0}, register_usage{0}, block_size_limit{0}, shared_usage{0}
{}

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
void KernelFormat::build(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                         const int register_usage_in, const int shared_usage_in,
                         const int block_subdivision, const GpuDetails &gpu,
                         const std::string &kernel_name_in) {

  // Compute the block dimension and the launch grid dimension.
  block_size_limit = lb_max_threads_per_block;
  shared_usage = shared_usage_in;
  if (kernel_name_in.size() > 0LLU) {
    kernel_name = kernel_name_in;
  }
  block_dimension = (block_size_limit / block_subdivision / warp_size_int) * warp_size_int;
  grid_dimension = block_subdivision * lb_min_blocks_per_smp * gpu.getSMPCount();

  // Get the register usage (this needs refinement and is unreliable with current cudart function
  // calls).
  register_usage = register_usage_in;
  const std::vector<int> register_break_points = {  0, 40, 48, 56, 64, 72, 80, 128, 256 };
  const std::vector<int> register_warp_counts  = { 48, 40, 36, 32, 28, 24, 16,   8 };
  const int register_bin = findBin(register_break_points, register_usage, ExceptionResponse::WARN);
  register_usage = register_warp_counts[register_bin];
}

//-------------------------------------------------------------------------------------------------
void KernelFormat::build(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                         const int register_usage_in, const int shared_usage_in,
                         const GpuDetails &gpu, const std::string &kernel_name_in) {
  build(lb_max_threads_per_block, lb_min_blocks_per_smp, register_usage_in, shared_usage_in, 1,
        gpu, kernel_name_in);
}

//-------------------------------------------------------------------------------------------------
KernelManager::KernelManager(const GpuDetails &gpu_in) :
    valence_kernel_de_dims{}, valence_kernel_dfsm_dims{}, valence_kernel_dfsa_dims{},
    valence_kernel_dfesm_dims{}, valence_kernel_dfesa_dims{}, valence_kernel_fe_dims{},
    valence_kernel_ffsm_dims{}, valence_kernel_ffwm_dims{}, valence_kernel_ffsa_dims{},
    valence_kernel_ffwa_dims{}, valence_kernel_ffewm_dims{}, valence_kernel_ffesm_dims{},
    valence_kernel_ffewa_dims{}, valence_kernel_ffesa_dims{}, nonbond_kernel_de_dims{},
    nonbond_kernel_dfs_dims{}, nonbond_kernel_dfes_dims{}, nonbond_kernel_fe_dims{},
    nonbond_kernel_ffs_dims{}, nonbond_kernel_ffw_dims{}, nonbond_kernel_ffes_dims{},
    nonbond_kernel_ffew_dims{}, reduction_kernel_gt_dims{}, reduction_kernel_sc_dims{},
    reduction_kernel_ar_dims{}, gpu{gpu_in}
{}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getValenceKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const ForceAccumulationMethod acc_meth,
                                         const VwuGoal purpose) const {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose);
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
                                           const ForceAccumulationMethod acc_meth) const {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getNonbondedKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getReductionKernelDims(const ReductionStage process) const {
  switch (process) {
  case ReductionStage::GATHER:
    return reduction_kernel_gt_dims.getLaunchParameters();
  case ReductionStage::SCATTER:
    return reduction_kernel_sc_dims.getLaunchParameters();
  case ReductionStage::ALL_REDUCE:
    return reduction_kernel_ar_dims.getLaunchParameters();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const GpuDetails& KernelManager::getGpu() const {
  return gpu;
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogValenceKernel(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const ForceAccumulationMethod acc_meth,
                                         const VwuGoal purpose, const int lb_max_threads_per_block,
                                         const int lb_max_blocks_per_smp, const int register_count,
                                         const int shared_usage, const int block_subdivision) {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + "  not found in the kernel map.",
          "KernelManager", "getValenceKernelDims");
  }
  k_dictionary[k_key].build(lb_max_threads_per_block, lb_max_blocks_per_smp, register_count,
                            shared_usage, block_subdivision, gpu);
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogNonbondedKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const ForceAccumulationMethod acc_meth,
                                           const int lb_max_threads_per_block,
                                           const int lb_max_blocks_per_smp,
                                           const int register_count, const int shared_usage,
                                           const int block_subdivision) {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + "  not found in the kernel map.",
          "KernelManager", "getNonbondedKernelDims");
  }
  k_dictionary[k_key].build(lb_max_threads_per_block, lb_max_blocks_per_smp, register_count,
                            shared_usage, block_subdivision, gpu);
}

//-------------------------------------------------------------------------------------------------
void KernelManager::setReductionKernelAttributes(const ReductionStage process,
                                                 const int thread_limit, const int register_count,
                                                 const int shared_usage, const int block_mult) {
  switch (process) {
  case ReductionStage::GATHER:
    reduction_kernel_gt_dims.build(register_count, thread_limit, shared_usage, block_mult,
                                   thread_limit, gpu);
  case ReductionStage::SCATTER:
    reduction_kernel_sc_dims.build(register_count, thread_limit, shared_usage, block_mult,
                                   thread_limit, gpu);
  case ReductionStage::ALL_REDUCE:
    reduction_kernel_ar_dims.build(register_count, thread_limit, shared_usage, block_mult,
                                   thread_limit, gpu);
  }
}

//-------------------------------------------------------------------------------------------------
void KernelManager::printLaunchParameters(const std::string &k_key) const {
  if (k_key.size() == 0) {
    return;
  }
  if (strcmpCased(k_key, "all", CaseSensitivity::NO)) {
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
std::string valenceKernelKey(const PrecisionModel prec, const EvaluateForce eval_force,
                             const EvaluateEnergy eval_nrg, const ForceAccumulationMethod acc_meth,
                             const VwuGoal purpose) {
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
    case ForceAccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case ForceAccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case ForceAccumulationMethod::AUTOMATIC:
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
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string nonbondedKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const EvaluateForce eval_force, const EvaluateEnergy eval_nrg,
                               const ForceAccumulationMethod acc_meth) {
  std::string k_key("nonb_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "s";
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
    case ForceAccumulationMethod::SPLIT:
      k_key += "s";
      break;
    case ForceAccumulationMethod::WHOLE:
      k_key += "w";
      break;
    case ForceAccumulationMethod::AUTOMATIC:
      break;
    }
  }
  return k_key;
}

} // namespace card
} // namespace omni
