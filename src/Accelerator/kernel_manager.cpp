// -*-c++-*-
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#ifdef OMNI_USE_HPC
#include "Math/hpc_reduction.cuh"
#include "Potential/hpc_nonbonded_potential.cuh"
#include "Potential/hpc_valence_potential.cuh"
#endif
#include "kernel_manager.h"

namespace omni {
namespace card {

using constants::ExceptionResponse;
#ifdef OMNI_USE_HPC
using energy::queryValenceKernelRequirements;
using energy::queryNonbondedKernelRequirements;
using math::queryReductionKernelRequirements;
#endif
using math::findBin;
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

#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
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
    valence_block_multiplier{valenceBlockMultiplier(gpu_in, poly_ag)},
    nonbond_block_multiplier{nonbondedBlockMultiplier(gpu_in, poly_ag.getUnitCellType())},
    reduction_block_multiplier{reductionBlockMultiplier(gpu_in, poly_ag)},
    k_dictionary{}
{
#ifdef OMNI_USE_HPC
  // Valence kernel entries
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::NO, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kdsValenceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       "kdsValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kdsValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       "kdsValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kdsValenceForceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::NO, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kfsValenceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       "kfsValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                       "kfValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kfsValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                       "kfValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                       "kfValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       "kfsValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                       "kfValenceForceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       "kfsValenceForceEnergyAccumulation");

  // Non-bonded kernel entries
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                           EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                           "ktgdNonbondedEnergy");
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, ForceAccumulationMethod::SPLIT,
                           "ktgdsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                           "ktgdsNonbondedForceEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                           EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                           "ktgfNonbondedEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, ForceAccumulationMethod::SPLIT,
                           "ktgfsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, ForceAccumulationMethod::WHOLE,
                           "ktgfsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, ForceAccumulationMethod::SPLIT,
                           "ktgfsNonbondedForceEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, ForceAccumulationMethod::WHOLE,
                           "ktgfsNonbondedForceEnergy");
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }

  // Reduction kernel entries
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, "kdgtConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, "kdscConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, "kdrsConjGrad");
  catalogReductionKernel(PrecisionModel::DOUBLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, "kdrdConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::GATHER, "kfgtConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::SCATTER, "kfscConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::RESCALE, "kfrsConjGrad");
  catalogReductionKernel(PrecisionModel::SINGLE, ReductionGoal::CONJUGATE_GRADIENT,
                         ReductionStage::ALL_REDUCE, "kfrdConjGrad");
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogValenceKernel(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const ForceAccumulationMethod acc_meth,
                                         const VwuGoal purpose, const std::string &kernel_name) {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogValenceKernel");
  }
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  const cudaFuncAttributes attr = queryValenceKernelRequirements(prec, eval_force, eval_nrg,
                                                                 acc_meth, purpose);
  k_dictionary[k_key] = KernelFormat(attr, valence_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogNonbondedKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const ForceAccumulationMethod acc_meth,
                                           const std::string &kernel_name) {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogNonbondedKernel");
  }
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  const cudaFuncAttributes attr = queryNonbondedKernelRequirements(prec, kind, eval_force,
                                                                   eval_nrg, acc_meth);
  k_dictionary[k_key] = KernelFormat(attr, nonbond_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogReductionKernel(const PrecisionModel prec, const ReductionGoal purpose,
                                           const ReductionStage process,
                                           const std::string &kernel_name) {
  const std::string k_key = reductionKernelKey(prec, purpose, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Reduction kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogReductionKernel");
  }
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  const cudaFuncAttributes attr = queryReductionKernelRequirements(prec, purpose, process);
  k_dictionary[k_key] = KernelFormat(attr, reduction_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

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
int valenceBlockMultiplier(const GpuDetails &gpu, const AtomGraphSynthesis &poly_ag) {
#ifdef OMNI_USE_HPC
  return 1;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
int nonbondedBlockMultiplier(const GpuDetails &gpu, const UnitCellType unit_cell) {
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  switch (unit_cell) {
  case UnitCellType::NONE:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 2 : 3;
  }
#  else
  // Other vendors are not known to make GPUs that have special requirements
  switch (unit_cell) {
  case UnitCellType::NONE:
    return 5;
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
int reductionBlockMultiplier(const GpuDetails &gpu, const AtomGraphSynthesis &poly_ag) {
#ifdef OMNI_USE_HPC
  return 1;
#else
  return 1;
#endif
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
    k_key += "cg";
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

} // namespace card
} // namespace omni
