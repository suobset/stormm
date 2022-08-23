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
using math::queryReductionKernelRequirements;
using math::optReductionKernelSubdivision;
using structure::queryVirtualSiteKernelRequirements;
using synthesis::optValenceKernelSubdivision;
using synthesis::optVirtualSiteKernelSubdivision;
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
    nonbond_block_multiplier{nonbondedBlockMultiplier(gpu_in, poly_ag.getUnitCellType())},
    gbradii_block_multiplier{gbRadiiBlockMultiplier(gpu_in)},
    reduction_block_multiplier{reductionBlockMultiplier()},
    virtual_site_block_multiplier{virtualSiteBlockMultiplier()},
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
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::NO, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_d_div, "kdsValenceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       valence_d_div, "kdsValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_d_div, "kdsValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       valence_d_div, "kdsValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_d_div, "kdsValenceForceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::NO, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_fxe_div, "kfsValenceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       valence_ffx_div, "kfsValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                       valence_ffx_div, "kfValenceAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_ffx_div, "kfsValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                       AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                       valence_ffx_div, "kfValenceForceAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                       valence_ffe_div, "kfValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                       valence_ffe_div, "kfsValenceEnergyAtomUpdate");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                       valence_ffe_div, "kfValenceForceEnergyAccumulation");
  catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                       AccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                       valence_ffe_div, "kfsValenceForceEnergyAccumulation");

  // Non-bonded kernel entries
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                           EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                           "ktgdNonbondedEnergy");
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                           "ktgdsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::DOUBLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                           "ktgdsNonbondedForceEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::NO,
                           EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                           "ktgfNonbondedEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                           "ktgfsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::NO, AccumulationMethod::WHOLE,
                           "ktgfsNonbondedForce");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                           "ktgfsNonbondedForceEnergy");
    catalogNonbondedKernel(PrecisionModel::SINGLE, NbwuKind::TILE_GROUPS, EvaluateForce::YES,
                           EvaluateEnergy::YES, AccumulationMethod::WHOLE,
                           "ktgfsNonbondedForceEnergy");
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
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
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogValenceKernel(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
                                         const VwuGoal purpose, const int subdivision,
                                         const std::string &kernel_name) {
  const std::string k_key = valenceKernelKey(prec, eval_force, eval_nrg, acc_meth, purpose);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Valence kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogValenceKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryValenceKernelRequirements(prec, eval_force, eval_nrg,
                                                                 acc_meth, purpose);
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
                                           const std::string &kernel_name) {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogNonbondedKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryNonbondedKernelRequirements(prec, kind, eval_force,
                                                                   eval_nrg, acc_meth);
  k_dictionary[k_key] = KernelFormat(attr, nonbond_block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void KernelManager::catalogBornRadiiKernel(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth,
                                           const std::string &kernel_name) {
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Born radii kernel identifier " + k_key + " already exists in the kernel map.",
          "KernelManager", "catalogBornRadiiKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryBornRadiiKernelRequirements(prec, kind, acc_meth);
  k_dictionary[k_key] = KernelFormat(attr, gbradii_block_multiplier, 1, gpu, kernel_name);
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
  k_dictionary[k_key] = KernelFormat(attr, virtual_site_block_multiplier, subdivision, gpu,
                                     kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getValenceKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                         const EvaluateEnergy eval_nrg,
                                         const AccumulationMethod acc_meth,
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
                                           const AccumulationMethod acc_meth) const {
  const std::string k_key = nonbondedKernelKey(prec, kind, eval_force, eval_nrg, acc_meth);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Non-bonded kernel identifier " + k_key + " was not found in the kernel map.",
          "KernelManager", "getNonbondedKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getBornRadiiKernelDims(const PrecisionModel prec, const NbwuKind kind,
                                           const AccumulationMethod acc_meth) const {
  const std::string k_key = bornRadiiKernelKey(prec, kind, acc_meth);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Born radii computation kernel identifier " + k_key + " was not found in the kernel "
          "map.", "KernelManager", "getBornRadiiKernelDims");
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
int nonbondedBlockMultiplier(const GpuDetails &gpu, const UnitCellType unit_cell) {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
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
int gbRadiiBlockMultiplier(const GpuDetails &gpu) {
#ifdef STORMM_USE_HPC
  return (gpu.getArchMajor() == 7 && gpu.getArchMinor() >= 5) ? 4 : 5;
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
int virtualSiteBlockMultiplier() {
#ifdef STORMM_USE_HPC
  return 4;
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
std::string valenceKernelKey(const PrecisionModel prec, const EvaluateForce eval_force,
                             const EvaluateEnergy eval_nrg, const AccumulationMethod acc_meth,
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
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string nonbondedKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const EvaluateForce eval_force, const EvaluateEnergy eval_nrg,
                               const AccumulationMethod acc_meth) {
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
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string bornRadiiKernelKey(const PrecisionModel prec, const NbwuKind kind,
                               const AccumulationMethod acc_meth) {
  std::string k_key("gbrd_");
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
  
} // namespace card
} // namespace stormm
