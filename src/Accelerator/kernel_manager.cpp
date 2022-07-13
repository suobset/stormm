// -*-c++-*-
#include "Constants/behavior.h"
#include "Math/vector_ops.h"
#include "kernel_manager.h"

namespace omni {
namespace card {

using constants::ExceptionResponse;
using math::findBin;
  
//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat() :
    block_dimension{0}, grid_dimension{0}, register_usage{0}, block_size_limit{0}, shared_usage{0}
{}

//-------------------------------------------------------------------------------------------------
int2 KernelFormat::getLaunchParameters() const {
  return { block_dimension, grid_dimension };
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
void KernelFormat::build(const int register_usage_in, const int block_size_limit_in,
                         const int shared_usage_in, const int block_dimension_in,
                         const GpuDetails &gpu) {
  register_usage = register_usage_in;
  block_size_limit = block_size_limit_in;
  shared_usage = shared_usage_in;
  const int register_file_size = gpu.getRegistersPerSMP();
  const int gpu_block_size_limit = gpu.getMaxThreadsPerBlock();
  const int smp_thread_limit = gpu.getMaxThreadsPerSMP();
  block_dimension = (block_dimension_in == 0) ? block_size_limit : block_dimension_in;
#ifdef OMNI_USE_HPC
#  ifdef OMNI_USE_CUDA
  const std::vector<int> register_break_points = {  0, 40, 48, 56, 64, 72, 80, 128, 256 };
  const std::vector<int> register_warp_counts  = { 48, 40, 36, 32, 28, 24, 16,   8 };
  const int register_bin = findBin(register_break_points, register_usage, ExceptionResponse::WARN);
  register_usage = register_break_points[register_bin];
#  endif
#endif
  grid_dimension = register_file_size / (block_dimension * register_usage);
}

//-------------------------------------------------------------------------------------------------
void KernelFormat::build(const int register_usage_in, const int block_size_limit_in,
                         const int shared_usage_in, const GpuDetails &gpu) {
  build(register_usage_in, block_size_limit_in, shared_usage_in, 0, gpu);
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
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          return valence_kernel_dfesa_dims.getLaunchParameters();
        case VwuGoal::MOVE_PARTICLES:
          return valence_kernel_dfesm_dims.getLaunchParameters();
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          return valence_kernel_dfsa_dims.getLaunchParameters();
        case VwuGoal::MOVE_PARTICLES:
          return valence_kernel_dfsm_dims.getLaunchParameters();
        }
        break;
      }
    case EvaluateForce::NO:
      return valence_kernel_de_dims.getLaunchParameters();
    }
    break;
  case PrecisionModel::SINGLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (acc_meth) {
      case ForceAccumulationMethod::SPLIT:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return valence_kernel_ffesa_dims.getLaunchParameters();
          case VwuGoal::MOVE_PARTICLES:
            return valence_kernel_ffesm_dims.getLaunchParameters();
          }
          break;
        case EvaluateEnergy::NO:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return valence_kernel_ffsa_dims.getLaunchParameters();
          case VwuGoal::MOVE_PARTICLES:
            return valence_kernel_ffsm_dims.getLaunchParameters();
          }
          break;
        }
        break;
      case ForceAccumulationMethod::WHOLE:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return valence_kernel_ffewa_dims.getLaunchParameters();
          case VwuGoal::MOVE_PARTICLES:
            return valence_kernel_ffewm_dims.getLaunchParameters();
          }
          break;
        case EvaluateEnergy::NO:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return valence_kernel_ffwa_dims.getLaunchParameters();
          case VwuGoal::MOVE_PARTICLES:
            return valence_kernel_ffwm_dims.getLaunchParameters();
          }
          break;
        }
        break;
      case ForceAccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    case EvaluateForce::NO:
      return valence_kernel_fe_dims.getLaunchParameters();
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getNonbondedKernelDims(const PrecisionModel prec,
                                           const EvaluateForce eval_force,
                                           const EvaluateEnergy eval_nrg,
                                           const ForceAccumulationMethod acc_meth) const {
  switch (prec) {
  case PrecisionModel::SINGLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          return nonbond_kernel_ffes_dims.getLaunchParameters();
        case ForceAccumulationMethod::WHOLE:
          return nonbond_kernel_ffew_dims.getLaunchParameters();
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          return nonbond_kernel_ffs_dims.getLaunchParameters();
        case ForceAccumulationMethod::WHOLE:
          return nonbond_kernel_ffw_dims.getLaunchParameters();
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      return nonbond_kernel_fe_dims.getLaunchParameters();
    }
    break;
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        return nonbond_kernel_dfes_dims.getLaunchParameters();
      case EvaluateEnergy::NO:
        return nonbond_kernel_dfs_dims.getLaunchParameters();
      }
    case EvaluateForce::NO:
      return nonbond_kernel_de_dims.getLaunchParameters();
    }
    break;
  }
  __builtin_unreachable();
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
void KernelManager::setValenceKernelAttributes(const PrecisionModel prec,
                                               const EvaluateForce eval_force,
                                               const EvaluateEnergy eval_nrg,
                                               const ForceAccumulationMethod acc_meth,
                                               const VwuGoal purpose, const int thread_limit,
                                               const int register_count, const int shared_usage) {
  switch (prec) {
  case PrecisionModel::SINGLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffsa_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffsm_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          }
          break;
        case ForceAccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffewa_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffewm_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          }
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }        
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffsa_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffsm_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          }
          break;
        case ForceAccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffwa_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffwm_dims.build(register_count, thread_limit, shared_usage, gpu);
            break;
          }
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      valence_kernel_fe_dims.build(register_count, thread_limit, shared_usage, gpu);
      break;
    }
    break;
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          valence_kernel_dfsa_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case VwuGoal::MOVE_PARTICLES:
          valence_kernel_dfsm_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          valence_kernel_dfsa_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case VwuGoal::MOVE_PARTICLES:
          valence_kernel_dfsm_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      valence_kernel_de_dims.build(register_count, thread_limit, shared_usage, gpu);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void KernelManager::setNonbondedKernelAttributes(const PrecisionModel prec,
                                                 const EvaluateForce eval_force,
                                                 const EvaluateEnergy eval_nrg,
                                                 const ForceAccumulationMethod acc_meth,
                                                 const NbwuKind kind, const int thread_limit,
                                                 const int register_count,
                                                 const int shared_usage) {
  switch (prec) {
  case PrecisionModel::SINGLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          nonbond_kernel_ffes_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case ForceAccumulationMethod::WHOLE:
          nonbond_kernel_ffew_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          nonbond_kernel_ffs_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case ForceAccumulationMethod::WHOLE:
          nonbond_kernel_ffw_dims.build(register_count, thread_limit, shared_usage, gpu);
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      nonbond_kernel_fe_dims.build(register_count, thread_limit, shared_usage, gpu);
      break;
    }
    break;
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        nonbond_kernel_dfes_dims.build(register_count, thread_limit, shared_usage, gpu);
        break;
      case EvaluateEnergy::NO:
        nonbond_kernel_dfs_dims.build(register_count, thread_limit, shared_usage, gpu);
        break;
      }
      break;
    case EvaluateForce::NO:
      nonbond_kernel_de_dims.build(register_count, thread_limit, shared_usage, gpu);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void KernelManager::setReductionKernelAttributes(const ReductionStage process,
                                                 const int thread_limit, const int register_count,
                                                 const int shared_usage) {
  switch (process) {
  case ReductionStage::GATHER:
    reduction_kernel_gt_dims.build(register_count, thread_limit, shared_usage, gpu);
  case ReductionStage::SCATTER:
    reduction_kernel_sc_dims.build(register_count, thread_limit, shared_usage, gpu);
  case ReductionStage::ALL_REDUCE:
    reduction_kernel_ar_dims.build(register_count, thread_limit, shared_usage, gpu);
  }
}

} // namespace card
} // namespace omni
