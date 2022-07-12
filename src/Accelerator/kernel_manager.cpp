// -*-c++-*-
#include "kernel_manager.h"

namespace omni {
namespace card {

//-------------------------------------------------------------------------------------------------
KernelManager::KernelManager() :
    valence_kernel_de_dims{1, 1, 0, 0}, valence_kernel_dfsm_dims{1, 1, 0, 0},
    valence_kernel_dfsa_dims{1, 1, 0, 0}, valence_kernel_dfesm_dims{1, 1, 0, 0},
    valence_kernel_dfesa_dims{1, 1, 0, 0}, valence_kernel_fe_dims{1, 1, 0, 0},
    valence_kernel_ffsm_dims{1, 1, 0, 0}, valence_kernel_ffwm_dims{1, 1, 0, 0},
    valence_kernel_ffsa_dims{1, 1, 0, 0}, valence_kernel_ffwa_dims{1, 1, 0, 0},
    valence_kernel_ffewm_dims{1, 1, 0, 0}, valence_kernel_ffesm_dims{1, 1, 0, 0},
    valence_kernel_ffewa_dims{1, 1, 0, 0}, valence_kernel_ffesa_dims{1, 1, 0, 0},
    nonbond_kernel_de_dims{1, 1, 0, 0}, nonbond_kernel_dfs_dims{1, 1, 0, 0},
    nonbond_kernel_dfes_dims{1, 1, 0, 0}, nonbond_kernel_fe_dims{1, 1, 0, 0},
    nonbond_kernel_ffs_dims{1, 1, 0, 0}, nonbond_kernel_ffw_dims{1, 1, 0, 0},
    nonbond_kernel_ffes_dims{1, 1, 0, 0}, nonbond_kernel_ffew_dims{1, 1, 0, 0},
    reduction_kernel_gt_dims{1, 1, 0, 0}, reduction_kernel_sc_dims{1, 1, 0, 0},
    reduction_kernel_ar_dims{1, 1, 0, 0}
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
          return { getSelectedBlockDim(valence_kernel_dfesa_dims),
                   getSelectedGridDim(valence_kernel_dfesa_dims) };
        case VwuGoal::MOVE_PARTICLES:
          return { getSelectedBlockDim(valence_kernel_dfesm_dims),
                   getSelectedGridDim(valence_kernel_dfesm_dims) };
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          return { getSelectedBlockDim(valence_kernel_dfsa_dims),
                   getSelectedGridDim(valence_kernel_dfsa_dims) };
        case VwuGoal::MOVE_PARTICLES:
          return { getSelectedBlockDim(valence_kernel_dfsm_dims),
                   getSelectedGridDim(valence_kernel_dfsm_dims) };
        }
        break;
      }
    case EvaluateForce::NO:
      return { getSelectedBlockDim(valence_kernel_de_dims),
               getSelectedGridDim(valence_kernel_de_dims) };
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
            return { getSelectedBlockDim(valence_kernel_ffesa_dims),
                     getSelectedGridDim(valence_kernel_ffesa_dims) };
          case VwuGoal::MOVE_PARTICLES:
            return { getSelectedBlockDim(valence_kernel_ffesm_dims),
                     getSelectedGridDim(valence_kernel_ffesm_dims) };
          }
          break;
        case EvaluateEnergy::NO:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return { getSelectedBlockDim(valence_kernel_ffsa_dims),
                     getSelectedGridDim(valence_kernel_ffsa_dims) };
          case VwuGoal::MOVE_PARTICLES:
            return { getSelectedBlockDim(valence_kernel_ffsm_dims),
                     getSelectedGridDim(valence_kernel_ffsm_dims) };
          }
          break;
        }
        break;
      case ForceAccumulationMethod::WHOLE:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return { getSelectedBlockDim(valence_kernel_ffewa_dims),
                     getSelectedGridDim(valence_kernel_ffewa_dims) };
          case VwuGoal::MOVE_PARTICLES:
            return { getSelectedBlockDim(valence_kernel_ffewm_dims),
                     getSelectedGridDim(valence_kernel_ffewm_dims) };
          }
          break;
        case EvaluateEnergy::NO:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            return { getSelectedBlockDim(valence_kernel_ffwa_dims),
                     getSelectedGridDim(valence_kernel_ffwa_dims) };
          case VwuGoal::MOVE_PARTICLES:
            return { getSelectedBlockDim(valence_kernel_ffwm_dims),
                     getSelectedGridDim(valence_kernel_ffwm_dims) };
          }
          break;
        }
        break;
      case ForceAccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    case EvaluateForce::NO:
      return { getSelectedBlockDim(valence_kernel_fe_dims),
               getSelectedGridDim(valence_kernel_fe_dims) };
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
          return { getSelectedBlockDim(nonbond_kernel_ffes_dims),
                   getSelectedGridDim(nonbond_kernel_ffes_dims) };
        case ForceAccumulationMethod::WHOLE:
          return { getSelectedBlockDim(nonbond_kernel_ffew_dims),
                   getSelectedGridDim(nonbond_kernel_ffew_dims) };
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          return { getSelectedBlockDim(nonbond_kernel_ffs_dims),
                   getSelectedGridDim(nonbond_kernel_ffs_dims) };
        case ForceAccumulationMethod::WHOLE:
          return { getSelectedBlockDim(nonbond_kernel_ffw_dims),
                   getSelectedGridDim(nonbond_kernel_ffw_dims) };
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      return { getSelectedBlockDim(nonbond_kernel_fe_dims),
               getSelectedGridDim(nonbond_kernel_fe_dims) };
    }
    break;
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        return { getSelectedBlockDim(nonbond_kernel_dfes_dims),
                 getSelectedGridDim(nonbond_kernel_dfes_dims) };
      case EvaluateEnergy::NO:
        return { getSelectedBlockDim(nonbond_kernel_dfs_dims),
                 getSelectedGridDim(nonbond_kernel_dfs_dims) };
      }
    case EvaluateForce::NO:
      return { getSelectedBlockDim(nonbond_kernel_de_dims),
               getSelectedGridDim(nonbond_kernel_de_dims) };
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::getReductionKernelDims(const ReductionStage process) const {
  switch (process) {
  case ReductionStage::GATHER:
    return { reduction_kernel_ar_dims.x, reduction_kernel_gt_dims.y };
  case ReductionStage::SCATTER:
    return { reduction_kernel_ar_dims.x, reduction_kernel_sc_dims.y };
  case ReductionStage::ALL_REDUCE:
    return { reduction_kernel_ar_dims.x, reduction_kernel_ar_dims.y };
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
            valence_kernel_ffsa_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffsa_dims.w = shared_usage;
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffsm_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffsm_dims.w = shared_usage;
            break;
          }
          break;
        case ForceAccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffewa_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffewa_dims.w = shared_usage;
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffewm_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffewm_dims.w = shared_usage;
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
            valence_kernel_ffsa_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffsa_dims.w = shared_usage;
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffsm_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffsm_dims.w = shared_usage;
            break;
          }
          break;
        case ForceAccumulationMethod::WHOLE:
          switch (purpose) {
          case VwuGoal::ACCUMULATE:
            valence_kernel_ffwa_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffwa_dims.w = shared_usage;
            break;
          case VwuGoal::MOVE_PARTICLES:
            valence_kernel_ffwm_dims.z = ((register_count << 16) | thread_limit);
            valence_kernel_ffwm_dims.w = shared_usage;
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
      valence_kernel_fe_dims.z = ((register_count << 16) | thread_limit);
      valence_kernel_fe_dims.z = shared_usage;
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
          valence_kernel_dfsa_dims.z = ((register_count << 16) | thread_limit);
          valence_kernel_dfsa_dims.w = shared_usage;
          break;
        case VwuGoal::MOVE_PARTICLES:
          valence_kernel_dfsm_dims.z = ((register_count << 16) | thread_limit);
          valence_kernel_dfsm_dims.w = shared_usage;
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          valence_kernel_dfsa_dims.z = ((register_count << 16) | thread_limit);
          valence_kernel_dfsa_dims.w = shared_usage;
          break;
        case VwuGoal::MOVE_PARTICLES:
          valence_kernel_dfsm_dims.z = ((register_count << 16) | thread_limit);
          valence_kernel_dfsm_dims.w = shared_usage;
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      valence_kernel_de_dims.z = ((register_count << 16) | thread_limit);
      valence_kernel_de_dims.w = shared_usage;
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
          nonbond_kernel_ffes_dims.z = ((register_count << 16) | thread_limit);
          nonbond_kernel_ffes_dims.w = shared_usage;
          break;
        case ForceAccumulationMethod::WHOLE:
          nonbond_kernel_ffew_dims.z = ((register_count << 16) | thread_limit);
          nonbond_kernel_ffew_dims.w = shared_usage;
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (acc_meth) {
        case ForceAccumulationMethod::SPLIT:
          nonbond_kernel_ffs_dims.z = ((register_count << 16) | thread_limit);
          nonbond_kernel_ffs_dims.w = shared_usage;
          break;
        case ForceAccumulationMethod::WHOLE:
          nonbond_kernel_ffw_dims.z = ((register_count << 16) | thread_limit);
          nonbond_kernel_ffw_dims.w = shared_usage;
          break;
        case ForceAccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      nonbond_kernel_fe_dims.z = ((register_count << 16) | thread_limit);
      nonbond_kernel_fe_dims.w = shared_usage;
      break;
    }
    break;
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        nonbond_kernel_dfes_dims.z = ((register_count << 16) | thread_limit);
        nonbond_kernel_dfes_dims.w = shared_usage;
        break;
      case EvaluateEnergy::NO:
        nonbond_kernel_dfs_dims.z = ((register_count << 16) | thread_limit);
        nonbond_kernel_dfs_dims.w = shared_usage;
        break;
      }
      break;
    case EvaluateForce::NO:
      nonbond_kernel_de_dims.z = ((register_count << 16) | thread_limit);
      nonbond_kernel_de_dims.w = shared_usage;
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void KernelManager::setReductionKernelAttributes(const ReductionStage process,
                                                 const int thread_limit, const int register_count,
                                                 const int shared_usage) {

}

//-------------------------------------------------------------------------------------------------
int KernelManager::getSelectedBlockDim(const int4 kval) const {
  return kval.x;
}

//-------------------------------------------------------------------------------------------------
int KernelManager::getSelectedGridDim(const int4 kval) const {
  return kval.y;
}

//-------------------------------------------------------------------------------------------------
int KernelManager::getMaximumThreadsPerBlock(const int4 kval) const {
  return (kval.z & 0xffff);
}

//-------------------------------------------------------------------------------------------------
int KernelManager::getRegistersPerThread(const int4 kval) const {
  return ((kval.z >> 16) & 0xffff);
}

//-------------------------------------------------------------------------------------------------
int KernelManager::getMaximumSharedMemoryPerBlock(const int4 kval) const {
  return kval.w;
}

//-------------------------------------------------------------------------------------------------
int2 KernelManager::setLaunchDims(const int4 kernel_layout, const GpuDetails &gpu) {
  int block_dim = 1;
  int grid_dim = 1;
  const int register_file_size = gpu.getRegistersPerSMP();
  const int registers_per_thread = ((kernel_layout.z >> 16) & 0xffff);
  const int gpu_block_size_limit = gpu.getMaxThreadsPerBlock();
  const int compilation_block_size_limit = (kernel_layout.z & 0xffff);
  const int smp_thread_limit = gpu.getMaxThreadsPerSMP();
  
  return { block_dim, grid_dim };
}

} // namespace card
} // namespace omni
