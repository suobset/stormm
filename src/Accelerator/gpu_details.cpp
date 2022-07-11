// -*-c++-*-
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "Reporting/error_format.h"
#include "gpu_details.h"

#ifndef OMNI_USE_HPC
omni::card::GpuDetails null_gpu;
#endif

namespace omni {
namespace card {

using numerics::getForceAccumulationMethodName;

//-------------------------------------------------------------------------------------------------
GpuDetails::GpuDetails() :
    available{false}, supported{false}, arch_major{0}, arch_minor{0}, smp_count{0}, card_ram{0},
    max_threads_per_block{0}, max_threads_per_smp{0}, max_blocks_per_smp{0},
    max_shared_per_block{0}, max_shared_per_smp{0}, registers_per_block{0}, registers_per_smp{0},
    card_name{std::string("blank_gpu")}
{}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getAvailability() const {
  return available;
}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getGpuSupported() const {
  return supported;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMajor() const {
  return arch_major;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMinor() const {
  return arch_minor;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getSMPCount() const {
  return smp_count;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getCardRam() const {
  return card_ram;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerBlock() const {
  return max_threads_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerSMP() const {
  return max_threads_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxBlocksPerSMP() const {
  return max_blocks_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerBlock() const {
  return max_shared_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerSMP() const {
  return max_shared_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerBlock() const {
  return registers_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerSMP() const {
  return registers_per_smp;
}

//-------------------------------------------------------------------------------------------------
std::string GpuDetails::getCardName() const {
  return card_name;
}

//-------------------------------------------------------------------------------------------------
GpuLaunch::GpuLaunch() :
  valence_kernel_des_dims{1, 1}, valence_kernel_dfs_dims{1, 1}, valence_kernel_dfes_dims{1, 1},
  valence_kernel_fes_dims{1, 1}, valence_kernel_ffs_dims{1, 1}, valence_kernel_ffes_dims{1, 1},
  valence_kernel_few_dims{1, 1}, valence_kernel_ffw_dims{1, 1}, valence_kernel_ffew_dims{1, 1},
  nonbond_kernel_de_dims{1, 1}, nonbond_kernel_df_dims{1, 1}, nonbond_kernel_dfe_dims{1, 1},
  reduction_kernel_dims{1, 1}
{}

//-------------------------------------------------------------------------------------------------
int2 GpuLaunch::getValenceKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                     const EvaluateEnergy eval_nrg,
                                     const ForceAccumulationMethod acc_meth) {
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        return { valence_kernel_dfes_dims.getSelectedBlockDim(),
                 valence_kernel_dfes_dims.getSelectedGridDim() };
      case EvaluateEnergy::NO:
        return { valence_kernel_dfs_dims.getSelectedBlockDim(),
                 valence_kernel_dfs_dims.getSelectedGridDim() };
      }
    case EvaluateForce::NO:
      return { valence_kernel_des_dims.getSelectedBlockDim(),
               valence_kernel_des_dims.getSelectedGridDim() };
    }
    break;
  case PrecisionModel::SINGLE:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (acc_meth) {
      case ForceAccumulationMethod::SPLIT:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          return { valence_kernel_ffes_dims.getSelectedBlockDim(),
                   valence_kernel_ffes_dims.getSelectedGridDim() };
        case EvaluateEnergy::NO:
          return { valence_kernel_ffs_dims.getSelectedBlockDim(),
                   valence_kernel_ffs_dims.getSelectedGridDim() };
        }
        break;
      case ForceAccumulationMethod::WHOLE:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          return { valence_kernel_ffew_dims.getSelectedBlockDim(),
                   valence_kernel_ffew_dims.getSelectedGridDim() };
        case EvaluateEnergy::NO:
          return { valence_kernel_ffw_dims.getSelectedBlockDim(),
                   valence_kernel_ffw_dims.getSelectedGridDim() };
        }
        break;
      case ForceAccumulationMethod::AUTOMATIC:
        rtErr("Force accumulation cannot be specified as " +
              getForceAccumulationMethodName(ForceAccumulationMethod::AUTOMATIC) + " in kernel "
              "launch dimension dertermination.", "GpuLaunch", "getValenceKernelDims");
      }
      break;
    case EvaluateForce::NO:
      return { valence_kernel_fe_dims.getSelectedBlockDim(),
               valence_kernel_fe_dims.getSelectedGridDim() };
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int2 GpuLaunch::getNonbondedKernelDims(const PrecisionModel prec, const EvaluateForce eval_force,
                                       const EvaluateEnergy eval_nrg,
                                       const ForceAccumulationMethod acc_meth) const {
}

//-------------------------------------------------------------------------------------------------
int2 GpuLaunch::getReductionKernelDims() const {
  return { reduction_kernel_dims.x, reduction_kernel_dims.y };
}

//-------------------------------------------------------------------------------------------------
int GpuLaunch::getSelectedBlockDim(const int4 kval) const {
  return kval.x;
}

//-------------------------------------------------------------------------------------------------
int GpuLaunch::getSelectedGridDim(const int4 kval) const {
  return kval.y;
}

//-------------------------------------------------------------------------------------------------
int GpuLaunch::getMaximumThreadsPerBlock(const int4 kval) const {
  return (kval.z & 0xffff);
}

//-------------------------------------------------------------------------------------------------
int GpuLaunch::getRegistersPerThread(const int4 kval) const {
  return ((kval.z >> 16) & 0xffff);
}

//-------------------------------------------------------------------------------------------------
int GpuLaunch::getMaximumSharedMemoryPerBlock(const int4 kval) const {
  return kval.w;
}

//-------------------------------------------------------------------------------------------------
int4 GpuLaunch::setLaunchDims(const int registers_per_thread, const int kernel_max_threads,
                              const int shared_memory_per_block, const GpuDetails &gpu) {

}
 
} // namespace card
} // namespace omni
