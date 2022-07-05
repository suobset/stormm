// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Math/rounding.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_valence_potential.cuh"

namespace omni {
namespace energy {

using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using data_types::int95_t;
using math::roundUp;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
using numerics::max_llint_accumulation_f;
using numerics::max_llint_accumulation_ll;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::near_to_one_f;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using symbols::inverse_twopi_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::VwuAbstractMap;
using synthesis::VwuGoal;
using synthesis::vwu_abstract_length;
using topology::TorsionKind;
  
#include "accumulation.cui"
#include "Math/rounding.cui"

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ double3 crossProduct(const double3 a, const double3 b) {
  return { (a.y * b.z) - (a.z * b.y), (a.z * b.x) - (a.x * b.z), (a.x * b.y) - (a.y * b.x) };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ float3 crossProduct(const float3 a, const float3 b) {
  return { (a.y * b.z) - (a.z * b.y), (a.z * b.x) - (a.x * b.z), (a.x * b.y) - (a.y * b.x) };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ float angleVerification(const float costheta, const float3 crabbc,
                                                   const float3 crbccd, const float3 bc,
                                                   const float3 scr) {
  if (fabsf(costheta) >= near_to_one_f) {

    // The floating-point representation of costheta is numerically ill-conditioned.  Compute the
    // distance from atom I to the plane of atoms J, K, and L to get the angle by the arcsin of an
    // extremely acute angle.
    const float mg_crabbc = 1.0f / sqrtf((crabbc.x * crabbc.x) + (crabbc.y * crabbc.y) +
                                         (crabbc.z * crabbc.z));
    const float mg_crbccd = 1.0f / sqrtf((crbccd.x * crbccd.x) + (crbccd.y * crbccd.y) +
                                         (crbccd.z * crbccd.z));
    const float nx_abbc = crabbc.x * mg_crabbc;
    const float ny_abbc = crabbc.y * mg_crabbc;
    const float nz_abbc = crabbc.z * mg_crabbc;
    const float nx_bccd = crbccd.x * mg_crbccd;
    const float ny_bccd = crbccd.y * mg_crbccd;
    const float nz_bccd = crbccd.z * mg_crbccd;
    float rdx = nx_bccd - nx_abbc;
    float rdy = ny_bccd - ny_abbc;
    float rdz = nz_bccd - nz_abbc;
    float rs = sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    if (fabsf(rs) > 1.0f) {
      rdx = nx_bccd + nx_abbc;
      rdy = ny_bccd + ny_abbc;
      rdz = nz_bccd + nz_abbc;
      rs = pi_f - sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    }
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0f) ? rs : -rs;
  }
  else {
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0f) ?
            acosf(costheta) : -acosf(costheta);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
double3 restraintDelta(const double2 init_k, const double2 final_k, const double4 init_r,
                      const double4 final_r, const double2 mixwt, const double dr) {
  const double r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const double r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const double r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const double r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const double k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const double k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  double dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return { keq, dl, du };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
float3 restraintDelta(const float2 init_k, const float2 final_k, const float4 init_r,
                      const float4 final_r, const float2 mixwt, const float dr) {
  const float r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const float r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const float r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const float r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const float k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const float k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  float dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return { keq, dl, du };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
double2 computeRestraintMixtureD(const int step_number, const int init_step,
                                 const int final_step) {
  if (step_number < init_step) {

    // If the restraint has not yet engaged, neither its initial or final values have any weight
    return { (double)(0.0), (double)(0.0) };
  }
  else if (init_step == final_step) {

    // The step count is far enough along that the restraint has been engaged, and it is constant.
    // Only the initial value matters.
    return { (double)(1.0), (double)(0.0) };
  }
  else if (step_number < final_step) {
    const double wslide = (double)(step_number - init_step) / (double)(final_step - init_step);

    // The difference between the initial and final steps is nonzero.  The mixture is a linear
    // combination of the two end points.
    return { (double)(1.0) - wslide, wslide };
  }

  // The step number has advanced beyond the point at which the restraint is mature.
  return { (double)(0.0), (double)(1.0) };
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
float2 computeRestraintMixtureF(const int step_number, const int init_step, const int final_step) {
  if (step_number < init_step) {

    // If the restraint has not yet engaged, neither its initial or final values have any weight
    return { (float)(0.0), (float)(0.0) };
  }
  else if (init_step == final_step) {

    // The step count is far enough along that the restraint has been engaged, and it is constant.
    // Only the initial value matters.
    return { (float)(1.0), (float)(0.0) };
  }
  else if (step_number < final_step) {
    const float wslide = (float)(step_number - init_step) / (float)(final_step - init_step);

    // The difference between the initial and final steps is nonzero.  The mixture is a linear
    // combination of the two end points.
    return { (float)(1.0) - wslide, wslide };
  }

  // The step number has advanced beyond the point at which the restraint is mature.
  return { (float)(0.0), (float)(1.0) };
}

// Single-precision floating point definitions
#define TCALC float
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define SIN_FUNC  sinf
#  define ABS_FUNC  fabsf
#  define MIX_FUNC  computeRestraintMixtureF
#  define CHECK_COSARG

#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      if (__CUDA_ARCH__ == 610)
#        define VALENCE_KERNEL_THREAD_COUNT 512
#      elif (__CUDA_ARCH__ == 600) || (__CUDA_ARCH__ == 700)
#        define VALENCE_KERNEL_THREAD_COUNT 1024
#      else
#        define VALENCE_KERNEL_THREAD_COUNT large_block_size
#      endif
#      define KERNEL_NAME kfsValenceForceAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME  
#      define UPDATE_ATOMS
#        define KERNEL_NAME kfsValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        if (__CUDA_ARCH__ == 610)
#          define VALENCE_KERNEL_THREAD_COUNT 448
#        else
#          define VALENCE_KERNEL_THREAD_COUNT 896
#        endif
#        define KERNEL_NAME kfsValenceForceEnergyAccumulation
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        define UPDATE_ATOMS
#          define KERNEL_NAME kfsValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef UPDATE_ATOMS
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef SPLIT_FORCE_ACCUMULATION
#    if (__CUDA_ARCH__ == 610)
#      define VALENCE_KERNEL_THREAD_COUNT medium_block_size
#    else
#      define VALENCE_KERNEL_THREAD_COUNT large_block_size
#    endif
#    define KERNEL_NAME kfValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define UPDATE_ATOMS
#      define KERNEL_NAME kfValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      if (__CUDA_ARCH__ == 610)
#        define VALENCE_KERNEL_THREAD_COUNT 448
#      else
#        define VALENCE_KERNEL_THREAD_COUNT 896
#      endif
#      define KERNEL_NAME kfValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kfValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    if (__CUDA_ARCH__ == 610)
#      define VALENCE_KERNEL_THREAD_COUNT medium_block_size
#    else
#      define VALENCE_KERNEL_THREAD_COUNT large_block_size
#    endif
#    define KERNEL_NAME kfValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef  COMPUTE_ENERGY

// Clear single-precision floating point definitions
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#  undef CHECK_COSARG
#undef TCALC

// Double-precision floating point definitions
#define TCALC double
#  if (__CUDA_ARCH__ == 610)
#    define VALENCE_KERNEL_THREAD_COUNT small_block_size
#  else
#    define VALENCE_KERNEL_THREAD_COUNT medium_block_size
#  endif  
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define ACOS_FUNC acos
#  define COS_FUNC  cos
#  define SIN_FUNC  sin
#  define ABS_FUNC  fabs
#  define MIX_FUNC  computeRestraintMixtureD

#  define COMPUTE_FORCE
#    define KERNEL_NAME kdValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define UPDATE_ATOMS
#      define KERNEL_NAME kdValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kdValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear double-precision floating point definitions
#  undef VALENCE_KERNEL_THREAD_COUNT
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void valenceKernelSetup() {
  cudaFuncSetSharedMemConfig(kfValenceAtomUpdate,              cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kfValenceEnergyAtomUpdate,        cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kfValenceForceAccumulation,       cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kfValenceEnergyAccumulation,      cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kfValenceForceEnergyAccumulation, cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kdValenceAtomUpdate,              cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kdValenceEnergyAtomUpdate,        cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kdValenceForceAccumulation,       cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kdValenceEnergyAccumulation,      cudaSharedMemBankSizeEightByte);
  cudaFuncSetSharedMemConfig(kdValenceForceEnergyAccumulation, cudaSharedMemBankSizeEightByte);
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceDp(const SyValenceKit<double> &poly_vk,
                            const SyRestraintKit<double, double2, double4> &poly_rk,
                            MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                            ScoreCardWriter *scw, CacheResourceKit<double> *gmem_r,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const GpuDetails &gpu) {
  const int blocks_multiplier = (gpu.getArchMajor() == 6 && gpu.getArchMinor() == 1) ? 2 : 1;
  const int nblocks = gpu.getSMPCount() * blocks_multiplier;
  const int nthreads = gpu.getMaxThreadsPerBlock() / (2 * blocks_multiplier);
  switch (purpose) {
  case VwuGoal::ACCUMULATE:

    // When the goal is to accumulate energies, forces, or both, the force accumulation method
    // is set to use int64 data.  A 95-bit method that splits the accumulation with overflow into
    // a secondary 32-bit int may be added, and likewise become the sole option for
    // double-precision computations.
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdValenceForceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                                *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdValenceForceAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          *gmem_r);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdValenceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         *scw, *gmem_r);
      break;
    }
    break;
  case VwuGoal::MOVE_PARTICLES:

    // When the goal is to move particles, evaluating the force is obligatory, but the manner in
    // which forces are accumulated is still important.  Whether to accumulate energies while
    // evaluating forces and moving the particles remains a consideration in choosing the proper
    // kernel.
    switch (eval_energy) {
    case EvaluateEnergy::YES:
      kdValenceEnergyAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                       *gmem_r);
      break;
    case EvaluateEnergy::NO:
      kdValenceAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceDp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                            PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
  const SyRestraintKit<double,
                       double2, double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit(tier);
  MMControlKit<double> ctrl = mmctrl->dpData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<double> gmem_r = tb_space->dpData(tier); 
  launchValenceDp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_r, eval_force, eval_energy,
                  purpose, gpu);
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceSp(const SyValenceKit<float> &poly_vk,
                            const SyRestraintKit<float, float2, float4> &poly_rk,
                            MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                            ScoreCardWriter *scw, CacheResourceKit<float> *gmem_r,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const ForceAccumulationMethod force_sum,
                            const GpuDetails &gpu) {

  // Determine the number of blocks to use, based on the architecture.  Other aspects of the
  // valence work units themselves will be adjusted based on the architecture in order to avoid
  // over-stuffing the L1 memory on the cards.  The L1 space is growing with newer generations,
  // but a lowest common denominator of 128kB will be assumed in most aspects of the code.
  const int major_arch = gpu.getArchMajor();
  const int minor_arch = gpu.getArchMinor();
  const int blocks_multiplier = (major_arch == 6 && minor_arch == 1) ? 2 : 1;
  const int nblocks = gpu.getSMPCount() * blocks_multiplier;
  
  // Have ready a maximum number of threads per block, and a fallback number if the register
  // pressure might be too high to employ the full thread complement.
  const int max_threads   = gpu.getMaxThreadsPerBlock() / blocks_multiplier;
  const int trim_threads  = roundUp<int>((max_threads * 7) / 8, twice_warp_size_int);
  const int lower_threads = roundUp<int>((max_threads * 3) / 4, twice_warp_size_int);
  int nthreads;
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    switch (eval_energy) {
    case EvaluateEnergy::YES:

      // The combination of energy and force evaluation accumulates register pressure
      switch (eval_force) {
      case EvaluateForce::YES:
        nthreads = trim_threads;
        break;
      case EvaluateForce::NO:
        nthreads = max_threads;
        break;
      }
      break;
    case EvaluateEnergy::NO:
      nthreads = max_threads;
      break;
    }
    
    // When the goal is to accumulate energies, forces, or both, the force accumulation method
    // becomes a critical detail when choosing the kernel.
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (force_sum) {
      case ForceAccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceForceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceForceAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             *gmem_r);
          break;
        }
        break;
      case ForceAccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceForceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceForceAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *gmem_r);
          break;
        }
        break;
      case ForceAccumulationMethod::AUTOMATIC:
        if (poly_psw->frc_bits <= 23) {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfsValenceForceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl,
                                                                     *poly_psw, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfsValenceForceAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                               *gmem_r);
            break;
          }
        }
        else {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfValenceForceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfValenceForceAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              *gmem_r);
            break;
          }
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      kfValenceEnergyAccumulation<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         *scw, *gmem_r);
      break;
    }
    break;
  case VwuGoal::MOVE_PARTICLES:
    switch (force_sum) {
    case ForceAccumulationMethod::SPLIT:
      nthreads = (major_arch == 7 && minor_arch == 0) ? trim_threads : max_threads;
      break;
    case ForceAccumulationMethod::WHOLE:
      nthreads = max_threads;
      break;
    }

    // When the goal is to move particles, evaluating the force is obligatory, but the manner in
    // which forces are accumulated is still important.  Whether to accumulate energies while
    // evaluating forces and moving the particles remains a consideration in choosing the proper
    // kernel.
    switch (force_sum) {
    case ForceAccumulationMethod::SPLIT:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kfsValenceEnergyAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                          *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kfsValenceAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::WHOLE:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kfValenceEnergyAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                         *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kfValenceAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::AUTOMATIC:
      if (poly_psw->frc_bits <= 23) {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceEnergyAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      else {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceEnergyAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceAtomUpdate<<<nblocks, nthreads>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceSp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                            PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const ForceAccumulationMethod force_sum,
                            const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
  const SyRestraintKit<float,
                       float2, float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit(tier);
  MMControlKit<float> ctrl = mmctrl->spData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<float> gmem_r = tb_space->spData(tier);
  launchValenceSp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_r, eval_force, eval_energy,
                  purpose, force_sum, gpu);
}
  
} // namespace energy
} // namespace omni
