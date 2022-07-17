// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Accelerator/gpu_details.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Math/rounding.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_valence_potential.cuh"

namespace omni {
namespace energy {

using card::GpuDetails;
using card::KernelManager;
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
using numerics::ForceAccumulationMethod;
using numerics::max_int_accumulation;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
using numerics::max_llint_accumulation_f;
using numerics::PrecisionModel;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::near_to_one_f;
using symbols::near_to_one_lf;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using symbols::inverse_twopi_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::ReductionStage;
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
__device__ __forceinline__ double angleVerification(const double costheta, const double3 crabbc,
                                                    const double3 crbccd, const double3 bc,
                                                    const double3 scr) {
  if (fabs(costheta) >= near_to_one_lf) {

    // The double-precision arccosine function is also vulnerable to numerical instability near
    // zero, so planar dihedral angles can still generate divergent forces on the order of 3.0e-7
    // kcal/mol-A.  Correct this with a similar strategy to the single-precision case.
    const double mg_crabbc = 1.0 / sqrt((crabbc.x * crabbc.x) + (crabbc.y * crabbc.y) +
                                        (crabbc.z * crabbc.z));
    const double mg_crbccd = 1.0 / sqrt((crbccd.x * crbccd.x) + (crbccd.y * crbccd.y) +
                                        (crbccd.z * crbccd.z));
    const double nx_abbc = crabbc.x * mg_crabbc;
    const double ny_abbc = crabbc.y * mg_crabbc;
    const double nz_abbc = crabbc.z * mg_crabbc;
    const double nx_bccd = crbccd.x * mg_crbccd;
    const double ny_bccd = crbccd.y * mg_crbccd;
    const double nz_bccd = crbccd.z * mg_crbccd;
    double rdx = nx_bccd - nx_abbc;
    double rdy = ny_bccd - ny_abbc;
    double rdz = nz_bccd - nz_abbc;
    double rs = sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    if (fabs(rs) > 1.0) {
      rdx = nx_bccd + nx_abbc;
      rdy = ny_bccd + ny_abbc;
      rdz = nz_bccd + nz_abbc;
      rs = pi - sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    }
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0) ? rs : -rs;
  }
  else {
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0) ?
            acos(costheta) : -acos(costheta);
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
#  if (__CUDA_ARCH__ == 610)
#    define VALENCE_BLOCK_MULTIPLICITY 2
#  else
#    define VALENCE_BLOCK_MULTIPLICITY 1
#  endif
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define LLCONV_FUNC __float2ll_rn
#  define SPLITCONV_FUNC convertSplitFixedPrecision
#  define SPLIT_TYPE int2
#  define SQRT_FUNC sqrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define SIN_FUNC  sinf
#  define ABS_FUNC  fabsf
#  define MIX_FUNC  computeRestraintMixtureF
#  define TCALC_IS_SINGLE
  
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      if (__CUDA_ARCH__ == 610)
#        define VALENCE_KERNEL_THREAD_COUNT medium_block_size
#      elif (__CUDA_ARCH__ == 600) || (__CUDA_ARCH__ == 700)
#        define VALENCE_KERNEL_THREAD_COUNT large_block_size
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
#  undef VALENCE_BLOCK_MULTIPLICITY
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SPLITCONV_FUNC
#  undef SPLIT_TYPE
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#  undef TCALC_IS_SINGLE
#undef TCALC

// Double-precision floating point definitions
#define TCALC double
#  if (__CUDA_ARCH__ == 610)
#    define VALENCE_KERNEL_THREAD_COUNT small_block_size
#    define VALENCE_BLOCK_MULTIPLICITY  2
#  else
#    define VALENCE_KERNEL_THREAD_COUNT medium_block_size
#    define VALENCE_BLOCK_MULTIPLICITY  1
#  endif  
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
#  define LLCONV_FUNC __double2ll_rn
#  define SPLITCONV_FUNC convertSplitFixedPrecision95
#  define SPLIT_TYPE int95_t
#  define SQRT_FUNC sqrt
#  define ACOS_FUNC acos
#  define COS_FUNC  cos
#  define SIN_FUNC  sin
#  define ABS_FUNC  fabs
#  define MIX_FUNC  computeRestraintMixtureD
#  define SPLIT_FORCE_ACCUMULATION

#  define COMPUTE_FORCE
#    define KERNEL_NAME kdsValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define UPDATE_ATOMS
#      define KERNEL_NAME kdsValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kdsValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdsValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdsValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef  COMPUTE_ENERGY

// Clear double-precision floating point definitions
#  undef VALENCE_KERNEL_THREAD_COUNT
#  undef VALENCE_BLOCK_MULTIPLICITY
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SPLITCONV_FUNC
#  undef SPLIT_TYPE
#  undef SQRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#  undef SPLIT_FORCE_ACCUMULATION
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void valenceKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kfValenceAtomUpdate, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfValenceAtomUpdate __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfValenceEnergyAtomUpdate, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfValenceEnergyAtomUpdate __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfValenceForceAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfValenceForceAccumulation __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfValenceEnergyAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfValenceEnergyAccumulation __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfValenceForceEnergyAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfValenceForceEnergyAccumulation __shared__ memory bank size to eight "
          "bytes.", "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdsValenceAtomUpdate, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdsValenceAtomUpdate __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdsValenceEnergyAtomUpdate, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdsValenceEnergyAtomUpdate __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdsValenceForceAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdsValenceForceAccumulation __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdsValenceEnergyAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdsValenceEnergyAccumulation __shared__ memory bank size to eight bytes.",
          "valenceKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdsValenceForceEnergyAccumulation, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdsValenceForceEnergyAccumulation __shared__ memory bank size to eight "
          "bytes.", "valenceKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
int2 testValenceKernelSubdivision(const int max_threads, const int smp_count, const int vwu_size,
                                  const int vwu_count) {
  const std::vector<int> block_size_options = { 128, 160, 192, 224, 256, 288, 320, 384, 448, 512,
                                                768, 896, 1024 };
  const int n_options = block_size_options.size();
  const double dvwu_count   = static_cast<double>(vwu_count);
  const double dmax_threads = static_cast<double>(max_threads);
  double best_eff = 0.0;
  int best_block_size = max_threads;
  int best_block_mult = 1;
  for (size_t i = 0; i < n_options; i++) {
    const int tmp_bdim = block_size_options[i];

    // In some cases, the valence work unit size can be slightly greater than the thread count
    // per block, i.e. 1024 atoms and 896 threads in some of the kernels that compute forces and
    // energies.  However, if the number of atoms is approaching the number of threads per block,
    // efficiency will already be fairly high.
    if (vwu_size > tmp_bdim || tmp_bdim > max_threads) {
      continue;
    }
    const int tmp_mult    = (max_threads / tmp_bdim);
    const int tmp_blocks  = tmp_mult * smp_count;
    const double tmp_batches = static_cast<double>((vwu_count + tmp_blocks - 1) / tmp_blocks);
    const double sm_eff  = static_cast<double>(tmp_bdim * tmp_mult) / dmax_threads;

    // The efficiency within an individual valence work unit is estimated as the number of tasks
    // (five task estimates per atom) divided by the number of threads available to work on them,
    // cycling until each task can be done by one of the threads.  If there are more tasks than
    // threads, this implies threads do additional batches of work, and the idle threads in the
    // final batch are less an dless inefficiency as the workload grows.  The more batches a single
    // thread performs, the less memory reading ratio is also carries, so the threads' batch count
    // is given a minor bonus in the estimated score.
    const int vwu_batches = ((5 * vwu_size) + tmp_bdim - 1) / tmp_bdim;
    const double vwu_eff = (static_cast<double>(5 * vwu_size) /
                            static_cast<double>(vwu_batches * tmp_bdim)) *
                           static_cast<double>(vwu_batches + 5) /
                           static_cast<double>(vwu_batches + 6);
    const double net_eff = sm_eff * vwu_eff * dvwu_count / (tmp_batches *
                                                            static_cast<double>(tmp_blocks));
    if (net_eff >= best_eff) {
      best_eff = net_eff;
      best_block_size = tmp_bdim;
      best_block_mult = tmp_mult;
    }
  }
  return { best_block_mult, best_block_size };
}

//-------------------------------------------------------------------------------------------------
extern void queryValenceKernelRequirements(KernelManager *wisdom,
                                           const AtomGraphSynthesis &poly_ag) {

  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes attr;
  const GpuDetails wgpu = wisdom->getGpu();
  const int arch_block_multiplier = (wgpu.getArchMajor() == 6 && wgpu.getArchMinor() == 1) ? 2 : 1;
  const int vwu_size  = poly_ag.getValenceWorkUnitSize();
  const int vwu_count = poly_ag.getValenceWorkUnitCount();
  if (cudaFuncGetAttributes(&attr, kfsValenceForceAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfsValenceAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfsValenceForceEnergyAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceForceEnergyAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfsValenceEnergyAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfsValenceEnergyAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfValenceForceAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfValenceAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfValenceAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfValenceForceEnergyAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfValenceForceEnergyAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfValenceEnergyAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfValenceEnergyAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::WHOLE, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kfValenceEnergyAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kfValenceEnergyAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::SINGLE, EvaluateForce::NO, EvaluateEnergy::YES,
                               ForceAccumulationMethod::WHOLE, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kdsValenceForceAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kdsValenceForceAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kdsValenceAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kdsValenceAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::NO,
                               ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kdsValenceForceEnergyAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kdsValenceForceEnergyAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kdsValenceEnergyAtomUpdate) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kdsValenceEnergyAtomUpdate.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::YES, EvaluateEnergy::YES,
                               ForceAccumulationMethod::SPLIT, VwuGoal::MOVE_PARTICLES,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
  if (cudaFuncGetAttributes(&attr, kdsValenceEnergyAccumulation) != cudaSuccess) {
    rtErr("Error obtaining attributes for kernel kdsValenceEnergyAccumulation.",
          "queryValenceKernelRequirements");
  }
  wisdom->catalogValenceKernel(PrecisionModel::DOUBLE, EvaluateForce::NO, EvaluateEnergy::YES,
                               ForceAccumulationMethod::SPLIT, VwuGoal::ACCUMULATE,
                               attr.maxThreadsPerBlock, arch_block_multiplier, attr.numRegs,
                               attr.sharedSizeBytes);
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceDp(const SyValenceKit<double> &poly_vk,
                            const SyRestraintKit<double, double2, double4> &poly_rk,
                            MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                            ScoreCardWriter *scw, CacheResourceKit<double> *gmem_r,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const KernelManager &launcher) {
  const int2 bt = launcher.getValenceKernelDims(PrecisionModel::DOUBLE, eval_force, eval_energy,
                                                ForceAccumulationMethod::SPLIT, purpose);
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
        kdsValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdsValenceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                   *gmem_r);
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
      kdsValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                 *gmem_r);
      break;
    case EvaluateEnergy::NO:
      kdsValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceDp(const AtomGraphSynthesis &poly_ag, MolecularMechanicsControls *mmctrl,
                            PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, CacheResource *tb_space,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
  const SyRestraintKit<double,
                       double2, double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit(tier);
  MMControlKit<double> ctrl = mmctrl->dpData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<double> gmem_r = tb_space->dpData(tier); 
  launchValenceDp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_r, eval_force, eval_energy,
                  purpose, launcher);
}

//-------------------------------------------------------------------------------------------------
extern void launchValenceSp(const SyValenceKit<float> &poly_vk,
                            const SyRestraintKit<float, float2, float4> &poly_rk,
                            MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                            ScoreCardWriter *scw, CacheResourceKit<float> *gmem_r,
                            const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                            const VwuGoal purpose, const ForceAccumulationMethod force_sum,
                            const KernelManager &launcher) {
  const int2 bt = launcher.getValenceKernelDims(PrecisionModel::SINGLE, eval_force, eval_energy,
                                                ForceAccumulationMethod::SPLIT, purpose);
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    
    // When the goal is to accumulate energies, forces, or both, the force accumulation method
    // becomes a critical detail when choosing the kernel.
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (force_sum) {
      case ForceAccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
        break;
      case ForceAccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
        break;
      case ForceAccumulationMethod::AUTOMATIC:
        if (poly_psw->frc_bits <= 23) {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfsValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfsValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        *gmem_r);
            break;
          }
        }
        else {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            kfValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             *scw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            kfValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                       *gmem_r);
            break;
          }
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      kfValenceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                  *gmem_r);
      break;
    }
    break;
  case VwuGoal::MOVE_PARTICLES:

    // When the goal is to move particles, evaluating the force is obligatory, but the manner in
    // which forces are accumulated is still important.  Whether to accumulate energies while
    // evaluating forces and moving the particles remains a consideration in choosing the proper
    // kernel.
    switch (force_sum) {
    case ForceAccumulationMethod::SPLIT:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kfsValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                   *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kfsValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::WHOLE:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kfValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                  *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kfValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case ForceAccumulationMethod::AUTOMATIC:
      if (poly_psw->frc_bits <= 23) {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfsValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                     *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfsValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
      }
      else {
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kfValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                    *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kfValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
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
                            const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
  const SyRestraintKit<float,
                       float2, float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit(tier);
  MMControlKit<float> ctrl = mmctrl->spData(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  CacheResourceKit<float> gmem_r = tb_space->spData(tier);
  launchValenceSp(poly_vk, poly_rk, &ctrl, &poly_psw, &scw, &gmem_r, eval_force, eval_energy,
                  purpose, force_sum, launcher);
}
  
} // namespace energy
} // namespace omni
