// -*-c++-*-
#include "Accelerator/ptx_macros.h"
#include "Constants/scaling.h"
#include "hpc_reduction.cuh"

namespace omni {
namespace math {
  
//-------------------------------------------------------------------------------------------------
// Perform an accumulation over the conjugate gradient work units' relevant atoms to obtain the
// sum of squared gradients (gg) and the evolution of the gradient (dgg).
//
// Overloaded:
//   - Operate on standard fixed-precision data in long long int format
//   - Operate on extended fixed-precision data in { long long int, int } format
//
// Arguments:
//   gg_collector:     __shared__ memory resource for storing gg accumulation per warp
//   dgg_collector:    __shared__ memory resource for storing dgg accumulation per warp
//   start_pos:        First atom affected by this work unit in the following global arrays
//   end_pos:          Upper limit of atoms affected by this work unit
//   xfrc:             Current forces on particles in the Cartesian X direction
//   yfrc:             Current forces on particles in the Cartesian Y direction
//   zfrc:             Current forces on particles in the Cartesian Z direction
//   [x,y,z]frc_ovrf:  Overflow in X, Y, or Z forces
//   xprv:             Prior iteration forces on particles in the Cartesian X direction
//   yprv:             Prior iteration forces on particles in the Cartesian Y direction
//   zprv:             Prior iteration forces on particles in the Cartesian Z direction
//   [x,y,z]prv_ovrf:  Overflow in prior iteration X, Y, or Z forces
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
void conjGradCoreGather(double* gg_collector, double* dgg_collector, const int start_pos,
                        const int end_pos, const double inv_frc_scale, const llint* xfrc,
                        const llint* yfrc, const llint* zfrc, const llint* xprv, const llint* yprv,
                        const llint* zprv) {
  double gg  = 0.0;
  double dgg = 0.0;
  for (int tpos = start_pos + threadIdx.x; tpos < end_pos; tpos += blockDim.x) {
    const double dpx = (double)(xprv[tpos]) * inv_frc_scale;
    const double dpy = (double)(yprv[tpos]) * inv_frc_scale;
    const double dpz = (double)(zprv[tpos]) * inv_frc_scale;
    gg += (dpx * dpx) + (dpy * dpy) + (dpz * dpz);
    const double dfx = (double)(xfrc[tpos]) * inv_frc_scale;
    const double dfy = (double)(yfrc[tpos]) * inv_frc_scale;
    const double dfz = (double)(zfrc[tpos]) * inv_frc_scale;
    const double ddx = dfx - dpx;
    const double ddy = dfy - dpy;
    const double ddz = dfz - dpz;
    dgg += (ddx * dfx) + (ddy * dfy) + (ddz * dfx);
  }
  WARP_REDUCE_DOWN(gg);
  WARP_REDUCE_DOWN(dgg);
  if (lane_idx == 0) {
    const size_t warp_idx = (threadIdx.x >> warp_bits);
    gg_collector[warp_idx]  =  gg;
    dgg_collector[warp_idx] = dgg;
  }  
}

__device__ __forceinline__
void conjGradCoreGather(double* gg_collector, double* dgg_collector, const int start_pos,
                        const int end_pos, const double inv_frc_scale, const llint* xfrc,
                        const int* xfrc_ovrf, const llint* yfrc, const int* yfrc_ovrf,
                        const llint* zfrc, const int* zfrc_ovrf, const llint* xprv,
                        const int* xprv_ovrf, const llint* yprv, const int* yprv_ovrf,
                        const llint* zprv, const int* zprv_ovrf) {
  double gg  = 0.0;
  double dgg = 0.0;
  for (int tpos = start_pos + threadIdx.x; tpos < end_pos; tpos += blockDim.x) {
    const double dpx = (((double)(xprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(xprv[tpos])) * inv_frc_scale;
    const double dpy = (((double)(yprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(yprv[tpos])) * inv_frc_scale;
    const double dpz = (((double)(zprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(zprv[tpos])) * inv_frc_scale;
    gg += (dpx * dpx) + (dpy * dpy) + (dpz * dpz);
    const double dfx = (((double)(xfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(xfrc[tpos])) * inv_frc_scale;
    const double dfy = (((double)(yfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(yfrc[tpos])) * inv_frc_scale;
    const double dfz = (((double)(zfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(zfrc[tpos])) * inv_frc_scale;
    dgg += (ddx * dfx) + (ddy * dfy) + (ddz * dfx);
  }
  WARP_REDUCE_DOWN(gg);
  WARP_REDUCE_DOWN(dgg);
  if (lane_idx == 0) {
    const size_t warp_idx = (threadIdx.x >> warp_bits);
    gg_collector[warp_idx]  =  gg;
    dgg_collector[warp_idx] = dgg;
  }
}

//-------------------------------------------------------------------------------------------------
// Distribute the results of the conjugate gradient calculation to normalize the forces and update
// prior holdings for the next cycle.
//
// Arguments:
//   gam:              Contribution factor applied to prior conjugate gradient results
//   start_pos:        First atom affected by this work unit in the following global arrays
//   end_pos:          Upper limit of atoms affected by this work unit
//   xfrc:             Current forces on particles in the Cartesian X direction
//   yfrc:             Current forces on particles in the Cartesian Y direction
//   zfrc:             Current forces on particles in the Cartesian Z direction
//   [x,y,z]frc_ovrf:  Overflow in X, Y, or Z forces
//   xprv:             Prior iteration forces on particles in the Cartesian X direction
//   yprv:             Prior iteration forces on particles in the Cartesian Y direction
//   zprv:             Prior iteration forces on particles in the Cartesian Z direction
//   [x,y,z]prv_ovrf:  Overflow in prior iteration X, Y, or Z forces
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
void conjGradScatter(const double gam, const int start_pos, const int end_pos, const llint* xfrc,
                     const llint* yfrc, const llint* zfrc, const llint* xprv, const llint* yprv,
                     const llint* zprv) {
  for (int tpos = atom_start_pos + threadIdx.x; tpos < atom_end_pos; tpos += blockDim.x) {

    // Because each work unit handles an exclusive subset of the atoms, the next __syncthreads()
    // command will be effective to ensure that the global memory manipulations handled here
    // are safe for read / write / read access.
    const llint ifx = xfrc[tpos];
    const llint ify = yfrc[tpos];
    const llint ifz = zfrc[tpos];
    xprv[tpos] = ifx;
    yprv[tpos] = ify;
    zprv[tpos] = ifz;
    const llint cg_x = __double2ll_rn((double)(ifx) + (gam * (double)(x_cg_temp[tpos])));
    const llint cg_y = __double2ll_rn((double)(ify) + (gam * (double)(y_cg_temp[tpos])));
    const llint cg_z = __double2ll_rn((double)(ifz) + (gam * (double)(z_cg_temp[tpos])));
    x_cg_temp[tpos] = cg_x;
    y_cg_temp[tpos] = cg_y;
    z_cg_temp[tpos] = cg_z;
    xfrc[tpos] = cg_x;
    yfrc[tpos] = cg_y;
    zfrc[tpos] = cg_z;
  }
}

__device__ __forceinline__
void conjGradScatter(const double gam, const int start_pos, const int end_pos, const llint* xfrc,
                     const int* xfrc_ovrf, const llint* yfrc, const int* yfrc_ovrf,
                     const llint* zfrc, const int* zfrc_ovrf, const llint* xprv,
                     const int* xprv_ovrf, const llint* yprv, const int* yprv_ovrf,
                     const llint* zprv, const int* zprv_ovrf, const llint* x_cg_temp,
                     const int* x_cg_temp_ovrf, const llint* y_cg_temp, const int* y_cg_temp_ovrf,
                     const llint* z_cg_temp, const int* z_cg_temp_ovrf) {
  for (int tpos = atom_start_pos + threadIdx.x; tpos < atom_end_pos; tpos += blockDim.x) {
    const llint ifx = xfrc[tpos];
    const llint ify = yfrc[tpos];
    const llint ifz = zfrc[tpos];
    xprv[tpos] = ifx;
    yprv[tpos] = ify;
    zprv[tpos] = ifz;
    const int ifx_ovrf = xfrc_ovrf[tpos];
    const int ify_ovrf = yfrc_ovrf[tpos];
    const int ifz_ovrf = zfrc_ovrf[tpos];
    xprv_ovrf[tpos] = ifx_ovrf;
    yprv_ovrf[tpos] = ify_ovrf;
    zprv_ovrf[tpos] = ifz_ovrf;
    const double  fx_part = ((double)(ifx_ovrf) * max_llint_accumulation) + (double)(ifx);
    const double cgx_part = ((double)(x_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(x_cg_temp[tpos]);
    const double  fy_part = ((double)(ify_ovrf) * max_llint_accumulation) + (double)(ify);
    const double cgy_part = ((double)(y_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(y_cg_temp[tpos]);
    const double  fz_part = ((double)(ifz_ovrf) * max_llint_accumulation) + (double)(ifz);
    const double cgz_part = ((double)(z_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(z_cg_temp[tpos]);
    const int95_t cg_x = convertSplitFixedPrecision95(fx_part + cgx_part);
    const int95_t cg_y = convertSplitFixedPrecision95(fy_part + cgy_part);
    const int95_t cg_z = convertSplitFixedPrecision95(fz_part + cgz_part);
    x_cg_temp[tpos] = cg_x.x;
    y_cg_temp[tpos] = cg_y.x;
    z_cg_temp[tpos] = cg_z.x;
    xfrc[tpos] = cg_x.x;
    yfrc[tpos] = cg_y.x;
    zfrc[tpos] = cg_z.x;
    x_cg_temp_ovrf[tpos] = cg_x.y;
    y_cg_temp_ovrf[tpos] = cg_y.y;
    z_cg_temp_ovrf[tpos] = cg_z.y;
    xfrc_ovrf[tpos] = cg_x.y;
    yfrc_ovrf[tpos] = cg_y.y;
    zfrc_ovrf[tpos] = cg_z.y;
  }
}

// Single-precision floating point conjugate gradient definitions.  The single-precision form
// still does its accumulation in double-precision, but does not store its data in the extended
// fixed-precision format which the double-precision forms of the kernels read and write.
#define KGATHER_NAME     kfgtConjGrad
#define KSCATTER_NAME    kfscConjGrad
#define KALLREDUCE_NAME  kfrdConjGrad
#include "conjugate_gradient.cui"
#undef KGATHER_NAME
#undef KSCATTER_NAME
#undef KALLREDUCE_NAME

// Double-precision floating point conjugate gradient definitions
#define KGATHER_NAME     kdgtConjGrad
#define KSCATTER_NAME    kdscConjGrad
#define KALLREDUCE_NAME  kdrdConjGrad
#include "conjugate_gradient.cui"
#undef KGATHER_NAME
#undef KSCATTER_NAME
#undef KALLREDUCE_NAME

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryReductionKernelRequirements(const PrecisionModel prec,
                                                           const ReductionGoal purpose,
                                                           const ReductionStage process) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (purpose) {
    case ReductionGoal::NORMALIZE:
    case ReductionGoal::CENTER_ON_ZERO:
      break;
    case ReductionGoal::CONJUGATE_GRADIENT:
      switch (process) {
      case ReductionStage::GATHER:
        if (cudaFuncGetAttributes(&result, kdgtConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdgtConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::SCATTER:
        if (cudaFuncGetAttributes(&result, kdscConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdscConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::ALL_REDUCE:
        if (cudaFuncGetAttributes(&result, kdrdConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdrdConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (purpose) {
    case ReductionGoal::NORMALIZE:
    case ReductionGoal::CENTER_ON_ZERO:
      break;
    case ReductionGoal::CONJUGATE_GRADIENT:
      switch (process) {
      case ReductionStage::GATHER:
        if (cudaFuncGetAttributes(&result, kfgtConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfgtConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::SCATTER:
        if (cudaFuncGetAttributes(&result, kfscConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfscConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::ALL_REDUCE:
        if (cudaFuncGetAttributes(&result, kfrdConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfrdConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      }
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradientDp(const ReductionKit redk, ConjGradSubstrate cgsbs,
                                      MMControlKit<double> *ctrl, const KernelManager &launcher) {

  // All conjugate gradient kernels take the same launch parameters.
  const int2 bt = launcher.getReductionKernelDims(PrecisionModel::DOUBLE,
                                                  ReductionGoal::CONJUGATE_GRADIENT,
                                                  ReductionStage::ALL_REDUCE);
  switch (redk.rps) {
  case RdwuPerSystem::ONE:
    kdrdConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    break;
  case RdwuPerSystem::MULTIPLE:
    kdgtConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    kdscConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradientDp(const AtomGraphSynthesis poly_ag,
                                      PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                                      MMControlKit<double> *ctrl, const KernelManager &launcher) {
  ReductionKit redk(poly_ag, HybridTargetLevel::DEVICE);
  ConjGradSubstrate cgsbs(poly_ps, rbg, HybridTargetLevel::DEVICE);
  MMControlKit<double> ctrl = mmctrl->dpData();
  launchConjugateGradientDp(redk, cgsbs, &ctrl, launcher);
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradientSp(const ReductionKit redk, ConjGradSubstrate cgsbs,
                                      MMControlKit<float> *ctrl, const KernelManager &launcher) {

  // All conjugate gradient kernels take the same launch parameters.
  const int2 bt = launcher.getReductionKernelDims(PrecisionModel::SINGLE,
                                                  ReductionGoal::CONJUGATE_GRADIENT,
                                                  ReductionStage::ALL_REDUCE);
  switch (redk.rps) {
  case RdwuPerSystem::ONE:
    kfrdConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    break;
  case RdwuPerSystem::MULTIPLE:
    kfgtConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    kfscConjGrad<<<bt.x, bt.y>>>(redk, cgsbs, *ctrl);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradientSp(const AtomGraphSynthesis poly_ag,
                                      PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                                      MolecularMechanicsControls *mmctrl,
                                      const KernelManager &launcher) {
  ReductionKit redk(poly_ag, HybridTargetLevel::DEVICE);
  ConjGradSubstrate cgsbs(poly_ps, rbg, HybridTargetLevel::DEVICE);
  MMControlKit<float> ctrl = mmctrl->spData();
  launchConjugateGradientDp(redk, cgsbs, &ctrl, launcher);
}

} // namespace synthesis
} // namespace omni
