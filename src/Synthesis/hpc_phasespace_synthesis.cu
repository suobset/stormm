// -*-c++-*-
#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Reporting/error_format.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_phasespace_synthesis.h"
#include "hpc_phasespace_synthesis.cuh"

namespace stormm {
namespace synthesis {

using data_types::int95_t;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::force_scale_nonoverflow_bits;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
  
#include "Numerics/accumulation.cui"
#include "Math/rounding.cui"
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSystemTransfer(PsSynthesisWriter destination, PsSynthesisWriter source, const int low_index,
                const int high_index, const TrajectoryKind material) {
  __shared__ int sread_start, sread_last_system, sread_last_offset;

  // Pre-compute some constants for transferring positions
  const int mtrx_stride = ((9 + warp_bits_mask_int) / warp_size_int) * warp_size_int;
  const int mtrx_read_start = low_index * mtrx_stride;
  const int mtrx_read_end   = (high_index + 1) * mtrx_stride;
  const int bdim_stride = ((6 + warp_bits_mask_int) / warp_size_int) * warp_size_int;
  const int bdim_read_start = low_index * bdim_stride;
  const int bdim_read_end   = (high_index + 1) * bdim_stride;

  // Assume that the bounds arrays for each system have been copied to the device, where they are
  // visible with the device view pointers.
  if (threadIdx.x == 0) {
    sread_start = source.atom_starts[low_index];
  }
  else if (threadIdx.x == warp_size_int) {
    sread_last_system = source.atom_starts[high_index - 1];
  }
  else if (threadIdx.x == 2 * warp_size_int) {
    sread_last_offset = source.atom_counts[high_index - 1];
  }
  __syncthreads();
  const int atom_read_start = sread_start;
  const int atom_read_end   = sread_last_system + sread_last_offset;

  // The amount of total information transfer is now known.  Assess the number of available warps,
  // assign each an index, and determine a way to subdivide the workload among the available warps.
  int atom_work_units, mtrx_work_units, bdim_work_units, box_work_units, box_work_per_system;
  atom_work_units = (atom_read_end - atom_read_start + warp_bits_mask_int) / warp_size_int;
  switch (material) {
  case TrajectoryKind::POSITIONS:
    {
      mtrx_work_units = mtrx_stride / warp_size_int;
      bdim_work_units = bdim_stride / warp_size_int;
      box_work_per_system = (5 * mtrx_work_units) + (2 * bdim_work_units);
      box_work_units = box_work_per_system * (high_index - low_index);
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    atom_work_units *= 3;
    box_work_per_system = 0;
    box_work_units = 0;
    break;
  }
  const float total_work_units = (float)(atom_work_units + box_work_units);
  const float total_warps = (float)((blockDim.x / warp_size_int) * gridDim.x);
  int tmp_box_warps = (float)(box_work_units) / total_work_units * (float)total_warps;

  // This setup will not fail on any known GPU, but if the GPU were exceptionally tiny (fewer than
  // seven warps per block and only one block available) it would fail.  But, the GPUs are only
  // getting bigger.
  const int box_warps = (box_work_per_system > 0) ?
                        ((tmp_box_warps + (tmp_box_warps == 0 && box_work_units > 0) +
                          box_work_per_system - 1) / box_work_per_system) * box_work_per_system :
                        0;
  const int atom_warps = total_warps - box_warps;
  
  // Warps now serve atom or box data transfer.
  const int warp_idx = (threadIdx.x >> warp_bits) + ((blockDim.x >> warp_bits) * blockIdx.x);
  const int tgx = (threadIdx.x & warp_bits_mask_int);
  switch (material) {
  case TrajectoryKind::POSITIONS:
    {
      if (warp_idx < atom_warps) {
        int read_pos = atom_read_start + (warp_idx * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          destination.xcrd[read_pos] = source.xcrd[read_pos];
          destination.ycrd[read_pos] = source.ycrd[read_pos];
          destination.zcrd[read_pos] = source.zcrd[read_pos];
          destination.xprv[read_pos] = source.xprv[read_pos];
          destination.yprv[read_pos] = source.yprv[read_pos];
          destination.zprv[read_pos] = source.zprv[read_pos];
          if (source.gpos_bits > globalpos_scale_nonoverflow_bits) {
            destination.xcrd_ovrf[read_pos] = source.xcrd_ovrf[read_pos];
            destination.ycrd_ovrf[read_pos] = source.ycrd_ovrf[read_pos];
            destination.zcrd_ovrf[read_pos] = source.zcrd_ovrf[read_pos];
            destination.xprv_ovrf[read_pos] = source.xprv_ovrf[read_pos];
            destination.yprv_ovrf[read_pos] = source.yprv_ovrf[read_pos];
            destination.zprv_ovrf[read_pos] = source.zprv_ovrf[read_pos];
          }
          read_pos += atom_warps * warp_size_int;
        }
      }
      else {
        const int warp_nature = (warp_idx - atom_warps) % box_work_per_system;
        const int box_warp_sets = box_warps / box_work_per_system;
        const int mtrx_warps = mtrx_work_units * box_warp_sets;
        const int bdim_warps = bdim_work_units * box_warp_sets;
        const int box_warp_idx = warp_idx - atom_warps;
        if (warp_nature < mtrx_warps) {
          int read_pos = mtrx_read_start + (box_warp_idx * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.boxvecs[read_pos] = source.boxvecs[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 2 * mtrx_warps) {
          int read_pos = mtrx_read_start + ((box_warp_idx - mtrx_warps) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.umat[read_pos] = source.umat[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 3 * mtrx_work_units) {
          int read_pos = mtrx_read_start +
                         ((box_warp_idx - (2 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.invu[read_pos] = source.invu[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 4 * mtrx_work_units) {
          int read_pos = mtrx_read_start +
                         ((box_warp_idx - (3 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.sp_umat[read_pos] = source.sp_umat[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 5 * mtrx_work_units) {
          int read_pos = mtrx_read_start +
                         ((box_warp_idx - (4 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.sp_invu[read_pos] = source.sp_invu[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < (5 * mtrx_work_units) + bdim_work_units) {
          int read_pos = bdim_read_start +
                         ((box_warp_idx - (5 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < bdim_read_end) {
            destination.boxdims[read_pos] = source.boxdims[read_pos];
            read_pos += bdim_warps * warp_size_int;
          }
        }
        else {
          int read_pos = bdim_read_start +
                         ((box_warp_idx - (5 * mtrx_warps) - bdim_warps) * warp_size_int) + tgx;
          while (read_pos < bdim_read_end) {
            destination.sp_boxdims[read_pos] = source.sp_boxdims[read_pos];
            read_pos += bdim_warps * warp_size_int;
          }
        }
      }
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    {
      const int pdim_warps = atom_warps / 3;
      if (warp_idx < pdim_warps) {
        int read_pos = atom_read_start + (warp_idx * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.xvel[read_pos] = source.xvel[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.xvel_ovrf[read_pos] = source.xvel_ovrf[read_pos];              
            }
          }
          else {
            destination.xfrc[read_pos] = source.xfrc[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.xfrc_ovrf[read_pos] = source.xfrc_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 2 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - pdim_warps) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.yvel[read_pos] = source.yvel[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.yvel_ovrf[read_pos] = source.yvel_ovrf[read_pos];              
            }
          }
          else {
            destination.yfrc[read_pos] = source.yfrc[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.yfrc_ovrf[read_pos] = source.yfrc_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 3 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - (2 * pdim_warps)) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.zvel[read_pos] = source.zvel[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.zvel_ovrf[read_pos] = source.zvel_ovrf[read_pos];              
            }
          }
          else {
            destination.zfrc[read_pos] = source.zfrc[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.zfrc_ovrf[read_pos] = source.zfrc_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
    }
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void systemTransfer(PsSynthesisWriter *destination, PsSynthesisWriter *source,
                           const TrajectoryKind kind, const int low_index,
                           const int high_index, const GpuDetails &gpu) {
  kSystemTransfer<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*destination, *source,
                                                                      low_index, high_index, kind);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyInitializeForces(PsSynthesisWriter psyw, const int index) {
  int minpos, maxpos;
  if (index < 0) {
    minpos = 0;
    maxpos = psyw.atom_starts[psyw.system_count - 1] + psyw.atom_counts[psyw.system_count - 1];
  }
  else {
    minpos = psyw.atom_starts[index];
    maxpos = minpos + psyw.atom_counts[index];
  }
  for (int pos = minpos + threadIdx.x + (blockDim.x * blockIdx.x);
       pos < maxpos; pos += gridDim.x * blockDim.x) {
    psyw.xfrc[pos] = 0LL;
    psyw.yfrc[pos] = 0LL;
    psyw.zfrc[pos] = 0LL;
    if (psyw.frc_bits > force_scale_nonoverflow_bits) {
      psyw.xfrc_ovrf[pos] = 0;
      psyw.yfrc_ovrf[pos] = 0;
      psyw.zfrc_ovrf[pos] = 0;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void psyInitializeForces(PsSynthesisWriter *psyw, const int index, const GpuDetails &gpu) {
  kPsyInitializeForces<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*psyw, index);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyPrimeConjugateGradient(PsSynthesisWriter psyw) {
  const int minpos = threadIdx.x + (blockDim.x * blockIdx.x);
  const int last_system = psyw.system_count - 1;
  const int maxpos = psyw.atom_starts[last_system] + psyw.atom_counts[last_system];
  for (int pos = minpos; pos < maxpos; pos += blockDim.x * gridDim.x) {
    psyw.xprv[pos] = psyw.xfrc[pos];
    psyw.yprv[pos] = psyw.yfrc[pos];
    psyw.zprv[pos] = psyw.zfrc[pos];
    psyw.xvel[pos] = 0LL;
    psyw.yvel[pos] = 0LL;
    psyw.zvel[pos] = 0LL;
    if (psyw.frc_bits > force_scale_nonoverflow_bits) {
      psyw.xprv_ovrf[pos] = psyw.xfrc_ovrf[pos];
      psyw.yprv_ovrf[pos] = psyw.yfrc_ovrf[pos];
      psyw.zprv_ovrf[pos] = psyw.zfrc_ovrf[pos];
      psyw.xvel_ovrf[pos] = 0;
      psyw.yvel_ovrf[pos] = 0;
      psyw.zvel_ovrf[pos] = 0;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void psyPrimeConjugateGradient(PsSynthesisWriter *psyw, const GpuDetails &gpu) {
  kPsyPrimeConjugateGradient<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*psyw);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyCopySystem(PhaseSpaceWriter psw, const PsSynthesisReader poly_psw, const int index) {
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t atom_offset = poly_psw.atom_starts[index];
  pos = splitFPToReal(psw.xcrd, pos,  0, &poly_psw.xcrd[atom_offset],
                      &poly_psw.xcrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.ycrd, pos,  1, &poly_psw.ycrd[atom_offset],
                      &poly_psw.ycrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.zcrd, pos,  2, &poly_psw.zcrd[atom_offset],
                      &poly_psw.zcrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.xprv, pos,  3, &poly_psw.xprv[atom_offset],
                      &poly_psw.xprv_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.yprv, pos,  4, &poly_psw.yprv[atom_offset],
                      &poly_psw.yprv_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.zprv, pos,  5, &poly_psw.zprv[atom_offset],
                      &poly_psw.zprv_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.xnxt, pos,  6, &poly_psw.xnxt[atom_offset],
                      &poly_psw.xnxt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.ynxt, pos,  7, &poly_psw.ynxt[atom_offset],
                      &poly_psw.ynxt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.znxt, pos,  8, &poly_psw.znxt[atom_offset],
                      &poly_psw.znxt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.xvel, pos,  9, &poly_psw.xvel[atom_offset],
                      &poly_psw.xvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.yvel, pos, 10, &poly_psw.yvel[atom_offset],
                      &poly_psw.yvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.zvel, pos, 11, &poly_psw.zvel[atom_offset],
                      &poly_psw.zvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vxprv, pos, 12, &poly_psw.vxprv[atom_offset],
                      &poly_psw.vxprv_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vyprv, pos, 13, &poly_psw.vyprv[atom_offset],
                      &poly_psw.vyprv_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vzprv, pos, 14, &poly_psw.vzprv[atom_offset],
                      &poly_psw.vzprv_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vxnxt, pos, 15, &poly_psw.vxnxt[atom_offset],
                      &poly_psw.vxnxt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vynxt, pos, 16, &poly_psw.vynxt[atom_offset],
                      &poly_psw.vynxt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vznxt, pos, 17, &poly_psw.vznxt[atom_offset],
                      &poly_psw.vznxt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.xfrc, pos, 18, &poly_psw.xfrc[atom_offset],
                      &poly_psw.xfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.yfrc, pos, 19, &poly_psw.yfrc[atom_offset],
                      &poly_psw.yfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.zfrc, pos, 20, &poly_psw.zfrc[atom_offset],
                      &poly_psw.zfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fxprv, pos, 21, &poly_psw.fxprv[atom_offset],
                      &poly_psw.fxprv_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fyprv, pos, 22, &poly_psw.fyprv[atom_offset],
                      &poly_psw.fyprv_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fzprv, pos, 23, &poly_psw.fzprv[atom_offset],
                      &poly_psw.fzprv_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fxnxt, pos, 24, &poly_psw.fxnxt[atom_offset],
                      &poly_psw.fxnxt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fynxt, pos, 25, &poly_psw.fynxt[atom_offset],
                      &poly_psw.fynxt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fznxt, pos, 26, &poly_psw.fznxt[atom_offset],
                      &poly_psw.fznxt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);

  // Copy the transformation matrices and box dimensions
  const int mtrx_offset = devcRoundUp(9, warp_size_int) * index;
  const int bdim_offset = devcRoundUp(6, warp_size_int) * index;
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (warp_idx == 0) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.umat[npos] = poly_psw.umat[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 1) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.invu[npos] = poly_psw.invu[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 2) {
    int npos = lane_idx;
    while (npos < 6) {
      psw.boxdim[npos] = poly_psw.boxdims[bdim_offset + npos];
      npos += warp_size_int;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::extractSystem(PhaseSpaceWriter *psw, const int index,
                                        const GpuDetails &gpu,
                                        const HybridTargetLevel origin) const {
  const PsSynthesisReader poly_psr = this->data(origin);
  kPsyCopySystem<<<gpu.getSMPCount(), large_block_size>>>(*psw, poly_psr, index);
  cudaDeviceSynchronize();
}
  
} // namespace synthesis
} // namespace stormm
