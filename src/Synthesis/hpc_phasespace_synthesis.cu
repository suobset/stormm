// -*-c++-*-
#include <vector>
#ifdef OMNI_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "Constants/scaling.h"
#include "Reporting/error_format.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_phasespace_synthesis.h"
#include "hpc_phasespace_synthesis.cuh"

namespace omni {
namespace synthesis {

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
          }
          else {
            destination.xfrc[read_pos] = source.xfrc[read_pos];
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 2 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - pdim_warps) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.yvel[read_pos] = source.yvel[read_pos];
          }
          else {
            destination.yfrc[read_pos] = source.yfrc[read_pos];
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 3 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - (2 * pdim_warps)) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.zvel[read_pos] = source.zvel[read_pos];
          }
          else {
            destination.zfrc[read_pos] = source.zfrc[read_pos];
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
  }
}

//-------------------------------------------------------------------------------------------------
extern void psyInitializeForces(PsSynthesisWriter *psyw, const int index, const GpuDetails &gpu) {
  kPsyInitializeForces<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*psyw, index);
}

} // namespace synthesis
} // namespace omni
