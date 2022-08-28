// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/energy_abstracts.h"
#include "Synthesis/implicit_solvent_workspace.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_nonbonded_potential.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using data_types::int95_t;
using numerics::chooseAccumulationMethod;
using numerics::max_int_accumulation_f;
using numerics::max_int_accumulation_ll;
using numerics::max_llint_accumulation;
using numerics::max_llint_accumulation_f;
using synthesis::NbwuKind;
using synthesis::small_block_max_imports;
using synthesis::small_block_max_atoms;
using synthesis::tile_groups_wu_abstract_length;
using topology::ImplicitSolventModel;
  
//-------------------------------------------------------------------------------------------------
// Get the number atoms in a particular tile stretch.
//
// Arguments:
//   nbwu_map:  Details of the non-bonded work unit, condensed into a simple array of integers
//   pos:       Thread position in the list of atoms
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ int getTileSideAtomCount(const int* nbwu_map, const int pos) {
  const int key_idx  = pos / 4;
  const int key_slot = pos - (key_idx * 4);
  return ((nbwu_map[small_block_max_imports + 1 + key_idx] >> (8 * key_slot)) & 0xff);
}

#include "Numerics/accumulation.cui"
#include "Math/rounding.cui"

//-------------------------------------------------------------------------------------------------
// Load coordinates relating to atoms in a non-bonded tile suitable for isolated boundary
// conditions.
//
// Overloaded:
//   - Work with an array of long long integers appropriate for single-precision arithmetic
//   - Work with dual arrays of long long int and int types, appropriate for double-precision
//     arithmetic
//
// Arguments:
//   pos:              Position in the tile list (not the atom list)
//   import_count:     Number of groups of atoms imported to populate one side of one or more tiles
//   iter:             Number of passes made by this or related routines (incrementation of
//                     iter is essential to maintain the correct procession through all loads)
//   nbwu_map:         Non-bonded work unit details            
//   read_crd:         Array of coordinates to read from
//   write_crd:        Array of coordinates to write into
//   read_crd_ovrf:    Overflow buffers for coordinates to be read
//   write_crd_ovrf:   Overflow buffers for local copies of coordinates
//   sh_tile_cog:      Array holding mean values of the positions of each imported atom group (this
//                     later expedites computing the center of geometry for each complete tile)
//   gpos_scale:       Scaling factor for coordinates in the fixed-precision representation (this
//                     is needed only to place dummy atom coordinates for blank slots of a tile)
//-------------------------------------------------------------------------------------------------
__device__ int loadTileCoordinates(const int pos, const int iter, const int* nbwu_map,
                                   const llint* read_crd, llint* write_crd, float* sh_tile_cog,
                                   const float gpos_scale) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    float fval;
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        const llint ival = __ldcs(&read_crd[read_idx]);
        fval = (float)(ival);
        __stwb(&write_crd[write_idx], ival);
      }
      else {
        fval = (float)(0.0);
        __stwb(&write_crd[write_idx], (128 * rel_pos * tile_lane_idx) * gpos_scale);
      }
    }
    else {
      fval = (float)(0.0);
    }
    for (int i = half_tile_length; i > 0; i >>= 1) {
      fval += SHFL_DOWN(fval, i);
    }
    if (tile_lane_idx == 0 && rel_pos < import_count) {
      sh_tile_cog[rel_pos] = fval;
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int loadTileCoordinates(const int pos, const int iter, const int* nbwu_map,
                                   const llint* read_crd, llint* write_crd,
                                   const int* read_crd_ovrf, int* write_crd_ovrf,
                                   double* sh_tile_cog, const double gpos_scale) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    double fval;
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        const llint ival = __ldcs(&read_crd[read_idx]);
        fval = (double)(ival);
        __stwb(&write_crd[write_idx], ival);
        const int ival_ovrf = __ldcs(&read_crd_ovrf[read_idx]);
        fval += (double)(ival_ovrf) * max_llint_accumulation;
        __stwb(&write_crd_ovrf[write_idx], ival_ovrf);
      }
      else {
        fval = 0.0;
        const int95_t fake_val = doubleToInt95((128 * rel_pos * tile_lane_idx) * gpos_scale);
        __stwb(&write_crd[write_idx], fake_val.x);
        __stwb(&write_crd_ovrf[write_idx], fake_val.y);
      }
    }
    else {
      fval = 0.0;
    }
    for (int i = half_tile_length; i > 0; i >>= 1) {
      fval += SHFL_DOWN(fval, i);
    }
    if (tile_lane_idx == 0 && rel_pos < import_count) {
      sh_tile_cog[rel_pos] = fval;
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

//-------------------------------------------------------------------------------------------------
// Load scalar values (integral or real) from global memory into local arrays, on a tile-by-tile
// basis for non-bonded work units involving isolated systems with all-to-all interaction matrices.
//
// Overloaded:
//   - Copy the values directly
//   - Fold in a scalar multiple
//   - Fold in a scalar addition
//
// Parameter descriptors follow from loadTileCoordinates() above, with alterations:
//   read_array:   Generic array of (global) information to read from
//   write_array:  Generic (local) array of information to write 
//-------------------------------------------------------------------------------------------------
template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T* write_array) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]);
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T* write_array, T multiplier) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]) * multiplier;
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T increment, T* write_array) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]) + increment;
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

//-------------------------------------------------------------------------------------------------
// Write information about the atoms in tile groups back to global accumulators.  Relevant for
// systems with all-to-all interactions in isolated boundary conditions.
//
// Overloaded:
//   - Accept various combinations of single- or double-integer local accumulators to contribute
//     to the implied single- or double-integer global accumulators
//
// Arguments:
//   pos:                   Position in the tile list (not the atom list)
//   iter:                  Number of passes made by this or related routines (incrementation of
//                          iter is essential to maintain the correct procession through all loads)
//   nbwu_map:              Non-bonded work unit details
//   tile_prop:             Primary local accumulator for the tile-based computed property
//   tile_prop_ovrf:        Local overflow accumulator for the tile-based computed property
//   gbl_accumulator:       Primary (or, perhaps lone) global accumulator for the computed property
//   gbl_accumulator_ovrf:  Overflow global accumulator for the computed property
//-------------------------------------------------------------------------------------------------
__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const int* tile_prop, const int* tile_prop_ovrf,
                                      llint* gbl_accumulator) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        llint itp = tile_prop_ovrf[read_idx];
        itp *= max_int_accumulation_ll;
        itp += tile_prop[read_idx];
        atomicAdd((ullint*)&gbl_accumulator[write_idx], (ullint)(itp));
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const llint* tile_prop, llint* gbl_accumulator) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        atomicAdd((ullint*)&gbl_accumulator[write_idx], (ullint)(tile_prop[read_idx]));
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const llint* tile_prop, const int* tile_prop_ovrf,
                                      llint* gbl_accumulator, int* gbl_accumulator_ovrf) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        atomicSplit(tile_prop[read_idx], tile_prop_ovrf[read_idx], write_idx, gbl_accumulator,
                    gbl_accumulator_ovrf);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

// Single-precision non-bonded kernel floating point definitions
#define TCALC float
#  define TCALC2 float2
#  define TCALC_IS_SINGLE
#  if (__CUDA_ARCH__ >= 750) && (__CUDA_ARCH__ < 800)
#    define GBRADII_KERNEL_BLOCKS_MULTIPLIER 4
#    define NONBOND_KERNEL_BLOCKS_MULTIPLIER 4
#  else
#    define GBRADII_KERNEL_BLOCKS_MULTIPLIER 5
#    define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#  endif
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define LOG_FUNC  logf
#  define EXP_FUNC  expf
#  define TANH_FUNC tanhf
#  define SPLIT_FORCE_ACCUMULATION
#    define KERNEL_NAME ktgfsCalculateGBRadii
#      include "gbradii_tilegroups.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME ktgfsCalculateGBDerivatives
#      include "gbderivative_tilegroups.cui"
#    undef KERNEL_NAME
#  undef SPLIT_FORCE_ACCUMULATION
#  define KERNEL_NAME ktgfCalculateGBRadii
#    include "gbradii_tilegroups.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME ktgfCalculateGBDerivatives
#    include "gbderivative_tilegroups.cui"
#  undef KERNEL_NAME
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define COMPUTE_ENERGY
#        define KERNEL_NAME ktgfsNonbondedForceEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#      undef COMPUTE_ENERGY
#      define KERNEL_NAME ktgfsNonbondedForce
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef SPLIT_FORCE_ACCUMULATION
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgfNonbondedForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgfNonbondedForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgfNonbondedEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef LOG_FUNC
#  undef EXP_FUNC
#  undef TANH_FUNC
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  undef GBRADII_KERNEL_BLOCKS_MULTIPLIER
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#undef TCALC

// Double-precision non-bonded kernel floating point definitions
#define TCALC double
#  define TCALC2 double2
#  define SPLIT_FORCE_ACCUMULATION
#  define GBRADII_KERNEL_BLOCKS_MULTIPLIER 3
#  define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define LOG_FUNC  log
#  define EXP_FUNC  exp
#  define TANH_FUNC tanh
#  define KERNEL_NAME ktgdsCalculateGBRadii
#    include "gbradii_tilegroups.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME ktgdsCalculateGBDerivatives
#    include "gbderivative_tilegroups.cui"
#  undef KERNEL_NAME
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgdsNonbondedForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgdsNonbondedForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgdNonbondedEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef LOG_FUNC
#  undef EXP_FUNC
#  undef TANH_FUNC
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  undef SPLIT_FORCE_ACCUMULATION
#  undef TCALC2
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void nonbondedKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedForce, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedForce __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgfNonbondedForceEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgfNonbondedForceEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdsNonbondedForce, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdsNonbondedForce __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdNonbondedEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdNonbondedEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(ktgdsNonbondedForceEnergy, sms_eight) != cudaSuccess) {
    rtErr("Error setting ktgdsNonbondedForceEnergy __shared__ memory bank size to eight bytes.",
          "nonbondedKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryNonbondedKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                 const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                                 const AccumulationMethod acc_meth) {

  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes attr;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForceEnergy) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdsNonbondedForceEnergy.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case EvaluateEnergy::NO:
        if (cudaFuncGetAttributes(&attr, ktgdsNonbondedForce) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdsNonbondedForce.",
                "queryNonbondedKernelRequirements");
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      if (cudaFuncGetAttributes(&attr, ktgdNonbondedEnergy) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ktgdNonbondedEnergy.",
              "queryValenceKernelRequirements");
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (kind) {
      case NbwuKind::TILE_GROUPS:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfsNonbondedForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case AccumulationMethod::WHOLE:
            if (cudaFuncGetAttributes(&attr, ktgfNonbondedForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfNonbondedForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            if (cudaFuncGetAttributes(&attr, ktgfsNonbondedForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfsNonbondedForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case AccumulationMethod::WHOLE:
            if (cudaFuncGetAttributes(&attr, ktgfNonbondedForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfNonbondedForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:
        if (cudaFuncGetAttributes(&attr, ktgfNonbondedEnergy) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfNonbondedEnergy.",
                "queryValenceKernelRequirements");
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
      break;
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  }
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryBornRadiiKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                 const AccumulationMethod acc_meth) {
  cudaFuncAttributes attr;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBRadii) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ktgdCalculateGBRadii.",
              "queryNonbondedKernelRequirements");
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBRadii) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfsCalculateGBRadii.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case AccumulationMethod::WHOLE:
        if (cudaFuncGetAttributes(&attr, ktgfCalculateGBRadii) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfCalculateGBRadii.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  }
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryBornDerivativeKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                      const AccumulationMethod acc_meth) {
  cudaFuncAttributes attr;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBDerivatives) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel ktgdCalculateGBDerivatives.",
              "queryNonbondedKernelRequirements");
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBDerivatives) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfsCalculateGBDerivatives.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case AccumulationMethod::WHOLE:
        if (cudaFuncGetAttributes(&attr, ktgfCalculateGBDerivatives) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgfCalculateGBDerivatives.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
    break;
  }
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern void launchBornRadiiCalculation(const NbwuKind kind,
                                       const SyNonbondedKit<double, double2> &poly_nbk,
                                       MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                                       CacheResourceKit<double> *gmem_r,
                                       ISWorkspaceKit<double> *iswk, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    ktgdsCalculateGBRadii<<<bt.x, bt.y>>>(poly_nbk, *ctrl, *poly_psw, *gmem_r, *iswk);
    break;
  case NbwuKind::SUPERTILES:
    break;
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchBornRadiiCalculation(const NbwuKind kind,
                                       const SyNonbondedKit<float, float2> &poly_nbk,
                                       MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                                       CacheResourceKit<float> *gmem_r,
                                       ISWorkspaceKit<float> *iswk,
                                       const AccumulationMethod acc_meth, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      ktgfsCalculateGBRadii<<<bt.x, bt.y>>>(poly_nbk, *ctrl, *poly_psw, *gmem_r, *iswk);
      break;
    case AccumulationMethod::WHOLE:
      ktgfCalculateGBRadii<<<bt.x, bt.y>>>(poly_nbk, *ctrl, *poly_psw, *gmem_r, *iswk);
      break;
    case AccumulationMethod::AUTOMATIC:
      break;
    }
    break;
  case NbwuKind::SUPERTILES:
    break;
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchBornRadiiCalculation(const PrecisionModel prec,
                                       const AtomGraphSynthesis &poly_ag,
                                       MolecularMechanicsControls *mmctrl,
                                       PhaseSpaceSynthesis *poly_ps, CacheResource *tb_space,
                                       ImplicitSolventWorkspace *isw,
                                       const KernelManager &launcher,
                                       const AccumulationMethod acc_meth) {
  AccumulationMethod actual_acc_meth;
  switch (acc_meth) {
  case AccumulationMethod::SPLIT:
  case AccumulationMethod::WHOLE:
    actual_acc_meth = acc_meth;
    break;
  case AccumulationMethod::AUTOMATIC:
    actual_acc_meth = chooseAccumulationMethod(isw->getFixedPrecisionBits());
    break;
  }
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const NbwuKind kind = poly_ag.getNonbondedWorkType();
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const int2 bt = launcher.getBornRadiiKernelDims(prec, kind, actual_acc_meth);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(devc_tier);
      ISWorkspaceKit<double> iswk = isw->dpData(devc_tier);
      launchBornRadiiCalculation(kind, poly_nbk, &ctrl, &poly_psw, &gmem_r, &iswk, bt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(devc_tier);
      ISWorkspaceKit<float> iswk = isw->spData(devc_tier);
      launchBornRadiiCalculation(kind, poly_nbk, &ctrl, &poly_psw, &gmem_r, &iswk, actual_acc_meth,
                                 bt);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<double, double2> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<double> *ctrl,
                            PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                            CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        ktgdsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                  *gmem_r);
        break;
      case EvaluateEnergy::NO:
        ktgdsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
        break;
      }
      break;
    case EvaluateForce::NO:
      ktgdNonbondedEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *gmem_r);
      break;
    }
    break;
  case NbwuKind::SUPERTILES:
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<float, float2> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<float> *ctrl,
                            PsSynthesisWriter *poly_psw, ScoreCardWriter *scw,
                            CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy,
                            const AccumulationMethod force_sum, const int2 bt) {
  switch (kind) {
  case NbwuKind::TILE_GROUPS:
    switch (eval_force) {
    case EvaluateForce::YES:
      switch (force_sum) {
      case AccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                    *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          ktgfNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                   *gmem_r);
          break;
        case EvaluateEnergy::NO:
          ktgfNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        if (poly_psw->frc_bits <= 24) {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgfsNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                      *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgfsNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
            break;
          }
        }
        else {
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgfNonbondedForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                     *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgfNonbondedForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *gmem_r);
            break;
          }
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      ktgfNonbondedEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *gmem_r);
      break;
    }
    break;
  case NbwuKind::SUPERTILES:
  case NbwuKind::HONEYCOMB:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            ScoreCard *sc, CacheResource *tb_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy,
                            const AccumulationMethod force_sum,
                            const KernelManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const int2 bt = launcher.getNonbondedKernelDims(prec, nb_work_type, eval_force, eval_energy,
                                                  force_sum);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r,
                      eval_force, eval_energy, bt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &scw, &gmem_r,
                      eval_force, eval_energy, force_sum, bt);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            ScoreCard *sc, CacheResource *tb_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const KernelManager &launcher) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, sc, tb_space, eval_force, eval_energy,
                    AccumulationMethod::SPLIT, launcher);
  }
  else {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, sc, tb_space, eval_force, eval_energy,
                    AccumulationMethod::WHOLE, launcher);
  }
}

} // namespace energy
} // namespace stormm
