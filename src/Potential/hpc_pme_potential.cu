// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/synthesis_abstracts.h"
#include "hpc_pme_potential.h"

namespace stormm {
namespace energy {

using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using synthesis::SyNonbondedKit;
  
#include "Accelerator/syncwarp.cui"
#include "Numerics/accumulation.cui"

/// \brief A device function for evaluating a local exclusion mask based on a topology index delta.
///        This is equivalent to evaluateLocalMask() in the LocalExclusionMask library (see
///        local_exclusionmask.h).
///
__device__ __forceinline__ bool devcEvaluateLocalMask(int atom_i, int atom_j, ullint prof,
                                                      const uint2* secondary_ptr) {
  const int del_ij = atom_j - atom_i;
  switch (prof & lmask_mode_bitmask) {
  case lmask_mode_a:
    {
      return (abs(del_ij) <= lmask_long_local_span &&
              ((prof >> (lmask_long_local_span + del_ij)) & 0x1));
    }
    break;
  case lmask_mode_b:
    {
      const int abs_dij = abs(del_ij);
      if (abs_dij > lmask_b_max_reach) {
        return false;
      }
      else if (abs_dij <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }
      else if (del_ij > 0) {
        const int upper_shft = ((prof & lmask_b_upper_shft) >> lmask_b_upper_shft_pos);
        const int upper_mask_start = lmask_short_local_span + upper_shft;
        const int rel_ij = del_ij - upper_mask_start;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_upper_mask_pos)) & 0x1));
      }
      else {

        // The only remaining case is that del_ij < 0
        const int lower_shft = ((prof & lmask_b_lower_shft) >> lmask_b_lower_shft_pos);
        const int lower_mask_start = -lmask_short_local_span - lower_shft - lmask_short_extra_span;
        const int rel_ij = del_ij - lower_mask_start;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_lower_mask_pos)) & 0x1));
      }
    }  
    break;
  case lmask_mode_c:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // Forming the unsigned long long int on the r.h.s. and then converting it to a signed
      // short int will translate the bit string appropriately.
      const int alt_mask_shft = static_cast<short int>((prof & lmask_c_shft) >> lmask_c_shft_pos);

      // Run the shift in terms of the index atom
      const int rel_ij = del_ij - alt_mask_shft;
      return (rel_ij >= 0 && rel_ij < lmask_long_extra_span &&
              ((prof >> (rel_ij + lmask_c_alt_mask_pos)) & 0x1));
    }
    break;
  case lmask_mode_d:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // This is the best possible path.  Obtain the number of masks and loop over all of them.
      const size_t nmasks = ((prof & lmask_d_array_cnt) >> lmask_d_array_cnt_pos);
      const size_t start_idx = ((prof & lmask_d_array_idx) >> lmask_d_array_idx_pos);
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }  
      return false;
    }
    break;
  case lmask_mode_e:
    {
      // This is the best possible path and there is no local exclusion arrangement to test.  Loop
      // over all the masks.
      const size_t nmasks = ((prof & lmask_e_array_cnt) >> lmask_e_array_cnt_pos);
      const size_t start_idx = (prof & lmask_e_array_idx);
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }
      return false;
    }
    break;
  case lmask_mode_f:
    break;
  default:
    break;
  }
  __builtin_unreachable();
}

// Detect the chip cache size
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 800
#    define LARGE_CHIP_CACHE
#  endif
#endif

// Single-precision tile evaluation
#define TCALC float
#  define TCALC2 float2
#  define TCALC4 float4
#  define TCOORD float
#  define TCOORD4 float4
#  define TCOORD_IS_REAL
#  define TACC   int
#  define TCALC_IS_SINGLE
#  ifdef LARGE_CHIP_CACHE
#    define PMENB_BLOCK_MULTIPLICITY 4
#  else
#    define PMENB_BLOCK_MULTIPLICITY 3
#  endif

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define KERNEL_NAME kfPmePairsForceEnergyDualTinyNonClash
#              include "pme_potential.cui"
#            undef KERNEL_NAME
#          undef CLASH_FORGIVENESS
#          define KERNEL_NAME kfPmePairsForceEnergyDualTiny
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kfPmePairsForceEnergyDualNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceEnergyDual
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kfPmePairsForceEnergyTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceEnergyTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceEnergyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsForceEnergy
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kfPmePairsForceDualTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceDualTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceDualNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsForceDual
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsForceTinyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsForceTiny
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsForceNonClash
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef CLASH_FORGIVENESS
#    define KERNEL_NAME kfPmePairsForce
#      include "pme_potential.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kfPmePairsEnergyDualTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsEnergyDualTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsEnergyDualNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsEnergyDual
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kfPmePairsEnergyTinyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsEnergyTiny
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define KERNEL_NAME kfPmePairsEnergyNonClash
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef CLASH_FORGIVENESS
#    define KERNEL_NAME kfPmePairsEnergy
#      include "pme_potential.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY

#  undef TCALC2
#  undef TCALC4
#  undef TCOORD
#  undef TCOORD4
#  undef TCOORD_IS_REAL
#  undef TACC
#  undef TCALC_IS_SINGLE
#  undef PMENB_BLOCK_MULTIPLICITY
#undef TCALC
  
// Double-precision tile evaluation
#define TCALC double
#  define TCALC2 double2
#  define TCALC4 double4
#  define TCOORD double
#  define TCOORD4 double4
#  define TCOORD_IS_REAL
#  define TCOORD_IS_LONG
#  define TACC   llint
#  ifdef LARGE_CHIP_CACHE
#    define PMENB_BLOCK_MULTIPLICITY 3
#  else
#    define PMENB_BLOCK_MULTIPLICITY 2
#  endif

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define KERNEL_NAME kdPmePairsForceEnergyDualTinyNonClash
#              include "pme_potential.cui"
#            undef KERNEL_NAME
#          undef CLASH_FORGIVENESS
#          define KERNEL_NAME kdPmePairsForceEnergyDualTiny
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kdPmePairsForceEnergyDualNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceEnergyDual
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kdPmePairsForceEnergyTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceEnergyTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceEnergyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsForceEnergy
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kdPmePairsForceDualTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceDualTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceDualNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsForceDual
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsForceTinyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsForceTiny
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsForceNonClash
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef CLASH_FORGIVENESS
#    define KERNEL_NAME kdPmePairsForce
#      include "pme_potential.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define KERNEL_NAME kdPmePairsEnergyDualTinyNonClash
#            include "pme_potential.cui"
#          undef KERNEL_NAME
#        undef CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsEnergyDualTiny
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsEnergyDualNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsEnergyDual
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define KERNEL_NAME kdPmePairsEnergyTinyNonClash
#          include "pme_potential.cui"
#        undef KERNEL_NAME
#      undef CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsEnergyTiny
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define KERNEL_NAME kdPmePairsEnergyNonClash
#        include "pme_potential.cui"
#      undef KERNEL_NAME
#    undef CLASH_FORGIVENESS
#    define KERNEL_NAME kdPmePairsEnergy
#      include "pme_potential.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY

#  undef TCALC2
#  undef TCALC4
#  undef TCOORD
#  undef TCOORD4
#  undef TCOORD_IS_REAL
#  undef TCOORD_IS_LONG
#  undef TACC
#  undef PMENB_BLOCK_MULTIPLICITY
#undef TCALC

// Clear hardware-dependent definitions
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 800
#    undef LARGE_CHIP_CACHE
#  endif
#endif

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_CUDA
extern cudaFuncAttributes queryPmePairsKernelRequirements(const PrecisionModel prec,
                                                          const EvaluateForce eval_frc,
                                                          const EvaluateEnergy eval_nrg,
                                                          const NeighborListKind neighbor_list,
                                                          const TinyBoxPresence has_tiny_box,
                                                          const ClashResponse collision_handling) {

  // As with other kernel querying functions, the kernel manager calling this function will have
  // specifications of the GPU in use.  It is the overall thread occupancy and multiplicity of each
  // kernel that this function must return.
  cudaFuncAttributes result;
  cudaError_t cfa;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergy);
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForce);
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:

        // If the force is not being evaluated, the energy must be required
        switch (neighbor_list) {
        case NeighborListKind::DUAL:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyDualTiny);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyDual);
            break;
          }
          break;
        case NeighborListKind::MONO:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyTiny);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergy);
            break;
          }
          break;
        }
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergy);
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForce);
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:

        // If the force is not being evaluated, the energy must be required
        switch (neighbor_list) {
        case NeighborListKind::DUAL:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyDualTiny);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyDual);
            break;
          }
          break;
        case NeighborListKind::MONO:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyTiny);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergy);
            break;
          }
          break;
        }
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceEnergyNonClash);
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdPmePairsForceNonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:

        // If the force is not being evaluated, the energy must be required
        switch (neighbor_list) {
        case NeighborListKind::DUAL:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyDualTinyNonClash);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyDualNonClash);
            break;
          }
          break;
        case NeighborListKind::MONO:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyTinyNonClash);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kdPmePairsEnergyNonClash);
            break;
          }
          break;
        }
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceEnergyNonClash);
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfPmePairsForceNonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:

        // If the force is not being evaluated, the energy must be required
        switch (neighbor_list) {
        case NeighborListKind::DUAL:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyDualTinyNonClash);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyDualNonClash);
            break;
          }
          break;
        case NeighborListKind::MONO:
          switch (has_tiny_box) {
          case TinyBoxPresence::YES:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyTinyNonClash);
            break;
          case TinyBoxPresence::NO:
            cfa = cudaFuncGetAttributes(&result, kfPmePairsEnergyNonClash);
            break;
          }
          break;
        }
        break;
      }
      break;
    }
    break;
  }

  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel k");
    switch (prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "s";
      break;
    }
    error_message += "PmePairs";
    switch (eval_frc) {
    case EvaluateForce::YES:
      error_message += "Force";
      break;
    case EvaluateForce::NO:
      break;
    }
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      error_message += "Energy";
      break;
    case EvaluateEnergy::NO:
      break;
    }

    // Report the error
    rtErr(error_message, "queryPmePairsKernelRequirements");
  }
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPmePairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                           const PPIKit<double, double4> &nrg_tab, ScoreCardWriter *scw,
                           CellGridWriter<double, llint, double, double4> *cgw,
                           MMControlKit<double> *ctrl, const EvaluateForce eval_frc,
                           const EvaluateEnergy eval_nrg, const TinyBoxPresence has_tiny_box,
                           const int2 bt, const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPmePairs() will invoke DOUBLE
  // precision calculations, coordinates, and parameter sets.  The single neighbor list grid also
  // implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself will
  // be queried for a condition to see whether there is a "tiny" simulation box with 4 cells along
  // any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                     clash_ratio, *scw, *cgw, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                  clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                 clash_ratio, *scw, *cgw, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergy<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForce<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergy<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPmePairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                           const PPIKit<double, double4> &nrg_tab, ScoreCardWriter *scw,
                           CellGridWriter<double, llint, double, double4> *cgw_qq,
                           CellGridWriter<double, llint, double, double4> *cgw_lj,
                           MMControlKit<double> *ctrl, const EvaluateForce eval_frc,
                           const EvaluateEnergy eval_nrg, const TinyBoxPresence has_tiny_box,
                           const int2 bt, const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPmePairs() will invoke DOUBLE
  // precision calculations, coordinates, and parameter sets.  The single neighbor list grid also
  // implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself will
  // be queried for a condition to see whether there is a "tiny" simulation box with 4 cells along
  // any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *cgw_qq,
                                                          *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                         clash_distance, clash_ratio, *scw,
                                                         *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw_qq,
                                                      *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                     clash_ratio, *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw_qq, *cgw_lj,
                                                  *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq,
                                                 *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdPmePairsForceEnergyDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdPmePairsForceDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw_qq, *cgw_lj,
                                              *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdPmePairsEnergyDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq, *cgw_lj,
                                             *ctrl);
        break;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPmePairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                           const PPIKit<float, float4> &nrg_tab, ScoreCardWriter *scw,
                           CellGridWriter<float, int, float, float4> *cgw,
                           MMControlKit<float> *ctrl, const EvaluateForce eval_frc,
                           const EvaluateEnergy eval_nrg, const TinyBoxPresence has_tiny_box,
                           const int2 bt, const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPmePairs() will invoke SINGLE
  // precision calculations, coordinates, and parameter sets.  The single neighbor list grid also
  // implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself will
  // be queried for a condition to see whether there is a "tiny" simulation box with 4 cells along
  // any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                     clash_ratio, *scw, *cgw, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                  clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                 clash_ratio, *scw, *cgw, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergy<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForce<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergy<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPmePairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                           const PPIKit<float, float4> &nrg_tab, ScoreCardWriter *scw,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj,
                           MMControlKit<float> *ctrl, const EvaluateForce eval_frc,
                           const EvaluateEnergy eval_nrg, const TinyBoxPresence has_tiny_box,
                           const int2 bt, const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPmePairs() will invoke SINGLE
  // precision calculations, coordinates, and parameter sets.  The single neighbor list grid also
  // implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself will
  // be queried for a condition to see whether there is a "tiny" simulation box with 4 cells along
  // any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                          clash_distance, clash_ratio,
                                                          *cgw_qq, *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyDualTinyNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                         clash_distance, clash_ratio, *scw,
                                                         *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw_qq,
                                                      *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyDualNonClash<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, clash_distance,
                                                     clash_ratio, *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw_qq,
                                                  *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyDualTiny<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq,
                                                 *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfPmePairsForceEnergyDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfPmePairsForceDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *cgw_qq, *cgw_lj,
                                              *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfPmePairsEnergyDual<<<bt.x, bt.y>>>(poly_nbk, lemr, tlpn, nrg_tab, *scw, *cgw_qq, *cgw_lj,
                                             *ctrl);
        break;
      }
      break;
    }
  }
}

#if 0
//-------------------------------------------------------------------------------------------------
extern void launchPmePairs(const PrecisionModel prec, CellGrid<double, llint, double, double4> *cg,
                           const LocalExclusionMask &lem, const TileManager &tlmn,
                           const PPITable<
                           const LocalExclusionMaskReader &lemr, const TilePlan &tlpn,
                           const PPIKit<float4> &nrg_tab, ScoreCardWriter *scw,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj,
                           MMControlKit<float> *ctrl, const EvaluateForce eval_frc,
                           const EvaluateEnergy eval_nrg, const TinyBoxPresence has_tiny_box,
                           const int2 bt, const double clash_distance, const double clash_ratio) {
}
#endif
  
} // namespace energy
} // namespace stormm
