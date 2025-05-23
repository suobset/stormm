#include "coordinate_graft.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"

namespace stormm {
namespace trajectory {

using constants::PrecisionModel;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::force_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
void graftCoordinateXYZ(const TrajectoryKind kind, llint* xdest, llint* ydest, llint* zdest,
                        const double dest_scale, const llint* xorig, const llint* yorig,
                        const llint* zorig, const double orig_scale, const int2* atom_list,
                        const int nxfer, int* xdest_ovrf, int* ydest_ovrf, int* zdest_ovrf,
                        const int* xorig_ovrf, const int* yorig_ovrf, const int* zorig_ovrf) {

  // Both the origin and destination store their coordinates in a fixed-precision format.
  const int dest_bits = round(log2(dest_scale));
  const int orig_bits = round(log2(orig_scale));
  const int scl_bit_delta = dest_bits - orig_bits;
  const bool dest_ovrf = overflowRequired(kind, dest_bits, xdest_ovrf, ydest_ovrf, zdest_ovrf);
  const bool orig_ovrf = overflowRequired(kind, orig_bits, xdest_ovrf, ydest_ovrf, zdest_ovrf);
  if (dest_ovrf && orig_ovrf) {
    for (int i = 0; i < nxfer; i++) {
      const size_t dest_idx = atom_list[i].x;
      const size_t orig_idx = atom_list[i].y;
      const int95_t ox_tmp = { xorig[orig_idx], xorig_ovrf[orig_idx] };
      const int95_t oy_tmp = { yorig[orig_idx], yorig_ovrf[orig_idx] };
      const int95_t oz_tmp = { zorig[orig_idx], zorig_ovrf[orig_idx] };
      const int95_t dx_tmp = hostChangeFPBits(ox_tmp, orig_bits, dest_bits);
      const int95_t dy_tmp = hostChangeFPBits(oy_tmp, orig_bits, dest_bits);
      const int95_t dz_tmp = hostChangeFPBits(oz_tmp, orig_bits, dest_bits);
      xdest[dest_idx] = dx_tmp.x;
      ydest[dest_idx] = dy_tmp.x;
      zdest[dest_idx] = dz_tmp.x;
      xdest_ovrf[dest_idx] = dx_tmp.y;
      ydest_ovrf[dest_idx] = dy_tmp.y;
      zdest_ovrf[dest_idx] = dz_tmp.y;
    }
  }
  else if (dest_ovrf) {
    for (int i = 0; i < nxfer; i++) {
      const size_t dest_idx = atom_list[i].x;
      const size_t orig_idx = atom_list[i].y;
      const int95_t ox_tmp = { xorig[orig_idx], 0 };
      const int95_t oy_tmp = { yorig[orig_idx], 0 };
      const int95_t oz_tmp = { zorig[orig_idx], 0 };
      const int95_t dx_tmp = hostChangeFPBits(ox_tmp, orig_bits, dest_bits);
      const int95_t dy_tmp = hostChangeFPBits(oy_tmp, orig_bits, dest_bits);
      const int95_t dz_tmp = hostChangeFPBits(oz_tmp, orig_bits, dest_bits);
      xdest[dest_idx] = dx_tmp.x;
      ydest[dest_idx] = dy_tmp.x;
      zdest[dest_idx] = dz_tmp.x;
      xdest_ovrf[dest_idx] = dx_tmp.y;
      ydest_ovrf[dest_idx] = dy_tmp.y;
      zdest_ovrf[dest_idx] = dz_tmp.y;
    }
  }
  else if (orig_ovrf) {
    for (int i = 0; i < nxfer; i++) {
      const size_t dest_idx = atom_list[i].x;
      const size_t orig_idx = atom_list[i].y;
      const int95_t ox_tmp = { xorig[orig_idx], xorig_ovrf[orig_idx] };
      const int95_t oy_tmp = { yorig[orig_idx], yorig_ovrf[orig_idx] };
      const int95_t oz_tmp = { zorig[orig_idx], zorig_ovrf[orig_idx] };
      const int95_t dx_tmp = hostChangeFPBits(ox_tmp, orig_bits, dest_bits);
      const int95_t dy_tmp = hostChangeFPBits(oy_tmp, orig_bits, dest_bits);
      const int95_t dz_tmp = hostChangeFPBits(oz_tmp, orig_bits, dest_bits);
      xdest[dest_idx] = dx_tmp.x;
      ydest[dest_idx] = dy_tmp.x;
      zdest[dest_idx] = dz_tmp.x;
    }
  }
  else {

    // Redirect to the templated function using pointers to long long integer data, to avoid
    // repeating code.
    graftCoordinateXYZ(kind, xdest, ydest, zdest, dest_scale, xorig, yorig, zorig, orig_scale,
                       atom_list, nxfer);
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const CoordinateFrameReader &orig,
                const int2* atom_list, const int nxfer, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:
      graftCoordinateXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                         orig.xcrd, orig.ycrd, orig.zcrd, 1.0, atom_list, nxfer);
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // If the function passes through the above switch, the coordinate transfer involves some
  // traffic to and from the GPU.  Extract the appropriate pointers from the abstracts, which are
  // trusted to point to the correct tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, dest_tier, orig_tier);
  launchGraftCoordXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                      orig.xcrd, orig.ycrd, orig.zcrd, 1.0, atom_list, nxfer, double_type_index,
                      double_type_index);
  launchResolution(sync, dest_tier, orig_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const CoordinateFrame &orig, const Hybrid<int2> &atom_list,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateFrameWriter destw = dest->data();
        coordGraft(&destw, orig.data(), atom_list.data(), atom_list.size(), dest_tier,
                   orig_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateFrameWriter destw = dest->deviceViewToHostData();
        coordGraft(&destw, orig.data(orig_tier), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), dest_tier, orig_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), dest_tier, orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orig_tier), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), dest_tier, orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const PhaseSpaceReader &orig, const int2* atom_list,
                const int nxfer, const TrajectoryKind kind, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xcrd, orig.ycrd,
                           orig.zcrd, 1.0, atom_list, nxfer);
        break;
      case TrajectoryKind::VELOCITIES:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xvel, orig.yvel,
                           orig.zvel, 1.0, atom_list, nxfer);
        break;
      case TrajectoryKind::FORCES:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xfrc, orig.yfrc,
                           orig.zfrc, 1.0, atom_list, nxfer);
        break;
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // If the function passes through the above switch, the coordinate transfer involves some
  // traffic to and from the GPU.  Extract the appropriate pointers from the abstracts, which are
  // trusted to point to the correct tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, dest_tier, orig_tier);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xcrd, orig.ycrd,
                        orig.zcrd, 1.0, atom_list, nxfer, double_type_index, double_type_index);
    break;
  case TrajectoryKind::VELOCITIES:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xvel, orig.yvel,
                        orig.zvel, 1.0, atom_list, nxfer, double_type_index, double_type_index);
    break;
  case TrajectoryKind::FORCES:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xfrc, orig.yfrc,
                        orig.zfrc, 1.0, atom_list, nxfer, double_type_index, double_type_index);
    break;
  }
  launchResolution(sync, dest_tier, orig_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                const TrajectoryKind kind, const CoordinateCycle orientation,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    {
      const PhaseSpaceReader origr = orig.data(orientation, orig_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        {
          CoordinateFrameWriter destw = dest->data(dest_tier);
          coordGraft(&destw, origr, atom_list.data(), atom_list.size(), kind, dest_tier,
                     orig_tier, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateFrameWriter destw = dest->deviceViewToHostData();
          coordGraft(&destw, origr, atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(),
                     kind, dest_tier, orig_tier, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(orientation),
                   atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(), kind, dest_tier,
                   orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orientation, orig_tier),
                   atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(), kind, dest_tier,
                   orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                const TrajectoryKind kind, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  coordGraft(dest, orig, atom_list, kind, orig.getCyclePosition(), dest_tier, orig_tier, gpu,
             sync);
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const PhaseSpaceReader &orig, const int2* atom_list,
                const int nxfer, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  coordGraft(dest, orig, atom_list, nxfer, TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu,
             sync);
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpace &orig, const Hybrid<int2> &atom_list,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    {
      const PhaseSpaceReader origr = orig.data(orig_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        {
          CoordinateFrameWriter destw = dest->data(dest_tier);
          coordGraft(&destw, origr, atom_list.data(), atom_list.size(), TrajectoryKind::POSITIONS,
                     dest_tier, orig_tier, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateFrameWriter destw = dest->deviceViewToHostData();
          coordGraft(&destw, origr, atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(),
                     TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orig_tier), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const PsSynthesisReader &orig,
                const int2* atom_list, const int nxfer, TrajectoryKind kind,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xcrd, orig.ycrd,
                           orig.zcrd, orig.gpos_scale, atom_list, nxfer);
        break;
      case TrajectoryKind::VELOCITIES:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xvel, orig.yvel,
                           orig.zvel, orig.vel_scale, atom_list, nxfer);
        break;
      case TrajectoryKind::FORCES:
        graftCoordinateXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xfrc, orig.yfrc,
                           orig.zfrc, orig.frc_scale, atom_list, nxfer);
        break;
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // If the function passes through the above switch, the coordinate transfer involves some
  // traffic to and from the GPU.  Extract the appropriate pointers from the abstracts, which are
  // trusted to point to the correct tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, dest_tier, orig_tier);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xcrd, orig.ycrd,
                        orig.zcrd, orig.gpos_scale, atom_list, nxfer, double_type_index,
                        double_type_index);
    break;
  case TrajectoryKind::VELOCITIES:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xvel, orig.yvel,
                        orig.zvel, orig.vel_scale, atom_list, nxfer, double_type_index,
                        double_type_index);
    break;
  case TrajectoryKind::FORCES:
    launchGraftCoordXYZ(kind, dest->xcrd, dest->ycrd, dest->zcrd, 1.0, orig.xfrc, orig.yfrc,
                        orig.zfrc, orig.frc_scale, atom_list, nxfer, double_type_index,
                        double_type_index);
    break;
  }
  launchResolution(sync, dest_tier, orig_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list, const TrajectoryKind kind,
                const CoordinateCycle orientation, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    {
      const PsSynthesisReader origr = orig.data(orientation, orig_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        {
          CoordinateFrameWriter destw = dest->data(dest_tier);
          coordGraft(&destw, origr, atom_list.data(), atom_list.size(), kind, dest_tier,
                     orig_tier, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateFrameWriter destw = dest->deviceViewToHostData();
          coordGraft(&destw, origr, atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(),
                     kind, dest_tier, orig_tier, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(orientation),
                   atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(), kind, dest_tier,
                   orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orientation, orig_tier),
                   atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(), kind, dest_tier,
                   orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const PsSynthesisReader &orig, const int2* atom_list,
                const int nxfer, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  coordGraft(dest, orig, atom_list, nxfer, TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu,
             sync);
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list, const TrajectoryKind kind,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, const HpcKernelSync sync) {
  coordGraft(dest, orig, atom_list, kind, orig.getCyclePosition(), dest_tier, orig_tier, gpu,
             sync);
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const PhaseSpaceSynthesis &orig,
                const Hybrid<int2> &atom_list, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    {
      const PsSynthesisReader origr = orig.data(orig_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        {
          CoordinateFrameWriter destw = dest->data(dest_tier);
          coordGraft(&destw, origr, atom_list.data(), atom_list.size(), TrajectoryKind::POSITIONS,
                     dest_tier, orig_tier, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateFrameWriter destw = dest->deviceViewToHostData();
          coordGraft(&destw, origr, atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(),
                     TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orig_tier), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), TrajectoryKind::POSITIONS, dest_tier, orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrameWriter *dest, const CondensateReader &orig, const int2* atom_list,
                const int nxfer, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:
      switch (orig.mode) {
      case PrecisionModel::DOUBLE:
        graftCoordinateXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                           orig.xcrd, orig.ycrd, orig.zcrd, 1.0, atom_list, nxfer);
        break;
      case PrecisionModel::SINGLE:
        graftCoordinateXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                           orig.xcrd_sp, orig.ycrd_sp, orig.zcrd_sp, 1.0, atom_list, nxfer);
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // If the function passes through the above switch, the coordinate transfer involves some
  // traffic to and from the GPU.  Extract the appropriate pointers from the abstracts, which are
  // trusted to point to the correct tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, dest_tier, orig_tier);
  switch (orig.mode) {
  case PrecisionModel::DOUBLE:
    launchGraftCoordXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                        orig.xcrd, orig.ycrd, orig.zcrd, 1.0, atom_list, nxfer, double_type_index,
                        double_type_index);
    break;
  case PrecisionModel::SINGLE:
    launchGraftCoordXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                        orig.xcrd_sp, orig.ycrd_sp, orig.zcrd_sp, 1.0, atom_list, nxfer,
                        float_type_index, float_type_index);
    break;
  }
  launchResolution(sync, dest_tier, orig_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordGraft(CoordinateFrame *dest, const Condensate &orig, const Hybrid<int2> &atom_list,
                const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
                const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    {
      const CondensateReader origr = orig.data(orig_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        {
          CoordinateFrameWriter destw = dest->data(dest_tier);
          coordGraft(&destw, origr, atom_list.data(), atom_list.size(), dest_tier, orig_tier, gpu,
                     sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateFrameWriter destw = dest->deviceViewToHostData();
          coordGraft(&destw, origr, atom_list.data(HybridTargetLevel::DEVICE), atom_list.size(),
                     dest_tier, orig_tier, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter destw = dest->data(dest_tier);
      switch (orig_tier) {
      case HybridTargetLevel::HOST:
        coordGraft(&destw, orig.deviceViewToHostData(), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), dest_tier, orig_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordGraft(&destw, orig.data(orig_tier), atom_list.data(HybridTargetLevel::DEVICE),
                   atom_list.size(), dest_tier, orig_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }

}

} // namespace trajectory
} // namespace stormm
