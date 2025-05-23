// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {
  
//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void graftCoordinateXYZ(const TrajectoryKind kind, Tdest* xdest, Tdest* ydest, Tdest* zdest,
                        const double dest_scale, const Torig* xorig, const Torig* yorig,
                        const Torig* zorig, const double orig_scale, const int2* atom_list,
                        const int nxfer) {
  if (dest_scale < 1.01) {
    if (orig_scale < 1.01) {

      // Neither coordinate system uses a fixed-precision representation.
      for (int i = 0; i < nxfer; i++) {
        const size_t dest_idx = atom_list[i].x;
        const size_t orig_idx = atom_list[i].y;
        xdest[dest_idx] = xorig[orig_idx];
        ydest[dest_idx] = yorig[orig_idx];
        zdest[dest_idx] = zorig[orig_idx];
      }
    }
    else {

      // The origin's coordinates use a fixed-precision representation while the destination's
      // coordinates do not.
      const double inv_orig_scale = 1.0 / orig_scale;
      for (int i = 0; i < nxfer; i++) {
        const size_t dest_idx = atom_list[i].x;
        const size_t orig_idx = atom_list[i].y;
        xdest[dest_idx] = static_cast<Tdest>(xorig[orig_idx]) * inv_orig_scale;
        ydest[dest_idx] = static_cast<Tdest>(yorig[orig_idx]) * inv_orig_scale;
        zdest[dest_idx] = static_cast<Tdest>(zorig[orig_idx]) * inv_orig_scale;
      }
    }
  }
  else {
    const int dest_bits = round(log2(dest_scale));
    if (orig_scale < 1.01) {

      // The origin's coordinates are represented as real numbers while the destination's
      // coordinates have a fixed-precision representation.
      for (int i = 0; i < nxfer; i++) {
        const size_t dest_idx = atom_list[i].x;
        const size_t orig_idx = atom_list[i].y;
        xdest[dest_idx] = llround(static_cast<double>(xorig[orig_idx]) * dest_scale);
        ydest[dest_idx] = llround(static_cast<double>(yorig[orig_idx]) * dest_scale);
        zdest[dest_idx] = llround(static_cast<double>(zorig[orig_idx]) * dest_scale);
      }
    }
    else {

      // Both the origin and destination store their coordinates in a fixed-precision format.
      const int orig_bits = round(log2(orig_scale));
      const int scl_bit_delta = dest_bits - orig_bits;
      if (scl_bit_delta < 0) {
        const llint divisor = round(pow(2.0, -scl_bit_delta));
        for (int i = 0; i < nxfer; i++) {
          const size_t dest_idx = atom_list[i].x;
          const size_t orig_idx = atom_list[i].y;
          xdest[dest_idx] = xorig[orig_idx] / divisor;
          ydest[dest_idx] = yorig[orig_idx] / divisor;
          zdest[dest_idx] = zorig[orig_idx] / divisor;
        }
      }
      else if (scl_bit_delta == 0) {
        for (int i = 0; i < nxfer; i++) {
          const size_t dest_idx = atom_list[i].x;
          const size_t orig_idx = atom_list[i].y;
          xdest[dest_idx] = xorig[orig_idx];
          ydest[dest_idx] = yorig[orig_idx];
          zdest[dest_idx] = zorig[orig_idx];
        }
      }
      else {
        const llint multiplier = round(pow(2.0, scl_bit_delta));
        for (int i = 0; i < nxfer; i++) {
          const size_t dest_idx = atom_list[i].x;
          const size_t orig_idx = atom_list[i].y;
          xdest[dest_idx] = xorig[orig_idx] * multiplier;
          ydest[dest_idx] = yorig[orig_idx] * multiplier;
          zdest[dest_idx] = zorig[orig_idx] * multiplier;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordGraft(CoordinateFrameWriter *dest, const CoordinateSeriesReader<Torig> &orig,
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
  // Traffic to and from the GPU.  Extract the appropriate pointers from the abstracts, which are
  // trusted to point to the correct tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, dest_tier, orig_tier);
  launchGraftCoordXYZ(TrajectoryKind::POSITIONS, dest->xcrd, dest->ycrd, dest->zcrd, 1.0,
                      orig.xcrd, orig.ycrd, orig.zcrd, 1.0, atom_list, nxfer, double_type_index,
                      std::type_index(typeid(Torig)).hash_code());
  launchResolution(sync, dest_tier, orig_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordGraft(CoordinateFrame *dest, const CoordinateSeries<Torig> &orig,
                const Hybrid<int2> &atom_list, const HybridTargetLevel dest_tier,
                const HybridTargetLevel orig_tier, const GpuDetails &gpu,
                const HpcKernelSync sync) {
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

} // namespace trajectory
} // namespace stormm
