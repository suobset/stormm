#include "copyright.h"
#include "Reporting/error_format.h"
#include "coordinate_copy.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void coordCopyValidateAtomCounts(const int destination_atoms, const int origin_atoms) {
  if (origin_atoms != destination_atoms) {
    rtErr("The number of atoms in the destination frame (" + std::to_string(destination_atoms) +
          ") must equal the number of atoms in the origin frame (" + std::to_string(origin_atoms) +
          ").", "coordCopy");
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopyValidateFrameIndex(const int frame_index, const int frame_count) {
  if (frame_index > frame_count) {
    rtErr("Frame index " + std::to_string(frame_index) + " is invalid in a coordinate series of " +
          std::to_string(frame_count) + " frames.", "coordCopy");
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopyValidateSystemIndex(const int system_index, const int system_count) {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a synthesis of " +
          std::to_string(system_count) + " systems.", "coordCopy");
  }
}

//-------------------------------------------------------------------------------------------------
void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        const int index_dest, const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, const int index_orig, llint* boxvecs,
                        int* boxvec_ovrf) {
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_offset = xfrm_w * index_dest;
  const int dest_bdim_offset = bdim_w * index_dest;
  const int orig_xfrm_offset = xfrm_w * index_orig;
  const int orig_bdim_offset = bdim_w * index_orig;  
  for (int i = 0; i < 6; i++) {
    dest_boxdim[dest_bdim_offset + i] = orig_boxdim[orig_bdim_offset + i];
  }
  for (int i = 0; i < 9; i++) {
    dest_umat[dest_xfrm_offset + i] = orig_umat[orig_xfrm_offset + i];
    dest_invu[dest_xfrm_offset + i] = orig_invu[orig_xfrm_offset + i];
  }

  // Fill out the fixed-precision box vectors for PhaseSpaceSynthesis objects
  if (boxvecs != nullptr) {
    for (int i = 0; i < 9; i++) {
      const int95_t ival = hostDoubleToInt95(orig_invu[i]);
      boxvecs[dest_xfrm_offset + i]     = ival.x;
      boxvec_ovrf[dest_xfrm_offset + i] = ival.y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        int index_dest, const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, llint* boxvecs, int* boxvec_ovrf) {
  copyBoxInformation(dest_boxdim, dest_umat, dest_invu, index_dest, orig_boxdim, orig_umat,
                     orig_invu, 0, boxvecs, boxvec_ovrf);
}

//-------------------------------------------------------------------------------------------------
void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu, const int index_orig) {
  copyBoxInformation(dest_boxdim, dest_umat, dest_invu, 0, orig_boxdim, orig_umat, orig_invu,
                     index_orig);
}

//-------------------------------------------------------------------------------------------------
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, const double dest_scale, const llint* xorig,
                       const int* xorig_ovrf, const llint* yorig, const int* yorig_ovrf,
                       const llint* zorig, const int* zorig_ovrf, const double orig_scale,
                       const int natom) {
  const int dest_bits = round(log2(dest_scale));
  const int orig_bits = round(log2(orig_scale));
  for (int i = 0; i < natom; i++) {
    const int95_t xo = { xorig[i], xorig_ovrf[i] };
    const int95_t yo = { yorig[i], yorig_ovrf[i] };
    const int95_t zo = { zorig[i], zorig_ovrf[i] };
    const int95_t xn = hostChangeFPBits(xo, orig_bits, dest_bits);
    const int95_t yn = hostChangeFPBits(yo, orig_bits, dest_bits);
    const int95_t zn = hostChangeFPBits(zo, orig_bits, dest_bits);
    xdest[i] = xn.x;
    ydest[i] = yn.x;
    zdest[i] = zn.x;
    xdest_ovrf[i] = xn.y;
    ydest_ovrf[i] = yn.y;
    zdest_ovrf[i] = zn.y;
  }
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const CoordinateFrameReader &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (origin_tier) {
  case HybridTargetLevel::HOST:
    switch (destination_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                         origin.umat, origin.invu);
      copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xcrd,
                        origin.ycrd, origin.zcrd, origin.natom);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                                 double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                                 double_type_index, origin.natom, destination->umat,
                                 destination->invu, destination->boxdim, origin.umat, origin.invu,
                                 origin.boxdim, 0, 0, 0, 0, 0, 0, gpu);
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, origin.natom, destination->umat,
                               destination->invu, destination->boxdim, origin.umat, origin.invu,
                               origin.boxdim, 0, 0, 0, 0, 0, 0, gpu);
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateFrameWriter cfw = destination->data(destination_tier);
        coordCopy(&cfw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateFrameWriter cfw = destination->deviceViewToHostData();
        coordCopy(&cfw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter cfw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cfw, origin.deviceViewToHostData(), destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cfw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const PhaseSpaceReader &origin,
               const TrajectoryKind kind, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                           origin.boxdim, origin.umat, origin.invu);
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xcrd,
                          origin.ycrd, origin.zcrd, origin.natom);
        break;
      case TrajectoryKind::VELOCITIES:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xvel,
                          origin.yvel, origin.zvel, origin.natom);
        break;
      case TrajectoryKind::FORCES:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xfrc,
                          origin.yfrc, origin.zfrc, origin.natom);
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

  // If the function passed through the above switch, the copy involves some copying to and/or
  // from the device, any variant of which is accomplished by the code below.  Extract the
  // appropriate pointers from the given abstracts, which are trusted to point to the correct
  // tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, origin.natom, destination->umat,
                               destination->invu, destination->boxdim, origin.umat,
                               origin.invu, origin.boxdim, 0, 0, 0, 0, 0, 0, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                            double_type_index, origin.xvel, origin.yvel, origin.zvel, 1.0,
                            double_type_index, origin.natom, 0, 0, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                            double_type_index, origin.xfrc, origin.yfrc, origin.zfrc, 1.0,
                            double_type_index, origin.natom, 0, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const PhaseSpaceReader &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, TrajectoryKind::POSITIONS, destination_tier, origin_tier, gpu,
            sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin, const TrajectoryKind kind,
               const CoordinateCycle orientation, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateFrameWriter cfw = destination->data(destination_tier);
        coordCopy(&cfw, origin.data(orientation, origin_tier), kind, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateFrameWriter cfw = destination->deviceViewToHostData();
        coordCopy(&cfw, origin.data(orientation, origin_tier), kind, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter cfw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cfw, origin.deviceViewToHostData(orientation), kind, destination_tier,
                  origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cfw, origin.data(orientation, origin_tier), kind, destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, kind, origin.getCyclePosition(), destination_tier, origin_tier,
            gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, TrajectoryKind::POSITIONS, origin.getCyclePosition(),
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const PsSynthesisReader &origin,
               const int orig_atom_start, const int index_orig, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int atom_count = destination->natom;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                           origin.boxdims, origin.umat, origin.invu, index_orig);
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                          &origin.xcrd[orig_atom_start], &origin.xcrd_ovrf[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.ycrd_ovrf[orig_atom_start],
                          &origin.zcrd[orig_atom_start], &origin.zcrd_ovrf[orig_atom_start],
                          origin.gpos_scale, atom_count);
        break;
      case TrajectoryKind::VELOCITIES:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                          &origin.xvel[orig_atom_start], &origin.xvel_ovrf[orig_atom_start],
                          &origin.yvel[orig_atom_start], &origin.yvel_ovrf[orig_atom_start],
                          &origin.zvel[orig_atom_start], &origin.zvel_ovrf[orig_atom_start],
                          origin.vel_scale, atom_count);
        break;
      case TrajectoryKind::FORCES:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                          &origin.xfrc[orig_atom_start], &origin.xfrc_ovrf[orig_atom_start],
                          &origin.yfrc[orig_atom_start], &origin.yfrc_ovrf[orig_atom_start],
                          &origin.zfrc[orig_atom_start], &origin.zfrc_ovrf[orig_atom_start],
                          origin.frc_scale, atom_count);
        break;
      }

      // Return after completing the host-to-host copy
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

  // If the function passed through the above switch, the copy involves some copying to and/or
  // from the device, any variant of which is accomplished by the code below.  Extract the
  // appropriate pointers from the given abstracts, which are trusted to point to the correct
  // tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const size_t bdim_start = static_cast<size_t>(index_orig) * roundUp<size_t>(6, warp_size_zu);
  const size_t xfrm_start = static_cast<size_t>(index_orig) * roundUp<size_t>(9, warp_size_zu);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0, 0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd,
                               origin.xcrd_ovrf, origin.ycrd_ovrf, origin.zcrd_ovrf,
                               origin.gpos_scale, origin.gpos_bits, atom_count, destination->umat,
                               destination->invu, destination->boxdim, origin.umat, origin.invu,
                               origin.boxdims, 0, orig_atom_start, 0, xfrm_start, 0, bdim_start,
                               gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, 1.0, 0,
                            double_type_index, origin.xvel, origin.yvel, origin.zvel,
                            origin.xvel_ovrf, origin.yvel_ovrf, origin.zvel_ovrf, origin.vel_scale,
                            origin.vel_bits, atom_count, 0, orig_atom_start, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, 1.0, 0,
                            double_type_index, origin.xfrc, origin.yfrc, origin.zfrc,
                            origin.xfrc_ovrf, origin.yfrc_ovrf, origin.zfrc_ovrf, origin.frc_scale,
                            origin.frc_bits, atom_count, 0, orig_atom_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin,
               const int index_orig, const TrajectoryKind kind,
               const CoordinateCycle orientation, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount(index_orig));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateFrameWriter cfw = destination->data();
        coordCopy(&cfw, origin.data(orientation, origin_tier), orig_atom_start, index_orig, kind,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateFrameWriter cfw = destination->deviceViewToHostData();
        coordCopy(&cfw, origin.data(orientation, origin_tier), orig_atom_start, index_orig, kind,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter cfw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cfw, origin.deviceViewToHostData(orientation), orig_atom_start, index_orig,
                  kind, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cfw, origin.data(orientation, origin_tier), orig_atom_start, index_orig, kind,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin,
               const int index_orig, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, index_orig, kind, origin.getCyclePosition(), destination_tier,
            origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin,
               const int index_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, origin, index_orig, TrajectoryKind::POSITIONS, origin.getCyclePosition(),
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const CondensateReader &origin,
               const int orig_atom_start, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->natom;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                         origin.boxdims, origin.umat, origin.invu, index_orig);
      switch (origin.mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                          &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                          &origin.zcrd[orig_atom_start], natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                          &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                          &origin.zcrd_sp[orig_atom_start], natom);
        break;
      }

      // Return after completing the host-to-host copy
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

  // If the function passed through the above switch, the copy involves some copying to and/or
  // from the device, any variant of which is accomplished by the code below.  Extract the
  // appropriate pointers from the given abstracts, which are trusted to point to the correct
  // tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int orig_xfrm_start = roundUp(9, warp_size_int) * index_orig;
  const int orig_bdim_start = roundUp(6, warp_size_int) * index_orig;
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, natom, destination->umat, destination->invu,
                               destination->boxdim, origin.umat, origin.invu, origin.boxdims, 0,
                               orig_atom_start, 0, orig_xfrm_start, 0, orig_bdim_start, gpu);
    break;
  case PrecisionModel::SINGLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd_sp, origin.ycrd_sp, origin.zcrd_sp,
                               1.0, float_type_index, natom, destination->umat, destination->invu,
                               destination->boxdim, origin.umat, origin.invu, origin.boxdims, 0,
                               orig_atom_start, 0, orig_xfrm_start, 0, orig_bdim_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount(index_orig));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateFrameWriter cfw = destination->data(destination_tier);
        coordCopy(&cfw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateFrameWriter cfw = destination->deviceViewToHostData();
        coordCopy(&cfw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateFrameWriter cfw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cfw, origin.deviceViewToHostData(), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cfw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const TrajectoryKind kind,
               const CoordinateFrameReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                           origin.boxdim, origin.umat, origin.invu);
        copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xcrd,
                          origin.ycrd, origin.zcrd, origin.natom);
        break;
      case TrajectoryKind::VELOCITIES:
        copyCoordinateXYZ(destination->xvel, destination->yvel, destination->zvel, origin.xcrd,
                          origin.ycrd, origin.zcrd, origin.natom);
        break;
      case TrajectoryKind::FORCES:
        copyCoordinateXYZ(destination->xfrc, destination->yfrc, destination->zfrc, origin.xcrd,
                          origin.ycrd, origin.zcrd, origin.natom);
        break;
      }

      // Return after completing the host-to-host copy
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

  // If the function passed through the above switch, the copy involves some copying to and/or
  // from the device, any variant of which is accomplished by the code below.  Extract the
  // appropriate pointers from the given abstracts, which are trusted to point to the correct
  // tier of memory.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, origin.natom, destination->umat,
                               destination->invu, destination->boxdim, origin.umat, origin.invu,
                               origin.boxdim, 0, 0, 0, 0, 0, 0, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destination->xvel, destination->yvel, destination->zvel, 1.0,
                            double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                            double_type_index, origin.natom, 0, 0, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destination->xfrc, destination->yfrc, destination->zfrc, 1.0,
                            double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                            double_type_index, origin.natom, 0, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const CoordinateFrameReader &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, origin, destination_tier, origin_tier, gpu,
            sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PhaseSpaceWriter psw = destination->data(orientation, destination_tier);
        coordCopy(&psw, kind, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PhaseSpaceWriter psw = destination->deviceViewToHostData(orientation);
        coordCopy(&psw, kind, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PhaseSpaceWriter psw = destination->data(orientation, destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&psw, kind, origin.deviceViewToHostData(), destination_tier, origin_tier, gpu,
                  sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&psw, kind, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, kind, destination->getCyclePosition(), origin, destination_tier,
            origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:

      // Copy the box information
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                         origin.umat, origin.invu);
      copyBoxInformation(destination->boxdim_alt, destination->umat_alt, destination->invu_alt,
                         origin.boxdim_alt, origin.umat_alt, origin.invu_alt);

      // Copy present positions, velocities, and forces
      copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xcrd,
                        origin.ycrd, origin.zcrd, origin.natom);
      copyCoordinateXYZ(destination->xvel, destination->yvel, destination->zvel, origin.xvel,
                        origin.yvel, origin.zvel, origin.natom);
      copyCoordinateXYZ(destination->xfrc, destination->yfrc, destination->zfrc, origin.xfrc,
                        origin.yfrc, origin.zfrc, origin.natom);

      // Copy alternate positions, velocities, and forces
      copyCoordinateXYZ(destination->xalt, destination->yalt, destination->zalt, origin.xalt,
                        origin.yalt, origin.zalt, origin.natom);
      copyCoordinateXYZ(destination->vxalt, destination->vyalt, destination->vzalt, origin.vxalt,
                        origin.vyalt, origin.vzalt, origin.natom);
      copyCoordinateXYZ(destination->fxalt, destination->fyalt, destination->fzalt, origin.fxalt,
                        origin.fyalt, origin.fzalt, origin.natom);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchCopyCoordinates(destination, origin, gpu);
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchCopyCoordinates(destination, origin, gpu);
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const PhaseSpace &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PhaseSpaceWriter psw = destination->data(destination_tier);
        coordCopy(&psw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PhaseSpaceWriter psw = destination->deviceViewToHostData();
        coordCopy(&psw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PhaseSpaceWriter psw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&psw, origin.deviceViewToHostData(), destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&psw, origin.data(origin_tier), destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const PsSynthesisReader &origin,
               const int orig_atom_start, const int index_orig,
               const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->natom;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:

      // Copy the box information
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                         origin.boxdims, origin.umat, origin.invu, index_orig);
      copyBoxInformation(destination->boxdim_alt, destination->umat_alt, destination->invu_alt,
                         origin.alt_boxdims, origin.umat_alt, origin.invu_alt, index_orig);

      // Copy present positions, velocities, and forces
      copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                        &origin.xcrd[orig_atom_start], &origin.xcrd_ovrf[orig_atom_start],
                        &origin.ycrd[orig_atom_start], &origin.ycrd_ovrf[orig_atom_start],
                        &origin.zcrd[orig_atom_start], &origin.zcrd_ovrf[orig_atom_start],
                        origin.gpos_scale, natom);
      copyCoordinateXYZ(destination->xvel, destination->yvel, destination->zvel,
                        &origin.xvel[orig_atom_start], &origin.xvel_ovrf[orig_atom_start],
                        &origin.yvel[orig_atom_start], &origin.yvel_ovrf[orig_atom_start],
                        &origin.zvel[orig_atom_start], &origin.zvel_ovrf[orig_atom_start],
                        origin.vel_scale, natom);
      copyCoordinateXYZ(destination->xfrc, destination->yfrc, destination->zfrc,
                        &origin.xfrc[orig_atom_start], &origin.xfrc_ovrf[orig_atom_start],
                        &origin.yfrc[orig_atom_start], &origin.yfrc_ovrf[orig_atom_start],
                        &origin.zfrc[orig_atom_start], &origin.zfrc_ovrf[orig_atom_start],
                        origin.frc_scale, natom);

      // Copy alternate positions, velocities, and forces
      copyCoordinateXYZ(destination->xalt, destination->yalt, destination->zalt,
                        &origin.xalt[orig_atom_start], &origin.xalt_ovrf[orig_atom_start],
                        &origin.yalt[orig_atom_start], &origin.yalt_ovrf[orig_atom_start],
                        &origin.zalt[orig_atom_start], &origin.zalt_ovrf[orig_atom_start],
                        origin.gpos_scale, natom);
      copyCoordinateXYZ(destination->vxalt, destination->vyalt, destination->vzalt,
                        &origin.vxalt[orig_atom_start], &origin.vxalt_ovrf[orig_atom_start],
                        &origin.vyalt[orig_atom_start], &origin.vyalt_ovrf[orig_atom_start],
                        &origin.vzalt[orig_atom_start], &origin.vzalt_ovrf[orig_atom_start],
                        origin.vel_scale, natom);
      copyCoordinateXYZ(destination->vxalt, destination->vyalt, destination->vzalt,
                        &origin.vxalt[orig_atom_start], &origin.vxalt_ovrf[orig_atom_start],
                        &origin.vyalt[orig_atom_start], &origin.vyalt_ovrf[orig_atom_start],
                        &origin.vzalt[orig_atom_start], &origin.vzalt_ovrf[orig_atom_start],
                        origin.frc_scale, natom);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchCopyCoordinates(destination, origin, orig_atom_start, index_orig, gpu);
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchCopyCoordinates(destination, origin, orig_atom_start, index_orig, gpu);
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const PhaseSpaceSynthesis &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount(index_orig));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (destination_tier) {
    case HybridTargetLevel::HOST:
      {
        PhaseSpaceWriter psw = destination->data(destination_tier);
        coordCopy(&psw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PhaseSpaceWriter psw = destination->deviceViewToHostData();
        coordCopy(&psw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PhaseSpaceWriter psw = destination->data(destination_tier);
      switch (destination_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&psw, origin.deviceViewToHostData(), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&psw, origin.data(origin_tier), orig_atom_start, index_orig, destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const TrajectoryKind kind,
               const CondensateReader &origin, const int orig_atom_start, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->natom;
  double *xdest, *ydest, *zdest;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    xdest = destination->xcrd;
    ydest = destination->ycrd;
    zdest = destination->zcrd;
    break;
  case TrajectoryKind::VELOCITIES:
    xdest = destination->xvel;
    ydest = destination->yvel;
    zdest = destination->zvel;
    break;
  case TrajectoryKind::FORCES:
    xdest = destination->xfrc;
    ydest = destination->yfrc;
    zdest = destination->zfrc;
    break;
  }
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                           origin.boxdims, origin.umat, origin.invu, index_orig);
        break;
      case TrajectoryKind::VELOCITIES:
      case TrajectoryKind::FORCES:
        break;
      }
      switch (origin.mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(xdest, ydest, zdest, &origin.xcrd[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                          destination->natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(xdest, ydest, zdest, &origin.xcrd_sp[orig_atom_start],
                          &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start],
                          destination->natom);
        break;
      }

      // Return after completing the host-to-host copy case
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

  // If the function passed through the above switch, the copy involves some copying to and/or
  // from the device, any variant of which is accomplished by the code below.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int orig_xfrm_start = roundUp(9, warp_size_int) * index_orig;
  const int orig_bdim_start = roundUp(6, warp_size_int) * index_orig;
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      launchCopyCoordinateXYZBox(xdest, ydest, zdest, 1.0, double_type_index, origin.xcrd,
                                 origin.ycrd, origin.zcrd, 1.0, double_type_index, natom,
                                 destination->umat, destination->invu, destination->boxdim,
                                 origin.umat, origin.invu, origin.boxdims, 0, orig_atom_start,
                                 0, orig_xfrm_start, 0, orig_bdim_start, gpu);
      break;
    case TrajectoryKind::VELOCITIES:
    case TrajectoryKind::FORCES:
      launchCopyCoordinateXYZ(xdest, ydest, zdest, 1.0, double_type_index, origin.xcrd,
                              origin.ycrd, origin.zcrd, 1.0, double_type_index, natom, 0,
                              orig_atom_start, gpu);
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      launchCopyCoordinateXYZBox(xdest, ydest, zdest, 1.0, double_type_index, origin.xcrd_sp,
                                 origin.ycrd_sp, origin.zcrd_sp, 1.0, float_type_index, natom,
                                 destination->umat, destination->invu, destination->boxdim,
                                 origin.umat, origin.invu, origin.boxdims, 0, orig_atom_start,
                                 0, orig_xfrm_start, 0, orig_bdim_start, gpu);
      break;
    case TrajectoryKind::VELOCITIES:
    case TrajectoryKind::FORCES:
      launchCopyCoordinateXYZ(xdest, ydest, zdest, 1.0, double_type_index, origin.xcrd_sp,
                              origin.ycrd_sp, origin.zcrd_sp, 1.0, float_type_index, natom, 0,
                              orig_atom_start, gpu);
      break;
    }
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const CondensateReader &origin,
               const int orig_atom_start, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, origin, orig_atom_start, index_orig,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount(index_orig));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PhaseSpaceWriter psw = destination->data(orientation, destination_tier);
        coordCopy(&psw, kind, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PhaseSpaceWriter psw = destination->deviceViewToHostData(orientation);
        coordCopy(&psw, kind, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PhaseSpaceWriter psw = destination->data(orientation, destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&psw, kind, origin.deviceViewToHostData(), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&psw, kind, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }    
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind, const Condensate &origin,
               const int index_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, kind, destination->getCyclePosition(), origin, index_orig,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
            index_orig, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrameReader &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      destination->import(origin, index_dest, orientation, kind);

      // Return after completing the host-to-host transfer
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

  // Funnel all device-dependent transfers past the above switch here.  A rare unrolling of an
  // enumeration is used to avoid excessive code replication and other messy effects.
#ifdef STORMM_USE_HPC
  PsSynthesisWriter destr = (destination_tier == HybridTargetLevel::HOST) ?
                            destination->deviceViewToHostData() :
                            destination->data(destination_tier);
  llint *xptr, *yptr, *zptr, *iboxv_ptr;
  int *xptr_ovrf, *yptr_ovrf, *zptr_ovrf, *iboxv_ovrf_ptr;
  double *umat_ptr, *invu_ptr, *bdim_ptr;
  switch (orientation) {
  case CoordinateCycle::ALTERNATE:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      xptr = destr.xalt;
      xptr_ovrf = destr.xalt_ovrf;
      yptr = destr.yalt;
      yptr_ovrf = destr.yalt_ovrf;
      zptr = destr.zalt;
      zptr_ovrf = destr.zalt_ovrf;
      umat_ptr = destr.umat_alt;
      invu_ptr = destr.invu_alt;
      bdim_ptr = destr.alt_boxdims;
      iboxv_ptr = destr.alt_boxvecs;
      iboxv_ovrf_ptr = destr.alt_boxvec_ovrf;
      break;
    case TrajectoryKind::VELOCITIES:
      xptr = destr.vxalt;
      xptr_ovrf = destr.vxalt_ovrf;
      yptr = destr.vyalt;
      yptr_ovrf = destr.vyalt_ovrf;
      zptr = destr.vzalt;
      zptr_ovrf = destr.vzalt_ovrf;
      break;
    case TrajectoryKind::FORCES:
      xptr = destr.fxalt;
      xptr_ovrf = destr.fxalt_ovrf;
      yptr = destr.fyalt;
      yptr_ovrf = destr.fyalt_ovrf;
      zptr = destr.fzalt;
      zptr_ovrf = destr.fzalt_ovrf;
      break;
    }
    break;
  case CoordinateCycle::PRIMARY:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      xptr = destr.xcrd;
      xptr_ovrf = destr.xcrd_ovrf;
      yptr = destr.ycrd;
      yptr_ovrf = destr.ycrd_ovrf;
      zptr = destr.zcrd;
      zptr_ovrf = destr.zcrd_ovrf;
      umat_ptr = destr.umat;
      invu_ptr = destr.invu;
      bdim_ptr = destr.boxdims;
      iboxv_ptr = destr.boxvecs;
      iboxv_ovrf_ptr = destr.boxvec_ovrf;
      break;
    case TrajectoryKind::VELOCITIES:
      xptr = destr.xvel;
      xptr_ovrf = destr.xvel_ovrf;
      yptr = destr.yvel;
      yptr_ovrf = destr.yvel_ovrf;
      zptr = destr.zvel;
      zptr_ovrf = destr.zvel_ovrf;
      break;
    case TrajectoryKind::FORCES:
      xptr = destr.xfrc;
      xptr_ovrf = destr.xfrc_ovrf;
      yptr = destr.yfrc;
      yptr_ovrf = destr.yfrc_ovrf;
      zptr = destr.zfrc;
      zptr_ovrf = destr.zfrc_ovrf;
      break;
    }
    break;
  }
  const int dest_atom_offset = destination->getAtomOffset(index_dest);
  const int dest_xfrm_offset = roundUp(9, warp_size_int) * index_dest;
  const int dest_bdim_offset = roundUp(6, warp_size_int) * index_dest;
  launchPreparation(sync, destination_tier, origin_tier);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(xptr, yptr, zptr, xptr_ovrf, yptr_ovrf, zptr_ovrf, destr.gpos_scale,
                               destr.gpos_bits, origin.xcrd, origin.ycrd, origin.zcrd, 1.0, 0,
                               double_type_index, origin.natom, umat_ptr, invu_ptr, bdim_ptr,
                               origin.umat, origin.invu, origin.boxdim, iboxv_ptr, iboxv_ovrf_ptr,
                               dest_atom_offset, 0, dest_xfrm_offset, 0, dest_bdim_offset, 0, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(xptr, yptr, zptr, xptr_ovrf, yptr_ovrf, zptr_ovrf, destr.vel_scale,
                            destr.vel_bits, origin.xcrd, origin.ycrd, origin.zcrd, 1.0, 0,
                            double_type_index, origin.natom, dest_atom_offset, 0, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(xptr, yptr, zptr, xptr_ovrf, yptr_ovrf, zptr_ovrf, destr.frc_scale,
                            destr.frc_bits, origin.xcrd, origin.ycrd, origin.zcrd, 1.0, 0,
                            double_type_index, origin.natom, dest_atom_offset, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateFrameReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, kind, destination->getCyclePosition(), origin,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateFrameReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            origin, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {

  // No checks are needed here as the function to which this delegates will perform checks on both
  // the system in dex in the destination and the atom counts in both objects.
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      destination->import(origin, index_dest, orientation, kind);

      // Return after completing the host-to-host transfer
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

  // As in other situations, the HPC code funnels pasts the switch if device data is involved.
  // This call does not need to be hidden behind pre-processor pragmas as it is delegating work
  // to a variant of coordcopy() involving similar objects.  However, the CPU code will never
  // reach this point.
#ifdef STORMM_USE_HPC
  switch (origin_tier) {
  case HybridTargetLevel::HOST:
    coordCopy(destination, index_dest, kind, orientation, origin.deviceViewToHostData(),
              destination_tier, origin_tier, gpu, sync);
    break;
  case HybridTargetLevel::DEVICE:
    coordCopy(destination, index_dest, kind, orientation, origin.data(origin_tier),
              destination_tier, origin_tier, gpu, sync);
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateFrame &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, kind, destination->getCyclePosition(), origin,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateFrame &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            origin, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpaceReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      destination->import(origin, index_dest);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter destr = destination->deviceViewToHostData();
        launchCopyCoordinates(&destr, destination->getAtomOffset(index_dest), index_dest, origin,
                              gpu);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter destr = destination->data(destination_tier);
      launchCopyCoordinates(&destr, destination->getAtomOffset(index_dest), index_dest, origin,
                            gpu);
    }
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpace &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.getAtomCount());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      destination->import(origin, index_dest);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter destr = destination->deviceViewToHostData();
        launchCopyCoordinates(&destr, destination->getAtomOffset(index_dest), index_dest,
                              origin.data(origin_tier), gpu);
      }
      break;
#endif
    }
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter destr = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        launchCopyCoordinates(&destr, destination->getAtomOffset(index_dest), index_dest,
                              origin.deviceViewToHostData(), gpu);
        break;
      case HybridTargetLevel::DEVICE:
        launchCopyCoordinates(&destr, destination->getAtomOffset(index_dest), index_dest,
                              origin.data(origin_tier), gpu);
        break;
      }
    }
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int dest_atom_start, const int index_dest,
               const TrajectoryKind kind, const CondensateReader &origin,
               const int orig_atom_start, const int index_orig, const int natom,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  llint *xdest, *ydest, *zdest;
  int *xdest_ovrf, *ydest_ovrf, *zdest_ovrf;
  double dest_scale;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    xdest = destination->xcrd;
    ydest = destination->ycrd;
    zdest = destination->zcrd;
    xdest_ovrf = destination->xcrd_ovrf;
    ydest_ovrf = destination->ycrd_ovrf;
    zdest_ovrf = destination->zcrd_ovrf;
    dest_scale = destination->gpos_scale;
    break;
  case TrajectoryKind::VELOCITIES:
    xdest = destination->xvel;
    ydest = destination->yvel;
    zdest = destination->zvel;
    xdest_ovrf = destination->xvel_ovrf;
    ydest_ovrf = destination->yvel_ovrf;
    zdest_ovrf = destination->zvel_ovrf;
    dest_scale = destination->vel_scale;
    break;
  case TrajectoryKind::FORCES:
    xdest = destination->xfrc;
    ydest = destination->yfrc;
    zdest = destination->zfrc;
    xdest_ovrf = destination->xfrc_ovrf;
    ydest_ovrf = destination->yfrc_ovrf;
    zdest_ovrf = destination->zfrc_ovrf;
    dest_scale = destination->frc_scale;
    break;
  }

  // Perform the familiar nested switches over destination and origin data tiers, with fall-through
  // to an HPC kernel launch.
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdims, destination->umat, destination->invu, index_dest,
                           origin.boxdims, origin.umat, origin.invu, index_orig,
                           destination->boxvecs, destination->boxvec_ovrf);
        break;
      case TrajectoryKind::VELOCITIES:
      case TrajectoryKind::FORCES:
        break;
      }
      switch (origin.mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&xdest[dest_atom_start], &xdest_ovrf[dest_atom_start],
                          &ydest[dest_atom_start], &ydest_ovrf[dest_atom_start],
                          &zdest[dest_atom_start], &zdest_ovrf[dest_atom_start], dest_scale,
                          &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                          &origin.zcrd[orig_atom_start], natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&xdest[dest_atom_start], &xdest_ovrf[dest_atom_start],
                          &ydest[dest_atom_start], &ydest_ovrf[dest_atom_start],
                          &zdest[dest_atom_start], &zdest_ovrf[dest_atom_start], dest_scale,
                          &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                          &origin.zcrd_sp[orig_atom_start], natom);
        break;
      }

      // Return after completing the host-to-host copy.
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

  // Funnel any devicce-related transfers into this code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_start = xfrm_w * index_dest;
  const int orig_xfrm_start = xfrm_w * index_orig;
  const int dest_bdim_start = bdim_w * index_dest;
  const int orig_bdim_start = bdim_w * index_orig;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    switch (origin.mode) {
    case PrecisionModel::DOUBLE:
      launchCopyCoordinateXYZBox(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf, zdest_ovrf,
                                 dest_scale, destination->gpos_bits, origin.xcrd, origin.ycrd,
                                 origin.zcrd, 1.0, 0, double_type_index, natom, destination->umat,
                                 destination->invu, destination->boxdims, origin.umat, origin.invu,
                                 origin.boxdims, destination->boxvecs, destination->boxvec_ovrf,
                                 dest_atom_start, orig_atom_start, dest_xfrm_start,
                                 orig_xfrm_start, dest_bdim_start, orig_bdim_start, gpu);
      break;
    case PrecisionModel::SINGLE:
      launchCopyCoordinateXYZBox(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf, zdest_ovrf,
                                 dest_scale, destination->gpos_bits, origin.xcrd_sp,
                                 origin.ycrd_sp, origin.zcrd_sp, 1.0, 0, float_type_index, natom,
                                 destination->umat, destination->invu, destination->boxdims,
                                 origin.umat, origin.invu, origin.boxdims, destination->boxvecs,
                                 destination->boxvec_ovrf, dest_atom_start, orig_atom_start,
                                 dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                                 orig_bdim_start, gpu);
      break;
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    switch (origin.mode) {
    case PrecisionModel::DOUBLE:
      launchCopyCoordinateXYZ(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf, zdest_ovrf,
                              dest_scale, destination->gpos_bits, origin.xcrd, origin.ycrd,
                              origin.zcrd, 1.0, 0, double_type_index, natom, dest_atom_start,
                              orig_atom_start, gpu);
      break;
    case PrecisionModel::SINGLE:
      launchCopyCoordinateXYZ(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf, zdest_ovrf,
                              dest_scale, destination->gpos_bits, origin.xcrd_sp, origin.ycrd_sp,
                              origin.zcrd_sp, 1.0, 0, double_type_index, natom, dest_atom_start,
                              orig_atom_start, gpu);
      break;
    }
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int dest_atom_start, const int index_dest,
               const CondensateReader &origin, const int orig_atom_start, const int index_orig,
               const int natom, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, dest_atom_start, index_dest, TrajectoryKind::POSITIONS, origin,
            orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int dest_atom_start = destination->getAtomOffset(index_dest);
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  const int natom = origin.getAtomCount(index_orig);
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PsSynthesisWriter poly_psw = destination->data(orientation, destination_tier);
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter poly_psw = destination->deviceViewToHostData(orientation);
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter poly_psw = destination->data(orientation, destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.deviceViewToHostData(),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, kind, destination->getCyclePosition(), origin, index_orig,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const Condensate &origin,
               const int index_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            origin, index_orig, destination_tier, origin_tier, gpu, sync);
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int dest_atom_start, const int index_dest,
               const PsSynthesisReader &origin, const int orig_atom_start, const int index_orig,
               const int natom, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  launchPreparation(sync, destination_tier, origin_tier);
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      // Copy the present positions, velocities, and forces
      copyBoxInformation(destination->boxdims, destination->umat, destination->invu, index_dest,
                         origin.boxdims, origin.umat, origin.invu, index_orig,
                         destination->boxvecs, destination->boxvec_ovrf);
      copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                        &destination->xcrd_ovrf[dest_atom_start],
                        &destination->ycrd[dest_atom_start],
                        &destination->ycrd_ovrf[dest_atom_start],
                        &destination->zcrd[dest_atom_start],
                        &destination->zcrd_ovrf[dest_atom_start], destination->gpos_scale,
                        &origin.xcrd[orig_atom_start], &origin.xcrd_ovrf[orig_atom_start],
                        &origin.ycrd[orig_atom_start], &origin.ycrd_ovrf[orig_atom_start],
                        &origin.zcrd[orig_atom_start], &origin.zcrd_ovrf[orig_atom_start],
                        origin.gpos_scale, natom);
      copyCoordinateXYZ(&destination->xvel[dest_atom_start],
                        &destination->xvel_ovrf[dest_atom_start],
                        &destination->yvel[dest_atom_start],
                        &destination->yvel_ovrf[dest_atom_start],
                        &destination->zvel[dest_atom_start],
                        &destination->zvel_ovrf[dest_atom_start], destination->vel_scale,
                        &origin.xvel[orig_atom_start], &origin.xvel_ovrf[orig_atom_start],
                        &origin.yvel[orig_atom_start], &origin.yvel_ovrf[orig_atom_start],
                        &origin.zvel[orig_atom_start], &origin.zvel_ovrf[orig_atom_start],
                        origin.vel_scale, natom);
      copyCoordinateXYZ(&destination->xfrc[dest_atom_start],
                        &destination->xfrc_ovrf[dest_atom_start],
                        &destination->yfrc[dest_atom_start],
                        &destination->yfrc_ovrf[dest_atom_start],
                        &destination->zfrc[dest_atom_start],
                        &destination->zfrc_ovrf[dest_atom_start], destination->frc_scale,
                        &origin.xfrc[orig_atom_start], &origin.xfrc_ovrf[orig_atom_start],
                        &origin.yfrc[orig_atom_start], &origin.yfrc_ovrf[orig_atom_start],
                        &origin.zfrc[orig_atom_start], &origin.zfrc_ovrf[orig_atom_start],
                        origin.frc_scale, natom);

      // Copy the alternate positions, velocities, and forces
      copyBoxInformation(destination->alt_boxdims, destination->umat_alt, destination->invu_alt,
                         index_dest, origin.alt_boxdims, origin.umat_alt, origin.invu_alt,
                         index_orig, destination->alt_boxvecs, destination->alt_boxvec_ovrf);
      copyCoordinateXYZ(&destination->xalt[dest_atom_start],
                        &destination->xalt_ovrf[dest_atom_start],
                        &destination->yalt[dest_atom_start],
                        &destination->yalt_ovrf[dest_atom_start],
                        &destination->zalt[dest_atom_start],
                        &destination->zalt_ovrf[dest_atom_start], destination->gpos_scale,
                        &origin.xalt[orig_atom_start], &origin.xalt_ovrf[orig_atom_start],
                        &origin.yalt[orig_atom_start], &origin.yalt_ovrf[orig_atom_start],
                        &origin.zalt[orig_atom_start], &origin.zalt_ovrf[orig_atom_start],
                        origin.gpos_scale, natom);
      copyCoordinateXYZ(&destination->vxalt[dest_atom_start],
                        &destination->vxalt_ovrf[dest_atom_start],
                        &destination->vyalt[dest_atom_start],
                        &destination->vyalt_ovrf[dest_atom_start],
                        &destination->vzalt[dest_atom_start],
                        &destination->vzalt_ovrf[dest_atom_start], destination->vel_scale,
                        &origin.vxalt[orig_atom_start], &origin.vxalt_ovrf[orig_atom_start],
                        &origin.vyalt[orig_atom_start], &origin.vyalt_ovrf[orig_atom_start],
                        &origin.vzalt[orig_atom_start], &origin.vzalt_ovrf[orig_atom_start],
                        origin.vel_scale, natom);
      copyCoordinateXYZ(&destination->fxalt[dest_atom_start],
                        &destination->fxalt_ovrf[dest_atom_start],
                        &destination->fyalt[dest_atom_start],
                        &destination->fyalt_ovrf[dest_atom_start],
                        &destination->fzalt[dest_atom_start],
                        &destination->fzalt_ovrf[dest_atom_start], destination->frc_scale,
                        &origin.fxalt[orig_atom_start], &origin.fxalt_ovrf[orig_atom_start],
                        &origin.fyalt[orig_atom_start], &origin.fyalt_ovrf[orig_atom_start],
                        &origin.fzalt[orig_atom_start], &origin.fzalt_ovrf[orig_atom_start],
                        origin.frc_scale, natom);

      // Return after completing the host-to-host copy
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchCopyCoordinates(destination, dest_atom_start, index_dest, origin, orig_atom_start,
                            index_orig, gpu);
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchCopyCoordinates(destination, dest_atom_start, index_dest, origin, orig_atom_start,
                          index_orig, gpu);
    break;
#endif
  }
  launchResolution(sync, destination_tier, origin_tier);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int dest_atom_start = destination->getAtomOffset(index_dest);
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  const int natom = origin.getAtomCount(index_orig);
  coordCopyValidateAtomCounts(natom, destination->getAtomCount(index_dest));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PsSynthesisWriter poly_psw = destination->data(destination_tier);
        coordCopy(&poly_psw, dest_atom_start, index_dest, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter poly_psw = destination->deviceViewToHostData();
        coordCopy(&poly_psw, dest_atom_start, index_dest, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter poly_psw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&poly_psw, dest_atom_start, index_dest, origin.deviceViewToHostData(),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&poly_psw, dest_atom_start, index_dest, origin.data(origin_tier),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const CoordinateFrameReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdims, destination->umat, destination->invu, index_dest,
                         origin.boxdim, origin.umat, origin.invu);
      switch (destination->mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], origin.xcrd, origin.ycrd,
                          origin.zcrd, origin.natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                          &destination->ycrd_sp[dest_atom_start],
                          &destination->zcrd_sp[dest_atom_start], origin.xcrd, origin.ycrd,
                          origin.zcrd, origin.natom);
        break;
      }

      // Return after completing the host-to-host copy
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

  // Any copy operations involving device data will pass through the switch above and execute the
  // following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int dest_xfrm_start = roundUp(9, warp_size_int) * index_dest;
  const int dest_bdim_start = roundUp(6, warp_size_int) * index_dest;
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, origin.natom, destination->umat,
                               destination->invu, destination->boxdims, origin.umat, origin.invu,
                               origin.boxdim, dest_atom_start, 0, dest_xfrm_start, 0,
                               dest_bdim_start, 0, gpu);
    break;
  case PrecisionModel::SINGLE:
    launchCopyCoordinateXYZBox(destination->xcrd_sp, destination->ycrd_sp, destination->zcrd_sp,
                               1.0, float_type_index, origin.xcrd, origin.ycrd, origin.zcrd, 1.0,
                               double_type_index, origin.natom, destination->umat,
                               destination->invu, destination->boxdims, origin.umat, origin.invu,
                               origin.boxdim, dest_atom_start, 0, dest_xfrm_start, 0,
                               dest_bdim_start, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const CoordinateFrame &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  const int dest_atom_start = destination->getAtomOffset(index_dest);
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.getAtomCount());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter cdnsw = destination->data(destination_tier);
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter cdnsw = destination->deviceViewToHostData();
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter cdnsw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.deviceViewToHostData(),
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const PhaseSpaceReader &origin, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        {
          copyBoxInformation(destination->boxdims, destination->umat, destination->invu,
                             index_dest, origin.boxdim, origin.umat, origin.invu);
          const double* xorig = origin.xcrd;
          const double* yorig = origin.ycrd;
          const double* zorig = origin.zcrd;
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          }
        }
        break;
      case TrajectoryKind::VELOCITIES:
        {
          const double* xorig = origin.xvel;
          const double* yorig = origin.yvel;
          const double* zorig = origin.zvel;
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          }
        }
        break;
      case TrajectoryKind::FORCES:
        {
          const double* xorig = origin.xfrc;
          const double* yorig = origin.yfrc;
          const double* zorig = origin.zfrc;
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, yorig, zorig,
                              origin.natom);
            break;
          }
        }
        break;
      }

      // Return after completing the host-to-host copy
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

  // Any copy operations involving data on the device will pass through the above switch and
  // execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int dest_xfrm_start = roundUp(9, warp_size_int) * index_dest;
  const int dest_bdim_start = roundUp(6, warp_size_int) * index_dest;
  size_t ct_dest;
  void *dest_xptr, *dest_yptr, *dest_zptr;
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd);
    ct_dest = double_type_index;
    break;
  case PrecisionModel::SINGLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd_sp);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd_sp);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd_sp);
    ct_dest = float_type_index;
    break;
  }
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(dest_xptr, dest_yptr, dest_zptr, 1.0, ct_dest, origin.xcrd,
                               origin.ycrd, origin.zcrd, 1.0, double_type_index, origin.natom,
                               destination->umat, destination->invu, destination->boxdims,
                               origin.umat, origin.invu, origin.boxdim, dest_atom_start, 0,
                               dest_xfrm_start, 0, dest_bdim_start, 0, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(dest_xptr, dest_yptr, dest_zptr, 1.0, ct_dest, origin.xvel,
                            origin.yvel, origin.zvel, 1.0, double_type_index, origin.natom,
                            dest_atom_start, 0, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(dest_xptr, dest_yptr, dest_zptr, 1.0, ct_dest, origin.xfrc,
                            origin.yfrc, origin.zfrc, 1.0, double_type_index, origin.natom,
                            dest_atom_start, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const PhaseSpaceReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, dest_atom_start, index_dest, origin, TrajectoryKind::POSITIONS,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin,
               const TrajectoryKind kind, const CoordinateCycle orientation,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.getAtomCount());
  const int dest_atom_start = destination->getAtomOffset(index_dest);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter cdnsw = destination->data(destination_tier);
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier), kind,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter cdnsw = destination->deviceViewToHostData();
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier), kind,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter cdnsw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.deviceViewToHostData(orientation),
                  kind, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier),
                  kind, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin,
               const TrajectoryKind kind, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, origin, kind, origin.getCyclePosition());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, origin, TrajectoryKind::POSITIONS, origin.getCyclePosition());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const PsSynthesisReader &origin, const int orig_atom_start, const int index_orig,
               const int natom, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        {
          copyBoxInformation(destination->boxdims, destination->umat, destination->invu,
                             index_dest, origin.boxdims, origin.umat, origin.invu, index_orig);
          const llint* xorig = &origin.xcrd[orig_atom_start];
          const llint* yorig = &origin.ycrd[orig_atom_start];
          const llint* zorig = &origin.zcrd[orig_atom_start];
          const int* xorig_ovrf = &origin.xcrd_ovrf[orig_atom_start];
          const int* yorig_ovrf = &origin.ycrd_ovrf[orig_atom_start];
          const int* zorig_ovrf = &origin.zcrd_ovrf[orig_atom_start];
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.gpos_scale, natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.gpos_scale, natom);
            break;
          }
        }
        break;
      case TrajectoryKind::VELOCITIES:
        {
          const llint* xorig = &origin.xvel[orig_atom_start];
          const llint* yorig = &origin.yvel[orig_atom_start];
          const llint* zorig = &origin.zvel[orig_atom_start];
          const int* xorig_ovrf = &origin.xvel_ovrf[orig_atom_start];
          const int* yorig_ovrf = &origin.yvel_ovrf[orig_atom_start];
          const int* zorig_ovrf = &origin.zvel_ovrf[orig_atom_start];
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.vel_scale, natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.vel_scale, natom);
            break;
          }
        }
        break;
      case TrajectoryKind::FORCES:
        {
          const llint* xorig = &origin.xfrc[orig_atom_start];
          const llint* yorig = &origin.yfrc[orig_atom_start];
          const llint* zorig = &origin.zfrc[orig_atom_start];
          const int* xorig_ovrf = &origin.xfrc_ovrf[orig_atom_start];
          const int* yorig_ovrf = &origin.yfrc_ovrf[orig_atom_start];
          const int* zorig_ovrf = &origin.zfrc_ovrf[orig_atom_start];
          switch (destination->mode) {
          case PrecisionModel::DOUBLE:
            copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                              &destination->ycrd[dest_atom_start],
                              &destination->zcrd[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.frc_scale, natom);
            break;
          case PrecisionModel::SINGLE:
            copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                              &destination->ycrd_sp[dest_atom_start],
                              &destination->zcrd_sp[dest_atom_start], xorig, xorig_ovrf, yorig,
                              yorig_ovrf, zorig, zorig_ovrf, origin.frc_scale, natom);
            break;
          }
        }
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

  // Any copy procedures involving data on the GPU will fall through and then execute the code
  // below.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_start = xfrm_w * index_dest;
  const int dest_bdim_start = bdim_w * index_dest;
  const int orig_xfrm_start = xfrm_w * index_orig;
  const int orig_bdim_start = bdim_w * index_orig;
  size_t ct_dest;
  void *dest_xptr, *dest_yptr, *dest_zptr;
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd);
    ct_dest = double_type_index;
    break;
  case PrecisionModel::SINGLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd_sp);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd_sp);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd_sp);
    ct_dest = float_type_index;
    break;
  }
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(dest_xptr, dest_yptr, dest_zptr, 1.0, 0, ct_dest, origin.xcrd,
                               origin.ycrd, origin.zcrd, origin.xcrd_ovrf, origin.ycrd_ovrf,
                               origin.zcrd_ovrf, origin.gpos_scale, origin.gpos_bits, natom,
                               destination->umat, destination->invu, destination->boxdims,
                               origin.umat, origin.invu, origin.boxdims, dest_atom_start,
                               orig_atom_start, dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                               orig_bdim_start, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(dest_xptr, dest_yptr, dest_zptr, 1.0, 0, ct_dest, origin.xvel,
                            origin.yvel, origin.zvel, origin.xvel_ovrf, origin.yvel_ovrf,
                            origin.zvel_ovrf, origin.vel_scale, origin.vel_bits, natom,
                            dest_atom_start, orig_atom_start, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(dest_xptr, dest_yptr, dest_zptr, 1.0, 0, ct_dest, origin.xfrc,
                            origin.yfrc, origin.zfrc, origin.xfrc_ovrf, origin.yfrc_ovrf,
                            origin.zfrc_ovrf, origin.frc_scale, origin.frc_bits, natom,
                            dest_atom_start, orig_atom_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const PsSynthesisReader &origin, const int orig_atom_start, const int index_orig,
               const int natom, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, dest_atom_start, index_dest, origin, orig_atom_start, index_orig, natom,
            TrajectoryKind::POSITIONS, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
               const int index_orig, const TrajectoryKind kind, CoordinateCycle orientation,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int natom = origin.getAtomCount(index_orig);
  const size_t dest_atom_start = destination->getAtomOffset(index_dest);
  const size_t orig_atom_start = origin.getAtomOffset(index_dest);
  coordCopyValidateAtomCounts(natom, destination->getAtomCount(index_orig));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter cdnsw = destination->data(destination_tier);
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier),
                  orig_atom_start, index_orig, natom, kind, destination_tier, origin_tier, gpu,
                  sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter cdnsw = destination->deviceViewToHostData();
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier),
                  orig_atom_start, index_orig, natom, kind, destination_tier, origin_tier, gpu,
                  sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter cdnsw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.deviceViewToHostData(orientation),
                  orig_atom_start, index_orig, natom, kind, destination_tier, origin_tier, gpu,
                  sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(orientation, origin_tier),
                  orig_atom_start, index_orig, natom, kind, destination_tier, origin_tier, gpu,
                  sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
	       const int index_orig, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, origin, index_orig, kind, origin.getCyclePosition(),
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
               const int index_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, index_dest, origin, index_orig, TrajectoryKind::POSITIONS,
            origin.getCyclePosition(), destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int dest_atom_start, const int index_dest,
               const CondensateReader &origin, const int orig_atom_start, const int index_orig,
               const int natom, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:

      // Assume that the data concerns positions and copy whatever box information is present.
      copyBoxInformation(destination->boxdims, destination->umat, destination->invu, index_dest,
                         origin.boxdims, origin.umat, origin.invu, index_orig);
      switch (destination->mode) {
      case PrecisionModel::DOUBLE:
        switch (origin.mode) {
        case PrecisionModel::DOUBLE:
          copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                            &destination->ycrd[dest_atom_start],
                            &destination->zcrd[dest_atom_start], &origin.xcrd[orig_atom_start],
                            &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], natom);
          break;
        case PrecisionModel::SINGLE:
          copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                            &destination->ycrd[dest_atom_start],
                            &destination->zcrd[dest_atom_start], &origin.xcrd_sp[orig_atom_start],
                            &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start],
                            natom);
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (destination->mode) {
        case PrecisionModel::DOUBLE:
          copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                            &destination->ycrd_sp[dest_atom_start],
                            &destination->zcrd_sp[dest_atom_start], &origin.xcrd[orig_atom_start],
                            &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], natom);
          break;
        case PrecisionModel::SINGLE:
          copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                            &destination->ycrd_sp[dest_atom_start],
                            &destination->zcrd_sp[dest_atom_start],
                            &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                            &origin.zcrd_sp[orig_atom_start], natom);
          break;
        }
        break;
      }

      // Return after completing the host-to-host copy procedure.
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

  // Any oeprations concerning data on the GPU will fall through the switch above and execute the
  // following code.  Unroll the switch over the destination object's precision mode.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_start = xfrm_w * index_dest;
  const int dest_bdim_start = bdim_w * index_dest;
  const int orig_xfrm_start = xfrm_w * index_dest;
  const int orig_bdim_start = bdim_w * index_dest;
  size_t ct_dest;
  void *dest_xptr, *dest_yptr, *dest_zptr;
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd);
    ct_dest = double_type_index;
    break;
  case PrecisionModel::SINGLE:
    dest_xptr = reinterpret_cast<void*>(destination->xcrd_sp);
    dest_yptr = reinterpret_cast<void*>(destination->ycrd_sp);
    dest_zptr = reinterpret_cast<void*>(destination->zcrd_sp);
    ct_dest = float_type_index;
    break;
  }
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    launchCopyCoordinateXYZBox(dest_xptr, dest_yptr, dest_zptr, 1.0, ct_dest, origin.xcrd,
                               origin.ycrd, origin.zcrd, 1.0, double_type_index, natom,
                               destination->umat, destination->invu, destination->boxdims,
                               origin.umat, origin.invu, origin.boxdims, dest_atom_start,
                               orig_atom_start, dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                               orig_bdim_start, gpu);
    break;
  case PrecisionModel::SINGLE:
    launchCopyCoordinateXYZBox(dest_xptr, dest_yptr, dest_zptr, 1.0, ct_dest, origin.xcrd_sp,
                               origin.ycrd_sp, origin.zcrd_sp, 1.0, float_type_index, natom,
                               destination->umat, destination->invu, destination->boxdims,
                               origin.umat, origin.invu, origin.boxdims, dest_atom_start,
                               orig_atom_start, dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                               orig_bdim_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const Condensate &origin,
               const int index_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  const int natom = origin.getAtomCount(index_dest);
  const size_t dest_atom_start = destination->getAtomOffset(index_dest);
  const size_t orig_atom_start = origin.getAtomOffset(index_dest);
  coordCopyValidateAtomCounts(natom, destination->getAtomCount(index_dest));
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter cdnsw = destination->data(destination_tier);
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), orig_atom_start,
                  index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter cdnsw = destination->deviceViewToHostData();
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), orig_atom_start,
                  index_orig, natom, destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter cdnsw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.deviceViewToHostData(),
                  orig_atom_start, index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), orig_atom_start,
                  index_orig, natom, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

} // namespace trajectory
} // namespace stormm
