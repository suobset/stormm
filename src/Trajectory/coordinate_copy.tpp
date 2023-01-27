// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const int natom) {

  // Check the data types
  if (isFloatingPointScalarType<Tdest>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Tdest>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be filled without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  for (int i = 0; i < natom; i++) {
    xdest[i] = xorig[i];
    ydest[i] = yorig[i];
    zdest[i] = zorig[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       const int natom) {

  // The origin data must be of some floating point type
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be filled without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
  if (ct_dest == llint_type_index) {
    for (int i = 0; i < natom; i++) {
      xdest[i] = llround(xorig[i] * dest_scale);
      ydest[i] = llround(yorig[i] * dest_scale);
      zdest[i] = llround(zorig[i] * dest_scale);
    }
  }
  else if (ct_dest == int_type_index) {
    for (int i = 0; i < natom; i++) {
      xdest[i] = round(xorig[i] * dest_scale);
      ydest[i] = round(yorig[i] * dest_scale);
      zdest[i] = round(zorig[i] * dest_scale);
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      xdest[i] = xorig[i] * dest_scale;
      ydest[i] = yorig[i] * dest_scale;
      zdest[i] = zorig[i] * dest_scale;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const double orig_scale,
                       const int natom) {
  const double inv_orig_scale = 1.0 / orig_scale;
  for (int i = 0; i < natom; i++) {
    xdest[i] = static_cast<double>(xorig[i]) * inv_orig_scale;
    ydest[i] = static_cast<double>(yorig[i]) * inv_orig_scale;
    zdest[i] = static_cast<double>(zorig[i]) * inv_orig_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       const double orig_scale, const int natom, const HybridTargetLevel dest_tier,
                       const HybridTargetLevel orig_tier) {
  if (isSignedIntegralScalarType<Tdest>() && isSignedIntegralScalarType<Torig>()) {
    const int dest_bits = round(log2(dest_scale));
    const int orig_bits = round(log2(orig_scale));
    if (orig_bits < dest_bits) {
      const int shft_bits = dest_bits - orig_bits;
      for (int i = 0; i < natom; i++) {
        xdest[i] = xorig[i];
        ydest[i] = yorig[i];
        zdest[i] = zorig[i];
        xdest[i] <<= shft_bits;
        ydest[i] <<= shft_bits;
        zdest[i] <<= shft_bits;
      }
    }
    else if (orig_bits > dest_bits) {
      const int shft_bits = orig_bits - dest_bits;
      for (int i = 0; i < natom; i++) {
        xdest[i] = xorig[i];
        ydest[i] = yorig[i];
        zdest[i] = zorig[i];
        xdest[i] >>= shft_bits;
        ydest[i] >>= shft_bits;
        zdest[i] >>= shft_bits;
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        xdest[i] = xorig[i];
        ydest[i] = yorig[i];
        zdest[i] = zorig[i];
      }
    }
  }
  else {
    const double conv_factor = dest_scale / orig_scale;
    for (int i = 0; i < natom; i++) {
      xdest[i] = static_cast<double>(xorig[i]) * conv_factor;
      ydest[i] = static_cast<double>(yorig[i]) * conv_factor;
      zdest[i] = static_cast<double>(zorig[i]) * conv_factor;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const size_t dest_offset, const Torig* xorig, const Torig* yorig,
                       const Torig* zorig, const double orig_scale, const size_t orig_offset,
                       const int natom, const HybridTargetLevel dest_tier,
                       const HybridTargetLevel orig_tier, const GpuDetails &gpu) {
  switch (dest_tier) {
  case HybridTargetLevel::HOST:
    switch (orig_tier) {
    case HybridTargetLevel::HOST:

      // The host-to-host case can just be a call to the ordinary C++ function.  All other cases
      // will fire off the templated kernel. 
      copyCoordinateXYZ(&xdest[dest_offset], &ydest[dest_offset], &zdest[dest_offset], dest_scale,
                        &xorig[orig_offset], &yorig[orig_offset], &zorig[orig_offset], orig_scale,
                        natom);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
        const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
        void* vd_xdest = reinterpret_cast<void*>(xdest);
        void* vd_ydest = reinterpret_cast<void*>(ydest);
        void* vd_zdest = reinterpret_cast<void*>(zdest);
        void* vd_xorig = reinterpret_cast<void*>(xorig);
        void* vd_yorig = reinterpret_cast<void*>(yorig);
        void* vd_zorig = reinterpret_cast<void*>(zorig);
        launchCopyCoordinateXYZ(vd_xdest, vd_ydest, vd_zdest, dest_scale, dest_offset, ct_dest,
                                vd_xorig, vd_yorig, vd_zorig, orig_scale, orig_offset, ct_orig,
                                natom, gpu);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
      const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
      void* vd_xdest = reinterpret_cast<void*>(xdest);
      void* vd_ydest = reinterpret_cast<void*>(ydest);
      void* vd_zdest = reinterpret_cast<void*>(zdest);
      void* vd_xorig = reinterpret_cast<void*>(xorig);
      void* vd_yorig = reinterpret_cast<void*>(yorig);
      void* vd_zorig = reinterpret_cast<void*>(zorig);
      launchCopyCoordinateXYZ(vd_xdest, vd_ydest, vd_zdest, dest_scale, dest_offset, ct_dest,
                              vd_xorig, vd_yorig, vd_zorig, orig_scale, orig_offset, ct_orig,
                              natom, gpu);
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const	int* zorig_ovrf,
                       const double orig_scale, const int natom) {
  if (isFloatingPointScalarType<Tdest>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Tdest>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  const double inv_orig_scale = 1.0 / orig_scale;
  for (int i = 0; i < natom; i++) {
    xdest[i] = hostInt95ToDouble(xorig[i], xorig_ovrf[i]) * inv_orig_scale;
    ydest[i] = hostInt95ToDouble(yorig[i], yorig_ovrf[i]) * inv_orig_scale;
    zdest[i] = hostInt95ToDouble(zorig[i], zorig_ovrf[i]) * inv_orig_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const	int* zorig_ovrf,
                       const double orig_scale, const int natom) {
  const size_t ct = std::type_index(typeid(Tdest)).hash_code();
  if (ct == llint_type_index) {

    // Short of the 95-bit integer format, a long long integer representation is the only format
    // that has any possibility of exceeding the information content of double-precision floating
    // point numbers.  Use a full int95_t bit conversion procedure and then take the 64-bit
    // component, which must hold all of the information or the llint format is busted anyway.
    const int dest_bits = round(log2(dest_scale));
    const int orig_bits = round(log2(orig_scale));
    for (int i = 0; i < natom; i++) {
      const int95_t xo = { xorig[i], xorig_ovrf[i] };
      const int95_t yo = { yorig[i], yorig_ovrf[i] };
      const int95_t zo = { zorig[i], zorig_ovrf[i] };
      const int95_t xn = hostChangeFPBits(xo, orig_bits, dest_bits);
      const int95_t yn = hostChangeFPBits(xo, orig_bits, dest_bits);
      const int95_t zn = hostChangeFPBits(xo, orig_bits, dest_bits);
      xdest[i] = xn.x;
      ydest[i] = yn.x;
      zdest[i] = zn.x;
    }
  }
  else {
    const double conv_scale = dest_scale / orig_scale;
    for (int i = 0; i < natom; i++) {
      xdest[i] = hostInt95ToDouble(xorig[i], xorig_ovrf[i]) * conv_scale;
      ydest[i] = hostInt95ToDouble(yorig[i], yorig_ovrf[i]) * conv_scale;
      zdest[i] = hostInt95ToDouble(zorig[i], zorig_ovrf[i]) * conv_scale;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, const double dest_scale, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const int natom) {
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  for (int i = 0; i < natom; i++) {
    const int95_t xn = hostDoubleToInt95(xorig[i] * dest_scale);
    xdest[i] = xn.x;
    xdest_ovrf[i] = xn.y;
    const int95_t yn = hostDoubleToInt95(yorig[i] * dest_scale);
    ydest[i] = yn.x;
    ydest_ovrf[i] = yn.y;
    const int95_t zn = hostDoubleToInt95(zorig[i] * dest_scale);
    zdest[i] = zn.x;
    zdest_ovrf[i] = zn.y;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, const double dest_scale, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const double orig_scale,
                       const int natom) {
  const size_t ct = std::type_index(typeid(Torig)).hash_code();
  if (ct == llint_type_index) {
    const int orig_bits = round(log2(orig_scale));
    const int dest_bits = round(log2(dest_scale));
    
    // Again, handle the long long integer format with a special case as it is the only type that
    // might lose bits in a conversion via double-precision numbers.
    for (int i = 0; i < natom; i++) {
      const int95_t xo = { xorig[i], 0 };
      const int95_t yo = { yorig[i], 0 };
      const int95_t zo = { zorig[i], 0 };
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
  else {
    const double conv_factor = dest_scale / orig_scale;
    for (int i = 0; i < natom; i++) {
      const int95_t xn = hostDoubleToInt95(xorig[i] * conv_factor);
      xdest[i] = xn.x;
      xdest_ovrf[i] = xn.y;
      const int95_t yn = hostDoubleToInt95(yorig[i] * conv_factor);
      ydest[i] = yn.x;
      ydest_ovrf[i] = yn.y;
      const int95_t zn = hostDoubleToInt95(zorig[i] * conv_factor);
      zdest[i] = zn.x;
      zdest_ovrf[i] = zn.y;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateFrame *destination, const CoordinateSeries<T> &origin,
               int frame_orig) {
  origin.extractFrame(destination, frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const CoordinateSeries<T> &origin,
               int frame_orig) {
  origin.extractFrame(destination, TrajectoryKind::POSITIONS, CoordinateCycle::PRIMARY,
                      frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateSeries<T> &origin, int frame_orig) {
  origin.extractFrame(destination, kind, CoordinateCycle::PRIMARY, frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<T> &origin,
               int frame_orig) {
  origin.extractFrame(destination, kind, orientation, frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const CoordinateFrameReader &origin) {
  destination->import(origin, frame_dest);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const CoordinateFrame &origin) {
  destination->import(origin, frame_dest);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const PhaseSpaceReader &origin) {
  destination->import(origin, frame_dest);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const PhaseSpace &origin, const TrajectoryKind kind,
               const CoordinateCycle orientation) {
  destination->import(origin, frame_dest, kind, orientation);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, const size_t frame_dest,
               const CoordinateSeries<Torig> &origin, size_t frame_orig) {
  const size_t dest_atom_start = frame_dest * roundUp<size_t>(destination->natom, warp_size_zu);
  const size_t orig_atom_start = frame_orig * roundUp<size_t>(origin.natom, warp_size_zu);
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  copyCoordinateXYZ<Tdest, Torig>(&destination->xcrd[dest_atom_start],
                                  &destination->ycrd[dest_atom_start],
                                  &destination->zcrd[dest_atom_start], destination->gpos_scale,
                                  &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                                  &origin.zcrd[orig_atom_start], origin.gpos_scale, origin.natom);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, const int frame_dest,
               const CoordinateSeries<Torig> &origin, int frame_orig) {
  coordCopy<Tdest, Torig>(destination->data(), frame_dest, origin.data(), frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, const size_t frame_dest,
               const PsSynthesisReader &origin, const size_t index_orig,
               const TrajectoryKind kind) {
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->natom;
  coordCopyValidateAtomCounts(natom, origin.atom_counts[index_orig]);
  const size_t dest_atom_start = roundUp<size_t>(destination->natom, warp_size_zu) * frame_dest;
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const size_t bdim_stride = roundUp(6, warp_size_int);
      const size_t xfrm_stride = roundUp(9, warp_size_int);
      const size_t dest_bdim_start = frame_dest * bdim_stride;
      const size_t dest_xfrm_start = frame_dest * xfrm_stride;
      const size_t orig_bdim_start = index_orig * bdim_stride;
      const size_t orig_xfrm_start = index_orig * xfrm_stride;
      copyBoxInformation(&destination->boxdims[dest_bdim_start],
                         &destination->umat[dest_xfrm_start], &destination->invu[dest_xfrm_start],
                         &origin.boxdims[orig_bdim_start], &origin.umat[orig_xfrm_start],
                         &origin.invu[orig_xfrm_start]);
      copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                        &destination->zcrd[dest_atom_start], destination->gpos_scale,
                        &origin.xcrd[orig_atom_start], &origin.xcrd_ovrf[orig_atom_start],
                        &origin.ycrd[orig_atom_start], &origin.ycrd_ovrf[orig_atom_start],
                        &origin.zcrd[orig_atom_start], &origin.zcrd_ovrf[orig_atom_start],
                        origin.gpos_scale, natom);
    }
    break;
  case TrajectoryKind::VELOCITIES:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], destination->gpos_scale,
                      &origin.xvel[orig_atom_start], &origin.xvel_ovrf[orig_atom_start],
                      &origin.yvel[orig_atom_start], &origin.yvel_ovrf[orig_atom_start],
                      &origin.zvel[orig_atom_start], &origin.zvel_ovrf[orig_atom_start],
                      origin.vel_scale, natom);
    break;
  case TrajectoryKind::FORCES:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], destination->gpos_scale,
                      &origin.xfrc[orig_atom_start], &origin.xfrc_ovrf[orig_atom_start],
                      &origin.yfrc[orig_atom_start], &origin.yfrc_ovrf[orig_atom_start],
                      &origin.zfrc[orig_atom_start], &origin.zfrc_ovrf[orig_atom_start],
                      origin.frc_scale, natom);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig, const TrajectoryKind kind,
               const CoordinateCycle orientation) {
  CoordinateSeriesWriter<T> csw = destination->data();
  coordCopy(&csw, frame_dest, origin.data(orientation), index_orig, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, size_t frame_dest,
               const CondensateReader &origin, size_t index_orig) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  coordCopyValidateAtomCounts(destination->natom, origin.atom_counts[index_orig]);
  const int natom = origin.atom_counts[index_orig];
  const size_t dest_atom_start = frame_dest * roundUp<size_t>(destination->natom, warp_size_zu);
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  const size_t dest_bdim_start = frame_dest * roundUp<size_t>(6, warp_size_int);
  const size_t dest_xfrm_start = frame_dest * roundUp<size_t>(9, warp_size_int);
  copyBoxDimensions(&destination->boxdims[dest_bdim_start], &destination->umat[dest_xfrm_start],
                    &destination->invu[dest_xfrm_start], origin, index_orig);
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], destination->gpos_scale,
                      &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                      &origin.zcrd[orig_atom_start], natom);
    break;
  case PrecisionModel::SINGLE:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], destination->gpos_scale,
                      &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                      &origin.zcrd_sp[orig_atom_start], natom);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const Condensate &origin, const int index_orig) {
  CoordinateSeriesWriter<T> csw = destination->data();
  coordCopy(&csw, frame_dest, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PsSynthesisWriter *destination, const int index_dest,
               const CoordinateSeriesReader<T> &origin, const int frame_orig) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, CoordinateCycle::PRIMARY, origin,
            frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PsSynthesisWriter *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateSeriesReader<T> &origin, const int frame_orig) {
  coordCopy(destination, index_dest, kind, CoordinateCycle::PRIMARY, origin, frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PsSynthesisWriter *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeriesReader<T> &origin,
               const int frame_orig) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  coordCopyValidateAtomCounts(destination->atom_counts[index_dest], origin.natom);
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t orig_atom_start = roundUp<size_t>(origin.natom, warp_size_zu) *
                                 static_cast<size_t>(frame_orig);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const size_t bdim_stride = roundUp(6, warp_size_int);
      const size_t xfrm_stride = roundUp(9, warp_size_int);
      const size_t dest_bdim_start = index_dest * bdim_stride;
      const size_t dest_xfrm_start = index_dest * xfrm_stride;
      const size_t orig_bdim_start = frame_orig * bdim_stride;
      const size_t orig_xfrm_start = frame_orig * xfrm_stride;
      copyBoxInformation(&destination->boxdims[dest_bdim_start],
                         &destination->umat[dest_xfrm_start], &destination->invu[dest_xfrm_start],
                         &origin.boxdim[orig_bdim_start], &origin.umat[orig_xfrm_start],
                         &origin.invu[orig_xfrm_start]);
      copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                        &destination->xcrd_ovrf[dest_atom_start],
                        &destination->ycrd[dest_atom_start],
                        &destination->ycrd_ovrf[dest_atom_start],
                        &destination->zcrd[dest_atom_start],
                        &destination->zcrd_ovrf[dest_atom_start],
                        destination->gpos_scale, &origin.xcrd[orig_atom_start],
                        &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                        origin.gpos_scale, origin.natom);
    }
    break;
  case TrajectoryKind::VELOCITIES:
    copyCoordinateXYZ(&destination->xvel[dest_atom_start],
                      &destination->xvel_ovrf[dest_atom_start],
                      &destination->yvel[dest_atom_start],
                      &destination->yvel_ovrf[dest_atom_start],
                      &destination->zvel[dest_atom_start],
                      &destination->zvel_ovrf[dest_atom_start],
                      destination->vel_scale, &origin.xcrd[orig_atom_start],
                      &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], 
                      origin.gpos_scale, origin.natom);
    break;
  case TrajectoryKind::FORCES:
    copyCoordinateXYZ(&destination->xfrc[dest_atom_start],
                      &destination->xfrc_ovrf[dest_atom_start],
                      &destination->yfrc[dest_atom_start],
                      &destination->yfrc_ovrf[dest_atom_start],
                      &destination->zfrc[dest_atom_start],
                      &destination->zfrc_ovrf[dest_atom_start],
                      destination->frc_scale, &origin.xcrd[orig_atom_start],
                      &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                      origin.gpos_scale, origin.natom);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateSeries<T> &origin, const int frame_orig) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, CoordinateCycle::PRIMARY, origin,
            frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateSeries<T> &origin, const int frame_orig) {
  coordCopy(destination, index_dest, kind, CoordinateCycle::PRIMARY, origin, frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<T> &origin,
               const int frame_orig) {
  PsSynthesisWriter poly_psw = destination->data(orientation);
  coordCopy(&poly_psw, index_dest, kind, orientation, origin.data(), frame_orig);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CondensateWriter *destination, int index_dest,
               const CoordinateSeriesReader<T> &origin, int frame_orig) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  coordCopyValidateAtomCounts(destination->atom_counts[index_dest], origin.natom);
  const int natom = destination->atom_counts[index_dest];
  const size_t orig_atom_start = static_cast<size_t>(frame_orig) *
                                 roundUp<size_t>(origin.natom, warp_size_zu);
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t orig_xfrm_start = static_cast<size_t>(frame_orig) * roundUp(9, warp_size_int);
  copyBoxInformation(destination, index_dest, &origin.umat[orig_xfrm_start],
                     &origin.invu[orig_xfrm_start]);
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], &origin.xcrd[orig_atom_start],
                      &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                      origin.gpos_scale, natom);
    break;
  case PrecisionModel::SINGLE:
    copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                      &destination->zcrd[dest_atom_start], &origin.xcrd_sp[orig_atom_start],
                      &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start],
                      origin.gpos_scale, natom);
    break;
  }
}

} // namespace trajectory
} // namespace stormm
