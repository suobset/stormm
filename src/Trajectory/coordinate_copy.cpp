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
                        const double* orig_boxdim, const double* orig_umat,
                        const double* orig_invu) {
  for (int i = 0; i < 6; i++) {
    dest_boxdim[i] = orig_boxdim[i];
  }
  for (int i = 0; i < 9; i++) {
    dest_umat[i] = orig_umat[i];
    dest_invu[i] = orig_invu[i];
  }
}

//-------------------------------------------------------------------------------------------------
void copyBoxInformation(double* dest_boxdim, double* dest_umat, double* dest_invu,
                        const CondensateReader &cdnsr, const int system_index) {

  // Extracting the box dimensions from a Condensate is slightly more involved, as there is no
  // explicit storage of the six-membered length and angle array.
  const size_t xfrm_start = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  extractBoxDimensions(&dest_boxdim[0], &dest_boxdim[1], &dest_boxdim[2], &dest_boxdim[3],
                       &dest_boxdim[4], &dest_boxdim[5], &cdnsr.invu[xfrm_start]);
  for (size_t i = 0; i < 9; i++) {
    dest_umat[i] = cdnsr.umat[xfrm_start + i];
    dest_invu[i] = cdnsr.invu[xfrm_start + i];
  }
}

//-------------------------------------------------------------------------------------------------
void copyBoxInformation(CondensateWriter *cdnsr, const int index_dest, const double* orig_umat,
                        const double* orig_invu) {
  const size_t xfrm_start = static_cast<size_t>(index_dest) * roundUp<size_t>(9, warp_size_zu);
  for (size_t i = 0; i < 9; i++) {
    cdnsr->umat[xfrm_start + i] = orig_umat[i];
    cdnsr->invu[xfrm_start + i] = orig_invu[i];
  }
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
void coordCopy(CoordinateFrameWriter *destination, const CoordinateFrameReader &origin) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                     origin.umat, origin.invu);
  copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd, origin.xcrd,
                    origin.ycrd, origin.zcrd, origin.natom);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const CoordinateFrame &origin) {
  CoordinateFrameWriter cfw = destination->data();
  coordCopy(&cfw, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const PhaseSpaceReader &origin,
               const TrajectoryKind kind) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                       origin.umat, origin.invu);
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
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpace &origin, const TrajectoryKind kind,
               const CoordinateCycle orientation) {
  CoordinateFrameWriter cfw = destination->data();
  coordCopy(&cfw, origin.data(orientation), kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const PsSynthesisReader &origin,
               const int index_orig, const TrajectoryKind kind) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const size_t bdim_start = static_cast<size_t>(index_orig) * roundUp<size_t>(6, warp_size_zu);
  const size_t xfrm_start = static_cast<size_t>(index_orig) * roundUp<size_t>(9, warp_size_zu);
  const size_t atom_start = origin.atom_starts[index_orig];
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                       &origin.boxdims[bdim_start], &origin.umat[xfrm_start],
                       &origin.invu[xfrm_start]);
    copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                      &origin.xcrd[atom_start], &origin.xcrd_ovrf[atom_start],
                      &origin.ycrd[atom_start], &origin.ycrd_ovrf[atom_start],
                      &origin.zcrd[atom_start], &origin.zcrd_ovrf[atom_start],
                      origin.gpos_scale, origin.atom_counts[index_orig]);
    break;
  case TrajectoryKind::VELOCITIES:
    copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                      &origin.xvel[atom_start], &origin.xvel_ovrf[atom_start],
                      &origin.yvel[atom_start], &origin.yvel_ovrf[atom_start],
                      &origin.zvel[atom_start], &origin.zvel_ovrf[atom_start],
                      origin.vel_scale, origin.atom_counts[index_orig]);
    break;
  case TrajectoryKind::FORCES:
    copyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                      &origin.xfrc[atom_start], &origin.xfrc_ovrf[atom_start],
                      &origin.yfrc[atom_start], &origin.yfrc_ovrf[atom_start],
                      &origin.zfrc[atom_start], &origin.zfrc_ovrf[atom_start],
                      origin.frc_scale, origin.atom_counts[index_orig]);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const PhaseSpaceSynthesis &origin,
               const int index_orig, const TrajectoryKind kind,
               const CoordinateCycle orientation) {
  CoordinateFrameWriter cfw = destination->data();
  coordCopy(&cfw, origin.data(orientation), index_orig, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrameWriter *destination, const CondensateReader &origin,
               const int index_orig) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  const int natom = origin.atom_counts[index_orig];
  coordCopyValidateAtomCounts(natom, destination->natom);
  copyBoxInformation(destination->boxdim, destination->umat, destination->invu,
                     origin, index_orig);
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
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CoordinateFrame *destination, const Condensate &origin, const int index_orig) {
  CoordinateFrameWriter cfw = destination->data();
  coordCopy(&cfw, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const TrajectoryKind kind,
               const CoordinateFrameReader &origin) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                       origin.umat, origin.invu);
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
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const CoordinateFrameReader &origin) {
  coordCopy(destination, TrajectoryKind::POSITIONS, origin);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const CoordinateFrame &origin) {
  PhaseSpaceWriter psw = destination->data();
  coordCopy(&psw, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind, const CoordinateFrame &origin) {
  PhaseSpaceWriter psw = destination->data();
  coordCopy(&psw, kind, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrame &origin) {
  PhaseSpaceWriter psw = destination->data(orientation);
  coordCopy(&psw, kind, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin) {
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin.boxdim,
                     origin.umat, origin.invu);

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
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const TrajectoryKind kind,
               const CondensateReader &origin, const int index_orig) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  coordCopyValidateAtomCounts(destination->natom, origin.atom_counts[index_orig]);
  const int natom = origin.atom_counts[index_orig];
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  double *xdest, *ydest, *zdest;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    copyBoxInformation(destination->boxdim, destination->umat, destination->invu, origin,
                       index_orig);
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
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    copyCoordinateXYZ(xdest, ydest, zdest, &origin.xcrd[orig_atom_start],
                      &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], natom);
    break;
  case PrecisionModel::SINGLE:
    copyCoordinateXYZ(xdest, ydest, zdest, &origin.xcrd_sp[orig_atom_start],
                      &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start], natom);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceWriter *destination, const CondensateReader &origin,
               const int index_orig) {
  coordCopy(destination, TrajectoryKind::POSITIONS, origin, index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const Condensate &origin, const int index_orig) {
  PhaseSpaceWriter psw = destination->data(orientation);
  coordCopy(&psw, kind, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind, const Condensate &origin,
               const int index_orig) {
  coordCopy(destination, kind, destination->getCyclePosition(), origin, index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const Condensate &origin, const int index_orig) {
  coordCopy(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
            index_orig);
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpace *destination, const PhaseSpace &origin) {
  PhaseSpaceWriter psw = destination->data();
  coordCopy(&psw, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateFrameReader &origin) {
  destination->import(origin, index_dest);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateFrameReader &origin) {
  destination->import(origin, index_dest, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrameReader &origin) {
  destination->import(origin, index_dest, orientation, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateFrame &origin) {
  destination->import(origin, index_dest);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateFrame &origin) {
  destination->import(origin, index_dest, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateFrame &origin) {
  destination->import(origin, index_dest, orientation, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpaceReader &origin) {
  destination->import(origin, index_dest);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpace &origin) {
  destination->import(origin, index_dest);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int index_dest, const TrajectoryKind kind,
               const CondensateReader &origin, const int index_orig) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = origin.atom_counts[index_orig];
  coordCopyValidateAtomCounts(natom, destination->atom_counts[index_dest]);
  llint *xdest, *ydest, *zdest;
  int *xdest_ovrf, *ydest_ovrf, *zdest_ovrf;
  double dest_scale;
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const size_t padded_box  = roundUp(6, warp_size_int);
      const size_t padded_xfrm = roundUp(9, warp_size_int);
      const size_t dest_bdim_start = static_cast<size_t>(index_dest) * padded_box;
      const size_t dest_xfrm_start = static_cast<size_t>(index_dest) * padded_xfrm;
      copyBoxInformation(&destination->boxdims[dest_bdim_start],
                         &destination->umat[dest_xfrm_start],
                         &destination->invu[dest_xfrm_start], origin, index_orig);
      xdest = &destination->xcrd[dest_atom_start];
      ydest = &destination->ycrd[dest_atom_start];
      zdest = &destination->zcrd[dest_atom_start];
      xdest_ovrf = &destination->xcrd_ovrf[dest_atom_start];
      ydest_ovrf = &destination->ycrd_ovrf[dest_atom_start];
      zdest_ovrf = &destination->zcrd_ovrf[dest_atom_start];
      dest_scale = destination->gpos_scale;
    }
    break;
  case TrajectoryKind::VELOCITIES:
    xdest = &destination->xvel[dest_atom_start];
    ydest = &destination->yvel[dest_atom_start];
    zdest = &destination->zvel[dest_atom_start];
    xdest_ovrf = &destination->xvel_ovrf[dest_atom_start];
    ydest_ovrf = &destination->yvel_ovrf[dest_atom_start];
    zdest_ovrf = &destination->zvel_ovrf[dest_atom_start];
    dest_scale = destination->vel_scale;
    break;
  case TrajectoryKind::FORCES:
    xdest = &destination->xfrc[dest_atom_start];
    ydest = &destination->yfrc[dest_atom_start];
    zdest = &destination->zfrc[dest_atom_start];
    xdest_ovrf = &destination->xfrc_ovrf[dest_atom_start];
    ydest_ovrf = &destination->yfrc_ovrf[dest_atom_start];
    zdest_ovrf = &destination->zfrc_ovrf[dest_atom_start];
    dest_scale = destination->frc_scale;
    break;
  }
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    copyCoordinateXYZ(xdest, xdest_ovrf, ydest, ydest_ovrf, zdest, zdest_ovrf, dest_scale,
                      &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                      &origin.zcrd[orig_atom_start], natom);
    break;
  case PrecisionModel::SINGLE:
    copyCoordinateXYZ(xdest, xdest_ovrf, ydest, ydest_ovrf, zdest, zdest_ovrf, dest_scale,
                      &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                      &origin.zcrd_sp[orig_atom_start], natom);
    break;
  }  
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int index_dest,
               const CondensateReader &origin, const int index_orig) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, origin, index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const Condensate &origin, const int index_orig) {
  PsSynthesisWriter poly_psw = destination->data(orientation);
  coordCopy(&poly_psw, index_dest, kind, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const Condensate &origin, const int index_orig) {
  PsSynthesisWriter poly_psw = destination->data();
  coordCopy(&poly_psw, index_dest, kind, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const Condensate &origin,
               const int index_orig) {
  PsSynthesisWriter poly_psw = destination->data();
  coordCopy(&poly_psw, index_dest, TrajectoryKind::POSITIONS, origin.data(), index_orig);
}
  
//-------------------------------------------------------------------------------------------------
void coordCopy(PsSynthesisWriter *destination, const int index_dest,
               const PsSynthesisReader &origin, const int index_orig) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  coordCopyValidateAtomCounts(destination->atom_counts[index_dest],
                              origin.atom_counts[index_orig]);
  const size_t padded_box  = roundUp(6, warp_size_int);
  const size_t padded_xfrm = roundUp(9, warp_size_int);
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t dest_bdim_start = static_cast<size_t>(index_dest) * padded_box;
  const size_t dest_xfrm_start = static_cast<size_t>(index_dest) * padded_xfrm;
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  const size_t orig_bdim_start = static_cast<size_t>(index_orig) * padded_box;
  const size_t orig_xfrm_start = static_cast<size_t>(index_orig) * padded_xfrm;

  // Copy the present positions, velocities, and forces
  copyBoxInformation(&destination->boxdims[dest_bdim_start], &destination->umat[dest_xfrm_start],
                     &destination->invu[dest_xfrm_start], &origin.boxdims[orig_bdim_start],
                     &origin.umat[orig_xfrm_start], &origin.invu[orig_xfrm_start]);
  copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->xcrd_ovrf[dest_atom_start],
                    &destination->ycrd[dest_atom_start], &destination->ycrd_ovrf[dest_atom_start],
                    &destination->zcrd[dest_atom_start], &destination->zcrd_ovrf[dest_atom_start],
                    destination->gpos_scale, &origin.xcrd[orig_atom_start],
                    &origin.xcrd_ovrf[orig_atom_start], &origin.ycrd[orig_atom_start],
                    &origin.ycrd_ovrf[orig_atom_start], &origin.zcrd[orig_atom_start],
                    &origin.zcrd_ovrf[orig_atom_start], origin.gpos_scale,
                    origin.atom_counts[index_orig]);
  copyCoordinateXYZ(&destination->xvel[dest_atom_start], &destination->xvel_ovrf[dest_atom_start],
                    &destination->yvel[dest_atom_start], &destination->yvel_ovrf[dest_atom_start],
                    &destination->zvel[dest_atom_start], &destination->zvel_ovrf[dest_atom_start],
                    destination->vel_scale, &origin.xvel[orig_atom_start],
                    &origin.xvel_ovrf[orig_atom_start], &origin.yvel[orig_atom_start],
                    &origin.yvel_ovrf[orig_atom_start], &origin.zvel[orig_atom_start],
                    &origin.zvel_ovrf[orig_atom_start], origin.vel_scale,
                    origin.atom_counts[index_orig]);
  copyCoordinateXYZ(&destination->xfrc[dest_atom_start], &destination->xfrc_ovrf[dest_atom_start],
                    &destination->yfrc[dest_atom_start], &destination->yfrc_ovrf[dest_atom_start],
                    &destination->zfrc[dest_atom_start], &destination->zfrc_ovrf[dest_atom_start],
                    destination->frc_scale, &origin.xfrc[orig_atom_start],
                    &origin.xfrc_ovrf[orig_atom_start], &origin.yfrc[orig_atom_start],
                    &origin.yfrc_ovrf[orig_atom_start], &origin.zfrc[orig_atom_start],
                    &origin.zfrc_ovrf[orig_atom_start], origin.frc_scale,
                    origin.atom_counts[index_orig]);

  // Copy the alternate positions, velocities, and forces
  copyBoxInformation(&destination->alt_boxdims[dest_bdim_start],
                     &destination->umat_alt[dest_xfrm_start],
                     &destination->invu_alt[dest_xfrm_start], &origin.alt_boxdims[orig_bdim_start],
                     &origin.umat_alt[orig_xfrm_start], &origin.invu_alt[orig_xfrm_start]);
  copyCoordinateXYZ(&destination->xalt[dest_atom_start], &destination->xalt_ovrf[dest_atom_start],
                    &destination->yalt[dest_atom_start], &destination->yalt_ovrf[dest_atom_start],
                    &destination->zalt[dest_atom_start], &destination->zalt_ovrf[dest_atom_start],
                    destination->gpos_scale, &origin.xalt[orig_atom_start],
                    &origin.xalt_ovrf[orig_atom_start], &origin.yalt[orig_atom_start],
                    &origin.yalt_ovrf[orig_atom_start], &origin.zalt[orig_atom_start],
                    &origin.zalt_ovrf[orig_atom_start], origin.gpos_scale,
                    origin.atom_counts[index_orig]);
  copyCoordinateXYZ(&destination->vxalt[dest_atom_start],
                    &destination->vxalt_ovrf[dest_atom_start],
                    &destination->vyalt[dest_atom_start],
                    &destination->vyalt_ovrf[dest_atom_start],
                    &destination->vzalt[dest_atom_start],
                    &destination->vzalt_ovrf[dest_atom_start], destination->vel_scale,
                    &origin.vxalt[orig_atom_start], &origin.vxalt_ovrf[orig_atom_start],
                    &origin.vyalt[orig_atom_start], &origin.vyalt_ovrf[orig_atom_start],
                    &origin.vzalt[orig_atom_start], &origin.vzalt_ovrf[orig_atom_start],
                    origin.vel_scale, origin.atom_counts[index_orig]);
  copyCoordinateXYZ(&destination->fxalt[dest_atom_start],
                    &destination->fxalt_ovrf[dest_atom_start],
                    &destination->fyalt[dest_atom_start],
                    &destination->fyalt_ovrf[dest_atom_start],
                    &destination->fzalt[dest_atom_start],
                    &destination->fzalt_ovrf[dest_atom_start], destination->frc_scale,
                    &origin.fxalt[orig_atom_start], &origin.fxalt_ovrf[orig_atom_start],
                    &origin.fyalt[orig_atom_start], &origin.fyalt_ovrf[orig_atom_start],
                    &origin.fzalt[orig_atom_start], &origin.fzalt_ovrf[orig_atom_start],
                    origin.frc_scale, origin.atom_counts[index_orig]);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig) {
  PsSynthesisWriter poly_psw = destination->data();
  coordCopy(&poly_psw, index_dest, origin.data(), index_orig);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int index_dest,
               const CoordinateFrameReader &origin) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  const int natom = destination->atom_counts[index_dest];
  const size_t atom_start = destination->atom_counts[index_dest];
  coordCopyValidateAtomCounts(natom, origin.natom);
  copyBoxInformation(destination, index_dest, origin.umat, origin.invu);
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    copyCoordinateXYZ(&destination->xcrd[atom_start], &destination->ycrd[atom_start],
                      &destination->zcrd[atom_start], origin.xcrd, origin.ycrd, origin.zcrd,
                      natom);
    break;
  case PrecisionModel::SINGLE:
    copyCoordinateXYZ(&destination->xcrd_sp[atom_start], &destination->ycrd_sp[atom_start],
                      &destination->zcrd_sp[atom_start], origin.xcrd, origin.ycrd, origin.zcrd,
                      natom);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const CoordinateFrame &origin) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data());
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int index_dest,
               const PhaseSpaceReader &origin, const TrajectoryKind kind) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  const int natom = destination->atom_counts[index_dest];
  const size_t atom_start = destination->atom_starts[index_dest];
  coordCopyValidateAtomCounts(natom, origin.natom);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      copyBoxInformation(destination, index_dest, origin.umat, origin.invu);
      const double* xorig = origin.xcrd;
      const double* yorig = origin.ycrd;
      const double* zorig = origin.zcrd;
      switch (destination->mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&destination->xcrd[atom_start], &destination->ycrd[atom_start],
                          &destination->zcrd[atom_start], xorig, yorig, zorig, natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd_sp[atom_start], &destination->ycrd_sp[atom_start],
                          &destination->zcrd_sp[atom_start], xorig, yorig, zorig, natom);
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
        copyCoordinateXYZ(&destination->xcrd[atom_start], &destination->ycrd[atom_start],
                          &destination->zcrd[atom_start], xorig, yorig, zorig, natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd_sp[atom_start], &destination->ycrd_sp[atom_start],
                          &destination->zcrd_sp[atom_start], xorig, yorig, zorig, natom);
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
        copyCoordinateXYZ(&destination->xcrd[atom_start], &destination->ycrd[atom_start],
                          &destination->zcrd[atom_start], xorig, yorig, zorig, natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd_sp[atom_start], &destination->ycrd_sp[atom_start],
                          &destination->zcrd_sp[atom_start], xorig, yorig, zorig, natom);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int index_dest,
               const PhaseSpaceReader &origin) {
  coordCopy(destination, index_dest, origin, TrajectoryKind::POSITIONS);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin,
               const TrajectoryKind kind, const CoordinateCycle orientation) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(orientation), kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin,
               const TrajectoryKind kind) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(), kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpace &origin) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(), TrajectoryKind::POSITIONS);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, int index_dest, const PsSynthesisReader &origin,
               int index_orig, TrajectoryKind kind) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->atom_counts[index_dest];
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t orig_atom_start = origin.atom_starts[index_orig];
  coordCopyValidateAtomCounts(natom, origin.atom_counts[index_orig]);
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const size_t orig_xfrm_start = static_cast<size_t>(index_orig) *
                                     roundUp<size_t>(9, warp_size_zu);
      copyBoxInformation(destination, index_dest, &origin.umat[orig_xfrm_start],
                         &origin.invu[orig_xfrm_start]);
      const llint* xorig = &origin.xcrd[orig_atom_start];
      const llint* yorig = &origin.ycrd[orig_atom_start];
      const llint* zorig = &origin.zcrd[orig_atom_start];
      const int* xorig_ovrf = &origin.xcrd_ovrf[orig_atom_start];
      const int* yorig_ovrf = &origin.ycrd_ovrf[orig_atom_start];
      const int* zorig_ovrf = &origin.zcrd_ovrf[orig_atom_start];
      switch (destination->mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
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
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
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
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
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
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int index_dest,
               const PsSynthesisReader &origin, const int index_orig) {
  coordCopy(destination, index_dest, origin, index_orig, TrajectoryKind::POSITIONS);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
               const int index_orig, const TrajectoryKind kind, CoordinateCycle orientation) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(orientation), index_orig, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
	       const int index_orig, const TrajectoryKind kind) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(), index_orig, kind);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const PhaseSpaceSynthesis &origin,
               const int index_orig) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(), index_orig, TrajectoryKind::POSITIONS);
}

//-------------------------------------------------------------------------------------------------
void coordCopy(CondensateWriter *destination, const int index_dest, const CondensateReader &origin,
	       const int index_orig) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->atom_counts[index_dest];
  const size_t dest_atom_start = destination->atom_starts[index_dest];
  const size_t orig_atom_start = origin.atom_starts[index_dest];
  coordCopyValidateAtomCounts(natom, origin.atom_counts[index_orig]);
  const size_t orig_xfrm_start = static_cast<size_t>(index_orig) *
                                 roundUp<size_t>(9, warp_size_zu);
  copyBoxInformation(destination, index_dest, &origin.umat[orig_xfrm_start],
                     &origin.invu[orig_xfrm_start]);
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    switch (origin.mode) {
    case PrecisionModel::DOUBLE:
      copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                        &destination->zcrd[dest_atom_start], &origin.xcrd[orig_atom_start],
                        &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], natom);
      break;
    case PrecisionModel::SINGLE:
      copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                        &destination->zcrd[dest_atom_start], &origin.xcrd_sp[orig_atom_start],
                        &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start], natom);
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
                        &destination->zcrd_sp[dest_atom_start], &origin.xcrd_sp[orig_atom_start],
                        &origin.ycrd_sp[orig_atom_start], &origin.zcrd_sp[orig_atom_start], natom);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void coordCopy(Condensate *destination, const int index_dest, const Condensate &origin,
               const int index_orig) {
  CondensateWriter cdnsw = destination->data();
  coordCopy(&cdnsw, index_dest, origin.data(), index_orig);
}

} // namespace trajectory
} // namespace stormm
