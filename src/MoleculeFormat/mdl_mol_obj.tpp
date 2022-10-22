// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMolObj::impartCoordinates(const T* xcrd, const T* ycrd, const T* zcrd,
                                  const double scale_factor) {
  if (scale_factor != 1.0) {
    for (int i = 0; i < atom_count; i++) {
      coordinates[i].x = xcrd[i] * scale_factor;
      coordinates[i].y = ycrd[i] * scale_factor;
      coordinates[i].z = zcrd[i] * scale_factor;
    }
  }
  else {
    for (int i = 0; i < atom_count; i++) {
      coordinates[i].x = xcrd[i];
      coordinates[i].y = ycrd[i];
      coordinates[i].z = zcrd[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMolObj::impartCoordinates(const CoordinateSeriesReader<T> &csr, const int frame_index,
                                  const HybridTargetLevel tier) {
  if (frame_index < 0 || frame_index >= csr.nframe) {
    rtErr("Invalid frame index request " + std::to_string(frame_index) + " for a series of " +
          std::to_string(csr.nframe) + " frames.", "MdlMolObj", "impartCoordinates");
  }
  const size_t padded_natom = roundUp(csr.natom, warp_size_int);
  const size_t offset = static_cast<size_t>(padded_natom) * static_cast<size_t>(frame_index);
  impartCoordinates(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset], csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMolObj::impartCoordinates(const CoordinateSeriesWriter<T> &csw, const int frame_index,
                                  const HybridTargetLevel tier) {
  if (frame_index < 0 || frame_index >= csw.nframe) {
    rtErr("Invalid frame index request " + std::to_string(frame_index) + " for a series of " +
          std::to_string(csw.nframe) + " frames.", "MdlMolObj", "impartCoordinates");
  }
  const size_t padded_natom = roundUp(csw.natom, warp_size_int);
  const size_t offset = static_cast<size_t>(padded_natom) * static_cast<size_t>(frame_index);
  impartCoordinates(&csw.xcrd[offset], &csw.ycrd[offset], &csw.zcrd[offset], csw.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMolObj::impartCoordinates(const CoordinateSeries<T> &cs, const int frame_index,
                                  const HybridTargetLevel tier) {
  impartCoordinates(cs.data(tier), frame_index);
}
  
} // namespace structure
} // namespace stormm
