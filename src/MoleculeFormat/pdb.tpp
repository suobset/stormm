// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;

//-------------------------------------------------------------------------------------------------
template <typename T>
Pdb::Pdb(const AtomGraph *ag, const CoordinateSeries<T> *cs, const std::vector<int> &frames) :
    Pdb(ag->getAtomCount(), (frames.size() == 0) ? cs->getFrameCount() : frames.size())
{
  loadTopologicalData(ag);
  if (frames.size() == 0) {
    loadCoordinates(cs);
  }
  else {
    const int nframe = frames.size();
    std::vector<int2> frame_mapping(nframe);
    for (int i = 0; i < nframe; i++) {
      frame_mapping[i] = { frames[i], i };
    }
    loadCoordinates(cs, frame_mapping);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Pdb::Pdb(const AtomGraph &ag, const CoordinateSeries<T> &cs, const std::vector<int> &frames) :
    Pdb(ag.getSelfPointer(), cs.getSelfPointer(), frames)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Pdb::loadCoordinates(const T* xcrd_in, const T* ycrd_in, const T* zcrd_in,
                          const int model_index, const double gpos_scale) {
  validateModelIndex(model_index, "loadCoordinates");
  if (gpos_scale < 1.01) {
    for (int i = 0; i < atom_count; i++) {
      const int ip = (model_index * padded_atom_count) + i;
      x_coordinates[ip] = xcrd_in[i];
      y_coordinates[ip] = ycrd_in[i];
      z_coordinates[ip] = zcrd_in[i];
    }
  }
  else {
    const int inv_gpos_scale = 1.0 / gpos_scale;
    for (int i = 0; i < atom_count; i++) {
      const int ip = (model_index * padded_atom_count) + i;
      x_coordinates[ip] = xcrd_in[i] * inv_gpos_scale;
      y_coordinates[ip] = ycrd_in[i] * inv_gpos_scale;
      z_coordinates[ip] = zcrd_in[i] * inv_gpos_scale;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Pdb::loadCoordinates(const CoordinateSeries<T> *cs, const std::vector<int2> &frame_mapping) {
  const CoordinateSeriesReader<T> csr = cs->data();
  const size_t padded_frame_size = roundUp(csr.natom, warp_size_int);
  size_t nframe;
  if (frame_mapping.size() == 0) {
    nframe = std::min(model_count, csr.nframe);
    for (size_t i = 0; i < nframe; i++) {
      const size_t offset = padded_frame_size * i;
      loadCoordinates<T>(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset], i,
                         csr.gpos_scale);
    }
  }
  else {
    nframe = frame_mapping.size();
    for (size_t i = 0; i < nframe; i++) {
      if (frame_mapping[i].x >= csr.nframe) {
        rtErr("Frame index " + std::to_string(frame_mapping[i].x) + " is invalid for a "
              "CoordinateSeries with " + std::to_string(csr.nframe) + " frames.", "Pdb",
              "loadCoordinates");
      }
      const size_t offset = padded_frame_size * frame_mapping[i].x;
      loadCoordinates<T>(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset],
                         frame_mapping[i].y, csr.gpos_scale);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Pdb::loadCoordinates(const CoordinateSeries<T> &cs, const std::vector<int2> &frame_mapping) {
  loadCoordinates(cs.getSelfPointer(), frame_mapping);
}

} // namespace structure
} // namespace stormm
