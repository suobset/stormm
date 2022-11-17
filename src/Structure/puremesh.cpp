// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T> PureMeshReader PureMesh::data() const {
  return PureMeshReader<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}

//-------------------------------------------------------------------------------------------------
template <typename T> PureMeshWriter PureMesh::data() {
  return PureMeshWriter<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void PureMesh::allocate() {
  const size_t gsize = static_cast<size_t>(measurements.nx) *
                       static_cast<size_t>(measurements.ny) * static_cast<size_t>(measurements.nz);
  meshgrids.resize(8LLU * gsize);
  const size_t padded_gsize = roundUp(gsize, warp_size_zu);
  values.setPointer(meshgrids, 0LLU, padded_gsize);
  dx.setPointer(meshgrids, padded_gsize, padded_gsize);
  dy.setPointer(meshgrids, 2LLU * padded_gsize, padded_gsize);
  dz.setPointer(meshgrids, 3LLU * padded_gsize, padded_gsize);
  dxy.setPointer(meshgrids, 4LLU * padded_gsize, padded_gsize);
  dxz.setPointer(meshgrids, 5LLU * padded_gsize, padded_gsize);
  dyz.setPointer(meshgrids, 6LLU * padded_gsize, padded_gsize);
  dxyz.setPointer(meshgrids, 7LLU * padded_gsize, padded_gsize);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void PureMesh::initialize(const HybridTargetLevel tier) {
  const size_t all_gsize = 8LLU * roundUp( static_cast<size_t>(measurements.nx) * static_cast<>(measurements.ny) + static_cast<size_t>(measurements.nz), warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    memset(meshgrids.data(tier), all_gsize * sizeof(double), 0);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    cudaMemset(meshgrids.data(tier), all_gsize * sizeof(double), 0);
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T PureMesh::at(const size_t i, const size_t j, const size_t k, const HybridTargetLevel tier) {

}

} // namespace structure
} // namespace stormm
