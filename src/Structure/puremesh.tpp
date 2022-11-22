// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> Tcoord[] MeshParameters::getMeshOrigin() const {
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (isFloatingPointScalarType<Tcoord>()) {
    Tcoord result[3];
    result[0] = int95ToDouble(origin_x) * inverse_scale_factor;
    result[1] = int95ToDouble(origin_y) * inverse_scale_factor;
    result[2] = int95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else if (isFloatingPointHpcVectorType<Tcoord>()) {
    rtErr("The mesh coordinate origin is available as an HPC vector type through the "
          "getMeshOriginAsTuple() function.", "MeshParameter", "getMeshOrigin");
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.", "MeshParameter", "getMeshOrigin");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> Tcoord MeshParameters::getMeshOriginAsTuple() const {
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (isFloatingPointScalarType<Tcoord>()) {
    rtErr("The mesh coordinate origin is available as a 3-element array of the desired floating "
          "point type through the getMeshOrigin() function.", "MeshParameter",
          "getMeshOriginAsTuple");
  }
  else if (isFloatingPointHpcVectorType<Tcoord>()) {
    Tcoord result;
    result.x = int95ToDouble(origin_x) * inverse_scale_factor;
    result.y = int95ToDouble(origin_y) * inverse_scale_factor;
    result.z = int95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.  The fixed-precision representation is "
          "not available as a tuple.", "MeshParameter", "getMeshOriginAsTuple");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
Tcoord[] MeshParameters::getMeshElementVector(const UnitCellAxis dim) const {
  Tcoord result[3];
  const int icol = static_cast<int>(dim);
  for (int i = 0; i < 3; i++) {
    result[i] = element_invu[i + (3 * icol)];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
Tcoord[] MeshParameters::getMeshElementVector(const CartesianDimension dim) const {
  Tcoord result[3];
  const int icol = static_cast<int>(dim);
  for (int i = 0; i < 3; i++) {
    result[i] = element_invu[i + (3 * icol)];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> PureMeshReader<T> PureMesh<T>::data() const {
  return PureMeshReader<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}

//-------------------------------------------------------------------------------------------------
template <typename T> PureMeshWriter<T> PureMesh<T>::data() {
  return PureMeshWriter<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void PureMesh<T>::allocate() {
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
template <typename T> void PureMesh<T>::initialize(const HybridTargetLevel tier) {
  const size_t all_gsize = 8LLU * roundUp(static_cast<size_t>(measurements.nx) *
                                          static_cast<size_t>(measurements.ny) *
                                          static_cast<size_t>(measurements.nz), warp_size_zu);
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
T PureMesh<T>::at(const size_t i, const size_t j, const size_t k, const HybridTargetLevel tier) {

}

} // namespace structure
} // namespace stormm
