// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
PureMesh<T>::PureMesh(const AtomGraph *ag_in, const MeshParameters &measurements_in) :
    measurements{measurements_in},
    a_line_x{HybridKind::POINTER, "mesh_avector_x"},
    a_line_y{HybridKind::POINTER, "mesh_avector_y"},
    a_line_z{HybridKind::POINTER, "mesh_avector_z"},
    b_line_x{HybridKind::POINTER, "mesh_bvector_x"},
    b_line_y{HybridKind::POINTER, "mesh_bvector_y"},
    b_line_z{HybridKind::POINTER, "mesh_bvector_z"},
    c_line_x{HybridKind::POINTER, "mesh_cvector_x"},
    c_line_y{HybridKind::POINTER, "mesh_cvector_y"},
    c_line_z{HybridKind::POINTER, "mesh_cvector_z"},
    a_line_x_overflow{HybridKind::POINTER, "mesh_avector_x_ovrf"},
    a_line_y_overflow{HybridKind::POINTER, "mesh_avector_y_ovrf"},
    a_line_z_overflow{HybridKind::POINTER, "mesh_avector_z_ovrf"},
    b_line_x_overflow{HybridKind::POINTER, "mesh_bvector_x_ovrf"},
    b_line_y_overflow{HybridKind::POINTER, "mesh_bvector_y_ovrf"},
    b_line_z_overflow{HybridKind::POINTER, "mesh_bvector_z_ovrf"},
    c_line_x_overflow{HybridKind::POINTER, "mesh_cvector_x_ovrf"},
    c_line_y_overflow{HybridKind::POINTER, "mesh_cvector_y_ovrf"},
    c_line_z_overflow{HybridKind::POINTER, "mesh_cvector_z_ovrf"},
    coefficients{HybridKind::ARRAY, "mesh_tricubic_coef"},
    ag_pointer{const_cast<AtomGraph*>(ag_in)},
    frozen_atoms{HybridKind::ARRAY, "mesh_frozen_atoms"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"}
{}

//-------------------------------------------------------------------------------------------------
PureMesh::PureMesh(const AtomGraph *ag_in, const CoordinateFrame &cf,
                   const MeshParameters &measurements_in, const GpuDetails &gpu) :
    PureMesh(ag_in, measurements_in)
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
void PureMesh::allocate() {
  const int padded_na = roundUp(na, warp_size_int);
  const int padded_nb = roundUp(nb, warp_size_int);
  const int padded_nc = roundUp(nc, warp_size_int);
  llint_data.resize(3 * (padded_na + padded_nb + padded_nc));
  int_data.resize(3 * (padded_na + padded_nb + padded_nc));
  a_line_x.setPointer(&llint_data,             0, na);
  a_line_y.setPointer(&llint_data,     padded_na, na);
  a_line_z.setPointer(&llint_data, 2 * padded_na, na);
  a_line_x_overflow.setPointer(&int_data,             0, na);
  a_line_y_overflow.setPointer(&int_data,     padded_na, na);
  a_line_z_overflow.setPointer(&int_data, 2 * padded_na, na);
  int thus_far = 3 * padded_na;
  b_line_x.setPointer(&llint_data,                   thus_far, nb);
  b_line_y.setPointer(&llint_data,       padded_nb + thus_far, nb);
  b_line_z.setPointer(&llint_data, (2 * padded_nb) + thus_far, nb);
  b_line_x_overflow.setPointer(&int_data,                   thus_far, nb);
  b_line_y_overflow.setPointer(&int_data,       padded_nb + thus_far, nb);
  b_line_z_overflow.setPointer(&int_data, (2 * padded_nb) + thus_far, nb);
  thus_far += 3 * padded_nb;
  c_line_x.setPointer(&llint_data,                   thus_far, nc);
  c_line_y.setPointer(&llint_data,       padded_nc + thus_far, nc);
  c_line_z.setPointer(&llint_data, (2 * padded_nc) + thus_far, nc);
  c_line_x_overflow.setPointer(&int_data,                   thus_far, nc);
  c_line_y_overflow.setPointer(&int_data,       padded_nc + thus_far, nc);
  c_line_z_overflow.setPointer(&int_data, (2 * padded_nc) + thus_far, nc);
  const int nbits = static_cast<int>(sizeof(uint)) * 8;
  coefficients.resize(64LLU * static_cast<size_t>(na * nb * nc));
  frozen_atoms.resize((ag_pointer->getAtomCount() + nbits - 1) / nbits);
}

#if 0
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
#endif

} // namespace structure
} // namespace stormm
