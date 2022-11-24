// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const MeshParameters &measurements_in) :
    measurements{measurements_in}, kind{kind_in},
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
    neighbor_lists{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"}
{
  validateMeshKind();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const MeshParameters &measurements_in,
                                  const GpuDetails &gpu) :
    PureMesh(kind_in, field_in, ag_in, measurements_in)
{
  allocate();

  // Loop over all atoms and apply the potential
  switch (kind) {
  case GridDetail::OCCLUSION:
    break;
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }

  // Loop over all elements and create neighbor lists
  switch (kind) {
  case GridDetail::OCCLUSION:
  case GridDetail::NONBONDED_FIELD:
    break;
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const double buffer, const double spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(cf, buffer, std::vector<double>(3, spacing)), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const double buffer, const std::vector<double> &spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf, getMeasurements(cf, buffer, spacing),
                 gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(cf, mesh_bounds, std::vector<double>(3, spacing)), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(cf, mesh_bounds, spacing), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const double buffer, const double spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(cf, buffer, std::vector<double>(3, spacing)), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const double buffer,
                                  const std::vector<double> &spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf, getMeasurements(cf, buffer, spacing), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const std::vector<double> &mesh_bounds,
                                  const double spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(cf, mesh_bounds, std::vector<double>(3, spacing)), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf, getMeasurements(cf, mesh_bounds, spacing),
                 gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const double buffer, const double spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, -1, ag_in, cf,
                 getMeasurements(cf, buffer, std::vector<double>(3, spacing)), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const int scale_bits,
                                  const double buffer, const std::vector<double> &spacing,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, -1, ag_in, cf,
                 getMeasurements(cf, buffer, spacing), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const std::vector<double> &mesh_bounds,
                                  const double spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(cf, mesh_bounds, std::vector<double>(3, spacing)), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const int scale_bits, const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(cf, mesh_bounds, spacing), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const BackgroundMesh<T> &original) :
    measurements{original.measurements},
    kind{original.kind},
    field{original.field},
    a_line_x{original.a_line_x},
    a_line_y{original.a_line_y},
    a_line_z{original.a_line_z},
    b_line_x{original.b_line_x},
    b_line_y{original.b_line_y},
    b_line_z{original.b_line_z},
    c_line_x{original.c_line_x},
    c_line_y{original.c_line_y},
    c_line_z{original.c_line_z},
    a_line_x_overflow{original.a_line_x_overflow},
    a_line_y_overflow{original.a_line_y_overflow},
    a_line_z_overflow{original.a_line_z_overflow},
    b_line_x_overflow{original.b_line_x_overflow},
    b_line_y_overflow{original.b_line_y_overflow},
    b_line_z_overflow{original.b_line_z_overflow},
    c_line_x_overflow{original.c_line_x_overflow},
    c_line_y_overflow{original.c_line_y_overflow},
    c_line_z_overflow{original.c_line_z_overflow},
    coefficients{original.coefficients},
    ag_pointer{original.ag_pointer},
    frozen_atoms{original.frozen_atoms},
    neighbor_list{original.neighbor_list},
    neighbor_list_bounds{original.neighbor_list_bounds},
    int_data{original.int_data},
    llint_data{original.llint_data}
{
  rebase_pointers();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(BackgroundMesh<T> &&original) :
  measurements{std::move(original.measurements)},
    kind{original.kind},
    field{original.field},
    a_line_x{std::move(original.a_line_x)},
    a_line_y{std::move(original.a_line_y)},
    a_line_z{std::move(original.a_line_z)},
    b_line_x{std::move(original.b_line_x)},
    b_line_y{std::move(original.b_line_y)},
    b_line_z{std::move(original.b_line_z)},
    c_line_x{std::move(original.c_line_x)},
    c_line_y{std::move(original.c_line_y)},
    c_line_z{std::move(original.c_line_z)},
    a_line_x_overflow{std::move(original.a_line_x_overflow)},
    a_line_y_overflow{std::move(original.a_line_y_overflow)},
    a_line_z_overflow{std::move(original.a_line_z_overflow)},
    b_line_x_overflow{std::move(original.b_line_x_overflow)},
    b_line_y_overflow{std::move(original.b_line_y_overflow)},
    b_line_z_overflow{std::move(original.b_line_z_overflow)},
    c_line_x_overflow{std::move(original.c_line_x_overflow)},
    c_line_y_overflow{std::move(original.c_line_y_overflow)},
    c_line_z_overflow{std::move(original.c_line_z_overflow)},
    coefficients{std::move(original.coefficients)},
    ag_pointer{original.ag_pointer},
    frozen_atoms{std::move(original.frozen_atoms)},
    neighbor_list{std::move(original.neighbor_list)},
    neighbor_list_bounds{std::move(original.neighbor_list_bounds)},
    int_data{std::move(original.int_data)},
    llint_data{std::move(original.llint_data)}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::operator=(const BackgroundMesh<T> &other)  {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  measurements = other.measurements;
  kind = other.kind;
  field = other.field;
  a_line_x = other.a_line_x;
  a_line_y = other.a_line_y;
  a_line_z = other.a_line_z;
  b_line_x = other.b_line_x;
  b_line_y = other.b_line_y;
  b_line_z = other.b_line_z;
  c_line_x = other.c_line_x;
  c_line_y = other.c_line_y;
  c_line_z = other.c_line_z;
  a_line_x_overflow = other.a_line_x_overflow;
  a_line_y_overflow = other.a_line_y_overflow;
  a_line_z_overflow = other.a_line_z_overflow;
  b_line_x_overflow = other.b_line_x_overflow;
  b_line_y_overflow = other.b_line_y_overflow;
  b_line_z_overflow = other.b_line_z_overflow;
  c_line_x_overflow = other.c_line_x_overflow;
  c_line_y_overflow = other.c_line_y_overflow;
  c_line_z_overflow = other.c_line_z_overflow;
  coefficients = other.coefficients;
  ag_pointer = other.ag_pointer;
  frozen_atoms = other.frozen_atoms;
  neighbor_list = other.neighbor_list;
  neighbor_list_bounds = other.neighbor_list_bounds;
  int_data = other.int_data;
  llint_data = other.llint_data;

  // Repair pointers
  rebase_pointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::operator=(BackgroundMesh<T> &&other)  {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  measurements = std::move(other.measurements);
  kind = other.kind;
  field = other.field;
  a_line_x = std::move(other.a_line_x);
  a_line_y = std::move(other.a_line_y);
  a_line_z = std::move(other.a_line_z);
  b_line_x = std::move(other.b_line_x);
  b_line_y = std::move(other.b_line_y);
  b_line_z = std::move(other.b_line_z);
  c_line_x = std::move(other.c_line_x);
  c_line_y = std::move(other.c_line_y);
  c_line_z = std::move(other.c_line_z);
  a_line_x_overflow = std::move(other.a_line_x_overflow);
  a_line_y_overflow = std::move(other.a_line_y_overflow);
  a_line_z_overflow = std::move(other.a_line_z_overflow);
  b_line_x_overflow = std::move(other.b_line_x_overflow);
  b_line_y_overflow = std::move(other.b_line_y_overflow);
  b_line_z_overflow = std::move(other.b_line_z_overflow);
  c_line_x_overflow = std::move(other.c_line_x_overflow);
  c_line_y_overflow = std::move(other.c_line_y_overflow);
  c_line_z_overflow = std::move(other.c_line_z_overflow);
  coefficients = std::move(other.coefficients);
  ag_pointer = std::move(other.ag_pointer);
  frozen_atoms = std::move(other.frozen_atoms);
  neighbor_list = std::move(other.neighbor_list);
  neighbor_list_bounds = std::move(other.neighbor_list_bounds);
  int_data = std::move(other.int_data);
  llint_data = std::move(other.llint_data);

  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::allocate() {
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

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::rebase_pointers() {
  a_line_x.swapTarget(&llint_data);
  a_line_y.swapTarget(&llint_data);
  a_line_z.swapTarget(&llint_data);
  b_line_x.swapTarget(&llint_data);
  b_line_y.swapTarget(&llint_data);
  b_line_z.swapTarget(&llint_data);
  c_line_x.swapTarget(&llint_data);
  c_line_y.swapTarget(&llint_data);
  c_line_z.swapTarget(&llint_data);
  a_line_x_overflow.swapTarget(&int_data);
  a_line_y_overflow.swapTarget(&int_data);
  a_line_z_overflow.swapTarget(&int_data);
  b_line_x_overflow.swapTarget(&int_data);
  b_line_y_overflow.swapTarget(&int_data);
  b_line_z_overflow.swapTarget(&int_data);
  c_line_x_overflow.swapTarget(&int_data);
  c_line_y_overflow.swapTarget(&int_data);
  c_line_z_overflow.swapTarget(&int_data);
}

#if 0
//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshReader<T> PureMesh<T>::data() const {
  return BackgroundMeshReader<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}

//-------------------------------------------------------------------------------------------------
template <typename T> BackgroundMeshWriter<T> PureMesh<T>::data() {
  return BackgroundMeshWriter<T>(measurements, values, dx, dy, dz, dxy, dxz, dyz, dxyz);
}
#endif

} // namespace structure
} // namespace stormm
