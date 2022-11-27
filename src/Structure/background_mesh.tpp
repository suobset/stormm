// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Txfrm, typename Tdata>
BackgroundMeshReader<Txfrm, Tdata>::
BackgroundMeshReader(const MeshParamAbstract<Txfrm> &dims_in, GridDetail kind_in,
                     NonbondedPotential field_in, const llint* avec_x_in, const llint* avec_y_in,
                     const llint* avec_z_in, const llint* bvec_x_in, const llint* bvec_y_in,
                     const llint* bvec_z_in, const llint* cvec_x_in, const llint* cvec_y_in,
                     const llint* cvec_z_in, const llint* avec_abs_x_in,
                     const llint* avec_abs_y_in, const llint* avec_abs_z_in,
                     const int* avec_x_ovrf_in, const int* avec_y_ovrf_in,
                     const int* avec_z_ovrf_in, const int* bvec_x_ovrf_in,
                     const int* bvec_y_ovrf_in, const int* bvec_z_ovrf_in,
                     const int* cvec_x_ovrf_in, const int* cvec_y_ovrf_in,
                     const int* cvec_z_ovrf_in, const int* avec_abs_x_ovrf_in,
                     const int* avec_abs_y_ovrf_in, const int* avec_abs_z_ovrf_in,
                     const Tdata* coeffs_in, const int* ngbr_in, const size_t* ngbr_bounds_in) :
    dims{dims_in}, kind{kind_in}, field{field_in}, avec_x{avec_x_in}, avec_y{avec_y_in},
    avec_z{avec_z_in}, bvec_x{bvec_x_in}, bvec_y{bvec_y_in}, bvec_z{bvec_z_in}, cvec_x{cvec_x_in},
    cvec_y{cvec_y_in}, cvec_z{cvec_z_in}, avec_abs_x{avec_abs_x_in}, avec_abs_y{avec_abs_y_in},
    avec_abs_z{avec_abs_z_in}, avec_x_ovrf{avec_x_ovrf_in}, avec_y_ovrf{avec_y_ovrf_in},
    avec_z_ovrf{avec_z_ovrf_in}, bvec_x_ovrf{bvec_x_ovrf_in}, bvec_y_ovrf{bvec_y_ovrf_in},
    bvec_z_ovrf{bvec_z_ovrf_in}, cvec_x_ovrf{cvec_x_ovrf_in}, cvec_y_ovrf{cvec_y_ovrf_in},
    cvec_z_ovrf{cvec_z_ovrf_in}, avec_abs_x_ovrf{avec_abs_x_ovrf_in},
    avec_abs_y_ovrf{avec_abs_y_ovrf_in}, avec_abs_z_ovrf{avec_abs_z_ovrf_in}, coeffs{coeffs_in},
    ngbr{ngbr_in}, ngbr_bounds{ngbr_bounds_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename Txfrm, typename Tdata>
BackgroundMeshWriter<Txfrm, Tdata>::
BackgroundMeshWriter(const MeshParamAbstract<Txfrm> &dims_in, GridDetail kind_in,
                     NonbondedPotential field_in, const llint* avec_x_in, const llint* avec_y_in,
                     const llint* avec_z_in, const llint* bvec_x_in, const llint* bvec_y_in,
                     const llint* bvec_z_in, const llint* cvec_x_in, const llint* cvec_y_in,
                     const llint* cvec_z_in, const llint* avec_abs_x_in,
                     const llint* avec_abs_y_in, const llint* avec_abs_z_in,
                     const int* avec_x_ovrf_in, const int* avec_y_ovrf_in,
                     const int* avec_z_ovrf_in, const int* bvec_x_ovrf_in,
                     const int* bvec_y_ovrf_in, const int* bvec_z_ovrf_in,
                     const int* cvec_x_ovrf_in, const int* cvec_y_ovrf_in,
                     const int* cvec_z_ovrf_in, const int* avec_abs_x_ovrf_in,
                     const int* avec_abs_y_ovrf_in, const int* avec_abs_z_ovrf_in,
                     Tdata* coeffs_in, int* ngbr_in, size_t* ngbr_bounds_in) :
    dims{dims_in}, kind{kind_in}, field{field_in}, avec_x{avec_x_in}, avec_y{avec_y_in},
    avec_z{avec_z_in}, bvec_x{bvec_x_in}, bvec_y{bvec_y_in}, bvec_z{bvec_z_in}, cvec_x{cvec_x_in},
    cvec_y{cvec_y_in}, cvec_z{cvec_z_in}, avec_abs_x{avec_abs_x_in}, avec_abs_y{avec_abs_y_in},
    avec_abs_z{avec_abs_z_in}, avec_x_ovrf{avec_x_ovrf_in}, avec_y_ovrf{avec_y_ovrf_in},
    avec_z_ovrf{avec_z_ovrf_in}, bvec_x_ovrf{bvec_x_ovrf_in}, bvec_y_ovrf{bvec_y_ovrf_in},
    bvec_z_ovrf{bvec_z_ovrf_in}, cvec_x_ovrf{cvec_x_ovrf_in}, cvec_y_ovrf{cvec_y_ovrf_in},
    cvec_z_ovrf{cvec_z_ovrf_in}, avec_abs_x_ovrf{avec_abs_x_ovrf_in},
    avec_abs_y_ovrf{avec_abs_y_ovrf_in}, avec_abs_z_ovrf{avec_abs_z_ovrf_in}, coeffs{coeffs_in},
    ngbr{ngbr_in}, ngbr_bounds{ngbr_bounds_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const int vdw_type_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const MeshParameters &measurements_in, const GpuDetails &gpu) :
    measurements{measurements_in}, kind{kind_in}, field{field_in}, probe_radius{probe_radius_in},
    a_line_x{HybridKind::POINTER, "mesh_avector_x"},
    a_line_y{HybridKind::POINTER, "mesh_avector_y"},
    a_line_z{HybridKind::POINTER, "mesh_avector_z"},
    b_line_x{HybridKind::POINTER, "mesh_bvector_x"},
    b_line_y{HybridKind::POINTER, "mesh_bvector_y"},
    b_line_z{HybridKind::POINTER, "mesh_bvector_z"},
    c_line_x{HybridKind::POINTER, "mesh_cvector_x"},
    c_line_y{HybridKind::POINTER, "mesh_cvector_y"},
    c_line_z{HybridKind::POINTER, "mesh_cvector_z"},
    a_abs_line_x{HybridKind::POINTER, "mesh_avector_abs_x"},
    a_abs_line_y{HybridKind::POINTER, "mesh_avector_abs_y"},
    a_abs_line_z{HybridKind::POINTER, "mesh_avector_abs_z"},
    a_line_x_overflow{HybridKind::POINTER, "mesh_avec_x_ovrf"},
    a_line_y_overflow{HybridKind::POINTER, "mesh_avec_y_ovrf"},
    a_line_z_overflow{HybridKind::POINTER, "mesh_avec_z_ovrf"},
    b_line_x_overflow{HybridKind::POINTER, "mesh_bvec_x_ovrf"},
    b_line_y_overflow{HybridKind::POINTER, "mesh_bvec_y_ovrf"},
    b_line_z_overflow{HybridKind::POINTER, "mesh_bvec_z_ovrf"},
    c_line_x_overflow{HybridKind::POINTER, "mesh_cvec_x_ovrf"},
    c_line_y_overflow{HybridKind::POINTER, "mesh_cvec_y_ovrf"},
    c_line_z_overflow{HybridKind::POINTER, "mesh_cvec_z_ovrf"},
    a_abs_line_x_overflow{HybridKind::POINTER, "mesh_avec_abs_x_ovrf"},
    a_abs_line_y_overflow{HybridKind::POINTER, "mesh_avec_abs_y_ovrf"},
    a_abs_line_z_overflow{HybridKind::POINTER, "mesh_avec_abs_z_ovrf"},
    coefficients{HybridKind::ARRAY, "mesh_tricubic_coef"},
    ag_pointer{const_cast<AtomGraph*>(ag_in)},
    cf_pointer{const_cast<CoordinateFrame*>(cf_in)},
    frozen_atoms{HybridKind::ARRAY, "mesh_frozen_atoms"},
    neighbor_list{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"}
{
  // Validate the mesh dimensions, then allocate memory
  validateMeshKind();
  validateScalingBits();
  allocate();

  // Loop over all atoms and apply the potential
  switch (kind) {
  case GridDetail::OCCLUSION:
    colorExclusionMesh(gpu);
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
                                  const CoordinateFrame *cf_in, const double buffer,
                                  const double spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame *cf_in, const double buffer,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, spacing, scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const double spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, spacing, scale_bits_in), gpu)
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
    a_abs_line_x{original.a_abs_line_x},
    a_abs_line_y{original.a_abs_line_y},
    a_abs_line_z{original.a_abs_line_z},
    a_line_x_overflow{original.a_line_x_overflow},
    a_line_y_overflow{original.a_line_y_overflow},
    a_line_z_overflow{original.a_line_z_overflow},
    b_line_x_overflow{original.b_line_x_overflow},
    b_line_y_overflow{original.b_line_y_overflow},
    b_line_z_overflow{original.b_line_z_overflow},
    c_line_x_overflow{original.c_line_x_overflow},
    c_line_y_overflow{original.c_line_y_overflow},
    c_line_z_overflow{original.c_line_z_overflow},
    a_abs_line_x_overflow{original.a_abs_line_x_overflow},
    a_abs_line_y_overflow{original.a_abs_line_y_overflow},
    a_abs_line_z_overflow{original.a_abs_line_z_overflow},
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
    a_abs_line_x{std::move(original.a_abs_line_x)},
    a_abs_line_y{std::move(original.a_abs_line_y)},
    a_abs_line_z{std::move(original.a_abs_line_z)},
    a_line_x_overflow{std::move(original.a_line_x_overflow)},
    a_line_y_overflow{std::move(original.a_line_y_overflow)},
    a_line_z_overflow{std::move(original.a_line_z_overflow)},
    b_line_x_overflow{std::move(original.b_line_x_overflow)},
    b_line_y_overflow{std::move(original.b_line_y_overflow)},
    b_line_z_overflow{std::move(original.b_line_z_overflow)},
    c_line_x_overflow{std::move(original.c_line_x_overflow)},
    c_line_y_overflow{std::move(original.c_line_y_overflow)},
    c_line_z_overflow{std::move(original.c_line_z_overflow)},
    a_abs_line_x_overflow{std::move(original.a_abs_line_x_overflow)},
    a_abs_line_y_overflow{std::move(original.a_abs_line_y_overflow)},
    a_abs_line_z_overflow{std::move(original.a_abs_line_z_overflow)},
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
BackgroundMesh<T>& BackgroundMesh<T>::operator=(const BackgroundMesh<T> &other) {

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
  a_abs_line_x = other.a_abs_line_x;
  a_abs_line_y = other.a_abs_line_y;
  a_abs_line_z = other.a_abs_line_z;
  a_line_x_overflow = other.a_line_x_overflow;
  a_line_y_overflow = other.a_line_y_overflow;
  a_line_z_overflow = other.a_line_z_overflow;
  b_line_x_overflow = other.b_line_x_overflow;
  b_line_y_overflow = other.b_line_y_overflow;
  b_line_z_overflow = other.b_line_z_overflow;
  c_line_x_overflow = other.c_line_x_overflow;
  c_line_y_overflow = other.c_line_y_overflow;
  c_line_z_overflow = other.c_line_z_overflow;
  a_abs_line_x_overflow = other.a_abs_line_x_overflow;
  a_abs_line_y_overflow = other.a_abs_line_y_overflow;
  a_abs_line_z_overflow = other.a_abs_line_z_overflow;
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
BackgroundMesh<T>& BackgroundMesh<T>::operator=(BackgroundMesh<T> &&other)  {

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
  a_abs_line_x = std::move(other.a_abs_line_x);
  a_abs_line_y = std::move(other.a_abs_line_y);
  a_abs_line_z = std::move(other.a_abs_line_z);
  a_line_x_overflow = std::move(other.a_line_x_overflow);
  a_line_y_overflow = std::move(other.a_line_y_overflow);
  a_line_z_overflow = std::move(other.a_line_z_overflow);
  b_line_x_overflow = std::move(other.b_line_x_overflow);
  b_line_y_overflow = std::move(other.b_line_y_overflow);
  b_line_z_overflow = std::move(other.b_line_z_overflow);
  c_line_x_overflow = std::move(other.c_line_x_overflow);
  c_line_y_overflow = std::move(other.c_line_y_overflow);
  c_line_z_overflow = std::move(other.c_line_z_overflow);
  a_abs_line_x_overflow = std::move(other.a_abs_line_x_overflow);
  a_abs_line_y_overflow = std::move(other.a_abs_line_y_overflow);
  a_abs_line_z_overflow = std::move(other.a_abs_line_z_overflow);
  coefficients = std::move(other.coefficients);
  ag_pointer = other.ag_pointer;
  frozen_atoms = std::move(other.frozen_atoms);
  neighbor_list = std::move(other.neighbor_list);
  neighbor_list_bounds = std::move(other.neighbor_list_bounds);
  int_data = std::move(other.int_data);
  llint_data = std::move(other.llint_data);

  // As usual, no pointer repair is needed for the move assignment operator (or move constructor).
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T> const AtomGraph* BackgroundMesh<T>::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
template <typename T> const CoordinateFrame* BackgroundMesh<T>::getCoordinatePointer() const {
  return cf_pointer;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshReader<double, T> BackgroundMesh<T>::dpData(const HybridTargetLevel tier) const {
  return BackgroundMeshReader<double, T>(measurements.dpData(), kind, field, a_line_x.data(tier),
                                         a_line_y.data(tier), a_line_z.data(tier),
                                         b_line_x.data(tier), b_line_y.data(tier),
                                         b_line_z.data(tier), c_line_x.data(tier),
                                         c_line_y.data(tier), c_line_z.data(tier),
                                         a_abs_line_x.data(tier), a_abs_line_y.data(tier),
                                         a_abs_line_z.data(tier), a_line_x_overflow.data(tier),
                                         a_line_y_overflow.data(tier),
                                         a_line_z_overflow.data(tier),
                                         b_line_x_overflow.data(tier),
                                         b_line_y_overflow.data(tier),
                                         b_line_z_overflow.data(tier),
                                         c_line_x_overflow.data(tier),
                                         c_line_y_overflow.data(tier),
                                         c_line_z_overflow.data(tier),
                                         a_abs_line_x_overflow.data(tier),
                                         a_abs_line_y_overflow.data(tier),
                                         a_abs_line_z_overflow.data(tier), coefficients.data(tier),
                                         neighbor_list.data(tier),
                                         neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<double, T> BackgroundMesh<T>::dpData(const HybridTargetLevel tier) {
  return BackgroundMeshReader<double, T>(measurements.dpData(), kind, field, a_line_x.data(tier),
                                         a_line_y.data(tier), a_line_z.data(tier),
                                         b_line_x.data(tier), b_line_y.data(tier),
                                         b_line_z.data(tier), c_line_x.data(tier),
                                         c_line_y.data(tier), c_line_z.data(tier),
                                         a_abs_line_x.data(tier), a_abs_line_y.data(tier),
                                         a_abs_line_z.data(tier), a_line_x_overflow.data(tier),
                                         a_line_y_overflow.data(tier),
                                         a_line_z_overflow.data(tier),
                                         b_line_x_overflow.data(tier),
                                         b_line_y_overflow.data(tier),
                                         b_line_z_overflow.data(tier),
                                         c_line_x_overflow.data(tier),
                                         c_line_y_overflow.data(tier),
                                         c_line_z_overflow.data(tier),
                                         a_abs_line_x_overflow.data(tier),
                                         a_abs_line_y_overflow.data(tier),
                                         a_abs_line_z_overflow.data(tier), coefficients.data(tier),
                                         neighbor_list.data(tier),
                                         neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshReader<float, T> BackgroundMesh<T>::spData(const HybridTargetLevel tier) const {
  return BackgroundMeshReader<float, T>(measurements.spData(), kind, field, a_line_x.data(tier),
                                        a_line_y.data(tier), a_line_z.data(tier),
                                        b_line_x.data(tier), b_line_y.data(tier),
                                        b_line_z.data(tier), c_line_x.data(tier),
                                        c_line_y.data(tier), c_line_z.data(tier),
                                        a_abs_line_x.data(tier), a_abs_line_y.data(tier),
                                        a_abs_line_z.data(tier), a_line_x_overflow.data(tier),
                                        a_line_y_overflow.data(tier), a_line_z_overflow.data(tier),
                                        b_line_x_overflow.data(tier), b_line_y_overflow.data(tier),
                                        b_line_z_overflow.data(tier), c_line_x_overflow.data(tier),
                                        c_line_y_overflow.data(tier), c_line_z_overflow.data(tier),
                                        a_abs_line_x_overflow.data(tier),
                                        a_abs_line_y_overflow.data(tier),
                                        a_abs_line_z_overflow.data(tier), coefficients.data(tier),
                                        neighbor_list.data(tier), neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<float, T> BackgroundMesh<T>::spData(const HybridTargetLevel tier) {
  return BackgroundMeshReader<float, T>(measurements.spData(), kind, field, a_line_x.data(tier),
                                        a_line_y.data(tier), a_line_z.data(tier),
                                        b_line_x.data(tier), b_line_y.data(tier),
                                        b_line_z.data(tier), c_line_x.data(tier),
                                        c_line_y.data(tier), c_line_z.data(tier),
                                        a_abs_line_x.data(tier), a_abs_line_y.data(tier),
                                        a_abs_line_z.data(tier), a_line_x_overflow.data(tier),
                                        a_line_y_overflow.data(tier), a_line_z_overflow.data(tier),
                                        b_line_x_overflow.data(tier), b_line_y_overflow.data(tier),
                                        b_line_z_overflow.data(tier), c_line_x_overflow.data(tier),
                                        c_line_y_overflow.data(tier), c_line_z_overflow.data(tier),
                                        a_abs_line_x_overflow.data(tier),
                                        a_abs_line_y_overflow.data(tier),
                                        a_abs_line_z_overflow.data(tier), coefficients.data(tier),
                                        neighbor_list.data(tier), neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::setTopologyPointer(const AtomGraph *ag_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::setCoordinatePointer(const CoordinateFrame *cf_in) {
  cf_pointer = const_cast<CoordinateFrame*>(cf_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                          const double padding, const double spacing,
                                          const int scale_bits_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
  cf_pointer = const_cast<CoordinateFrame*>(cf_in);
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(ag_pointer, cf_pointer, padding, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                          const double padding,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
  cf_pointer = const_cast<CoordinateFrame*>(cf_in);
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(ag_pointer, cf_pointer, padding, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const double spacing,
                                          const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(ag_pointer, cf_pointer, padding, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(ag_pointer, cf_pointer, padding, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const double spacing, const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(mesh_bounds, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(mesh_bounds, spacing, actual_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::colorExclusionMesh(const GpuDetails &gpu) {

  // Use the HPC kernel to color the mesh if a GPU is available
#ifdef STORMM_USE_HPC
  if (gpu != null_gpu) {
    coefficients.download();
    return;
  }
#endif

  // Color the mesh on the CPU
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const MeshParamAbstract<double> mps = measurements.dpData();
  const CoordinateFrameReader cfr = cf_pointer->data();
  
  // Initialize the mesh in CPU memory.
  const int n_elem = mps.na * mps.nb * mps.nc;
  T* coeff_ptr = coefficients.data();
  for (int pos = 0; pos < n_elem; pos++) {
    coeff_ptr[pos] = 0LLU;
  }

  // Replicate the mesh's grid vectors.  Subtract the origin coordinates from the "b" and "c"
  // vectors so that a[i] + b[j] + c[k] = the Cartesian coordinates of any grid point (i,j,k).

  std::vector<double> avx(mps.na + 1), avy(mps.na + 1), avz(mps.na + 1);
  std::vector<double> bvx(mps.nb + 1), bvy(mps.nb + 1), bvz(mps.nb + 1);
  std::vector<double> cvx(mps.nc + 1), cvy(mps.nc + 1), cvz(mps.nc + 1);
  int95ToDouble(&avx, &avy, &avz, a_abs_line_x.readHost(), a_abs_line_x_overflow.readHost(),
                a_abs_line_y.readHost(), a_abs_line_y_overflow.readHost(),
                a_abs_line_z.readHost(), a_abs_line_z_overflow.readHost(), mps.na, mps.inv_scale);
  int95ToDouble(&bvx, &bvy, &bvz, b_line_x.readHost(), b_line_x_overflow.readHost(),
                b_line_y.readHost(), b_line_y_overflow.readHost(), b_line_z.readHost(),
                b_line_z_overflow.readHost(), mps.nb, mps.inv_scale);
  int95ToDouble(&cvx, &cvy, &cvz, c_line_x.readHost(), c_line_x_overflow.readHost(),
                c_line_y.readHost(), c_line_y_overflow.readHost(), c_line_z.readHost(),
                c_line_z_overflow.readHost(), mps.nc, mps.inv_scale);
  const double dorig_x = int95ToDouble(mps.orig_x);
  const double dorig_y = int95ToDouble(mps.orig_y);
  const double dorig_z = int95ToDouble(mps.orig_z);

  // Compute the increment within one element as the 16 x 16 x 16 mesh gets traced
  const double de_ax = 0.0625 * (avx[1] - avx[0]);
  const double de_ay = 0.0625 * (avy[1] - avy[0]);
  const double de_az = 0.0625 * (avz[1] - avz[0]);
  const double de_bx = 0.0625 * (bvx[1] - bvx[0]);
  const double de_by = 0.0625 * (bvy[1] - bvy[0]);
  const double de_bz = 0.0625 * (bvz[1] - bvz[0]);
  const double de_cx = 0.0625 * (cvx[1] - cvx[0]);
  const double de_cy = 0.0625 * (cvy[1] - cvy[0]);
  const double de_cz = 0.0625 * (cvz[1] - cvz[0]);

  // Prepare a buffer to hold subgrid results
  std::vector<ullint> cube_buffer(64, 0LLU);

  // Loop over all atoms in the system
  const int lj_idx_offset = nbk.n_lj_types + 1;
  for (int pos = 0; pos < nbk.natom; pos++) {
    const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
    const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                         1.0 / 6.0);
    const double color_radius = atom_radius + probe_radius;
    const double color_radius_sq = color_radius * color_radius;
    const double atom_x = cfr.xcrd[pos];
    const double atom_y = cfr.ycrd[pos];
    const double atom_z = cfr.zcrd[pos];
    const int ixmin = std::max(static_cast<int>(floor((atom_x - color_radius - dorig_x) /
                                                      mps.widths[0])), 0);
    const int iymin = std::max(static_cast<int>(floor((atom_y - color_radius - dorig_y) /
                                                      mps.widths[1])), 0);
    const int izmin = std::max(static_cast<int>(floor((atom_z - color_radius - dorig_z) /
                                                      mps.widths[2])), 0);
    const int ixmax = std::min(static_cast<int>(ceil((atom_x + color_radius - dorig_x) /
                                                     mps.widths[0])), mps.na - 1);
    const int iymax = std::min(static_cast<int>(ceil((atom_y + color_radius - dorig_y) /
                                                     mps.widths[1])), mps.nb - 1);
    const int izmax = std::min(static_cast<int>(ceil((atom_z + color_radius - dorig_z) /
                                                     mps.widths[2])), mps.nc - 1);
    for (int i = ixmin; i < ixmax; i++) {
      for (int j = iymin; j < iymax; j++) {
        for (int k = izmin; k < izmax; k++) {

          // Color the buffer for this atom and this element
          const double base_x = avx[ixmin] + bvx[iymin] + cvx[izmin];
          const double base_y = avy[ixmin] + bvy[iymin] + cvy[izmin];
          const double base_z = avz[ixmin] + bvz[iymin] + cvz[izmin];
          for (size_t m = 0LLU; m < 64LLU; m++) {
            cube_buffer[m] = 0LLU;
          }
          int grid_i = 0;
          for (double di = 0.5; di < 16.0; di += 1.0) {
            int grid_j = 0;
            for (double dj = 0.5; dj < 16.0; dj += 1.0) {
              int grid_k = 0;
              for (double dk = 0.5; dk < 16.0; dk += 1.0) {
                const double xpt = base_x + (di * de_ax) + (dj * de_bx) + (dk * de_cx);
                const double ypt = base_y + (di * de_ay) + (dj * de_by) + (dk * de_cy);
                const double zpt = base_z + (di * de_az) + (dj * de_bz) + (dk * de_cz);
                const double dx = xpt - atom_x;
                const double dy = ypt - atom_y;
                const double dz = zpt - atom_z;
                if ((dx * dx) + (dy * dy) + (dz * dz) < color_radius_sq) {
                  const int cubelet_i = grid_i / 4;
                  const int cubelet_j = grid_j / 4;
                  const int cubelet_k = grid_k / 4;
                  const int cubelet_idx = (((cubelet_i * 4) + cubelet_j) * 4) + cubelet_k;
                  const int bit_i = grid_i - (4 * cubelet_i);
                  const int bit_j = grid_j - (4 * cubelet_j);
                  const int bit_k = grid_k - (4 * cubelet_k);
                  const int bit_idx = (((bit_i * 4) + bit_j) * 4) + bit_k;
                  cube_buffer[cubelet_idx] |= (0x1LLU << bit_idx);
                }
                grid_k++;
              }
              grid_j++;
            }
            grid_i++;
          }

          // Accumulate the atom's mapping onto the grid
          const size_t coef_base_idx = 64LLU *
                                       static_cast<size_t>((((i * mps.nb) + j) * mps.nc) + k);
          for (size_t m = 0LLU; m < 64LLU; m++) {

            // This is necessary due to the templated nature of the BackgroundMesh object.  In
            // practice, the only way to get here is to have the coefficients array already be of
            // type ullint.
            ullint tcoef = coeff_ptr[coef_base_idx + m];
            tcoef |= cube_buffer[m];
            coeff_ptr[coef_base_idx + m] = tcoef;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::mapElectrostatics(const GpuDetails &gpu) {

}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                                                  const double padding, const double spacing,
                                                  const int scale_bits_in) const {
  return getMeasurements(ag, cf, padding, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                                                  const std::vector<double> &mesh_bounds,
                                                  const double spacing,
                                                  const int scale_bits_in) const {
  return getMeasurements(ag, cf, mesh_bounds, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                                                  const double padding,
                                                  const std::vector<double> &spacing,
                                                  const int scale_bits_in) const {
  if (cf == nullptr || ag == nullptr) {
    rtErr("When creating a mesh based on a zone surrounding a system of interest, the topology "
          "and coordinates must be defined.  Supply the mesh boundaries explicitly to avoid "
          "relying on pre-defined coordinates, or submit the system to the mesh prior to calling "
          "this function.", "BackgroundMesh", "getMeasurements");
  }
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const CoordinateFrameReader cfr = cf->data();
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bool points_unset = true;
  const int lj_idx_offset = nbk.n_lj_types + 1;
  for (int pos = 0; pos < nbk.natom; pos++) {
    if (ag->getAtomMobility(pos)) {
      const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
      const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                           1.0 / 6.0);
      if (points_unset) {
        xmin = cfr.xcrd[pos] - atom_radius;
        xmax = cfr.xcrd[pos] + atom_radius;
        ymin = cfr.ycrd[pos] - atom_radius;
        ymax = cfr.ycrd[pos] + atom_radius;
        zmin = cfr.zcrd[pos] - atom_radius;
        zmax = cfr.zcrd[pos] + atom_radius;
        points_unset = false;
      }
      else {
        xmin = std::min(xmin, cfr.xcrd[pos] - atom_radius);
        xmax = std::max(xmax, cfr.xcrd[pos] + atom_radius);
        ymin = std::min(ymin, cfr.ycrd[pos] - atom_radius);
        ymax = std::max(ymax, cfr.ycrd[pos] + atom_radius);
        zmin = std::min(zmin, cfr.zcrd[pos] - atom_radius);
        zmax = std::max(zmax, cfr.zcrd[pos] + atom_radius);
      }
    }
  }
  xmin -= padding;
  xmax += padding;
  ymin -= padding;
  ymax += padding;
  zmin -= padding;
  zmax += padding;
  const std::vector<double> limits = { xmin, ymin, zmin, xmax, ymax, zmax };
  return getMeasurements(limits, spacing, scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const std::vector<double> &mesh_bounds,
                                                  const std::vector<double> &spacing,
                                                  const int scale_bits_in) const {
  if (mesh_bounds.size() != 6LLU) {
    rtErr("An array of six elements, the minimum X, Y, and Z Cartesian coordinates followed by "
          "the maximum coordinates, is required.  " + std::to_string(mesh_bounds.size()) +
          " elements were provided.", "BackgroundMesh", "getMeasurements");
  }
  if (spacing.size() != 3LLU && spacing.size() != 9LLU) {
    rtErr("An array of three elements, the length, width, and height (Cartesian X, Y, and Z "
          "dimensions) of a rectilinear mesh element, or nine elements defining the bounding "
          "vectors of a triclininc element, is required.  " +
          std::to_string(spacing.size()) + " elements were provided.", "BackgroundMesh",
          "getMeasurements");
  }
  const std::vector<double> mesh_limits = {
    std::min(mesh_bounds[0], mesh_bounds[3]), std::max(mesh_bounds[0], mesh_bounds[3]),
    std::min(mesh_bounds[1], mesh_bounds[4]), std::max(mesh_bounds[1], mesh_bounds[4]),
    std::min(mesh_bounds[2], mesh_bounds[5]), std::max(mesh_bounds[2], mesh_bounds[5]) };  
  const int pna = ceil((mesh_limits[3] - mesh_limits[0]) / spacing[0]);
  const int pnb = ceil((mesh_limits[4] - mesh_limits[1]) / spacing[1]);
  const int pnc = ceil((mesh_limits[5] - mesh_limits[2]) / spacing[2]);
  const double dnx = static_cast<double>(pna) * spacing[0];
  const double dny = static_cast<double>(pnb) * spacing[1];
  const double dnz = static_cast<double>(pnc) * spacing[2];
  const double overshoot_x = 0.5 * (dnx - (mesh_limits[3] - mesh_limits[0]));
  const double overshoot_y = 0.5 * (dny - (mesh_limits[4] - mesh_limits[1]));
  const double overshoot_z = 0.5 * (dnz - (mesh_limits[5] - mesh_limits[2]));
  MeshParameters result(pna, pnb, pnc, mesh_limits[0] - overshoot_x, mesh_limits[1] - overshoot_y,
                        mesh_limits[2] - overshoot_z, spacing, scale_bits_in);  
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::allocate() {
  MeshParamAbstract<double> mps = measurements.dpData();
  const int padded_na = roundUp(mps.na, warp_size_int);
  const int padded_nb = roundUp(mps.nb, warp_size_int);
  const int padded_nc = roundUp(mps.nc, warp_size_int);
  llint_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  int_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  a_line_x.setPointer(&llint_data,             0, mps.na);
  a_line_y.setPointer(&llint_data,     padded_na, mps.na);
  a_line_z.setPointer(&llint_data, 2 * padded_na, mps.na);
  a_line_x_overflow.setPointer(&int_data,             0, mps.na);
  a_line_y_overflow.setPointer(&int_data,     padded_na, mps.na);
  a_line_z_overflow.setPointer(&int_data, 2 * padded_na, mps.na);
  int thus_far = 3 * padded_na;
  b_line_x.setPointer(&llint_data,                   thus_far, mps.nb);
  b_line_y.setPointer(&llint_data,       padded_nb + thus_far, mps.nb);
  b_line_z.setPointer(&llint_data, (2 * padded_nb) + thus_far, mps.nb);
  b_line_x_overflow.setPointer(&int_data,                   thus_far, mps.nb);
  b_line_y_overflow.setPointer(&int_data,       padded_nb + thus_far, mps.nb);
  b_line_z_overflow.setPointer(&int_data, (2 * padded_nb) + thus_far, mps.nb);
  thus_far += 3 * padded_nb;
  c_line_x.setPointer(&llint_data,                   thus_far, mps.nc);
  c_line_y.setPointer(&llint_data,       padded_nc + thus_far, mps.nc);
  c_line_z.setPointer(&llint_data, (2 * padded_nc) + thus_far, mps.nc);
  c_line_x_overflow.setPointer(&int_data,                   thus_far, mps.nc);
  c_line_y_overflow.setPointer(&int_data,       padded_nc + thus_far, mps.nc);
  c_line_z_overflow.setPointer(&int_data, (2 * padded_nc) + thus_far, mps.nc);
  thus_far += 3 * padded_nc;
  a_abs_line_x.setPointer(&llint_data,                   thus_far, mps.na);
  a_abs_line_y.setPointer(&llint_data,       padded_na + thus_far, mps.na);
  a_abs_line_z.setPointer(&llint_data, (2 * padded_na) + thus_far, mps.na);
  a_abs_line_x_overflow.setPointer(&int_data,                   thus_far, mps.na);
  a_abs_line_y_overflow.setPointer(&int_data,       padded_na + thus_far, mps.na);
  a_abs_line_z_overflow.setPointer(&int_data, (2 * padded_na) + thus_far, mps.na);
  const int nbits = static_cast<int>(sizeof(uint)) * 8;
  coefficients.resize(64LLU * static_cast<size_t>(mps.na * mps.nb * mps.nc));
  if (ag_pointer != nullptr) {
    frozen_atoms.resize((ag_pointer->getAtomCount() + nbits - 1) / nbits);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::rebase_pointers() {
  a_line_x.swapTarget(&llint_data);
  a_line_y.swapTarget(&llint_data);
  a_line_z.swapTarget(&llint_data);
  b_line_x.swapTarget(&llint_data);
  b_line_y.swapTarget(&llint_data);
  b_line_z.swapTarget(&llint_data);
  c_line_x.swapTarget(&llint_data);
  c_line_y.swapTarget(&llint_data);
  c_line_z.swapTarget(&llint_data);
  a_abs_line_x.swapTarget(&llint_data);
  a_abs_line_y.swapTarget(&llint_data);
  a_abs_line_z.swapTarget(&llint_data);
  a_line_x_overflow.swapTarget(&int_data);
  a_line_y_overflow.swapTarget(&int_data);
  a_line_z_overflow.swapTarget(&int_data);
  b_line_x_overflow.swapTarget(&int_data);
  b_line_y_overflow.swapTarget(&int_data);
  b_line_z_overflow.swapTarget(&int_data);
  c_line_x_overflow.swapTarget(&int_data);
  c_line_y_overflow.swapTarget(&int_data);
  c_line_z_overflow.swapTarget(&int_data);
  a_abs_line_x_overflow.swapTarget(&int_data);
  a_abs_line_y_overflow.swapTarget(&int_data);
  a_abs_line_z_overflow.swapTarget(&int_data);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::validateMeshKind() const {
  switch (kind) {
  case GridDetail::OCCLUSION:    
    if (std::type_index(typeid(T)).hash_code() != ullint_type_index) {
      if (isScalarType<T>()) {
        rtErr("An occlusion mask requires " + getStormmScalarTypeName<ullint>() + " data type, "
              "not " + getStormmScalarTypeName<T>() + ".", "BackgroundMesh", "validateMeshKind");
      }
      else {
        rtErr("An occlusion mask requires " + getStormmScalarTypeName<ullint>() + " data type.",
              "BackgroundMesh", "validateMeshKind");
      }
    }
    break;
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    if (isFloatingPointScalarType<T>() == false) {
      rtErr("An non-bonded potential field requires " + getStormmScalarTypeName<float>() + " or " +
            getStormmScalarTypeName<double>() + " data type, not " + getStormmScalarTypeName<T>() +
            ".", "BackgroundMesh", "validateMeshKind");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::validateScalingBits() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();  
  if (ct == ullint_type_index || ct == float_type_index) {
    const std::vector<double> mesh_origin = measurements.getMeshOrigin<double>();
    std::vector<double> side_a = measurements.getMeshElementVector<double>(UnitCellAxis::A);
    std::vector<double> side_b = measurements.getMeshElementVector<double>(UnitCellAxis::B);
    std::vector<double> side_c = measurements.getMeshElementVector<double>(UnitCellAxis::C);
    elementwiseMultiply<double>(&side_a, measurements.getAxisElementCount(UnitCellAxis::A));
    elementwiseMultiply<double>(&side_b, measurements.getAxisElementCount(UnitCellAxis::B));
    elementwiseMultiply<double>(&side_c, measurements.getAxisElementCount(UnitCellAxis::C));
    double max_log_bound = 0.0;
    for (int i = 0; i < 2; i++) {
      const double dix = static_cast<double>(i) * side_a[0];
      const double diy = static_cast<double>(i) * side_a[1];
      const double diz = static_cast<double>(i) * side_a[2];
      for (int j = 0; j < 2; j++) {
        const double djx = static_cast<double>(j) * side_b[0];
        const double djy = static_cast<double>(j) * side_b[1];
        const double djz = static_cast<double>(j) * side_b[2];
        for (int k = 0; k < 2; k++) {
          const double dkx = static_cast<double>(k) * side_c[0];
          const double dky = static_cast<double>(k) * side_c[1];
          const double dkz = static_cast<double>(k) * side_c[2];
          const std::vector<double> mesh_corner = { mesh_origin[0] + dix + djx + dkx,
                                                    mesh_origin[1] + diy + djy + dky,
                                                    mesh_origin[2] + diz + djz + dkz };
          for (int m = 0; m < 3; m++) {
            if (fabs(mesh_corner[m]) > 0.0) {
              max_log_bound = std::max(max_log_bound, log2(fabs(mesh_corner[m])));
            }
          }
        }
      }
    }
    if (max_log_bound + measurements.getScalingBits() >= 63) {
      rtErr("Occlusion meshes, and non-bonded field meshes using single-precision coefficients, "
            "must keep all elements in the range +/- 32768.0.", "BackgroundMesh",
            "validateScalingBits");
    }
    if (measurements.getScalingBits() > mesh_nonoverflow_bits) {
      rtErr("Occlusion meshes, and non-bonded field meshes using single-precision coefficients, "
            "are assumed to not need positional representations in the extended fixed-precision "
            "model (the overflow bits).  A maximum of " + std::to_string(mesh_nonoverflow_bits) +
            " bits can be used in the represenation (this will rival or exceed the precision of "
            " 64-bit floating point representations for most of space).  A setting of " +
            std::to_string(measurements.getScalingBits()) + " is not acceptable.",
            "BackgroundMesh", "validateScalingBits");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::computeMeshAxisCoordinates() {

}

} // namespace structure
} // namespace stormm
