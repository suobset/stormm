// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Txfrm, typename Tdata>
BackgroundMeshReader<Txfrm, Tdata>::
BackgroundMeshReader(const MeshParamKit<Txfrm> &dims_in, GridDetail kind_in,
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
BackgroundMeshWriter(const MeshParamKit<Txfrm> &dims_in, GridDetail kind_in,
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
                                  const MeshParameters &measurements_in) :
    measurements{measurements_in}, kind{kind_in}, field{field_in}, probe_radius{0.0},
    well_depth{0.0}, mixing_protocol{VdwCombiningRule::LORENTZ_BERTHELOT},
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
    ag_pointer{nullptr},
    cf_pointer{nullptr},
    frozen_atoms{HybridKind::ARRAY, "mesh_frozen_atoms"},
    neighbor_list{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"}
{
  validateMeshKind();
  allocate();
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const MeshParameters &measurements_in, const GpuDetails &gpu) :
    BackgroundMesh(kind_in, field_in, measurements_in)
{
  // Set the system
  setTopologyPointer(ag_in);
  setCoordinatePointer(cf_in);
  
  // Validate the mesh dimensions, then allocate memory
  validateScalingBits();
  allocate();
  computeMeshAxisCoordinates();

  // Map the potential and atomic near-neighbor interactions
  setProbeRadius(probe_radius_in);
  setWellDepth(well_depth_in);
  setCombiningRule(mixing_protocol_in);
  computeField(gpu);
  computeNeighborLists(gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                 mixing_protocol_in,
                 getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                 mixing_protocol_in, getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in),
                 gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                 mixing_protocol_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                 mixing_protocol_in,
                 getMeasurements(ag_in, cf_in, mesh_bounds, spacing, scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in,  ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
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
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
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
  return BackgroundMeshWriter<double, T>(measurements.dpData(), kind, field, a_line_x.data(tier),
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
  setMeshParameters(ag_in, cf_in, padding, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                          const double padding,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  ag_pointer = (ag_in == nullptr) ? ag_pointer : const_cast<AtomGraph*>(ag_in);
  cf_pointer = (cf_in == nullptr) ? cf_pointer : const_cast<CoordinateFrame*>(cf_in);
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(ag_pointer, cf_pointer, padding, spacing, actual_scale_bits);
  validateScalingBits();
  allocate();
  computeMeshAxisCoordinates();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const double spacing,
                                          const int scale_bits_in) {
  setMeshParameters(ag_pointer, cf_pointer, padding, std::vector<double>(3, spacing),
                    scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  setMeshParameters(ag_pointer, cf_pointer, padding, spacing, scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const double spacing, const int scale_bits_in) {
  setMeshParameters(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(mesh_bounds, spacing, actual_scale_bits);
  validateScalingBits();
  allocate();
  computeMeshAxisCoordinates();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setProbeRadius(const double probe_radius_in) {
  if (probe_radius_in < 0.0) {
    rtErr("A probe radius of " + realToString(probe_radius, 7, 4, NumberFormat::STANDARD_REAL) +
          " Angstroms is invalid.", "BackgroundMesh", "setProbeRadius");
  }
  probe_radius = probe_radius_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setWellDepth(const double well_depth_in) {
  if (well_depth_in < 0.0) {
    rtErr("A negative well depth of " +
          realToString(probe_radius, 7, 4, NumberFormat::STANDARD_REAL) + " kcal/mol is invalid.  "
          "Use positive numbers to define the depth of a potential energy minimum.",
          "BackgroundMesh", "setWellDepth");
  }
  well_depth = well_depth_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setCombiningRule(const VdwCombiningRule mixing_protocol_in) {
  mixing_protocol = mixing_protocol_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::computeField(const GpuDetails &gpu) {

  // Loop over all atoms and apply the potential
  switch (kind) {
  case GridDetail::OCCLUSION:
    colorExclusionMesh(gpu);
    break;
  case GridDetail::NONBONDED_FIELD:
    mapPureNonbondedPotential(gpu);
    break;
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::computeNeighborLists(const GpuDetails &gpu) {
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
      continue;
    }
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
    std::min(mesh_bounds[0], mesh_bounds[3]), std::min(mesh_bounds[1], mesh_bounds[4]),
    std::min(mesh_bounds[2], mesh_bounds[5]), std::max(mesh_bounds[0], mesh_bounds[3]),
    std::max(mesh_bounds[1], mesh_bounds[4]), std::max(mesh_bounds[2], mesh_bounds[5]) };  
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
  MeshParamKit<double> mps = measurements.dpData();
  const int padded_na = roundUp(mps.na + 1, warp_size_int);
  const int padded_nb = roundUp(mps.nb + 1, warp_size_int);
  const int padded_nc = roundUp(mps.nc + 1, warp_size_int);
  llint_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  int_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  a_line_x.setPointer(&llint_data,             0, mps.na + 1);
  a_line_y.setPointer(&llint_data,     padded_na, mps.na + 1);
  a_line_z.setPointer(&llint_data, 2 * padded_na, mps.na + 1);
  a_line_x_overflow.setPointer(&int_data,             0, mps.na + 1);
  a_line_y_overflow.setPointer(&int_data,     padded_na, mps.na + 1);
  a_line_z_overflow.setPointer(&int_data, 2 * padded_na, mps.na + 1);
  int thus_far = 3 * padded_na;
  b_line_x.setPointer(&llint_data,                   thus_far, mps.nb + 1);
  b_line_y.setPointer(&llint_data,       padded_nb + thus_far, mps.nb + 1);
  b_line_z.setPointer(&llint_data, (2 * padded_nb) + thus_far, mps.nb + 1);
  b_line_x_overflow.setPointer(&int_data,                   thus_far, mps.nb + 1);
  b_line_y_overflow.setPointer(&int_data,       padded_nb + thus_far, mps.nb + 1);
  b_line_z_overflow.setPointer(&int_data, (2 * padded_nb) + thus_far, mps.nb + 1);
  thus_far += 3 * padded_nb;
  c_line_x.setPointer(&llint_data,                   thus_far, mps.nc + 1);
  c_line_y.setPointer(&llint_data,       padded_nc + thus_far, mps.nc + 1);
  c_line_z.setPointer(&llint_data, (2 * padded_nc) + thus_far, mps.nc + 1);
  c_line_x_overflow.setPointer(&int_data,                   thus_far, mps.nc + 1);
  c_line_y_overflow.setPointer(&int_data,       padded_nc + thus_far, mps.nc + 1);
  c_line_z_overflow.setPointer(&int_data, (2 * padded_nc) + thus_far, mps.nc + 1);
  thus_far += 3 * padded_nc;
  a_abs_line_x.setPointer(&llint_data,                   thus_far, mps.na + 1);
  a_abs_line_y.setPointer(&llint_data,       padded_na + thus_far, mps.na + 1);
  a_abs_line_z.setPointer(&llint_data, (2 * padded_na) + thus_far, mps.na + 1);
  a_abs_line_x_overflow.setPointer(&int_data,                   thus_far, mps.na + 1);
  a_abs_line_y_overflow.setPointer(&int_data,       padded_na + thus_far, mps.na + 1);
  a_abs_line_z_overflow.setPointer(&int_data, (2 * padded_na) + thus_far, mps.na + 1);
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
  const MeshParamKit<double> mps = measurements.dpData();
  fixedPrecisionGrid(&a_line_x, &a_line_x_overflow, { 0LL, 0 }, mps.fp_invu[0]);
  fixedPrecisionGrid(&a_line_y, &a_line_y_overflow, { 0LL, 0 }, mps.fp_invu[1]);
  fixedPrecisionGrid(&a_line_z, &a_line_z_overflow, { 0LL, 0 }, mps.fp_invu[2]);
  fixedPrecisionGrid(&b_line_x, &b_line_x_overflow, { 0LL, 0 }, mps.fp_invu[3]);
  fixedPrecisionGrid(&b_line_y, &b_line_y_overflow, { 0LL, 0 }, mps.fp_invu[4]);
  fixedPrecisionGrid(&b_line_z, &b_line_z_overflow, { 0LL, 0 }, mps.fp_invu[5]);
  fixedPrecisionGrid(&c_line_x, &c_line_x_overflow, { 0LL, 0 }, mps.fp_invu[6]);
  fixedPrecisionGrid(&c_line_y, &c_line_y_overflow, { 0LL, 0 }, mps.fp_invu[7]);
  fixedPrecisionGrid(&c_line_z, &c_line_z_overflow, { 0LL, 0 }, mps.fp_invu[8]);
  fixedPrecisionGrid(&a_abs_line_x, &a_abs_line_x_overflow, mps.orig_x, mps.fp_invu[0]);
  fixedPrecisionGrid(&a_abs_line_y, &a_abs_line_y_overflow, mps.orig_y, mps.fp_invu[1]);
  fixedPrecisionGrid(&a_abs_line_z, &a_abs_line_z_overflow, mps.orig_z, mps.fp_invu[2]);
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
  const std::vector<bool> mobile_atom = ag_pointer->getAtomMobility();
  const MeshParamKit<double> mps = measurements.dpData();
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
  hostInt95ToDouble(&avx, &avy, &avz, a_abs_line_x.readHost(), a_abs_line_x_overflow.readHost(),
                    a_abs_line_y.readHost(), a_abs_line_y_overflow.readHost(),
                    a_abs_line_z.readHost(), a_abs_line_z_overflow.readHost(), mps.inv_scale);
  hostInt95ToDouble(&bvx, &bvy, &bvz, b_line_x.readHost(), b_line_x_overflow.readHost(),
                    b_line_y.readHost(), b_line_y_overflow.readHost(), b_line_z.readHost(),
                    b_line_z_overflow.readHost(), mps.inv_scale);
  hostInt95ToDouble(&cvx, &cvy, &cvz, c_line_x.readHost(), c_line_x_overflow.readHost(),
                    c_line_y.readHost(), c_line_y_overflow.readHost(), c_line_z.readHost(),
                    c_line_z_overflow.readHost(), mps.inv_scale);
  const double dorig_x = hostInt95ToDouble(mps.orig_x);
  const double dorig_y = hostInt95ToDouble(mps.orig_y);
  const double dorig_z = hostInt95ToDouble(mps.orig_z);

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
    if (mobile_atom[pos]) {
      continue;
    }
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
                  const int cubelet_idx = (((cubelet_k * 4) + cubelet_j) * 4) + cubelet_i;
                  const int bit_i = grid_i - (4 * cubelet_i);
                  const int bit_j = grid_j - (4 * cubelet_j);
                  const int bit_k = grid_k - (4 * cubelet_k);
                  const int bit_idx = (((bit_k * 4) + bit_j) * 4) + bit_i;
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
                                       static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i);
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
void BackgroundMesh<T>::computeCoefficients(const std::vector<double> &u_grid,
                                            const std::vector<double> &dudx_grid,
                                            const std::vector<double> &dudy_grid,
                                            const std::vector<double> &dudz_grid,
                                            const std::vector<double> &dudxy_grid,
                                            const std::vector<double> &dudxz_grid,
                                            const std::vector<double> &dudyz_grid,
                                            const std::vector<double> &dudxyz_grid,
                                            const std::vector<double> &dudxx_grid,
                                            const std::vector<double> &dudyy_grid,
                                            const std::vector<double> &dudzz_grid) {
  
  // Get the measurements abstract locally--better than passing a reference of the same size
  // and then having to de-reference it many times.
  const MeshParamKit<double> mps = measurements.dpData();

  // Compute the weights matrix
  std::vector<double> tc_weights = getTricubicMatrix().readHost();  

  // Lay out the mesh cell bounds
  std::vector<double> tc_bounds;
  const UnitCellType unit_cell = measurements.getMeshCellType();
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    tc_bounds.resize(6);
    tc_bounds[3] = hostInt95ToDouble(mps.fp_invu[0]);
    tc_bounds[4] = hostInt95ToDouble(mps.fp_invu[4]);
    tc_bounds[5] = hostInt95ToDouble(mps.fp_invu[8]);
    break;
  case UnitCellType::TRICLINIC:
    tc_bounds.resize(12);
    for (int i = 0; i < 9; i++) {
      tc_bounds[i + 3] = hostInt95ToDouble(mps.fp_invu[i]);
    }
    break;
  }

  // Set constants and pointers to mesh lines for more rapid access
  const size_t ngrid_a = mps.na + 1;
  const size_t ngrid_b = mps.nb + 1;
  const size_t ngrid_c = mps.nc + 1;
  const llint* bvec_x_ptr = b_line_x.data();
  const llint* bvec_y_ptr = b_line_y.data();
  const llint* bvec_z_ptr = b_line_z.data();
  const llint* cvec_x_ptr = c_line_x.data();
  const llint* cvec_y_ptr = c_line_y.data();
  const llint* cvec_z_ptr = c_line_z.data();
  const llint* bvec_x_overflow_ptr = b_line_x.data();
  const llint* bvec_y_overflow_ptr = b_line_y.data();
  const llint* bvec_z_overflow_ptr = b_line_z.data();
  const llint* cvec_x_overflow_ptr = c_line_x.data();
  const llint* cvec_y_overflow_ptr = c_line_y.data();
  const llint* cvec_z_overflow_ptr = c_line_z.data();
  T* coeff_ptr = coefficients.data();
  
  // Compute coefficeints
  std::vector<double> u_elem(8), dudx_elem(8), dudy_elem(8), dudz_elem(8);
  std::vector<double> dudxy_elem(8), dudxz_elem(8), dudyz_elem(8), dudxyz_elem(8);
  std::vector<double> dudxx_elem(8), dudyy_elem(8), dudzz_elem(8);
  for (int i = 0; i < mps.na; i++) {
    const int95_t mesh_ax = { a_abs_line_x.readHost(i), a_abs_line_x_overflow.readHost(i) };
    const int95_t mesh_ay = { a_abs_line_y.readHost(i), a_abs_line_y_overflow.readHost(i) };
    const int95_t mesh_az = { a_abs_line_z.readHost(i), a_abs_line_z_overflow.readHost(i) };
    for (int j = 0; j < mps.nb; j++) {
      const int95_t mesh_abx = hostSplitFPSum(mesh_ax, bvec_x_ptr[j], bvec_x_overflow_ptr[j]);
      const int95_t mesh_aby = hostSplitFPSum(mesh_ay, bvec_y_ptr[j], bvec_y_overflow_ptr[j]);
      const int95_t mesh_abz = hostSplitFPSum(mesh_az, bvec_z_ptr[j], bvec_z_overflow_ptr[j]);
      for (int k = 0; k < mps.nc; k++) {
        const int95_t mesh_abcx = hostSplitFPSum(mesh_abx, cvec_x_ptr[k], cvec_x_overflow_ptr[k]);
        const int95_t mesh_abcy = hostSplitFPSum(mesh_aby, cvec_y_ptr[k], cvec_y_overflow_ptr[k]);
        const int95_t mesh_abcz = hostSplitFPSum(mesh_abz, cvec_z_ptr[k], cvec_z_overflow_ptr[k]);

        // Compose the input vectors
        for (int ci = 0; ci < 2; ci++) {
          const size_t ici = i + ci;
          for (int cj = 0; cj < 2; cj++) {
            const size_t jcj = j + cj;
            for (int ck = 0; ck < 2; ck++) {
              const size_t nv = (((ck * 2) + cj) * 2) + ci;
              const size_t kck = k + ck;
              const size_t nijk = (((kck * ngrid_b) + jcj) * ngrid_a) + ici;
              u_elem[nv] = u_grid[nijk];
              dudx_elem[nv] = dudx_grid[nijk];
              dudy_elem[nv] = dudy_grid[nijk];
              dudz_elem[nv] = dudz_grid[nijk];
              dudxy_elem[nv] = dudxy_grid[nijk];
              dudxz_elem[nv] = dudxz_grid[nijk];
              dudyz_elem[nv] = dudyz_grid[nijk];
              dudxyz_elem[nv] = dudxyz_grid[nijk];
              if (unit_cell == UnitCellType::TRICLINIC) {
                dudxx_elem[nv] = dudxx_grid[nijk];
                dudyy_elem[nv] = dudyy_grid[nijk];
                dudzz_elem[nv] = dudzz_grid[nijk];
              }
            }
          }
        }

        // Compose the mesh element.  Complete the bounds array by adding the origin, then
        // compute the tricubic coefficients.
        tc_bounds[0] = hostInt95ToDouble(mesh_abcx);
        tc_bounds[1] = hostInt95ToDouble(mesh_abcy);
        tc_bounds[2] = hostInt95ToDouble(mesh_abcz);
        const TricubicCell<double> tc_elem(tc_weights, tc_bounds, u_elem, dudx_elem, dudy_elem,
                                           dudz_elem, dudxy_elem, dudxz_elem, dudyz_elem,
                                           dudxyz_elem, dudxx_elem, dudyy_elem, dudzz_elem);
        const std::vector<double> tc_coeffs = tc_elem.getCoefficients();
        const size_t tc_offset = static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i) * 64LLU;
        for (size_t cpos = 0; cpos < 64; cpos++) {
          coeff_ptr[tc_offset + cpos] = tc_coeffs[cpos];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::mapPureNonbondedPotential(const GpuDetails &gpu,
                                                  const std::vector<double> &eps_table,
                                                  const std::vector<double> &sigma_table) {

  // Check for problematic input cases
  switch (field) {
  case NonbondedPotential::ELECTROSTATIC:
    break;
  case NonbondedPotential::VAN_DER_WAALS:
    if (mixing_protocol == VdwCombiningRule::NBFIX &&
        (static_cast<int>(eps_table.size())   != ag_pointer->getAtomTypeCount() ||
         static_cast<int>(sigma_table.size()) != ag_pointer->getAtomTypeCount())) {
      rtErr("In order to map a van-der Waals potential using the " +
            getEnumerationName(mixing_protocol) + " method, tables of the particle-particle well "
            "depth and interaction sigma values must be provided for the probe interacting with "
            "all types present in the receptor's topology.  Number of atom types in the "
            "topology: " + std::to_string(ag_pointer->getAtomTypeCount()) + ".  Number of atom "
            "types present in each table provided: " + std::to_string(eps_table.size()) +
            " (epsilon) and " + std::to_string(sigma_table.size()) + " (sigma).", "BackgroundMesh",
            "mapPureNonbondedPotential");
    }
    break;
  case NonbondedPotential::CLASH:

    // This would have to be an occlusion mask and would not have an associated potential field.
    break;
  }
  
  // Use the HPC kernel to color the mesh if a GPU is available
#ifdef STORMM_USE_HPC
  if (gpu != null_gpu) {
    
    coefficients.download();
    return;
  }
#endif

  // Color the mesh on the CPU
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const std::vector<bool> mobile_atom = ag_pointer->getAtomMobility();
  const MeshParamKit<double> mps = measurements.dpData();
  const CoordinateFrameReader cfr = cf_pointer->data();
  
  // Initialize the mesh in CPU memory.
  const int n_elem = mps.na * mps.nb * mps.nc;
  T* coeff_ptr = coefficients.data();
  const T zero = 0.0;
  for (int pos = 0; pos < n_elem; pos++) {
    coeff_ptr[pos] = zero;
  }

  // Set pointers to mesh lines for more rapid access
  const llint* bvec_x_ptr = b_line_x.data();
  const llint* bvec_y_ptr = b_line_y.data();
  const llint* bvec_z_ptr = b_line_z.data();
  const llint* cvec_x_ptr = c_line_x.data();
  const llint* cvec_y_ptr = c_line_y.data();
  const llint* cvec_z_ptr = c_line_z.data();
  const int* bvec_x_overflow_ptr = b_line_x_overflow.data();
  const int* bvec_y_overflow_ptr = b_line_y_overflow.data();
  const int* bvec_z_overflow_ptr = b_line_z_overflow.data();
  const int* cvec_x_overflow_ptr = c_line_x_overflow.data();
  const int* cvec_y_overflow_ptr = c_line_y_overflow.data();
  const int* cvec_z_overflow_ptr = c_line_z_overflow.data();
  
  // Create a series of eight three-dimensional grids that will yield the mesh coefficients for
  // each element.  This allocates a small amount of additional memory, relative to the mesh that
  // will be produced, but reduces the pre-computations for getting those elements by a factor
  // approaching eight.  Since this is all being done on the CPU, use std::vector<double> objects
  // and compute the coefficients in double-precision (64-bit), even if they will be used in
  // 32-bit representations.  Here, "grid" will refer to the three-dimensional series of regular
  // points on which the electrostatic potential U, dU/dx, dU/dy, ..., d2U/dxdy, ..., d3U/dxdydz
  // are computed, while "mesh" will refer to the collection of coefficients for each corresponding
  // element.
  const size_t ngrid_a = mps.na + 1;
  const size_t ngrid_b = mps.nb + 1;
  const size_t ngrid_c = mps.nc + 1;
  const size_t ngabc = ngrid_a * ngrid_b * ngrid_c;
  std::vector<double> u_grid(ngabc, 0.0), dudx_grid(ngabc, 0.0), dudy_grid(ngabc, 0.0);
  std::vector<double> dudz_grid(ngabc, 0.0), dudxy_grid(ngabc, 0.0), dudxz_grid(ngabc, 0.0);
  std::vector<double> dudyz_grid(ngabc, 0.0), dudxyz_grid(ngabc, 0.0);
  std::vector<double> dudxx_grid, dudyy_grid, dudzz_grid;
  const UnitCellType unit_cell = measurements.getMeshCellType();
  if (unit_cell == UnitCellType::TRICLINIC) {
    dudxx_grid.resize(ngabc, 0.0);
    dudyy_grid.resize(ngabc, 0.0);
    dudzz_grid.resize(ngabc, 0.0);
  }
  for (int pos = 0; pos < nbk.natom; pos++) {
    if (mobile_atom[pos]) {
      continue;
    }
    
    // Compute the particle position in the mesh's native fixed-precision format
    const int95_t atom_x = hostDoubleToInt95(cfr.xcrd[pos] * mps.scale);
    const int95_t atom_y = hostDoubleToInt95(cfr.ycrd[pos] * mps.scale);
    const int95_t atom_z = hostDoubleToInt95(cfr.zcrd[pos] * mps.scale);

    // Get the particle's charge.  The mesh coefficients will be computed in kcal/mol-A^n,
    // n = { 0, 1, 2, 3 } depending on the degree of the derivative.
    const double atom_q = nbk.charge[pos] * nbk.coulomb_constant;
    const int atom_lj_idx = nbk.lj_idx[pos];
    const double atom_lja = nbk.lja_coeff[atom_lj_idx * (nbk.n_lj_types + 1)];
    const double atom_ljb = nbk.ljb_coeff[atom_lj_idx * (nbk.n_lj_types + 1)];
    const double atom_sigma = (atom_lja >= 1.0e-6 && atom_ljb >= 1.0e-6) ?
                              pow(atom_lja / atom_ljb, 1.0 / 6.0) : 0.0;
    const double atom_eps = (atom_sigma > 0.0) ? 0.25 * atom_ljb / pow(atom_sigma, 6.0) :
                                                 0.0;
    double pair_eps, pair_sigma;
    switch (mixing_protocol) {
    case VdwCombiningRule::LORENTZ_BERTHELOT:
      pair_eps = (atom_eps * well_depth > 1.0e-6) ? sqrt(atom_eps * well_depth) : 0.0;
      pair_sigma = 0.5 * (atom_sigma + probe_radius);
      break;
    case VdwCombiningRule::GEOMETRIC:
      pair_eps = (atom_eps * well_depth > 1.0e-6) ? sqrt(atom_eps * well_depth) : 0.0;
      pair_sigma = 0.5 * (atom_sigma + probe_radius);
      break;
    case VdwCombiningRule::NBFIX:
      pair_eps   = eps_table[atom_lj_idx];
      pair_sigma = sigma_table[atom_lj_idx];
      break;
    }
    const double pair_lja = 4.0 * pair_eps * pow(pair_sigma, 12.0);
    const double pair_ljb = 4.0 * pair_eps * pow(pair_sigma, 6.0);
    
    // Loop over all grid points
    for (size_t i = 0; i < ngrid_a; i++) {
      const int95_t mesh_ax = { a_abs_line_x.readHost(i), a_abs_line_x_overflow.readHost(i) };
      const int95_t mesh_ay = { a_abs_line_y.readHost(i), a_abs_line_y_overflow.readHost(i) };
      const int95_t mesh_az = { a_abs_line_z.readHost(i), a_abs_line_z_overflow.readHost(i) };
      for (size_t j = 0; j < ngrid_b; j++) {
        const int95_t mesh_abx = hostSplitFPSum(mesh_ax, bvec_x_ptr[j], bvec_x_overflow_ptr[j]);
        const int95_t mesh_aby = hostSplitFPSum(mesh_ay, bvec_y_ptr[j], bvec_y_overflow_ptr[j]);
        const int95_t mesh_abz = hostSplitFPSum(mesh_az, bvec_z_ptr[j], bvec_z_overflow_ptr[j]);
        for (size_t k = 0; k < ngrid_c; k++) {
          const int95_t mesh_abcx = hostSplitFPSum(mesh_abx, cvec_x_ptr[k],
                                                   cvec_x_overflow_ptr[k]);
          const int95_t mesh_abcy = hostSplitFPSum(mesh_aby, cvec_y_ptr[k],
                                                   cvec_y_overflow_ptr[k]);
          const int95_t mesh_abcz = hostSplitFPSum(mesh_abz, cvec_z_ptr[k],
                                                   cvec_z_overflow_ptr[k]);

          // Compute the displacements using the mesh's fixed-precision representation, then
          // immediately convert to double for real-valued computations.
          const int95_t fp_disp_x = hostSplitFPSum(mesh_abcx, -atom_x.x, -atom_x.y);
          const int95_t fp_disp_y = hostSplitFPSum(mesh_abcy, -atom_y.x, -atom_y.y);
          const int95_t fp_disp_z = hostSplitFPSum(mesh_abcz, -atom_z.x, -atom_z.y);
          const double disp_x = hostInt95ToDouble(fp_disp_x) * mps.inv_scale;
          const double disp_y = hostInt95ToDouble(fp_disp_y) * mps.inv_scale;
          const double disp_z = hostInt95ToDouble(fp_disp_z) * mps.inv_scale;
          const double r2 = (disp_x * disp_x) + (disp_y * disp_y) + (disp_z * disp_z);
          const double r  = sqrt(r2);
          const double invr2 = 1.0 / r2;

          // Compute one of a variety of potentials
          double u, du_dx, du_dy, du_dz, du_dxy, du_dxz, du_dyz, du_dxyz, du_dxx, du_dyy, du_dzz;
          switch (field) {
          case NonbondedPotential::ELECTROSTATIC:
            u = atom_q / r;
            du_dx = -u * disp_x * invr2;
            du_dy = -u * disp_y * invr2;
            du_dz = -u * disp_z * invr2;
            du_dxy = -3.0 * du_dx * disp_y * invr2;
            du_dxz = -3.0 * du_dx * disp_z * invr2;
            du_dyz = -3.0 * du_dy * disp_z * invr2;
            du_dxyz = -5.0 * du_dxy * disp_z * invr2;
            if (unit_cell == UnitCellType::TRICLINIC) {
              du_dxx = -3.0 * du_dx * disp_x * invr2;
              du_dyy = -3.0 * du_dy * disp_y * invr2;
              du_dzz = -3.0 * du_dz * disp_z * invr2;
            }
            break;
          case NonbondedPotential::VAN_DER_WAALS:
            {
              const double invr6 = invr2 * invr2 * invr2;
              const double invr8 = invr6 * invr2;
              const double invr10 = invr8 * invr2;
              const double d_factor = ((6.0 * pair_ljb) - (12.0 * pair_lja * invr6)) * invr8;
              const double d2_factor = ((168.0 * pair_lja * invr6) - (48.0 * pair_ljb)) * invr10;
              const double d3_factor = ((480.0 * pair_ljb) - (2688.0 * pair_lja * invr6)) *
                                       invr6 * invr6;
              u = ((pair_lja * invr6) - pair_ljb) * invr6;
              du_dx = d_factor * disp_x;
              du_dy = d_factor * disp_y;
              du_dz = d_factor * disp_z;
              du_dxy = d2_factor * disp_x * disp_y;
              du_dxz = d2_factor * disp_x * disp_z;
              du_dyz = d2_factor * disp_y * disp_z;
              du_dxyz = d3_factor * disp_x * disp_y * disp_z;
              if (unit_cell == UnitCellType::TRICLINIC) {
                du_dxx = d2_factor * disp_x * disp_x;
                du_dyy = d2_factor * disp_y * disp_y;
                du_dzz = d2_factor * disp_z * disp_z;
              }
            }
            break;
          case NonbondedPotential::CLASH:
            break;
          }

          // Log the results
          const size_t nijk = (((k * ngrid_b) + j) * ngrid_a) + i;
          u_grid[nijk] += u;
          dudx_grid[nijk] += du_dx;
          dudy_grid[nijk] += du_dy;
          dudz_grid[nijk] += du_dz;
          dudxy_grid[nijk] += du_dxy;
          dudxz_grid[nijk] += du_dxz;
          dudyz_grid[nijk] += du_dyz;
          dudxyz_grid[nijk] += du_dxyz;
          if (unit_cell == UnitCellType::TRICLINIC) {
            dudxx_grid[nijk] += du_dxx;
            dudyy_grid[nijk] += du_dyy;
            dudzz_grid[nijk] += du_dzz;
          }
        }
      }
    }
  }

  // Convert the computed grid data into a functional mesh
  computeCoefficients(u_grid, dudx_grid, dudy_grid, dudz_grid, dudxy_grid, dudxz_grid,
                      dudyz_grid, dudxyz_grid);
}

//-------------------------------------------------------------------------------------------------
template <typename Txfrm, typename Tdata, typename Tcoord, typename Tprop>
std::vector<Txfrm> interpolate(const BackgroundMeshReader<Txfrm, Tdata> &bgmr, const Tcoord* xcrd,
                               const Tcoord* ycrd, const Tcoord* zcrd, const Tprop* prop_a,
                               const Tprop* prop_b, const int natom, const int* xcrd_ovrf,
                               const int* ycrd_ovrf, const int* zcrd_ovrf,
                               const int coord_scaling_bits) {

  // Check whether the grid works in a format that would restrict the fixed-precision math to just
  // the first 64 bits of the representation.
  const size_t data_ct = std::type_index(typeid(Tdata)).hash_code();
  const bool mesh_data_is_double = (data_ct == double_type_index);  
  const bool coord_is_integral = isSignedIntegralScalarType<Tcoord>();
  const bool coordinate_overflow_active = (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr &&
                                           zcrd_ovrf != nullptr);
  std::vector<Txfrm> result;

  // Loop over all particles
  for (int pos = 0; pos < natom; pos++) {
  
    // Determine the grid bin using a transformation into mesh space followed by a local check and,
    // if necessary, a correction.  The precision of the mesh coefficients will determine the
    // precision of the fixed-precision coordinate representation, following the standard set forth
    // in validateScalingBits() above.
    if (mesh_data_is_double) {
      int95_t ixcrd, iycrd, izcrd;
      if (coord_is_integral) {

        // The coordinates may not have the same precision as the positions expected by the mesh.
        // Harmonize the representations.
        if (coord_scaling_bits == bgmr.dims.scale_bits) {
          ixcrd.x = xcrd[pos];
          iycrd.x = ycrd[pos];
          izcrd.x = zcrd[pos];
          if (coordinate_overflow_active) {
            ixcrd.y = xcrd[pos];
            iycrd.y = ycrd[pos];
            izcrd.y = zcrd[pos];          
          }
        }
        else {
          int95_t orep_xcrd, orep_ycrd, orep_zcrd;
          if (coordinate_overflow_active) {
            orep_xcrd = { xcrd[pos], xcrd_ovrf[pos] };
            orep_ycrd = { ycrd[pos], ycrd_ovrf[pos] };
            orep_zcrd = { zcrd[pos], zcrd_ovrf[pos] };
          }
          else {
            orep_xcrd = { xcrd[pos], 0 };
            orep_ycrd = { ycrd[pos], 0 };
            orep_zcrd = { zcrd[pos], 0 };
          }
          ixcrd = changeFPBits(orep_xcrd, coord_scaling_bits, bgmr.dims.scale_bits);
          iycrd = changeFPBits(orep_ycrd, coord_scaling_bits, bgmr.dims.scale_bits);
          izcrd = changeFPBits(orep_zcrd, coord_scaling_bits, bgmr.dims.scale_bits);
        }
      }
      else {
        ixcrd = hostDoubleToInt95(xcrd[pos]);
        iycrd = hostDoubleToInt95(ycrd[pos]);
        izcrd = hostDoubleToInt95(zcrd[pos]);
      }

      // Obtain the coordinates of the atom, relative to the mesh origin, in the precision of the
      // transformation matrices.  Estimate the appropriate mesh element, then test by computing
      // the location of the atom relative to the origin of this element.
      const int95_t ipt_rel_x = hostSplitFPSum(ixcrd, -bgmr.dims.orig_x.x, -bgmr.dims.orig_x.y);
      const int95_t ipt_rel_y = hostSplitFPSum(iycrd, -bgmr.dims.orig_y.x, -bgmr.dims.orig_y.y);
      const int95_t ipt_rel_z = hostSplitFPSum(izcrd, -bgmr.dims.orig_z.x, -bgmr.dims.orig_z.y);
      const Txfrm pt_rel_x = hostInt95ToDouble(ipt_rel_x);
      const Txfrm pt_rel_y = hostInt95ToDouble(ipt_rel_y);
      const Txfrm pt_rel_z = hostInt95ToDouble(ipt_rel_z);
      const Txfrm pt_grid_a = (bgmr.dims.umat[0] * pt_rel_x) + (bgmr.dims.umat[3] * pt_rel_y) +
                              (bgmr.dims.umat[6] * pt_rel_z);
      const Txfrm pt_grid_b = (bgmr.dims.umat[1] * pt_rel_x) + (bgmr.dims.umat[4] * pt_rel_y) +
                              (bgmr.dims.umat[7] * pt_rel_z);
      const Txfrm pt_grid_c = (bgmr.dims.umat[2] * pt_rel_x) + (bgmr.dims.umat[5] * pt_rel_y) +
                              (bgmr.dims.umat[8] * pt_rel_z);
      int cell_a = floor(pt_grid_a);
      int cell_b = floor(pt_grid_b);
      int cell_c = floor(pt_grid_c);

      // If the initial estimate is already off the grid, bail out
      if (cell_a < 0 || cell_a >= bgmr.dims.na || cell_b < 0 || cell_b >= bgmr.dims.nb ||
          cell_c < 0 || cell_c >= bgmr.dims.nc) {
        result[pos] = 0.0;
        continue;
      }
      int95_t element_origin_x = { bgmr.avec_abs_x[cell_a], bgmr.avec_abs_x_ovrf[cell_a] };
      int95_t element_origin_y = { bgmr.avec_abs_y[cell_a], bgmr.avec_abs_y_ovrf[cell_a] };
      int95_t element_origin_z = { bgmr.avec_abs_z[cell_a], bgmr.avec_abs_z_ovrf[cell_a] };
      element_origin_x = hostSplitFPSum(element_origin_x, bgmr.bvec_x[cell_b],
                                        bgmr.bvec_x_ovrf[cell_b]);
      element_origin_y = hostSplitFPSum(element_origin_y, bgmr.bvec_y[cell_b],
                                        bgmr.bvec_y_ovrf[cell_b]);
      element_origin_z = hostSplitFPSum(element_origin_z, bgmr.bvec_z[cell_b],
                                        bgmr.bvec_z_ovrf[cell_b]);
      element_origin_x = hostSplitFPSum(element_origin_x, bgmr.cvec_x[cell_b],
                                        bgmr.cvec_x_ovrf[cell_b]);
      element_origin_y = hostSplitFPSum(element_origin_y, bgmr.cvec_y[cell_b],
                                        bgmr.cvec_y_ovrf[cell_b]);
      element_origin_z = hostSplitFPSum(element_origin_z, bgmr.cvec_z[cell_b],
                                        bgmr.cvec_z_ovrf[cell_b]);
      const int95_t idisp_x = hostSplitFPSum(ixcrd, -element_origin_x.x, -element_origin_x.y);
      const int95_t idisp_y = hostSplitFPSum(iycrd, -element_origin_y.x, -element_origin_y.y);
      const int95_t idisp_z = hostSplitFPSum(izcrd, -element_origin_z.x, -element_origin_z.y);
      const Txfrm disp_x = hostInt95ToDouble(idisp_x);
      const Txfrm disp_y = hostInt95ToDouble(idisp_y);
      const Txfrm disp_z = hostInt95ToDouble(idisp_z);
      Txfrm test_a = (bgmr.dims.umat[0] * disp_x) + (bgmr.dims.umat[3] * disp_y) +
                     (bgmr.dims.umat[6] * disp_z);
      Txfrm test_b = (bgmr.dims.umat[1] * disp_x) + (bgmr.dims.umat[4] * disp_y) +
                     (bgmr.dims.umat[7] * disp_z);
      Txfrm test_c = (bgmr.dims.umat[2] * disp_x) + (bgmr.dims.umat[5] * disp_y) +
                     (bgmr.dims.umat[8] * disp_z);
      cell_a += (test_a < 0.0) - (test_a >= 1.0);
      cell_b += (test_b < 0.0) - (test_b >= 1.0);
      cell_c += (test_c < 0.0) - (test_c >= 1.0);

    }
    else {
    }
  
    switch (bgmr.kind) {
    case GridDetail::OCCLUSION:
      break;
    case GridDetail::NONBONDED_FIELD:
    case GridDetail::NONBONDED_ATOMIC:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Txfrm, typename Tdata, typename Tcoord, typename Tprop>
Txfrm interpolate(const BackgroundMeshReader<Txfrm, Tdata> &bgmr, const std::vector<Tcoord> &xcrd,
                  const std::vector<Tcoord> &ycrd, const std::vector<Tcoord> &zcrd,
                  const std::vector<Tprop> &prop_a, const std::vector<Tprop> &prop_b,
                  const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                  const int coord_scaling_bits) {

}

} // namespace structure
} // namespace stormm
