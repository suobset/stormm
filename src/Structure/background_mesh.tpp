// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const int vdw_type_in,
                                  const AtomGraph *ag_in, const MeshParameters &measurements_in) :
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
    neighbor_list{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"}
{
  validateMeshKind();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const int vdw_type_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const MeshParameters &measurements_in, const GpuDetails &gpu) :
    BackgroundMesh(kind_in, field_in, probe_radius_in, vdw_type_in, ag_in, measurements_in)
{
  allocate();

  // Loop over all atoms and apply the potential
  switch (kind) {
  case GridDetail::OCCLUSION:
    colorExclusionMesh(cf, gpu);
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
                                  const CoordinateFrame &cf, const double buffer,
                                  const double spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf, const double buffer,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, spacing, scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const int vdw_type_in, const AtomGraph *ag_in,
                                  const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, vdw_type_in, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, spacing, scale_bits_in), gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, field_in, 0.0, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, buffer, spacing, scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds,
                                  const double spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, std::vector<double>(3, spacing),
                                 scale_bits_in), gpu)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const GpuDetails &gpu) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, probe_radius_in, -1, ag_in, cf,
                 getMeasurements(ag_in, cf, mesh_bounds, spacing, scale_bits_in), gpu)
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
template <typename T>
BackgroundMeshReader<double, T> BackgroundMesh<T>::dpData(const HybridTargetLevel tier) const {
  return BackgroundMeshReader<double, T>(measurements.dpData(), kind, field, a_line_x.data(tier),
                                         a_line_y.data(tier), a_line_z.data(tier),
                                         b_line_x.data(tier), b_line_y.data(tier),
                                         b_line_z.data(tier), c_line_x.data(tier),
                                         c_line_y.data(tier), c_line_z.data(tier),
                                         a_line_x_overflow.data(tier),
                                         a_line_y_overflow.data(tier),
                                         a_line_z_overflow.data(tier),
                                         b_line_x_overflow.data(tier),
                                         b_line_y_overflow.data(tier),
                                         b_line_z_overflow.data(tier),
                                         c_line_x_overflow.data(tier),
                                         c_line_y_overflow.data(tier),
                                         c_line_z_overflow.data(tier), coefficients.data(tier),
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
                                         a_line_x_overflow.data(tier),
                                         a_line_y_overflow.data(tier),
                                         a_line_z_overflow.data(tier),
                                         b_line_x_overflow.data(tier),
                                         b_line_y_overflow.data(tier),
                                         b_line_z_overflow.data(tier),
                                         c_line_x_overflow.data(tier),
                                         c_line_y_overflow.data(tier),
                                         c_line_z_overflow.data(tier), coefficients.data(tier),
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
                                        a_line_x_overflow.data(tier), a_line_y_overflow.data(tier),
                                        a_line_z_overflow.data(tier), b_line_x_overflow.data(tier),
                                        b_line_y_overflow.data(tier), b_line_z_overflow.data(tier),
                                        c_line_x_overflow.data(tier), c_line_y_overflow.data(tier),
                                        c_line_z_overflow.data(tier), coefficients.data(tier),
                                        neighbor_list.data(tier),
                                        neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<float, T> BackgroundMesh<T>::spData(const HybridTargetLevel tier) {
  return BackgroundMeshReader<float, T>(measurements.spData(), kind, field, a_line_x.data(tier),
                                        a_line_y.data(tier), a_line_z.data(tier),
                                        b_line_x.data(tier), b_line_y.data(tier),
                                        b_line_z.data(tier), c_line_x.data(tier),
                                        c_line_y.data(tier), c_line_z.data(tier),
                                        a_line_x_overflow.data(tier), a_line_y_overflow.data(tier),
                                        a_line_z_overflow.data(tier), b_line_x_overflow.data(tier),
                                        b_line_y_overflow.data(tier), b_line_z_overflow.data(tier),
                                        c_line_x_overflow.data(tier), c_line_y_overflow.data(tier),
                                        c_line_z_overflow.data(tier), coefficients.data(tier),
                                        neighbor_list.data(tier),
                                        neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::setTopologyPointer(const AtomGraph *ag_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::colorExclusionMesh(const CoordinateFrame &cf, const GpuDetails &gpu) {

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
  const CoordinateFrameReader cfr = cf.data();
  
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
  int95ToDouble(&avx, &avy, &avz, a_line_x, a_line_x_overflow, a_line_y, a_line_y_overflow,
                a_line_z, a_line_z_overflow, mps.na, mps.inv_scale);
  int95ToDouble(&bvx, &bvy, &bvz, b_line_x, b_line_x_overflow, b_line_y, b_line_y_overflow,
                b_line_z, b_line_z_overflow, mps.nb, mps.inv_scale);
  int95ToDouble(&cvx, &cvy, &cvz, c_line_x, c_line_x_overflow, c_line_y, c_line_y_overflow,
                c_line_z, c_line_z_overflow, mps.nc, mps.inv_scale);
  const double dorig_x = int95ToDouble(mps.orig_x);
  const double dorig_y = int95ToDouble(mps.orig_y);
  const double dorig_z = int95ToDouble(mps.orig_z);
  addScalarToVector(&bvx, -dorig_x);
  addScalarToVector(&bvy, -dorig_y);
  addScalarToVector(&bvz, -dorig_z);
  addScalarToVector(&cvx, -dorig_x);
  addScalarToVector(&cvy, -dorig_y);
  addScalarToVector(&cvz, -dorig_z);

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
            coefficients[coef_base_idx + m] |= cube_buffer[m];
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const double padding, const double spacing,
                                                  const int scale_bits_in) const {
  return getMeasurements(ag, cf, padding, std::vector<double>(3, spacing));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const std::vector<double> &mesh_bounds,
                                                  const double spacing,
                                                  const int scale_bits_in) const {
  return getMeasurements(ag, cf, mesh_bounds, std::vector<double>(3, spacing));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const double padding,
                                                  const std::vector<double> &spacing,
                                                  const int scale_bits_in) const {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const CoordinateFrameReader cfr = cf.data();
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
  return getMeasurements(ag, cf, limits, spacing, scale_bits_in);
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
  if (spacing.size() != 3LLU) {
    rtErr("An array of three elements, the length, width, and height (Cartesian X, Y, and Z "
          "dimensions) of a a rectilinear mesh element, is required.  " +
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
  llint_data.resize(3 * (padded_na + padded_nb + padded_nc));
  int_data.resize(3 * (padded_na + padded_nb + padded_nc));
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
  const int nbits = static_cast<int>(sizeof(uint)) * 8;
  coefficients.resize(64LLU * static_cast<size_t>(mps.na * mps.nb * mps.nc));
  frozen_atoms.resize((ag_pointer->getAtomCount() + nbits - 1) / nbits);
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

} // namespace structure
} // namespace stormm
