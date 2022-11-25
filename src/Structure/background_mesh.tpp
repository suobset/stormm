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
                                  const double probe_radius_in, const int vdw_type_in,
                                  const AtomGraph *ag_in, const CoordinateFrame &cf,
                                  const MeshParameters &measurements_in, const GpuDetails &gpu) :
  PureMesh(kind_in, field_in, probe_radius_in, vdw_type_in, ag_in, measurements_in)
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
                                  const std::vector<double> &spacing, const int scale_bits,
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
                                         neighblor_lists.data(tier),
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
                                         neighblor_lists.data(tier),
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
                                        neighblor_lists.data(tier),
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
                                        neighblor_lists.data(tier),
                                        neighbor_list_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
void BackgroundMesh<T>::colorExclusionMesh(const CoordinateFrame &cf) {

  // Use the HPC kernel to color the mesh if a GPU is available
#ifdef STORMM_USE_HPC
  if (gpu != null_gpu) {
    coefficients.download();
    return;
  }
#endif
  
  // Color the mesh on the CPU
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const MeshParameters mps = measurements.dpData();
  const int lj_idx_offset = nbk.n_lj_types + 1;
  const ullint* excl_mesh_ptr = exclusion_mesh.data();
  for (int pos = 0; pos < nbk.natom; pos++) {
    const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
    const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                         1.0 / 6.0);
    const double color_radius = atom_radius + probe_radius;
    const int ixmin = std::max(floor((xcrd[pos] - color_radius - mps.orig_x) /
                                     mps.orig), 0);
    const int iymin = std::max(floor((ycrd[pos] - color_radius - exclusion_mesh_y_origin) /
                                     exclusion_mesh_spacing), 0);
    const int izmin = std::max(floor((zcrd[pos] - color_radius - exclusion_mesh_z_origin) /
                                     exclusion_mesh_spacing), 0);
    const int ixmax = std::min(ceil((xcrd[pos] + color_radius - exclusion_mesh_x_origin) /
                                    exclusion_mesh_spacing), exclusion_mesh_nx);
    const int iymax = std::min(ceil((ycrd[pos] + color_radius - exclusion_mesh_y_origin) /
                                    exclusion_mesh_spacing), exclusion_mesh_ny);
    const int izmax = std::min(ceil((zcrd[pos] + color_radius - exclusion_mesh_z_origin) /
                                    exclusion_mesh_spacing), exclusion_mesh_nz);
    double dx = exclusion_mesh_x_origin + (exclusion_mesh_spacing * ixmin);
    const int mesh_supercube_x_dim = exclusion_mesh_nx / 16;
    const int mesh_supercube_y_dim = exclusion_mesh_ny / 16;
    const int mesh_supercube_z_dim = exclusion_mesh_nz / 16;
    const int mesh_cubelet_x_dim = exclusion_mesh_nx / 4;
    const int mesh_cubelet_y_dim = exclusion_mesh_ny / 4;
    const int mesh_cubelet_z_dim = exclusion_mesh_nz / 4;
    for (int i = ixmin; i < ixmax; i++) {
      double dy = exclusion_mesh_y_origin + (exclusion_mesh_spacing * iymin);
      for (int j = iymin; j < iymax; j++) {
        double dz = exclusion_mesh_z_origin + (exclusion_mesh_spacing * izmin);
        for (int k = izmin; k < izmax; k++) {
          if ((dx * dx) + (dy * dy) + (dz * dz)) {

            // Determine the supercube, then the cubelet to which this bit belongs, then
            // the number of the bit within the cubelet.
            const int supercube_ix = i / 16;
            const int supercube_iy = j / 16;
            const int supercube_iz = k / 16;
            const size_t supercube_idx = (((supercube_ix * mesh_supercube_y_dim) + supercube_iy) *
                                          mesh_supercube_z_dim) + supercube_iz;
            const int cubelet_ix = (i - (supercube_ix * 16)) / 4;
            const int cubelet_iy = (j - (supercube_iy * 16)) / 4;
            const int cubelet_iz = (k - (supercube_iz * 16)) / 4;
            const size_t cubelet_idx = (((cubelet_ix * 4) + cubelet_iy) * 4) + cubelet_iz;
            const int bit_ix = i - (supercube_ix * 16) - (cubelet_ix * 4);
            const int bit_iy = j - (supercube_iy * 16) - (cubelet_iy * 4);
            const int bit_iz = k - (supercube_iz * 16) - (cubelet_iz * 4);
            const int bit_idx = (((bit_ix * 4) + bit_iy) * 4) + bit_iz;
            excl_mesh_ptr[(supercube_idx * 64LLU) + cubelet_idx] |= (0x1 << bit_idx);
          }
          dz += exclusion_mesh_spacing;
        }
        dy += exclusion_mesh_spacing;
      }
      dx += exclusion_mesh_spacing;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const double padding, double spacing) const {
  return getMeasurements(ag, cf, padding, std::vector<double>(3, spacing));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const std::vector<double> &mesh_bounds,
                                                  double spacing) const {
  return getMeasurements(ag, cf, mesh_bounds, std::vector<double>(3, spacing));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters BackgroundMesh<T>::getMeasurements(const AtomGraph *ag, const CoordinateFrame &cf,
                                                  const double padding,
                                                  const std::vector<double> &spacing) const {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const CoordinateFrameReader cfr = cf.data();
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bool points_unset = true;
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
  return getMeasurements(ag, cf, limits, spacing);
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

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::validateMeshKind() {
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
  case GridDetail::NONBONDED_FIELD:
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
