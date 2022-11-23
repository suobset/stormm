// -*-c++-*-
#ifndef STORMM_PUREMESH_H
#define STORMM_PUREMESH_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "mesh_parameters.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using constants::CartesianDimension;
using constants::UnitCellAxis;
using topology::AtomGraph;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;
using trajectory::PhaseSpace;

/// \brief The templated, writeable abstract of a BackgroundMesh object.

  
/// \brief A workspace for constructing a pure potential mesh based on the frozen atoms of a
///        large molecule.  If the large molecule has nonrigid components, they must be excluded
///        from contributing to the grid.  In addition, any atoms up to 1:4 (connected by three
///        bonds or less) must also be excluded from the grid-based potential.  Computations on
///        these atoms will not be accurate off the grid, but since they are frozen the
///        consequences are mitigated.
template <typename T> class BackgroundMesh {
public:

  /// \brief The constructor takes all dimension parameters plus an indication of what type of
  ///        potential, the molecular system, and what mask of atoms is to be mapped.  Variants
  ///        include different ways to define the limits of the mesh.  If a GPU is available, it
  ///        will be used to compute the mesh.
  ///
  /// \param ag               System topology (markings in its mobile_atoms array will be used to
  ///                         determine which atoms to map)
  /// \param cf               Cartesian coordinates of all particles
  /// \param scale_bits       Number of bits after the decimal when taking coordinates in units of
  ///                         Angstroms into the fixed-precision format of the mesh framework
  /// \param buffer           Breadth around the molecule drawing an orthhombic region in which to
  ///                         map the mesh
  /// \param mesh_bounds      Boundaries of the mesh, a six-element vector describing the lower
  ///                         and upper Cartesian X, Y, and Z limits of an orthorhombic region
  ///                         in which to define a rectilinear mesh.
  /// \param spacing          Grid spacings for the mesh cells.  Provide either a single number or
  ///                         three dimensions for the Cartesian X, Y, and Z widths of rectilinear
  ///                         elements.
  /// \param measurements_in  A full description of the mesh parameters.  This provides a means to
  ///                         define non-orthorhombic meshes.
  /// \param gpu              Details of the available GPU
  /// \{
  BackgroundMesh(GridDetail kind_in = GridDetail::NONBONDED_FIELD,
                 const AtomGraph *ag_in = nullptr,
                 const MeshParameters &measurements_in = MeshParameters());
  
  BackgroundMesh(GridDetail kind_in, const AtomGraph *ag_in, const CoordinateFrame &cf,
                 const MeshParameters &measurements_in, const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, const AtomGraph *ag_in, const CoordinateFrame &cf,
                 int scale_bits, double buffer, double spacing, const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, const AtomGraph *ag_in, const CoordinateFrame &cf,
                 int scale_bits, double buffer, const std::vector<double> &spacing,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, const AtomGraph *ag_in, const CoordinateFrame &cf,
                 int scale_bits, const std::vector<double> &mesh_bounds, double spacing,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, const AtomGraph *ag_in, const CoordinateFrame &cf,
                 int scale_bits, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing, const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Copy and move constructors must be defined explicitly to accommodate POINTER-kind
  ///        Hybrid objects.
  ///
  /// \param original  The original object to copy or move
  /// \{
  BackgroundMesh(const BackgroundMesh<T> &original);
  BackgroundMesh(BackgroundMesh<T> &&original);
  /// \}

  /// \brief Copy and move assignment operators must be defined explicitly to accommodate
  ///        POINTER-kind Hybrid objects.
  ///
  /// \param other  Assignments are made based on this pre-existing object.
  /// \{
  BackgroundMesh& operator=(const BackgroundMesh<T> &original);
  BackgroundMesh& operator=(BackgroundMesh<T> &&original);
  /// \}
  
  /// \brief Get a const pointer to the topology responsible for creating this mesh.
  const AtomGraph* getTopologyPointer();
  
private:

  /// Measurements for the mesh, including various transformation matrices and fixed-precision
  /// representations.
  MeshParameters measurements;

  /// The detail level of the mesh--this carries implications for the allowed data types.  An
  /// "OCCLUSION" mesh must be of type ullint (64-bit unsigned integer, to encode 4 x 4 x 4
  /// cubelets of bitwise information for "occluded" or "non-occluded." Such meshes can be very
  /// dense, as every mesh element constitutes 4 x 4 x 4 cubelets, or a total of 512 bytes of
  /// information, as much as 20k mesh points on each side.  In contrast, a "NONBONDED_FIELD" type
  /// grid must have data type float or double to encode the appropriate non-bonded potential or
  /// influence field (each Lennard-Jones type and Generalized Born parameter combination will get
  /// its own field).  A "NONBONDED_ATOMIC" type mesh stores field data in floating-point scalar
  /// data as well as local neighbor lists in an additional array of integer data, with a bounds
  /// array.
  GridDetail kind;

  /// The type of non-bonded potential to map to the mesh.  This is irrelevant for "OCCLUSION"-kind
  /// meshes (the effective potential is a clash based on the Lennard-Jones sigma radii and a probe
  /// width).
  NonbondedPotential field;
  
  /// Fixed-precision Cartesian coordinates of the mesh grid points stepping along the lines
  /// [ i = 0...nx, 0, 0 ] (the "a" box vector), [ 0, j = 0...ny, 0 ] (the "b" box vector), and
  /// [ 0, 0, k = 0...nz ] (the "c" box vector).  Storing these values obviates the need to do
  /// expensive and complicated multiplications of fixed-precision numbers.
  /// \{
  Hybrid<llint> a_line_x;
  Hybrid<llint> a_line_y;
  Hybrid<llint> a_line_z;
  Hybrid<llint> b_line_x;
  Hybrid<llint> b_line_y;
  Hybrid<llint> b_line_z;
  Hybrid<llint> c_line_x;
  Hybrid<llint> c_line_y;
  Hybrid<llint> c_line_z;
  /// \}

  /// Overflow for fixed-precision Cartesian coordinates of the mesh grid points stepping along
  /// the a, b, and c box vectors.
  /// \{
  Hybrid<int> a_line_x_overflow;
  Hybrid<int> a_line_y_overflow;
  Hybrid<int> a_line_z_overflow;
  Hybrid<int> b_line_x_overflow;
  Hybrid<int> b_line_y_overflow;
  Hybrid<int> b_line_z_overflow;
  Hybrid<int> c_line_x_overflow;
  Hybrid<int> c_line_y_overflow;
  Hybrid<int> c_line_z_overflow;
  /// \}

  /// Coefficients for tricubic spline functions spanning each grid element.  Coefficients for
  /// element {Ea, Eb, Ec} traversing the "a", "b", and "c" axes, respectively, are held
  /// consecutively in a 64-element series with order {i, j, k} = (0, 0, 0), (0, 0, 1), ...,
  /// (0, 0, 3), (1, 0, 0), (1, 0, 1), ..., (1, 0, 3), ..., (3, 3, 3).  The 64-element series for
  /// element {Ea, Eb, Ec} will be found at 64 * ((((Ea * nb) + Eb) * nc) + Ec), where {na, nb, nc}
  /// are the mesh dimensions along each axis.
  Hybrid<T> coefficients;

  /// Pointer to the original topology
  AtomGraph *ag_pointer;

  /// Snapshot of the frozen atoms in the topology (in case the topology is modified later)
  Hybrid<uint> frozen_atoms;

  /// Element-wise neighbor listts for nearby frozen atoms
  Hybrid<int> neighbor_list;

  /// Bounds array for the neighbor list atoms array above
  Hybrid<int> neighbor_list_bounds;
  
  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> int_data;

  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> llint_data;

  /// \brief Allocate memory for the mesh coefficients, axis coordinates, and frozen atoms mask.
  void allocate();

  /// \brief Repair the POINTER-kind Hybrid objects in a newly copied or moved object.
  void rebase_pointers();
};

} // namespace structure
} // namespace stormm

#include "background_mesh.tpp"

#endif
