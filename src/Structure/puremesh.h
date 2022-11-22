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

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using constants::CartesianDimension;
using constants::UnitCellAxis;
using data_types::isFloatingPointScalarType;
using data_types::isFloatingPointHpcVectorType;
using topology::AtomGraph;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;
using trajectory::PhaseSpace;

/// \brief The default mesh fixed-precision scaling factor is higher than a typical simulation due
///        to the way the mesh occupies a confined region of space.
const int default_mesh_scaling_bits = 40;

/// \brief Encode the critical dimensions of a regular, rectilinear mesh.  The locations of mesh
///        points as well as spacings are stored as fixed-precision integers to ensure consistency
///        and high performance on architectures with deficient 64-bit floating point arithmetic.
class MeshParameters {
public:

  /// \brief The constructor takes formal arguments for all member variables.  Variants support
  ///        triclinic and orthorhombic meshes.
  /// \{
  MeshParameters(int na_in = 1, int nb_in = 1, int nc_in = 1, double origin_x_in = 0.0,
                 double origin_y_in = 0.0, double origin_z_in = 0.0,
                 int scale_bits_in = default_mesh_scaling_bits);

  MeshParameters(int na_in, int nb_in, int nc_in, double origin_x_in, double origin_y_in,
                 double origin_z_in, const std::vector<double> &element_vectors,
                 int scale_bits_in = default_mesh_scaling_bits);

  MeshParameters(int na_in, int nb_in, int nc_in, double origin_x_in, double origin_y_in,
                 double origin_z_in, double element_x, double element_y, double element_z,
                 int scale_bits_in = default_mesh_scaling_bits);
  /// \}

  /// \brief Get the number of points along one of the mesh axes
  ///
  /// Overloaded:
  ///   - Accept a unit cell axis ('a', 'b', or 'c')
  ///   - Accept a Cartesian dimension ('x', 'y', or 'z'), with x ~ a, y ~ b, z ~ c.  This variant
  ///     is specific to orthorhombic unit cells, but will be accepted in general cases for
  ///     convenience.
  ///
  /// \param dim  The dimension of interest
  /// \{
  int getAxisElementCount(UnitCellAxis dim) const;
  int getAxisElementCount(CartesianDimension dim) const;
  /// \}
  
  /// \brief Get the Cartesian origin of the mesh in floating-point numbers.
  ///
  /// Overloaded:
  ///   - Get a three-element vector of all three Cartesian coordinates
  ///   - Get a particular Cartesian coordinate of the origin
  ///
  /// \param dim  The specific Cartesian axis of interest
  /// \{
  template <typename Tcoord> std::vector<Tcoord> getMeshOrigin() const;
  double getMeshOrigin(CartesianDimension dim) const;
  /// \}

  /// \brief Get the Cartesian origin of the mesh as a tuple of floating-point numbers.
  template <typename Tcoord> Tcoord getMeshOriginAsTuple() const;

  /// \brief Get the Cartesian origin of the mesh in fixed-precision numbers.
  ///
  /// Overloaded:
  ///   - Get a three-element vector of all three Cartesian coordinates
  ///   - Get a particular Cartesian coordinate of the origin
  ///
  /// \param dim  The specific Cartesian axis of interest
  /// \{
  std::vector<int95_t> getMeshOriginAsFixedPrecision() const;
  int95_t getMeshOriginAsFixedPrecision(CartesianDimension dim) const;
  /// \}

  /// \brief Get the element vector along one of the unit cell axes.
  ///
  /// Overloaded:
  ///   - Accept a unit cell axis ('a', 'b', or 'c')
  ///   - Accept a Cartesian dimension ('x', 'y', or 'z'), with x ~ a, y ~ b, z ~ c.  This variant
  ///     is specific to orthorhombic unit cells, but will be accepted in general cases for
  ///     convenience.
  ///
  /// \param dim  The axis of interest
  /// \{
  template <typename Tcoord> std::vector<Tcoord> getMeshElementVector(UnitCellAxis dim) const;
  template <typename Tcoord>
  std::vector<Tcoord> getMeshElementVector(CartesianDimension dim) const;
  /// \}

  /// \brief Get the entire element space matrix in any format.  Real formats will have units of
  ///        inverse Angstroms.
  template <typename Tcoord> std::vector<Tcoord> getMeshTransform() const;

  /// \brief Get the inverse element transformation matrix in any format.  Real formats will have
  ///        units of inverse Angstroms.
  template <typename Tcoord> std::vector<Tcoord> getMeshInverseTransform() const;

  /// \brief Get the number of bits after the decimal in this mesh's fixed-precision coordinate
  ///        representations.
  int getScalingBits() const;
  
  /// \brief Get the scaling factor for this mesh's fixed-precision format
  double getScalingFactor() const;

  /// \brief Get the inverse scaling factor for this mesh's fixed-precision format
  double getInverseScalingFactor() const;

  /// \brief Get a vector of fixed-precision format coordinates of the line of grid points starting
  ///        at the origin and proceeding along one of the mesh axes.  One additional point is
  ///        provided to put an upper bound on the final element in whatever dimension.  There are
  ///        nine possible outputs: Cartesian X, Y, or Z coordinates of the mesh's "a", "b", or "c"
  ///        vectors.
  ///
  /// \param mesh_axis  The mesh axis of interest
  /// \param cart_axis  The Cartesian axis of interest
  /// \{
  std::vector<int95_t> getAxisCoordinates(UnitCellAxis mesh_axis,
                                          CartesianDimension cart_axis) const;
  /// \}

private:
  int na;                       ///< Mesh dimension along the unit cell "a" vector
  int nb;                       ///< Mesh dimension along the unit cell "b" vector
  int nc;                       ///< Mesh dimension along the unit cell "c" vector
  int95_t origin_x;             ///< Cartesian X origin of the mesh, expressed in fixed-precision
  int95_t origin_y;             ///< Cartesian Y origin of the mesh, expressed in fixed-precision
  int95_t origin_z;             ///< Cartesian Z origin of the mesh, expressed in fixed-precision
  int scale_bits;               ///< Number of bits after the decimal in positional fixed-precision
                                ///<   representations of grid coordinates and boundaries
  double scale_factor;          ///< Scaling factor to take Cartesian coordinates of grid points
                                ///<   into the fixed-precision representation
  double inverse_scale_factor;  ///< Scaling factor to take Cartesian coordinates of grid points
                                ///<   into the fixed-precision representation

  /// Inverse spacings along all three grid cell vectors, each given by by three consecutive
  /// elements of the array (Fortran order).  The vectors pertain to a single grid element and are
  /// given in units of Angstroms^(-1).  This matrix serves to take forces computed on the regular,
  /// rectilinear element space mesh back into real space.
  double element_umat[9];

  /// Single-precision variant of element_umat.  The fixed-precision representation of coordinates
  /// and mesh points makes it feasible to select the correct element with high fidelity even in
  /// in single-precision arithmetic.
  float sp_element_umat[9];

  /// Mesh spacing along all three grid cell vectors, each given by three consecutive elements of
  /// the array (Fortran order).  The vectors pertain to a single grid element and are given in
  /// units of Angstroms.  This also serves as the inverse transformation matrix for transforming
  /// element-space coordinates back into real space.
  double element_invu[9];

  /// Single-precision variant of element_invu.
  float sp_element_invu[9];
  
  /// The inverse element transformation matrix, represented in fixed precision.  The number of
  /// bits after the decimal in this representation should match that of the fixed precision
  /// coordinate representation.  Representing the grid origin and grid points in this manner
  /// ensures high-precision computations of the relative particle and mesh positions.
  int95_t fp_element_invu[9];
};

/// \brief A workspace for constructing a pure potential mesh based on the frozen atoms of a
///        large molecule.  If the large molecule has nonrigid components, they must be excluded
///        from contributing to the grid.  In addition, any atoms up to 1:4 (connected by three
///        bonds or less) must also be excluded from the grid-based potential.  Computations on
///        these atoms will not be accurate off the grid, but since they are frozen the
///        consequences are mitigated.
template <typename T> class PureMesh {
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
  PureMesh(const AtomGraph *ag, const CoordinateFrame &cf, const MeshParameters &measurements_in,
           const GpuDetails &gpu = null_gpu);

  PureMesh(const AtomGraph *ag, const CoordinateFrame &cf, int scale_bits, double buffer,
           double spacing, const GpuDetails &gpu = null_gpu);

  PureMesh(const AtomGraph *ag, const CoordinateFrame &cf, int scale_bits, double buffer,
           const std::vector<double> &spacing, const GpuDetails &gpu = null_gpu);

  PureMesh(const AtomGraph *ag, const CoordinateFrame &cf, int scale_bits,
           const std::vector<double> &mesh_bounds, double spacing,
           const GpuDetails &gpu = null_gpu);

  PureMesh(const AtomGraph *ag, const CoordinateFrame &cf, int scale_bits,
           const std::vector<double> &mesh_bounds, const std::vector<double> &spacing,
           const GpuDetails &gpu = null_gpu);
  /// \}
  
private:

  /// Measurements for the mesh, including various transformation matrices and fixed-precision
  /// representations.
  MeshParameters measurements;

  /// Fixed-precision Cartesian coordinates of the mesh grid points stepping along the lines
  /// [ i = 0...nx, 0, 0 ] (the "a" box vector), [ 0, j = 0...ny, 0 ] (the "b" box vector), and
  /// [ 0, 0, k = 0...nz ] (the "c" box vector).  Storing these values obviates the need to do
  /// expensive and complicated multiplications of fixed-precision numbers.
  /// \{
  Hybrid<llint> fp_a_line_x;
  Hybrid<llint> fp_a_line_y;
  Hybrid<llint> fp_a_line_z;
  Hybrid<llint> fp_b_line_x;
  Hybrid<llint> fp_b_line_y;
  Hybrid<llint> fp_b_line_z;
  Hybrid<llint> fp_c_line_x;
  Hybrid<llint> fp_c_line_y;
  Hybrid<llint> fp_c_line_z;
  /// \}

  /// Overflow for fixed-precision Cartesian coordinates of the mesh grid points stepping along
  /// the a, b, and c box vectors.
  /// \{
  Hybrid<int> fp_a_line_overflow_x;
  Hybrid<int> fp_a_line_overflow_y;
  Hybrid<int> fp_a_line_overflow_z;
  Hybrid<int> fp_b_line_overflow_x;
  Hybrid<int> fp_b_line_overflow_y;
  Hybrid<int> fp_b_line_overflow_z;
  Hybrid<int> fp_c_line_overflow_x;
  Hybrid<int> fp_c_line_overflow_y;
  Hybrid<int> fp_c_line_overflow_z;
  /// \}

  /// Coefficients for tricubic spline functions spanning each grid element.  Coefficients for
  /// element {Ea, Eb, Ec} traversing the "a", "b", and "c" axes, respectively, are held
  /// consecutively in a 64-element series with order {i, j, k} = (0, 0, 0), (0, 0, 1), ...,
  /// (0, 0, 3), (1, 0, 0), (1, 0, 1), ..., (1, 0, 3), ..., (3, 3, 3).  The 64-element series for
  /// element {Ea, Eb, Ec} will be found at 64 * ((((Ea * nb) + Eb) * nc) + Ec), where {na, nb, nc}
  /// are the mesh dimensions along each axis.
  Hybrid<T> coefficients;

  /// Pointer to the original topology
  const AtomGraph *ag_pointer;

  /// Snapshot of the frozen atoms in the topology (in case the topology is modified later)
  Hybrid<uint> frozen_atoms;

  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> int_data;

  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> llint_data;
};

} // namespace structure
} // namespace stormm

#endif
