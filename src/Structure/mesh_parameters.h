// -*-c++-*-
#ifndef STORMM_MESH_PARAMETERS_H
#define STORMM_MESH_PARAMETERS_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"

namespace stormm {
namespace structure {

using constants::CartesianDimension;
using constants::UnitCellAxis;
using data_types::isFloatingPointScalarType;
using data_types::isFloatingPointHpcVectorType;

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

  /// \brief With no const members or pointers, the default copy and move constructors as well as
  ///        copy and move assignment operators are valid.
  /// \{
  MeshParameters(const MeshParameters &original) = default;
  MeshParameters(MeshParameters &&original) = default;
  MeshParameters& operator=(const MeshParameters &other) = default;
  MeshParameters& operator=(MeshParameters &&other) = default;
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

} // namespace structure
} // namespace stormm

#include "mesh_parameters.tpp"

#endif
