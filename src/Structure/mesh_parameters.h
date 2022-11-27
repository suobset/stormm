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
using data_types::getStormmScalarTypeName;
using data_types::isFloatingPointScalarType;
using data_types::isFloatingPointHpcVectorType;

/// \brief The default mesh fixed-precision scaling factor is higher than a typical simulation due
///        to the way the mesh occupies a confined region of space.
constexpr int default_mesh_scaling_bits = 40;

/// \brief The maximum number of fixed-precision bits that can be used in certain situations when
///        particle position overflow bits will be assumed to not be in use.
constexpr int mesh_nonoverflow_bits = 48;

/// \brief The abstract of a MeshParameters object is read-only (modify the original to get a new
///        abstract if the dimensions change), but templated to prune the information present.
template <typename T> struct MeshParamAbstract {

  /// \brief The constructor takes arguments for all member variables.
  MeshParamAbstract(int na_in, int nb_in, int nc_in, int95_t orig_x_in, int95_t orig_y_in,
                    int95_t orig_z_in, T scale_in, T inv_scale_in, const T* umat_in,
                    const T* invu_in, const T* widths_in, const int95_t* fp_invu_in);

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  MeshParamAbstract(const MeshParamAbstract &original) = default;
  MeshParamAbstract(MeshParamAbstract &&original) = default;
  /// \}
  
  const int na;              ///< Number of mesh elements along the "a" (~x) axis
  const int nb;              ///< Number of mesh elements along the "b" (~y) axis
  const int nc;              ///< Number of mesh elements along the "c" (~z) axis
  const int95_t orig_x;      ///< Cartesian X origin of the mesh in fixed-precision format
  const int95_t orig_y;      ///< Cartesian Y origin of the mesh in fixed-precision format
  const int95_t orig_z;      ///< Cartesian Z origin of the mesh in fixed-precision format
  const T scale;             ///< Scaling factor for taking coordinates with units of Angstroms
                             ///<   into the fixed-percision format of the mesh
  const T inv_scale;         ///< Inverse scaling factor for taking fixed-precision coordinates
                             ///<   back into STORMM's internal units of Angstroms
  const T umat[9];           ///< Transformation matrix to take coordinates into element space
  const T invu[9];           ///< Inverse transformation matrix for each element, or the vectors
                             ///<   defining each side of the element.
  const T widths[3];         ///< Widths of the mesh element between faces defining vectors roughly
                             ///<   associated with the Cartesian X, Y, and Z axes (of these, the
                             ///<   X and Y axes do not exactly line up, but the normals between
                             ///<   faces defined by the element's "b" and "c" vectors and its "a"
                             ///<   and "c" vectors are close).
  const int95_t fp_invu[9];  ///< Fixed-precision inverse transformation matrix for each element
};
  
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
  template <typename T3> T3 getMeshOriginAsTuple() const;

  /// \brief Get the Cartesian origin of the mesh in fixed-precision numbers.
  ///
  /// Overloaded:
  ///   - Get a three-element vector of all three Cartesian coordinates
  ///   - Get a particular Cartesian coordinate of the origin
  ///
  /// \param dim  The specific Cartesian axis of interest
  /// \{
  std::vector<int95_t> getMeshOriginAsFP() const;
  int95_t getMeshOriginAsFP(CartesianDimension dim) const;
  /// \}

  /// \brief Get the element vector along one of the unit cell axes in floating-point numbers.
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

  /// \brief Get the element vector along one of the unit cell axes as a tuple of floating-point
  ///        numbers.
  ///
  /// Overloaded:
  ///   - Accept a unit cell axis ('a', 'b', or 'c')
  ///   - Accept a Cartesian dimension ('x', 'y', or 'z'), with x ~ a, y ~ b, z ~ c.  This variant
  ///     is specific to orthorhombic unit cells, but will be accepted in general cases for
  ///     convenience.
  ///
  /// \param dim  The axis of interest
  /// \{
  template <typename T3> T3 getMeshElementVectorAsTuple(UnitCellAxis dim) const;
  template <typename T3> T3 getMeshElementVectorAsTuple(CartesianDimension dim) const;
  /// \}
  
  /// \brief Get the element vector along one of the unit cell axes in fixed precision.
  ///
  /// Overloaded:
  ///   - Accept a unit cell axis ('a', 'b', or 'c')
  ///   - Accept a Cartesian dimension ('x', 'y', or 'z'), with x ~ a, y ~ b, z ~ c.  This variant
  ///     is specific to orthorhombic unit cells, but will be accepted in general cases for
  ///     convenience.
  ///
  /// \param dim  The axis of interest
  /// \{
  std::vector<int95_t> getMeshElementVectorAsFP(UnitCellAxis dim) const;
  std::vector<int95_t> getMeshElementVectorAsFP(CartesianDimension dim) const;
  /// \}

  /// \brief Get the entire element space matrix in any format.  Real formats will have units of
  ///        inverse Angstroms.
  template <typename Tcoord> std::vector<Tcoord> getMeshTransform() const;

  /// \brief Get the inverse element transformation matrix in real-number format, with units of
  ///        inverse Angstroms.
  template <typename Tcoord> std::vector<Tcoord> getMeshInverseTransform() const;

  /// \brief Get the inverse element transformation matrix in (authoritative) fixed-precision
  ///        format.
  std::vector<int95_t> getMeshInverseTransformAsFP() const;

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

  /// \brief Obtain a double-precision abstract for this object.
  MeshParamAbstract<double> dpData() const;

  /// \brief Obtain a single-precision abstract for this object.
  MeshParamAbstract<float> spData() const;

  /// \brief Set the number of mesh elements along the "a", "b", or "c" axes.
  ///
  /// Overloaded::
  ///   - Set one dimension or all three
  ///
  /// \param n_in       The dimension or dimensions to set
  /// \param mesh_axis  The mesh axis to define
  /// \{
  void setMeshDimension(int n_in, UnitCellAxis mesh_axis);
  void setMeshDimension(const std::vector<int> &n_in);
  /// \}
  
  /// \brief Set the origin's Cartesian X, Y, or Z coordinates.
  ///
  /// Overloaded::
  ///   - Provide a real-valued number in Angstroms
  ///   - Provide a fixed-precision value
  ///   - Set one coordinate or all three
  ///
  /// \param v          The origin coordinate
  /// \param cart_axis  The axis to set the origin along
  /// \{
  void setOrigin(double v, CartesianDimension cart_axis);
  void setOrigin(int95_t v, CartesianDimension cart_axis);
  void setOrigin(const std::vector<double> &v);
  void setOrigin(const std::vector<int95_t> &v);
  /// \}

  /// \brief Set the scaling bits.  This will also update the scaling factors.
  ///
  /// \param scale_bits_in  The number of bis after the decimal
  void setScalingBits(int scale_bits_in);

  /// \brief Define the basic element coordinates using three vectors in three-dimensional space.
  void defineElement(const std::vector<double> &element_vectors);
  
  /// \brief Validate the choice of mesh dimensions.
  void validateMeshDimensions() const;

  /// \brief Validate the mesh's fixed-precision representation.
  void validateFixedPrecisionBits() const;
  
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

  /// Single-precision variant of element_invu
  float sp_element_invu[9];

  /// Widths of the mesh element between faces defining vectors roughly associated with the
  /// Cartesian X, Y, and Z axes (of these, the X and Y axes do not exactly line up, but the
  /// normals between faces defined by the element's "b" and "c" vectors and its "a" and "c"
  /// vectors are close).
  double widths[3];

  /// Single-precision variant of the element widths
  float sp_widths[3];

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
