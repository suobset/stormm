// -*-c++-*-
#ifndef STORMM_TRICUBIC_CELL_H
#define STORMM_TRICUBIC_CELL_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "matrix_ops.h"

namespace stormm {
namespace math {

using card::Hybrid;
using constants::CartesianDimension;

/// \brief A list of the different boundary components that determine a tricubic spline, covering
///        the values and all derivatives at the boundaries of the grid element.
enum class TricubicBound {
  VALUE,   ///< Value of the function at a particular grid point
  DX,      ///< Derivative of the function along the Cartesian X direction
  DY,      ///< Derivative of the function along the Cartesian Y direction
  DZ,      ///< Derivative of the function along the Cartesian Z direction
  DXY,     ///< Cross-derivative of the function along the Cartesian X and Y directions
  DXZ,     ///< Cross-derivative of the function along the Cartesian X and Z directions
  DYZ,     ///< Cross-derivative of the function along the Cartesian Y and Z directions
  DXYZ     ///< Cross-derivative of the function along all three Cartesian directions
};

/// \brief Store the coefficients and bounding data for one regular grid cell of a
///        three-dimensional function to be interpolated by a tricubic spline.
template <typename T> class TricubicCell {
public:

  /// \brief The constructor can take nothing and simply initialize all values to zero, or accept
  ///        the tricubic weights matrix, a dimensions array, and details of the potential
  ///        function.
  ///
  /// \param weights_matrix  The inverse matrix of polynomial weights obtained from
  ///                        getTricubicMatrix() in this library
  /// \param bounds          Array containing the Cartesian x, y, and z coordinates of the grid
  ///                        cell origin plus its lengths along each axis (the grid cell is assumed
  ///                        to be orthorhombic).  This array can have four entries (in which case
  ///                        the final entry is assumed to be the isotropic length parameter) or
  ///                        six (for anisotropic cells).
  /// \{
  TricubicCell();

  TricubicCell(const std::vector<double> weights_matrix, const std::vector<double> &bounds,
               const std::vector<T> &f_in,
               const std::vector<T> &dx_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dy_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dz_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0),
                                               static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dxy_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dxz_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dyz_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0),
                                                static_cast<T>(0.0), static_cast<T>(0.0) },
               const std::vector<T> &dxyz_in = { static_cast<T>(0.0), static_cast<T>(0.0),
                                                 static_cast<T>(0.0), static_cast<T>(0.0),
                                                 static_cast<T>(0.0), static_cast<T>(0.0),
                                                 static_cast<T>(0.0), static_cast<T>(0.0) });
  /// \}

  /// \brief Retrieve one of the 64 coefficients Aijk for the tricubic spline.
  ///
  /// \param i  The ith coefficient relevant to the notation Aijk
  /// \param j  The jth coefficient relevant to the notation Aijk
  /// \param k  The kth coefficient relevant to the notation Aijk
  T getCoefficient(int i, int j, int k) const;

  /// \brief Set one of the 64 coefficients Aijk for the tricubic spline.  Parameter descriptions
  ///        follow from above, with the addition of:
  ///
  /// \param value  The value to apply
  void setCoefficient(T value, int i, int j, int k);

  /// \brief Get one of the data points from the boundary.  Parameter descriptions follow from
  ///        above, with the addition of:
  ///
  /// \param kind  The classification of the boundary condition as a derivative (or value)
  T getData(TricubicBound kind, int i, int j, int k) const;

  /// \brief Set one of the data items.  Parameter descriptions follow from above.
  void setData(T value, TricubicBound kind, int i, int j, int k);

  /// \brief Get the cell origin along one dimension.
  ///
  /// \param dim  The Cartesian dimension along which to return the origin coordinate
  T getCellOrigin(CartesianDimension dim) const;

  /// \brief Get the cell length along one dimension.
  ///
  /// \param dim  The Cartesian dimension along which to return the cell length
  T getCellLength(CartesianDimension dim) const;

  /// \brief Evaluate the function at a specific point in space.  This will take into account the
  ///        grid cell's origin and lengths to determine where in the grid cell the point of
  ///        interest lies.  If the point is outside the grid cell, produce an error.
  ///
  /// \param point_x  Cartesian X location of the point
  /// \param point_y  Cartesian Y location of the point
  /// \param point_z  Cartesian Z location of the point
  T evaluate(T x, T y, T z) const;
  
private:
  T coefficients[64];  ///< Solved coefficients of the tricubic spline that satisfies all boundary
                       ///<   conditions.

  // The following arrays store their series of values in "Fortran" order: (X0, Y0, Z0),
  // (X1, Y0, Z0), (X0, Y1, Z0), (X1, Y1, Z0), (X0, Y0, Z1), ...
  T f[8];      ///< Values of the function at the bounding grid points
  T dx[8];     ///< Cartesian X derivatives of the function at the bounding grid points
  T dy[8];     ///< Cartesian Y derivatives of the function at the bounding grid points
  T dz[8];     ///< Cartesian Z derivatives of the function at the bounding grid points
  T dxy[8];    ///< Cartesian X/Y cross-derivatives of the function at the bounding grid points
  T dxz[8];    ///< Cartesian X/Z cross-derivatives of the function at the bounding grid points
  T dyz[8];    ///< Cartesian Y/Z cross-derivatives of the function at the bounding grid points
  T dxyz[8];   ///< Cartesian X/Y/Z triple-derivatives of the function at the bounding grid points

  // The grid cell boundaries are stored in double precision for accuracy considerations.
  double origin_x;  ///< Cartesian X location of the grid cell origin
  double origin_y;  ///< Cartesian Y location of the grid cell origin
  double origin_z;  ///< Cartesian Z location of the grid cell origin
  double length_x;  ///< Cartesian X length of the grid cell
  double length_y;  ///< Cartesian Y length of the grid cell
  double length_z;  ///< Cartesian Z length of the grid cell
};

/// \brief Construct a matrix for grinding out tricubic spline coefficients given the eight vectors
///        found in the TricubicCell object.  The matrix is constructed in double-precision format
///        to ensure that all subsequent calculations are reasonably accurate, at least with
///        respect to this step.  Once the coefficients have been pre-computed they can be stored
///        and used in lower-precision formats.
Hybrid<double> getTricubicMatrix();

} // namespace math
} // namespace stormm

#include "tricubic_cell.tpp"

#endif
