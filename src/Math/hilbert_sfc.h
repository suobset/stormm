// -*-c++-*-
#ifndef STORMM_HILBERT_SFC_H
#define STORMM_HILBERT_SFC_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

using constants::CartesianDimension;
using constants::UnitCellAxis;
  
/// \brief A basic Hilbert cube which begins at the origin and ends on the Y axis.  Rotate to
///        obtain basic cubes that end on the Y or Z axes.
const int3 hilbert_cube_hll[] = { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 },
                                  { 1, 0, 1 }, { 1, 1, 1 }, { 1, 1, 0 }, { 1, 0, 0 } };

/// \brief A basic Hilbert cube which begins at the origin and ends at the far corner (displaced by
///        one along X, Y, and Z).
const int3 hilbert_cube_hhh[] = { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 },
                                  { 1, 0, 1 }, { 0, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 } };

/// \brief A basic Hilbert Square which lies in the XY plane, begins at the origin, and ends on the
///        X axis.  Rotate to obtain basic squres that lie in the XZ or YZ planes, begin at the
///        origin, and end on the Y or Z axes.
const int3 hilbert_square_hl[] = { { 0, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 1, 0, 0 } };
const int3 peano_square[] = { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 2, 0 },
                              { 1, 2, 0 }, { 1, 1, 0 }, { 1, 0, 0 },
                              { 2, 0, 0 }, { 2, 1, 0 }, { 2, 2, 0 } };
const int3 hileano_square[] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 },
                                { 0, 1, 0 }, { 0, 2, 0 }, { 1, 2, 0 },
                                { 2, 2, 0 }, { 2, 1, 0 }, { 2, 0, 0 } };
    
/// \}
  
/// \brief Draw a Hilbert space-filling curve.
class HilbertSFC {
public:

  /// \brief The constructor takes as inputs a trio of integers indicating the overall unitless
  ///        dimensions of the grid into which to inscribe the space-filling curve, as well as a
  ///        choice of methods for spanning the region of space. 
  HilbertSFC(int x_span_in, int y_span_in, int z_span_in,
             HilbertCurveMode method = HilbertCurveMode::OVERSPAN);

  /// \brief Get the dimensions of the region filled by the curve
  int3 getDimensions() const;

  /// \brief Get the dimension of the region filled by the curve along one axis.
  ///
  /// Overloaded:
  ///   - Provide the Cartesian axis
  ///   - Provide the unit cell axis
  ///
  /// \param dim  The dimension of interest
  /// \{
  int getDimension(CartesianDimension dim) const;

  int getDimension(UnitCellAxis dim) const;
  /// \}

  /// \brief Get the curve index to which a particular grid element corresponds.  The grid begins
  ///        at (0, 0, 0).  This function returns an error if the grid element is outside the
  ///        curve's valid region.
  ///
  /// \param grid_x  Unitless X position on the region filled by the curve
  /// \param grid_y  Unitless Y position on the region filled by the curve
  /// \param grid_z  Unitless Z position on the region filled by the curve
  int getCurveIndex(int grid_x, int grid_y, int grid_z) const;

  /// \brief Get the Hilbert space filling curve, expressed in unitless integers with components
  ///        for each Cartesian dimension.
  const std::vector<int3>& getCurve() const;
  
private:

  int x_span;                      ///< Span of the region of interest along the Cartesian X axis
  int y_span;                      ///< Span of the region of interest along the Cartesian Y axis
  int z_span;                      ///< Span of the region of interest along the Cartesian Z axis
  std::vector<int3> curve;         ///< The space-filling curve, with unitless X, Y, and Z
                                   ///<   coordinates given in the "x", "y", and "z" members of
                                   ///<   each tuple.
  std::vector<int> curve_indices;  ///< Indices to which grid points in the applicable region
                                   ///<   correspond on the list of points in the space-filling
                                   ///<   curve.  The space-filling curve visits Cartesian point
                                   ///<   (i, j, k) at the node indicated by
                                   ///<   curve_indicies[(((k * y_span) + j) * x_span) + i].
};

/// \brief Recursive function for propagating a Hilbert space-filling curve
///
/// \param s      The current span to consider (when recursion takes this to one, the array curve
///               is lengthened
/// \param curve  The developing curve
/// \param x      Current X location of the "cursor"
/// \param y      Current Y location of the "cursor"
/// \param z      Current Z location of the "cursor"
/// \param dx_a   Shift in the X cursor location based on movement along the X axis
/// \param dy_a   Shift in the Y cursor location based on movement along the X axis
/// \param dz_a   Shift in the Z cursor location based on movement along the X axis
/// \param dx_b   Shift in the X cursor location based on movement along the Y axis
/// \param dy_b   Shift in the Y cursor location based on movement along the Y axis
/// \param dz_b   Shift in the Z cursor location based on movement along the Y axis
/// \param dx_c   Shift in the X cursor location based on movement along the Z axis
/// \param dy_c   Shift in the Y cursor location based on movement along the Z axis
/// \param dz_c   Shift in the Z cursor location based on movement along the Z axis
void propagateHilbertCurve(int s, std::vector<int3> *curve, int x = 0, int y = 0, int z = 0,
                           int dx = 1, int dy = 0, int dz = 0, int dx2 = 0, int dy2 = 1,
                           int dz2 = 0, int dx3 = 0, int dy3 = 0, int dz3 = 1);

/// \brief Draw a three-dimensional Hilbert space-filling curve of the specified order.
///
/// \param order  The order of the curve, which will then cover a space 2^order indices in three
///               directions.  The maximum allowed value is 10, minimum 0.
std::vector<int3> drawHilbertSpaceCurve(int order);
  
  
} // namespace stmath
} // namespace stormm

#endif
