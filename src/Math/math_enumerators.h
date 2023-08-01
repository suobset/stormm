// -*-c++-*-
#ifndef STORMM_MATH_ENUMERATORS_H
#define STORMM_MATH_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace stmath {

/// \brief A list of the different boundary components that determine a tricubic spline, covering
///        the values and all derivatives at the boundaries of the grid element.
enum class FunctionLevel {
  VALUE,   ///< Value of the function at a particular grid point
  DX,      ///< Derivative of the function along the Cartesian X direction
  DY,      ///< Derivative of the function along the Cartesian Y direction
  DZ,      ///< Derivative of the function along the Cartesian Z direction
  DXX,     ///< Double derivative of the function along the Cartesian X direction
  DXY,     ///< Cross-derivative of the function along the Cartesian X and Y directions
  DXZ,     ///< Cross-derivative of the function along the Cartesian X and Z directions
  DYY,     ///< Double derivative of the function along the Cartesian Y direction
  DYZ,     ///< Cross-derivative of the function along the Cartesian Y and Z directions
  DZZ,     ///< Double derivative of the function along the Cartesian Z direction
  DXXX,    ///< Triple derivative of the function along the Cartesian X direction
  DXXY,    ///< Triple partial derivative of the function along the Cartesian X and Y directions
  DXXZ,    ///< Triple partial derivative of the function along the Cartesian X and Z directions
  DXYY,    ///< Triple partial derivative of the function along the Cartesian X and Y directions
  DXZZ,    ///< Triple partial derivative of the function along the Cartesian X and Z directions
  DXYZ,    ///< Triple partial derivative of the function along all three Cartesian directions
  DYYY,    ///< Triple derivative of the function along the Cartesian Y direction
  DYYZ,    ///< Triple partial derivative of the function along the Cartesian Y and Z directions
  DYZZ,    ///< Triple partial derivative of the function along the Cartesian Y and Z directions
  DZZZ     ///< Triple derivative of the function along the Cartesian Z direction
};

/// \brief For tricubic interpolation, there are many choices as to the 64 pieces of information
///        that will define the polynomial coefficients for any particular mesh element.  The
///        traditional approach taken by Francois Lekien and Jerrold Marsden favors maximal
///        continuity at the corners of the interpolated region of space, taking the values of the
///        function f(x,y,z), its first derivatives df/dx, df/dy, and df/dz, and its mixed partial
///        second and third derivatives d2f/dxy, d2f/dxz, d2f/dyz, and d3f/dxyz.  Despite employing
///        third derivatives, the resulting interpolated function is only C1 continuous--the first
///        derivative alone is guaranteed to be smooth throughout space.
///
///        As Lekien and Marsden showed, fitting local polynomial splines to the function values
///        and first derivatives at each mesh vertex is sufficient to guarantee C1 continuity
///        across mesh element boundaries.  Fitting to higher-order mixed partial derivatives
///        grants only partial C2 and C3 continuity.  d2f/dx2, d2f/dy2, and d2f/dz2, as well as
///        various third derivatives, are not guaranteed to be continuous if there is any degree of
///        approximation in the interpolant.  Furthermore, imposing these restrictions requires the
///        interpolated function to have C2 and C3 continuity of its own.  Finally, results suggest
///        the numerical accuracy of this interpolant is only strong in orthorhombic meshes.
///        Outside of this case, the higher-order derivatives appear to lead to significant
///        over-corrections in the interiors of mesh elements.
///
///        To avoid the limitations and improve the accuracy of tricubic interpolation in the
///        interiors of mesh elements, a second style of interpolant is provided which evaluates
///        additional values of the function f(x,y,z) at selected points in the interior of each
///        mesh element.
/// \{
enum class Interpolant {
  SMOOTHNESS,     ///< Evaluate f(x,y,z), its first derivatives, and various mixed partial
                  ///<   derivatives at mesh corners until sufficient criteria are found to fit the
                  ///<   necessary coefficients for the interpolant.
  FUNCTION_VALUE  ///< Evaluate f(x,y,z) and its first derivatives at each mesh vertex, then make
                  ///<   additional evaluations of f(x,y,z) at points in the interior of the
                  ///<   element to get the most accurate reproduction of f(x,y,z) throughout
                  ///<   space.  For tricubic interpolation, these additional f(x,y,z) evaluations
                  ///<   begin once C1 continuity is guaranteed.
};
/// \}

/// \brief Enumerate the ways to approach a limit.
enum class LimitApproach {
  BELOW,  ///< Approach the limit from below, with the function argument ascending
  ABOVE   ///< Approach the limit from above, with the function argument descending
};
 
/// \brief Get a human-readable string describing an enumeration of the provided type.  Various
///        overloads of this function serve enumerators across many libraries.
///
/// \param input  The enumeration to translate
/// \{
std::string getEnumerationName(FunctionLevel input);
std::string getEnumerationName(Interpolant input);
std::string getEnumerationName(LimitApproach input);
/// \}

} // namespace stmath
} // namespace stormm

#endif
