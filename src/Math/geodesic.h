// -*-c++-*-
#ifndef STORMM_GEODESIC_H
#define STORMM_GEODESIC_H

#include <vector>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "math_enumerators.h"
#include "matrix_ops.h"

namespace stormm {
namespace stmath {

using symbols::pi;
using symbols::twopi;
using symbols::tetrahedral_angle;

/// \brief Collect the triangles by subdividing polygons of a defined order found on the sphere
///        surface, then add new points to the sphere surface at the midpoint of each triangle.
///
/// \param midpt_llim  The lower limit of polygon midpoints, each of which will form one vertex of
///                    the nring triangles in the polygon
/// \param midpt_hlim  The upper limit of polygon midpoints
/// \param ring_llim   The lower limit of points on the sphere surface that can be used to define
///                    the polygons for forming triangles
/// \param ring_hlim   The upper limit of points on the sphere surface that can be used to define
///                    the polygons for forming triangles
/// \param nring       The order of the polygon from which to form triangles and then place new
///                    points
/// \param crit_rsq    The (squared) critical distance between points that can be used to identify
///                    whether they are consecutive when stepping around the polygon
/// \param result      The developing array of new points added at the midpoints of any triangles
///                    found.  Pre-allocated, modified by the function and returned.
/// \param nadd_pt     The number of points added thus far, an indicator of progress through the
///                    result array over successive calls to this function.  An update to nadd_pt
///                    will be the output of the function.
/// \param scaffold    Locations of points with which to form new triangles, the current collection
///                    of sphere surface points
template <typename T, typename T3>
int gatherTriangles(int midpt_llim, int midpt_hlim, int ring_llim, int ring_hlim, int nring,
                    T crit_rsq, std::vector<T3> *result, int nadd_pt,
                    const std::vector<double3> &scaffold);

/// \brief Produce a series of points on a unit sphere distributed as close to equidistant from one
///        another as possible.
///
/// Overloaded:
///   - Provide a mutable integer indicating the number of points which can be modified to
///     adapt to the actual size of the result
///   - Provide an immutable integer and take the actual number of returned points from the size
///     of the resulting vector
///
/// \param n          The requested number of points to distribute over the sphere surface.  This
///                   may be modified, specifically by the Deserno algorithm, if such a number of
///                   points could not be placed at regular intervals across the surface.
/// \{
template <typename T, typename T3>
std::vector<T3> surfaceDistribution(int *n, SpherePlacement method = SpherePlacement::POLYHEDRON);

template <typename T, typename T3>
std::vector<T3> surfaceDistribution(int n, SpherePlacement method = SpherePlacement::POLYHEDRON);
/// \}

} // namespace stmath
} // namespace stormm

#include "geodesic.tpp"

#endif
