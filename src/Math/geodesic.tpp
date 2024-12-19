// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T3>
int gatherTriangles(const int midpt_llim, const int midpt_hlim, const int ring_llim,
                    const int ring_hlim, const int nring,
                    const T crit_rsq, std::vector<T3> *result, int nadd_pt,
                    const std::vector<double3> &scaffold) {
  std::vector<int> plgn_points(nring);
  std::vector<int> ordr_points(nring);
  T3* result_ptr = result->data();
  int loc_nadd_pt = nadd_pt;
  for (int i = midpt_llim; i < midpt_hlim; i++) {
    plgn_points.resize(0);
    const double cn_x = scaffold[i].x;
    const double cn_y = scaffold[i].y;
    const double cn_z = scaffold[i].z;
    int npts = 0;
    for (int j = ring_llim; j < ring_hlim; j++) {
      const double dx = cn_x - scaffold[j].x;
      const double dy = cn_y - scaffold[j].y;
      const double dz = cn_z - scaffold[j].z;
      if ((dx * dx) + (dy * dy) + (dz * dz) < crit_rsq) {
        plgn_points.push_back(j);
      }
    }

    // Check that the ring has the expected number of points.
    if (static_cast<int>(plgn_points.size()) != nring) {
      rtWarn("With " + std::to_string(nring) + " points expected for a midpoint at { " +
             realToString(scaffold[i].x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(scaffold[i].y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(scaffold[i].z, 7, 4, NumberFormat::STANDARD_REAL) + "}, " +
             std::to_string(plgn_points.size()) + " were found within the critical radius of " +
             realToString(sqrt(crit_rsq), 7, 4, NumberFormat::STANDARD_REAL) + " Angstroms.",
             "gatherTriangles");
    }

    // Re-order the ring.
    std::vector<bool> taken(nring, false);
    taken[0] = true;
    ordr_points[0] = plgn_points[0];
    int curr_idx = plgn_points[0];
    for (int j = 1; j < nring; j++) {

      // The zeroth position is already taken and will therefore trigger the first unclaimed point
      // to qualify as a minimally distant point. 
      int min_loc = 0;
      double min_dr;
      for (int k = 0; k < nring; k++) {
        if (taken[k]) {
          continue;
        }
        const double dx = scaffold[plgn_points[k]].x - scaffold[curr_idx].x;
        const double dy = scaffold[plgn_points[k]].y - scaffold[curr_idx].y;
        const double dz = scaffold[plgn_points[k]].z - scaffold[curr_idx].z;
        const double dr = (dx * dx) + (dy * dy) + (dz * dz);
        if (min_loc == 0 || dr < min_dr) {
          min_dr = dr;
          min_loc = k;
        }
      }
      ordr_points[j] = plgn_points[min_loc];
      curr_idx = plgn_points[min_loc];
      taken[min_loc] = true;
    }

    // Loop over the ordered points and compute the centers of each triangle in the face.
    for (int j = 0; j < nring; j++) {
      const int jpt = ordr_points[j];
      const int jnpt = (j == nring - 1) ? ordr_points[0] : ordr_points[j + 1];
      const double tri_dx = cn_x + scaffold[jpt].x + scaffold[jnpt].x;
      const double tri_dy = cn_y + scaffold[jpt].y + scaffold[jnpt].y;
      const double tri_dz = cn_z + scaffold[jpt].z + scaffold[jnpt].z;
      const double inv_mag = 1.0 / sqrt((tri_dx * tri_dx) + (tri_dy * tri_dy) + (tri_dz * tri_dz));
      result_ptr[loc_nadd_pt].x = tri_dx * inv_mag;
      result_ptr[loc_nadd_pt].y = tri_dy * inv_mag;
      result_ptr[loc_nadd_pt].z = tri_dz * inv_mag;
      loc_nadd_pt++;
    }
  }

  // Update the number of added points
  return loc_nadd_pt;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T3>
std::vector<T3> surfaceDistribution(int *n, const SpherePlacement method) {
  std::vector<T3> result(*n);
  const int requested_n = *n;
  switch (method) {
  case SpherePlacement::POLYHEDRON:
    
    // Return for trivial cases
    if (requested_n <= 6 || requested_n == 8 || requested_n == 12) {
      result[0].x = 0.0;
      result[0].y = 0.0;
      result[0].z = 1.0;
    }
    if (requested_n == 2 || requested_n == 6 || requested_n == 12) {
      result[1].x = 0.0;
      result[1].y = 0.0;
      result[1].z = -1.0;
    }
    if (requested_n == 3) {

      // Equilaterial triangle in the Cartesian XZ plane
      result[1].x = sin(twopi / 3.0);
      result[1].y = 0.0;
      result[1].z = -0.5;
      result[2].x = -result[1].x;
      result[2].y = 0.0;
      result[2].z = -0.5;
    }
    else if (requested_n == 4 || requested_n == 8) {

      // Regular tetrahedron with one point at the sphere's north pole
      result[1].x = sin(tetrahedral_angle);
      result[1].y = 0.0;
      result[1].z = cos(tetrahedral_angle);
      result[2].x = result[1].x * cos(twopi / 3.0);
      result[2].y = result[1].x * sin(twopi / 3.0);
      result[2].z = result[1].z;
      result[3].x = result[2].x;
      result[3].y = -result[2].y;
      result[3].z = result[1].z;
      if (requested_n == 8) {

        // The cube is obtained by inverting and rotating the regular tetrahedron
        const std::vector<double3> tet = surfaceDistribution<double, double3>(4);
        const double3 z_axis = { 0.0, 0.0, 1.0 };
        const std::vector<double> u = rotationMatrixAboutVector(z_axis, pi / 3.0);
        for (int i = 0; i < 4; i++) {
          result[i].x = tet[i].x;
          result[i].y = tet[i].y;
          result[i].z = tet[i].z;
          result[i + 4].x = (u[0] * tet[i].x) + (u[3] * tet[i].y) + (u[6] * tet[i].z);
          result[i + 4].y = (u[1] * tet[i].x) + (u[4] * tet[i].y) + (u[7] * tet[i].z);
          result[i + 4].z = -((u[2] * tet[i].x) + (u[5] * tet[i].y) + (u[8] * tet[i].z));
        }
      }
    }
    else if (requested_n == 5) {

      // Make a tetragonal (Egyptian) pyramid such that the imbalance of the longer distance
      // between opposite corners of the base (opp_b), the distances between the base and the top
      // (top_b), and the distances between adjacent corners of the base (adj_b) is minimized
      // overall.
      //
      // T = (4 top_b) + (4 adj_b) + (2 opp_b)
      // opp_b = sqrt(2) adj_b
      // T = (4 top_b) + ((4 + 2 sqrt(2)) adj_b)
      //
      // All of this depends on an angle theta.
      // top_b = sqrt((1 - cos(theta))^2 + sin(theta)^2)
      // adj_b = sqrt(2) * sin(theta)
      // T = 4 * sqrt((1 - cos(theta))^2 + sin(theta)^2) + ((4 + 2 sqrt(2)) * sqrt(2) * sin(theta))
      //
      // This was evaluated numerically to theta = 1.8267045013301753.
      const double best_theta = 1.8267045013301753;
      result[1].x = sin(best_theta);
      result[1].y = 0.0;
      result[1].z = cos(best_theta);
      result[2].x = 0.0;
      result[2].y =  result[1].x;
      result[2].z =  result[1].z;
      result[3].x = -result[1].x;
      result[3].y = 0.0;
      result[3].z =  result[1].z;
      result[4].x = 0.0;
      result[4].y = -result[1].x;
      result[4].z =  result[1].z;
    }
    else if (requested_n == 6) {

      // Regular octahedron
      result[2].x = 1.0;
      result[2].y = 0.0;
      result[2].z = 0.0;
      result[3].x = 0.0;
      result[3].y = 1.0;
      result[3].z = 0.0;
      result[4].x = -1.0;
      result[4].y = 0.0;
      result[4].z = 0.0;
      result[5].x = 0.0;
      result[5].y = -1.0;
      result[5].z = 0.0;
    }
    else if (requested_n == 12) {

      // Regular icosahedron.  There are two bands of five points each, staggered relative to one
      // another if viewed down the Cartesian Z axis.  Compute the radius of the band by noting
      // that either band is 1.0 units from the origin and the upper band is as far from the first
      // point at (0, 0, 1) as any of its points are from their cognates in the second band at the
      // same distance from (0, 0, -1).  The band is thus defined by five line segments of unknown
      // length S, but they are isosceles triangles with a base of length S, and they connect to an
      // equilateral triangle of side length S.  The widest angle in the isosceles triangle is 72
      // degrees, that of the equilateral triangle 60 degrees.  The ratio of the short sides of the
      // isosceles triangle to its base is thus indicative of the ratio of the side lengths of the
      // base of a right triangle and its hypotenuse.
      const double isoc_short_side = 0.5 / sin(pi * 0.2);
      const double upper_band_gap = sqrt(1.0 - (isoc_short_side * isoc_short_side));

      // The upper_band_gap is the distance, in units of the length between two of the points in
      // the icosahedron, that the upper band sits below the very top of the icosahedron.  However,
      // it is not yet known what the actual radius of this upper band is, except in relation to
      // the distance between two of the points (isoc_short_side).  To find the distance in terms
      // of the radius of the sphere upon which all points are inscribed, consider the upper and
      // lower bans, separated by 2 x (G - upper_band_gap) and twisted by pi / 5 relative to one
      // another.  It is not known what G is, but it is the ratio of the edge length of the
      // icosahedron to the radius of the sphere.  The distance between a point on one and the same
      // point on the other will also be the edge length of the icosahedron.
      const double dx = isoc_short_side * (cos(pi * 0.2) - 1.0);
      const double dy = isoc_short_side * sin(pi * 0.2);
      const double dz = sqrt(1.0 - (dx * dx) - (dy * dy));
      const double isoc_edge_length = 2.0 / (dz + (2.0 * upper_band_gap));
      const double band_height = 1.0 - (isoc_edge_length * upper_band_gap);
      const double band_rad = isoc_edge_length * isoc_short_side;
      for (int i = 0; i < 5; i++) {
        const double di = i;
        result[2 + i].x = band_rad * sin(0.4 * pi * di);
        result[2 + i].y = band_rad * cos(0.4 * pi * di);
        result[2 + i].z = band_height;
        result[7 + i].x = band_rad * sin(0.4 * pi * (0.5 + di));
        result[7 + i].y = band_rad * cos(0.4 * pi * (0.5 + di));
        result[7 + i].z = -band_height;
      }
    }
    else if (requested_n == 20) {

      // Inscribe a dodecahedron in the icosahedron.
      const std::vector<double3> ico = surfaceDistribution<double, double3>(12);
      std::vector<int3> faces(20);
      for (int i = 0; i < 5; i++) {
        faces[i     ] = { 0, 2 + i, (i == 4) ? 2 : 3 + i };
        faces[i +  5] = { 1, 7 + i, (i == 4) ? 7 : 8 + i };
        faces[i + 10] = { 2 + i, (i == 4) ? 2 : 3 + i, 7 + i };
        faces[i + 15] = { 7 + i, (i == 0) ? 11 : 6 + i, 2 + i }; 
      }
      for (int i = 0; i < 20; i++) {
        const double dx = (ico[faces[i].x].x + ico[faces[i].y].x + ico[faces[i].z].x);
        const double dy = (ico[faces[i].x].y + ico[faces[i].y].y + ico[faces[i].z].y);
        const double dz = (ico[faces[i].x].z + ico[faces[i].y].z + ico[faces[i].z].z);
        const double inv_mag = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
        result[i].x = dx * inv_mag;
        result[i].y = dy * inv_mag;
        result[i].z = dz * inv_mag;
      }
    }
    else if (requested_n == 24) {

      // Truncate the octahedron.
      const std::vector<double3> octahedron = surfaceDistribution<double, double3>(6);
      std::vector<int2> edges(12);
      for (int i = 0; i < 4; i++) {
        edges[i] = { 0, i + 2 };
        edges[i + 4] = { 1, i + 2 };
        edges[i + 8] = { i + 2, (i == 3) ? 2 : i + 3 };
      }
      for (int i = 0; i < 12; i++) {
        const size_t atm_a = edges[i].x;
        const size_t atm_b = edges[i].y;
        result[i     ].x = ((2.0 * octahedron[atm_a].x) + octahedron[atm_b].x) / 3.0;
        result[i     ].y = ((2.0 * octahedron[atm_a].y) + octahedron[atm_b].y) / 3.0;
        result[i     ].z = ((2.0 * octahedron[atm_a].z) + octahedron[atm_b].z) / 3.0;
        result[i + 12].x = (octahedron[atm_a].x + (2.0 * octahedron[atm_b].x)) / 3.0;
        result[i + 12].y = (octahedron[atm_a].y + (2.0 * octahedron[atm_b].y)) / 3.0;
        result[i + 12].z = (octahedron[atm_a].z + (2.0 * octahedron[atm_b].z)) / 3.0;
      }
    }
    else if (requested_n == 32) {

      // Augment the icosahedron by taking the midpoints of the vertices of each face and then
      // projecting these points outward to the sphere surface.  This is equivalent to obtaining
      // the dodecahedron and adding it to the icosahedron.
      const std::vector<double3> ico = surfaceDistribution<double, double3>(12);
      for (int i = 0; i < 12; i++) {
        result[i].x = ico[i].x;
        result[i].y = ico[i].y;
        result[i].z = ico[i].z;
      }
      const std::vector<double3> dodeca = surfaceDistribution<double, double3>(20);
      for (int i = 0; i < 20; i++) {
        result[i + 12].x = dodeca[i].x;
        result[i + 12].y = dodeca[i].y;
        result[i + 12].z = dodeca[i].z;
      }
    }
    else if (requested_n == 60) {

      // Truncate the icosahedron.
      const std::vector<double3> icosahedron = surfaceDistribution<double, double3>(12);
      std::vector<int2> edges(30);
      for (int i = 0; i < 5; i++) {
        edges[     i] = { 0,  2 + i };
        edges[ 5 + i] = { 1,  7 + i };
        edges[10 + i] = { 2 + i, (i == 4) ? 2 : 3 + i };
        edges[15 + i] = { 7 + i, (i == 4) ? 7 : 8 + i };
        edges[20 + i] = { 2 + i, 7 + i };
        edges[25 + i] = { 2 + i, (i == 0) ? 11 : 6 + i };
      }
      std::vector<double3> tmp(60);
      for (int i = 0; i < 30; i++) {
        const double3 pa = icosahedron[edges[i].x];
        const double3 pb = icosahedron[edges[i].y];
        tmp[     i] = { (2.0 * pa.x) + pb.x, (2.0 * pa.y) + pb.y, (2.0 * pa.z) + pb.z };
        tmp[30 + i] = { (2.0 * pb.x) + pa.x, (2.0 * pb.y) + pa.y, (2.0 * pb.z) + pa.z };
      }
      const double inv_rad = 1.0 / sqrt((tmp[0].x * tmp[0].x) + (tmp[0].y * tmp[0].y) +
                                        (tmp[0].z * tmp[0].z));
      for (int i = 0; i < 60; i++) {
        result[i].x = tmp[i].x * inv_rad;
        result[i].y = tmp[i].y * inv_rad;
        result[i].z = tmp[i].z * inv_rad;
      }
    }
    else if (requested_n == 92) {

      // Make a geodesic of 60 points, then place additional verticies at the centers of each of
      // the twelve pentagons and twenty hexagons using the augmented icosahedron.  The truncated
      // and augmneted icosahedrons will complement each other without the need for further
      // rescaling or rotation.
      const std::vector<double3> soccer_ball = surfaceDistribution<double, double3>(60);
      for (int i = 0; i < 60; i++) {
        result[i].x = soccer_ball[i].x;
        result[i].y = soccer_ball[i].y;
        result[i].z = soccer_ball[i].z;
      }
      const std::vector<double3> aug_ico = surfaceDistribution<double, double3>(32);
      for (int i = 0; i < 32; i++) {
        result[60 + i].x = aug_ico[i].x;
        result[60 + i].y = aug_ico[i].y;
        result[60 + i].z = aug_ico[i].z;
      }
    }
    else if (requested_n == 272) {
      const std::vector<double3> aug_soccer = surfaceDistribution<double, double3>(92);

      // Find each pentagon and hexagon by scanning over the final 32 points and searching among
      // the first 60 for the closest partners.  A good metric of the distance between points is
      // the edge length of the icosahedron.  The edge lengths of the soccer ball (truncated
      // icosahedron) are one third of this, and the distances between the centers of the hexagons
      // to the surrounding six points will be slightly more than that as the center point gets
      // pushed out slightly to reach the surface of the unit sphere.  Take (the square of) 35% of
      // the edge length of the icosahedron as the measure to define each hexagon or pentagon.
      const double dx_ico = aug_soccer[61].x - aug_soccer[60].x;
      const double dy_ico = aug_soccer[61].y - aug_soccer[60].y;
      const double dz_ico = aug_soccer[61].z - aug_soccer[60].z;
      const double crit_rsq = ((dx_ico * dx_ico) + (dy_ico * dy_ico) + (dz_ico * dz_ico)) * 0.1225;

      // Points 60 up to 72 are the verticies of the icosahedron.  Each will collect a pentagon.
      int nadd_pt = 92;
      nadd_pt = gatherTriangles(60, 72, 0, 60, 5, crit_rsq, &result, nadd_pt, aug_soccer);

      // The final twenty points are the centers of faces of the icosahedron and will each collect
      // a hexagon.
      nadd_pt = gatherTriangles(72, 92, 0, 60, 6, crit_rsq, &result, nadd_pt, aug_soccer);

      // Add the original 92 points
      for (int i = 0; i < 92; i++) {
        result[i].x = aug_soccer[i].x;
        result[i].y = aug_soccer[i].y;
        result[i].z = aug_soccer[i].z;
      }
    }

    // Return the result of the polyhedral calculation
    if (requested_n <=   6 || requested_n ==   8 || requested_n ==  12 || requested_n ==  20 ||
        requested_n ==  24 || requested_n ==  32 || requested_n ==  60 || requested_n ==  92 ||
        requested_n == 272) {
      return result;
    }
    else {
      return surfaceDistribution<T, T3>(n, SpherePlacement::DIST_FIBONACCI);
    }
    break;
  case SpherePlacement::AREA_FIBONACCI:
    {
      // Distribute the points using the Fibonacci sphere algorithm
      const double phi = pi * (sqrt(5.0) - 1.0);
      const double dn_minus_one = requested_n - 1;
      for (int i = 0; i < requested_n; i++) {
        const double di = i;
        const double y = 1.0 - (di / dn_minus_one) * 2.0;
        const double rad_y = sqrt(1.0 - (y * y));
        const double theta = phi * di;
        const double x = cos(theta) * rad_y;
        const double z = sin(theta) * rad_y;
        result[i].x = x;
        result[i].y = y;
        result[i].z = z;
      }
    }
    break;
  case SpherePlacement::DIST_FIBONACCI:
    {
      const double golden_ratio = 0.5 * (1.0 + sqrt(5.0));
      const double dn = requested_n;
      for (int i = 0; i < requested_n; i++) {
        const double di = i;
        const double theta = twopi * di / golden_ratio;
        const double phi = acos(1.0 - 2.0 * (di + 0.5) / dn);
        result[i].x = cos(theta) * sin(phi);
        result[i].y = sin(theta) * sin(phi);
        result[i].z = cos(phi);
      }
    }
    break;
  case SpherePlacement::DESERNO:
    {
      int pt_idx = 0;
      const double dn = requested_n;
      const double alpha = (4.0 * pi / dn);
      const double m_nu = round(pi / sqrt(alpha));
      const int im_nu = m_nu;
      const double d_nu = pi / m_nu;
      const double d_phi = alpha / d_nu;
      for (int i = 0; i < im_nu; i++) {
        const double di = i;
        const double nu = pi * (di + 0.5) / m_nu;
        const double m_phi = round(twopi * sin(nu) / d_phi);
        const int im_phi = m_phi;
        for (int j = 0; j < im_phi; j++) {
          const double dj = j;
          const double phi = twopi * dj / m_phi;
          result[pt_idx].x = sin(nu) * cos(phi);
          result[pt_idx].y = sin(nu) * sin(phi);
          result[pt_idx].z = cos(nu);
          pt_idx++;
        }
      }
      result.resize(pt_idx);
      *n = pt_idx;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T3>
std::vector<T3> surfaceDistribution(const int n, const SpherePlacement method) {
  int loc_n = n;
  return surfaceDistribution<T, T3>(&loc_n, method);
}

} // namespace stmath
} // namespace stormm
