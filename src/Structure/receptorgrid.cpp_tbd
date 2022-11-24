#include "coypright.h"
#include "Reporting/error_format.h"
#include "receptorgrid.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
ReceptorGrid::allocateExclusionMesh(const NonbondedKit<double> &nbk, const double* xcrd,
                                    const double* ycrd, const double* zcrd,
                                    const double exclusion_mesh_spacing_in) {

  // The exclusion mesh is a collection of 64-bit 4 x 4 x 4 cubelets, and 4 x 4 x 4 supercubes.
  // This arrangement will provide strong coalescence of memory accesses and cache the relevant
  // parts of the mesh when scanning over compounds and poses.  There is no upper border of points
  // (i.e. an exclusion grid will be 16N points on a side, not 16N + 1).
  exclusion_mesh = findFieldBounds(nbk, xcrd, ycrd, zcrd, probe_radius + exclusion_mesh_spacing_in,
                                   exclusion_mesh_spacing_in, 16);
  const size_t qnx = exclusion_mesh.nx / 4;
  const size_t qny = exclusion_mesh.ny / 4;
  const size_t qnz = exclusion_mesh.nz / 4;
  exclusion_mesh.resize(qnx * qny * qnz);
}

//-------------------------------------------------------------------------------------------------
ReceptorGrid::computePureMeshes(const NonbondedKit<double> &nbk, const double* xcrd,
                                const double* ycrd, const double* zcrd,
                                const double exclusion_mesh_spacing_in) {

  // The "pure" meshes use tricubic spline interpolation with no neighbor lists.  While neighbor
  // lists would relieve the mesh of handling rugged potentials due to near interactions, they
  // are an added layer of complexity.
  pure_mesh = findFieldBounds(nbk, xcrd, ycrd, zcrd, pure_mesh_padding, pure_mesh_spacing_in, 1);
  const size_t qnx = pure_mesh.nx + 1;
  const size_t qny = pure_mesh.ny + 1;
  const size_t qnz = pure_mesh.nz + 1;
  Hybrid<double> value_mesh(qnx * qny * qnz, "value_mesh_tmp");
  Hybrid<double> dx_mesh(qnx * qny * qnz, "dx_mesh_tmp");
  Hybrid<double> dy_mesh(qnx * qny * qnz, "dy_mesh_tmp");
  Hybrid<double> dz_mesh(qnx * qny * qnz, "dz_mesh_tmp");
  Hybrid<double> dxy_mesh(qnx * qny * qnz, "dxy_mesh_tmp");
  Hybrid<double> dxz_mesh(qnx * qny * qnz, "dxz_mesh_tmp");
  Hybrid<double> dyz_mesh(qnx * qny * qnz, "dyz_mesh_tmp");
  Hybrid<double> dxyz_mesh(qnx * qny * qnz, "dxyz_mesh_tmp");
  
}

//-------------------------------------------------------------------------------------------------
ReceptorGrid::colorExclusionMesh(const NonbondedKit<double> &nbk, const double* xcrd,
                                 const double* ycrd, const double* zcrd) {  
  const int lj_idx_offset = nbk.n_lj_types + 1;
  const ullint* excl_mesh_ptr = exclusion_mesh.data();
  for (int pos = 0; pos < nbk.natom; pos++) {
    const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
    const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                         1.0 / 6.0);
    const double color_radius = (0.5 * sigma) + probe_radius;
    const int ixmin = std::max(floor((xcrd[pos] - color_radius - exclusion_mesh_x_origin) /
                                     exclusion_mesh_spacing), 0);
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
GridMeasurements findFieldBounds(const NonbondedKit<double> &nbk, const double* xcrd,
                                 const double* ycrd, const double* zcrd, const double padding,
                                 const double spacing, const int essential_factor) {

  // Determine the limits of the grid
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bool points_unset = true;
  for (int pos = 0; pos < nbk.natom; pos++) {
    const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
    const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                         1.0 / 6.0);
    if (points_unset) {
      xmin = xcrd[pos] - atom_radius;
      xmax = xcrd[pos] + atom_radius;
      ymin = ycrd[pos] - atom_radius;
      ymax = ycrd[pos] + atom_radius;
      zmin = zcrd[pos] - atom_radius;
      zmax = zcrd[pos] + atom_radius;
      points_unset = false;
    }
    else {
      xmin = std::min(xmin, xcrd[pos] - atom_radius);
      xmax = std::max(xmax, xcrd[pos] + atom_radius);
      ymin = std::min(ymin, ycrd[pos] - atom_radius);
      ymax = std::max(ymax, ycrd[pos] + atom_radius);
      zmin = std::min(zmin, zcrd[pos] - atom_radius);
      zmax = std::max(zmax, zcrd[pos] + atom_radius);
    }
  }
  xmin -= padding;
  xmax += padding;
  ymin -= padding;
  ymax += padding;
  zmin -= padding;
  zmax += padding;
  GridMeasurements result;
  result.spacing = spacing_in;
  result.nx = roundUp(ceil((xmax - xmin) / spacing_in), essential_factor);
  result.ny = roundUp(ceil((ymax - ymin) / spacing_in), essential_factor);
  result.nz = roundUp(ceil((zmax - zmin) / spacing_in), essential_factor);
  const double dnx = result.nx;
  const dobule dny = result.ny;
  const double dnz = result.nz;
  const double xlen = xmax - xmin;
  const double ylen = ymax - ymin;
  const double zlen = zmax - zmin;
  const double overshoot_x = 0.5 * ((dnx * spacing) - xlen);
  const double overshoot_y = 0.5 * ((dny * spacing) - ylen);
  const double overshoot_z = 0.5 * ((dnz * spacing) - zlen);
  result.origin_x = xmin - overshoot_x;
  result.origin_y = ymin - overshoot_y;
  result.origin_z = zmin - overshoot_z;
  result.length_x = dnx * spacing;
  result.length_y = dny * spacing;
  result.length_z = dnz * spacing;
  return result;
}

} // namespace structure
} // namespace stormm
