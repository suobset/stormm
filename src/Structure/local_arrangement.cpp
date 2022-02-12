#include <cmath>
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "local_arrangement.h"

namespace omni {
namespace structure {

using math::crossProduct;
using trajectory::getCoordinateFrameReader;

//-------------------------------------------------------------------------------------------------
void imageCoordinates(double *x, double *y, double *z, const double* umat, const double* invu,
                      const UnitCellType unit_cell) {
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    {
      double local_x = (*x) * umat[0];
      double local_y = (*y) * umat[4];
      double local_z = (*z) * umat[8];
      local_x = (local_x > 0.0) ? local_x - ceil(local_x) : local_x - floor(local_x);
      local_y = (local_y > 0.0) ? local_y - ceil(local_y) : local_y - floor(local_y);
      local_z = (local_z > 0.0) ? local_z - ceil(local_z) : local_z - floor(local_z);
      *x = local_x * invu[0];
      *y = local_y * invu[4];
      *z = local_z * invu[8];
    }
    break;
  case UnitCellType::TRICLINIC:
    {
      const double local_x = (*x);
      const double local_y = (*y);
      const double local_z = (*z);
      double ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      double ndy = (umat[1] * local_x) + (umat[4] * local_y) + (umat[7] * local_z);
      double ndz = (umat[2] * local_x) + (umat[5] * local_y) + (umat[8] * local_z);
      ndx = (ndx > 0.0) ? ndx - ceil(ndx) : ndx - floor(ndx);
      ndy = (ndy > 0.0) ? ndy - ceil(ndy) : ndy - floor(ndy);
      ndz = (ndz > 0.0) ? ndz - ceil(ndz) : ndz - floor(ndz);
      *x = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
      *y = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
      *z = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(double* x, double* y, double* z, const int length, const double* umat,
                      const double* invu, const UnitCellType unit_cell) {
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    for (int i = 0; i < length; i++) {
      double local_x = x[i] * umat[0];
      double local_y = y[i] * umat[4];
      double local_z = z[i] * umat[8];
      local_x = (local_x > 0.0) ? local_x - ceil(local_x) : local_x - floor(local_x);
      local_y = (local_y > 0.0) ? local_y - ceil(local_y) : local_y - floor(local_y);
      local_z = (local_z > 0.0) ? local_z - ceil(local_z) : local_z - floor(local_z);
      x[i] = local_x * invu[0];
      y[i] = local_y * invu[4];
      z[i] = local_z * invu[8];
    }
    break;
  case UnitCellType::TRICLINIC:
    for (int i = 0; i < length; i++) {
      const double local_x = x[i];
      const double local_y = y[i];
      const double local_z = z[i];
      double ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      double ndy = (umat[1] * local_x) + (umat[4] * local_y) + (umat[7] * local_z);
      double ndz = (umat[2] * local_x) + (umat[5] * local_y) + (umat[8] * local_z);
      ndx = (ndx > 0.0) ? ndx - ceil(ndx) : ndx - floor(ndx);
      ndy = (ndy > 0.0) ? ndy - ceil(ndy) : ndy - floor(ndy);
      ndz = (ndz > 0.0) ? ndz - ceil(ndz) : ndz - floor(ndz);
      x[i] = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
      y[i] = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
      z[i] = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(std::vector<double> *x, std::vector<double> *y, std::vector<double> *z,
                      const double* umat, const double* invu, const UnitCellType unit_cell) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(Hybrid<double> *x, Hybrid<double> *y, Hybrid<double> *z,
                      const double* umat, const double* invu,
                      const UnitCellType unit_cell) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrameReader &cfr) {
  double dx = cfr.xcrd[atom_j] - cfr.xcrd[atom_i];
  double dy = cfr.ycrd[atom_j] - cfr.ycrd[atom_i];
  double dz = cfr.zcrd[atom_j] - cfr.zcrd[atom_i];
  imageCoordinates(&dx, &dy, &dz, cfr.umat, cfr.invu, cfr.unit_cell);
  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrame &cf) {
  return distance(atom_i, atom_j, cf.data());
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpace &ps) {
  return distance(atom_i, atom_j, getCoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k,
             const CoordinateFrameReader &cfr) {

  // Image all three atoms and put the atom J at the origin
  double rix = cfr.xcrd[atom_i] - cfr.xcrd[atom_j];
  double riy = cfr.ycrd[atom_i] - cfr.ycrd[atom_j];
  double riz = cfr.zcrd[atom_i] - cfr.zcrd[atom_j];
  double rkx = cfr.xcrd[atom_k] - cfr.xcrd[atom_j];
  double rky = cfr.ycrd[atom_k] - cfr.ycrd[atom_j];
  double rkz = cfr.zcrd[atom_k] - cfr.zcrd[atom_j];
  imageCoordinates(&rix, &riy, &riz, cfr.umat, cfr.invu, cfr.unit_cell);
  imageCoordinates(&rkx, &rky, &rkz, cfr.umat, cfr.invu, cfr.unit_cell);

  // Compute the angle, in radians
  const double mgba = (rix * rix) + (riy * riy) + (riz * riz);
  const double mgbc = (rkx * rkx) + (rky * rky) + (rkz * rkz);
  const double invbabc = 1.0 / sqrt(mgba * mgbc);
  double costheta = ((rix * rkx) + (riy * rky) + (riz * rkz)) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  return acos(costheta);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const CoordinateFrame &cf) {
  return angle(atom_i, atom_j, atom_k, cf.data());
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const PhaseSpace &ps) {
  return angle(atom_i, atom_j, atom_k, getCoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l,
                      const CoordinateFrameReader &cfr) {

  // Image all four atoms and put atom K at the origin
  double rix = cfr.xcrd[atom_i] - cfr.xcrd[atom_k];
  double riy = cfr.ycrd[atom_i] - cfr.ycrd[atom_k];
  double riz = cfr.zcrd[atom_i] - cfr.zcrd[atom_k];
  double rjx = cfr.xcrd[atom_j] - cfr.xcrd[atom_k];
  double rjy = cfr.ycrd[atom_j] - cfr.ycrd[atom_k];
  double rjz = cfr.zcrd[atom_j] - cfr.zcrd[atom_k];
  double rlx = cfr.xcrd[atom_l] - cfr.xcrd[atom_k];
  double rly = cfr.ycrd[atom_l] - cfr.ycrd[atom_k];
  double rlz = cfr.zcrd[atom_l] - cfr.zcrd[atom_k];
  imageCoordinates(&rix, &riy, &riz, cfr.umat, cfr.invu, cfr.unit_cell);
  imageCoordinates(&rjx, &rjy, &rjz, cfr.umat, cfr.invu, cfr.unit_cell);
  imageCoordinates(&rlx, &rly, &rlz, cfr.umat, cfr.invu, cfr.unit_cell);

  // Compute the dihedral angle, in radians
  double ab[3], bc[3], cd[3];
  ab[0] = rjx - rix;
  ab[1] = rjy - riy;
  ab[2] = rjz - riz;
  bc[0] = -rjx;
  bc[1] = -rjy;
  bc[2] = -rjz;
  cd[0] = rlx;
  cd[1] = rly;
  cd[2] = rlz;

  // Compute cross products and then the angle between the planes
  double crabbc[3], crbccd[3], scr[3];
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                   (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  const double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                            -acos(costheta);
  return theta;
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const CoordinateFrame &cf) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, cf.data());
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const PhaseSpace &ps) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, getCoordinateFrameReader(ps));
}

} // namespace geometry
} // namespace omni
