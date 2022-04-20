// -*-c++-*-
namespace omni {
namespace structure {

using math::crossProduct;
using trajectory::CoordinateFrameWriter;
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcoord imageValue(const Tcoord x, const Tcalc range, const ImagingMethod style,
                 const Tcalc globalpos_scale_factor) {
  Tcalc x_frac;
  if (tcoord_is_sgnint) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_globalpos_scale_factor = 1.0 / globalpos_scale_factor;
    x_frac = static_cast<Tcalc>(x) * inv_globalpos_scale_factor / range;
  }
  else {
    x_frac = x / range;
  }
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc half = 0.5;
  switch (style) {
  case ImagingMethod::PRIMARY_UNIT_CELL:
    if (tcalc_is_double) {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floor(x_frac)) * range * globalpos_scale_factor);
      }
      else {
        return ((x_frac - floor(x_frac)) * range);
      }
    }
    else {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floorf(x_frac)) * range * globalpos_scale_factor);
      }
      else {
        return ((x_frac - floorf(x_frac)) * range);
      }
    }
  case ImagingMethod::MINIMUM_IMAGE:
    if (tcalc_is_double) {
      x_frac -= ((x_frac >= half) *  ceil(x_frac - half)) +
                ((x_frac < -half) * floor(x_frac + half));
    }
    else {
      x_frac -= ((x_frac >= half) *  ceilf(x_frac - half)) +
                ((x_frac < -half) * floorf(x_frac + half));
    }

    // The final subtraction covers the case of the re-imaged coordinate sitting right on 0.5
    if (tcoord_is_sgnint) {
      return llround((x_frac - (x_frac >= half)) * range * globalpos_scale_factor);
    }
    else {
      return (x_frac - (x_frac >= half)) * range;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void imageCoordinates(Tcoord *x, Tcoord *y, Tcoord *z, const Tcalc* umat, const Tcalc* invu,
                      const UnitCellType unit_cell, const ImagingMethod style) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc half = 0.5;
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    {
      Tcalc local_x = (*x) * umat[0];
      Tcalc local_y = (*y) * umat[4];
      Tcalc local_z = (*z) * umat[8];
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          local_x -= floor(local_x);
          local_y -= floor(local_y);
          local_z -= floor(local_z);
        }
        else {
          local_x -= floorf(local_x);
          local_y -= floorf(local_y);
          local_z -= floorf(local_z);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          local_x -= ((local_x >= 0.5) *  ceil(local_x - 0.5)) +
                     ((local_x < -0.5) * floor(local_x + 0.5));
          local_y -= ((local_y >= 0.5) *  ceil(local_y - 0.5)) +
                     ((local_y < -0.5) * floor(local_y + 0.5));
          local_z -= ((local_z >= 0.5) *  ceil(local_z - 0.5)) +
                     ((local_z < -0.5) * floor(local_z + 0.5));
          local_x -= (local_x >= 0.5);
          local_y -= (local_y >= 0.5);
          local_z -= (local_z >= 0.5);
        }
        else {
          local_x -= ((local_x >= 0.5) *  ceilf(local_x - 0.5)) +
                     ((local_x < -0.5) * floorf(local_x + 0.5));
          local_y -= ((local_y >= 0.5) *  ceilf(local_y - 0.5)) +
                     ((local_y < -0.5) * floorf(local_y + 0.5));
          local_z -= ((local_z >= 0.5) *  ceilf(local_z - 0.5)) +
                     ((local_z < -0.5) * floorf(local_z + 0.5));
          local_x -= (local_x >= 0.5);
          local_y -= (local_y >= 0.5);
          local_z -= (local_z >= 0.5);
        }
        break;
      }
      *x = local_x * invu[0];
      *y = local_y * invu[4];
      *z = local_z * invu[8];
    }
    break;
  case UnitCellType::TRICLINIC:
    {      
      const Tcoord local_x = (*x);
      const Tcoord local_y = (*y);
      const Tcoord local_z = (*z);
      double ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      double ndy = (umat[1] * local_x) + (umat[4] * local_y) + (umat[7] * local_z);
      double ndz = (umat[2] * local_x) + (umat[5] * local_y) + (umat[8] * local_z);
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        ndx -= floor(ndx);
        ndy -= floor(ndy);
        ndz -= floor(ndz);
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        ndx -= ((ndx >= 0.5) * ceil(ndx - 0.5)) + ((ndx < -0.5) * floor(ndx + 0.5));
        ndy -= ((ndy >= 0.5) * ceil(ndy - 0.5)) + ((ndy < -0.5) * floor(ndy + 0.5));
        ndz -= ((ndz >= 0.5) * ceil(ndz - 0.5)) + ((ndz < -0.5) * floor(ndz + 0.5));
        ndx -= (ndx >= 0.5);
        ndy -= (ndy >= 0.5);
        ndz -= (ndz >= 0.5);
        break;
      }
      *x = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
      *y = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
      *z = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(double* x, double* y, double* z, const int length, const double* umat,
                      const double* invu, const UnitCellType unit_cell,
                      const ImagingMethod style) {
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    for (int i = 0; i < length; i++) {
      double local_x = x[i] * umat[0];
      double local_y = y[i] * umat[4];
      double local_z = z[i] * umat[8];
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        local_x -= floor(local_x);
        local_y -= floor(local_y);
        local_z -= floor(local_z);
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        local_x -= ((local_x >= 0.5) *  ceil(local_x - 0.5)) +
                   ((local_x < -0.5) * floor(local_x + 0.5));
        local_y -= ((local_y >= 0.5) *  ceil(local_y - 0.5)) +
                   ((local_y < -0.5) * floor(local_y + 0.5));
        local_z -= ((local_z >= 0.5) *  ceil(local_z - 0.5)) +
                   ((local_z < -0.5) * floor(local_z + 0.5));
        local_x -= (local_x >= 0.5);
        local_y -= (local_y >= 0.5);
        local_z -= (local_z >= 0.5);
        break;
      }
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
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        ndx -= floor(ndx);
        ndy -= floor(ndy);
        ndz -= floor(ndz);
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        ndx -= ((ndx >= 0.5) * ceil(ndx - 0.5)) + ((ndx < -0.5) * floor(ndx + 0.5));
        ndy -= ((ndy >= 0.5) * ceil(ndy - 0.5)) + ((ndy < -0.5) * floor(ndy + 0.5));
        ndz -= ((ndz >= 0.5) * ceil(ndz - 0.5)) + ((ndz < -0.5) * floor(ndz + 0.5));
        ndx -= (ndx >= 0.5);
        ndy -= (ndy >= 0.5);
        ndz -= (ndz >= 0.5);
        break;
      }
      x[i] = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
      y[i] = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
      z[i] = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(std::vector<double> *x, std::vector<double> *y, std::vector<double> *z,
                      const double* umat, const double* invu, const UnitCellType unit_cell,
                      const ImagingMethod style) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(Hybrid<double> *x, Hybrid<double> *y, Hybrid<double> *z,
                      const double* umat, const double* invu, const UnitCellType unit_cell,
                      const ImagingMethod style) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(PhaseSpace *ps, const ImagingMethod style) {
  CoordinateFrameWriter cfw(ps);
  imageCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, cfw.umat, cfw.invu, cfw.unit_cell,
                   style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(CoordinateFrame *cf, const ImagingMethod style) {
  CoordinateFrameWriter cfw = cf->data();
  imageCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, cfw.umat, cfw.invu, cfw.unit_cell,
                   style);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const double* xcrd, const double* ycrd,
                const double* zcrd, const double* umat, const double* invu,
                const UnitCellType unit_cell) {
  double dx = xcrd[atom_j] - xcrd[atom_i];
  double dy = ycrd[atom_j] - ycrd[atom_i];
  double dz = zcrd[atom_j] - zcrd[atom_i];
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrameReader &cfr) {
  return distance(atom_i, atom_j, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrame &cf) {
  return distance(atom_i, atom_j, cf.data());
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpace &ps) {
  return distance(atom_i, atom_j, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const double* xcrd,
             const double* ycrd, const double* zcrd, const double* umat, const double* invu,
             const UnitCellType unit_cell) {

  // Image all three atoms and put the atom J at the origin
  double rix = xcrd[atom_i] - xcrd[atom_j];
  double riy = ycrd[atom_i] - ycrd[atom_j];
  double riz = zcrd[atom_i] - zcrd[atom_j];
  double rkx = xcrd[atom_k] - xcrd[atom_j];
  double rky = ycrd[atom_k] - ycrd[atom_j];
  double rkz = zcrd[atom_k] - zcrd[atom_j];
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth);
  imageCoordinates(&rkx, &rky, &rkz, umat, invu, unit_cell, imeth);

  // Compute the angle, in radians
  const double mgba = (rix * rix) + (riy * riy) + (riz * riz);
  const double mgbc = (rkx * rkx) + (rky * rky) + (rkz * rkz);
  const double invbabc = 1.0 / sqrt(mgba * mgbc);
  double costheta = ((rix * rkx) + (riy * rky) + (riz * rkz)) * invbabc;
  costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
  return acos(costheta);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k,
             const CoordinateFrameReader &cfr) {
  return angle(atom_i, atom_j, atom_k, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
               cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const CoordinateFrame &cf) {
  return angle(atom_i, atom_j, atom_k, cf.data());
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const PhaseSpace &ps) {
  return angle(atom_i, atom_j, atom_k, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const double* xcrd,
                      const double* ycrd, const double* zcrd, const double* umat,
                      const double* invu, const UnitCellType unit_cell) {

  // Image all four atoms and put atom K at the origin
  double rix = xcrd[atom_i] - xcrd[atom_k];
  double riy = ycrd[atom_i] - ycrd[atom_k];
  double riz = zcrd[atom_i] - zcrd[atom_k];
  double rjx = xcrd[atom_j] - xcrd[atom_k];
  double rjy = ycrd[atom_j] - ycrd[atom_k];
  double rjz = zcrd[atom_j] - zcrd[atom_k];
  double rlx = xcrd[atom_l] - xcrd[atom_k];
  double rly = ycrd[atom_l] - ycrd[atom_k];
  double rlz = zcrd[atom_l] - zcrd[atom_k];
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth);
  imageCoordinates(&rjx, &rjy, &rjz, umat, invu, unit_cell, imeth);
  imageCoordinates(&rlx, &rly, &rlz, umat, invu, unit_cell, imeth);

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
double dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l,
                      const CoordinateFrameReader &cfr) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                        cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const CoordinateFrame &cf) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, cf.data());
}

//-------------------------------------------------------------------------------------------------
double dihedral_angle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                      const PhaseSpace &ps) {
  return dihedral_angle(atom_i, atom_j, atom_k, atom_l, CoordinateFrameReader(ps));
}

} // namespace geometry
} // namespace omni
