// -*-c++-*-
namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcoord imageValue(const Tcoord x, const Tcalc range, const ImagingMethod style,
                  const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  Tcalc x_frac;
  if (tcoord_is_sgnint) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    x_frac = static_cast<Tcalc>(x) * inv_gpos_scale_factor / range;
  }
  else {
    x_frac = x / range;
  }
  switch (style) {
  case ImagingMethod::PRIMARY_UNIT_CELL:
    if (tcalc_is_double) {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floor(x_frac)) * range * gpos_scale_factor);
      }
      else {
        return ((x_frac - floor(x_frac)) * range);
      }
    }
    else {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floorf(x_frac)) * range * gpos_scale_factor);
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
      return llround((x_frac - (x_frac >= half)) * range * gpos_scale_factor);
    }
    else {
      return (x_frac - (x_frac >= half)) * range;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord *x, Tcoord *y, Tcoord *z, const Tcalc* umat, const Tcalc* invu,
                      const UnitCellType unit_cell, const ImagingMethod style,
                      const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        const Tcalc value_one = 1.0;
        const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
        local_x = static_cast<Tcalc>((*x) * inv_gpos_scale_factor) * umat[0];
        local_y = static_cast<Tcalc>((*y) * inv_gpos_scale_factor) * umat[4];
        local_z = static_cast<Tcalc>((*z) * inv_gpos_scale_factor) * umat[8];
      }
      else {
        local_x = (*x) * umat[0];
        local_y = (*y) * umat[4];
        local_z = (*z) * umat[8];
      }
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
          local_x -= ((local_x >= half) *  ceil(local_x - half)) +
                     ((local_x < -half) * floor(local_x + half));
          local_y -= ((local_y >= half) *  ceil(local_y - half)) +
                     ((local_y < -half) * floor(local_y + half));
          local_z -= ((local_z >= half) *  ceil(local_z - half)) +
                     ((local_z < -half) * floor(local_z + half));
        }
        else {
          local_x -= ((local_x >= half) *  ceilf(local_x - half)) +
                     ((local_x < -half) * floorf(local_x + half));
          local_y -= ((local_y >= half) *  ceilf(local_y - half)) +
                     ((local_y < -half) * floorf(local_y + half));
          local_z -= ((local_z >= half) *  ceilf(local_z - half)) +
                     ((local_z < -half) * floorf(local_z + half));
        }
        local_x -= (local_x >= half);
        local_y -= (local_y >= half);
        local_z -= (local_z >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        *x = llround(local_x * invu[0] * gpos_scale_factor);
        *y = llround(local_y * invu[4] * gpos_scale_factor);
        *z = llround(local_z * invu[8] * gpos_scale_factor);
      }
      else {
        *x = local_x * invu[0];
        *y = local_y * invu[4];
        *z = local_z * invu[8];
      }
    }
    break;
  case UnitCellType::TRICLINIC:
    {      
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        const Tcalc value_one = 1.0;
        const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
        local_x = static_cast<Tcalc>((*x) * inv_gpos_scale_factor);
        local_y = static_cast<Tcalc>((*y) * inv_gpos_scale_factor);
        local_z = static_cast<Tcalc>((*z) * inv_gpos_scale_factor);
      }
      else {
        local_x = (*x);
        local_y = (*y);
        local_z = (*z);
      }
      Tcalc ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      Tcalc ndy = (umat[1] * local_x) + (umat[4] * local_y) + (umat[7] * local_z);
      Tcalc ndz = (umat[2] * local_x) + (umat[5] * local_y) + (umat[8] * local_z);
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          ndx -= floor(ndx);
          ndy -= floor(ndy);
          ndz -= floor(ndz);
        }
        else {
          ndx -= floorf(ndx);
          ndy -= floorf(ndy);
          ndz -= floorf(ndz);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          ndx -= ((ndx >= half) * ceil(ndx - half)) + ((ndx < -half) * floor(ndx + half));
          ndy -= ((ndy >= half) * ceil(ndy - half)) + ((ndy < -half) * floor(ndy + half));
          ndz -= ((ndz >= half) * ceil(ndz - half)) + ((ndz < -half) * floor(ndz + half));
        }
        else {
          ndx -= ((ndx >= half) * ceilf(ndx - half)) + ((ndx < -half) * floorf(ndx + half));
          ndy -= ((ndy >= half) * ceilf(ndy - half)) + ((ndy < -half) * floorf(ndy + half));
          ndz -= ((ndz >= half) * ceilf(ndz - half)) + ((ndz < -half) * floorf(ndz + half));
        }
        ndx -= (ndx >= half);
        ndy -= (ndy >= half);
        ndz -= (ndz >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        *x = llround(((invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz)) * gpos_scale_factor);
        *y = llround(((invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz)) * gpos_scale_factor);
        *z = llround(((invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz)) * gpos_scale_factor);
      }
      else {
        *x = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
        *y = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
        *z = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord* x, Tcoord* y, Tcoord* z, const int length, const Tcalc* umat,
                      const Tcalc* invu, const UnitCellType unit_cell, const ImagingMethod style,
                      const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    for (int i = 0; i < length; i++) {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        local_x = static_cast<Tcalc>(x[i] * umat[0] * inv_gpos_scale_factor);
        local_y = static_cast<Tcalc>(y[i] * umat[4] * inv_gpos_scale_factor);
        local_z = static_cast<Tcalc>(z[i] * umat[8] * inv_gpos_scale_factor);
      }
      else {
        local_x = x[i] * umat[0];
        local_y = y[i] * umat[4];
        local_z = z[i] * umat[8];
      }
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
          local_x -= ((local_x >= half) *  ceil(local_x - half)) +
                     ((local_x < -half) * floor(local_x + half));
          local_y -= ((local_y >= half) *  ceil(local_y - half)) +
                     ((local_y < -half) * floor(local_y + half));
          local_z -= ((local_z >= half) *  ceil(local_z - half)) +
                     ((local_z < -half) * floor(local_z + half));
        }
        else {
          local_x -= ((local_x >= half) *  ceilf(local_x - half)) +
                     ((local_x < -half) * floorf(local_x + half));
          local_y -= ((local_y >= half) *  ceilf(local_y - half)) +
                     ((local_y < -half) * floorf(local_y + half));
          local_z -= ((local_z >= half) *  ceilf(local_z - half)) +
                     ((local_z < -half) * floorf(local_z + half));
        }
        local_x -= (local_x >= half);
        local_y -= (local_y >= half);
        local_z -= (local_z >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        x[i] = llround(local_x * invu[0] * gpos_scale_factor);
        y[i] = llround(local_y * invu[4] * gpos_scale_factor);
        z[i] = llround(local_z * invu[8] * gpos_scale_factor);
      }
      else {
        x[i] = local_x * invu[0];
        y[i] = local_y * invu[4];
        z[i] = local_z * invu[8];
      }
    }
    break;
  case UnitCellType::TRICLINIC:
    for (int i = 0; i < length; i++) {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        local_x = static_cast<Tcalc>(x[i] * inv_gpos_scale_factor);
        local_y = static_cast<Tcalc>(y[i] * inv_gpos_scale_factor);
        local_z = static_cast<Tcalc>(z[i] * inv_gpos_scale_factor);
      }
      else {
        local_x = x[i];
        local_y = y[i];
        local_z = z[i];
      }
      Tcalc ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      Tcalc ndy = (umat[1] * local_x) + (umat[4] * local_y) + (umat[7] * local_z);
      Tcalc ndz = (umat[2] * local_x) + (umat[5] * local_y) + (umat[8] * local_z);
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          ndx -= floor(ndx);
          ndy -= floor(ndy);
          ndz -= floor(ndz);
        }
        else {
          ndx -= floorf(ndx);
          ndy -= floorf(ndy);
          ndz -= floorf(ndz);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          ndx -= ((ndx >= half) * ceil(ndx - half)) + ((ndx < -half) * floor(ndx + half));
          ndy -= ((ndy >= half) * ceil(ndy - half)) + ((ndy < -half) * floor(ndy + half));
          ndz -= ((ndz >= half) * ceil(ndz - half)) + ((ndz < -half) * floor(ndz + half));
        }
        else {
          ndx -= ((ndx >= half) * ceilf(ndx - half)) + ((ndx < -half) * floorf(ndx + half));
          ndy -= ((ndy >= half) * ceilf(ndy - half)) + ((ndy < -half) * floorf(ndy + half));
          ndz -= ((ndz >= half) * ceilf(ndz - half)) + ((ndz < -half) * floorf(ndz + half));
        }
        ndx -= (ndx >= half);
        ndy -= (ndy >= half);
        ndz -= (ndz >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        x[i] = llround(((invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz)) * gpos_scale_factor);
        y[i] = llround(((invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz)) * gpos_scale_factor);
        z[i] = llround(((invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz)) * gpos_scale_factor);
      }
      else {
        x[i] = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
        y[i] = (invu[1] * ndx) + (invu[4] * ndy) + (invu[7] * ndz);
        z[i] = (invu[2] * ndx) + (invu[5] * ndy) + (invu[8] * ndz);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(std::vector<Tcoord> *x, std::vector<Tcoord> *y, std::vector<Tcoord> *z,
                      const Tcalc* umat, const Tcalc* invu, const UnitCellType unit_cell,
                      const ImagingMethod style, const Tcalc gpos_scale_factor) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style,
                   gpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Hybrid<Tcoord> *x, Hybrid<Tcoord> *y, Hybrid<Tcoord> *z,
                      const Tcalc* umat, const Tcalc* invu, const UnitCellType unit_cell,
                      const ImagingMethod style, const Tcalc gpos_scale_factor) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style,
                   gpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const Tcoord* xcrd, const Tcoord* ycrd,
               const Tcoord* zcrd, const Tcalc* umat, const Tcalc* invu,
               const UnitCellType unit_cell, const Tcalc gpos_scale_factor) {
  Tcoord dx = xcrd[atom_j] - xcrd[atom_i];
  Tcoord dy = ycrd[atom_j] - ycrd[atom_i];
  Tcoord dz = zcrd[atom_j] - zcrd[atom_i];
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE,
                   gpos_scale_factor);
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  if (tcoord_is_sgnint) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_gpos_factor = value_one / gpos_scale_factor;
    const Tcalc rdx = static_cast<Tcalc>(dx * inv_gpos_factor);
    const Tcalc rdy = static_cast<Tcalc>(dy * inv_gpos_factor);
    const Tcalc rdz = static_cast<Tcalc>(dz * inv_gpos_factor);
    return (tcalc_is_double) ? sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz)) :
                               sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));    
  }
  else {
    return (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                               sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k, const Tcoord* xcrd,
            const Tcoord* ycrd, const Tcoord* zcrd, const Tcalc* umat, const Tcalc* invu,
            const UnitCellType unit_cell, const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;

  // Image all three atoms and put the atom J at the origin
  Tcalc rix, riy, riz, rkx, rky, rkz;
  if (tcoord_is_sgnint) {
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    rix = static_cast<Tcalc>(xcrd[atom_i] - xcrd[atom_j]) * inv_gpos_scale_factor;
    riy = static_cast<Tcalc>(ycrd[atom_i] - ycrd[atom_j]) * inv_gpos_scale_factor;
    riz = static_cast<Tcalc>(zcrd[atom_i] - zcrd[atom_j]) * inv_gpos_scale_factor;
    rkx = static_cast<Tcalc>(xcrd[atom_k] - xcrd[atom_j]) * inv_gpos_scale_factor;
    rky = static_cast<Tcalc>(ycrd[atom_k] - ycrd[atom_j]) * inv_gpos_scale_factor;
    rkz = static_cast<Tcalc>(zcrd[atom_k] - zcrd[atom_j]) * inv_gpos_scale_factor;
  }
  else {
    rix = xcrd[atom_i] - xcrd[atom_j];
    riy = ycrd[atom_i] - ycrd[atom_j];
    riz = zcrd[atom_i] - zcrd[atom_j];
    rkx = xcrd[atom_k] - xcrd[atom_j];
    rky = ycrd[atom_k] - ycrd[atom_j];
    rkz = zcrd[atom_k] - zcrd[atom_j];
  }
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rkx, &rky, &rkz, umat, invu, unit_cell, imeth, gpos_scale_factor);

  // Compute the angle, in radians
  const Tcalc mgba = (rix * rix) + (riy * riy) + (riz * riz);
  const Tcalc mgbc = (rkx * rkx) + (rky * rky) + (rkz * rkz);
  const Tcalc invbabc = (tcalc_is_double) ? value_one / sqrt(mgba * mgbc) :
                                            value_one / sqrtf(mgba * mgbc);
  Tcalc costheta = ((rix * rkx) + (riy * rky) + (riz * rkz)) * invbabc;
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  return (tcalc_is_double) ? acos(costheta) : acosf(costheta);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc dihedral_angle(int atom_i, int atom_j, int atom_k, int atom_l, const Tcoord* xcrd,
                     const Tcoord* ycrd, const Tcoord* zcrd, const Tcalc* umat,
                     const Tcalc* invu, const UnitCellType unit_cell,
                     const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;

  // Image all four atoms and put atom K at the origin
  Tcalc rix, riy, riz, rjx, rjy, rjz, rlx, rly, rlz;
  if (tcoord_is_sgnint) {
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    rix = static_cast<Tcalc>(xcrd[atom_i] - xcrd[atom_k]) * inv_gpos_scale_factor;
    riy = static_cast<Tcalc>(ycrd[atom_i] - ycrd[atom_k]) * inv_gpos_scale_factor;
    riz = static_cast<Tcalc>(zcrd[atom_i] - zcrd[atom_k]) * inv_gpos_scale_factor;
    rjx = static_cast<Tcalc>(xcrd[atom_j] - xcrd[atom_k]) * inv_gpos_scale_factor;
    rjy = static_cast<Tcalc>(ycrd[atom_j] - ycrd[atom_k]) * inv_gpos_scale_factor;
    rjz = static_cast<Tcalc>(zcrd[atom_j] - zcrd[atom_k]) * inv_gpos_scale_factor;
    rlx = static_cast<Tcalc>(xcrd[atom_l] - xcrd[atom_k]) * inv_gpos_scale_factor;
    rly = static_cast<Tcalc>(ycrd[atom_l] - ycrd[atom_k]) * inv_gpos_scale_factor;
    rlz = static_cast<Tcalc>(zcrd[atom_l] - zcrd[atom_k]) * inv_gpos_scale_factor;
  }
  else {
    rix = xcrd[atom_i] - xcrd[atom_k];
    riy = ycrd[atom_i] - ycrd[atom_k];
    riz = zcrd[atom_i] - zcrd[atom_k];
    rjx = xcrd[atom_j] - xcrd[atom_k];
    rjy = ycrd[atom_j] - ycrd[atom_k];
    rjz = zcrd[atom_j] - zcrd[atom_k];
    rlx = xcrd[atom_l] - xcrd[atom_k];
    rly = ycrd[atom_l] - ycrd[atom_k];
    rlz = zcrd[atom_l] - zcrd[atom_k];
  }
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rjx, &rjy, &rjz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rlx, &rly, &rlz, umat, invu, unit_cell, imeth, gpos_scale_factor);

  // Compute the dihedral angle, in radians
  Tcalc ab[3], bc[3], cd[3];
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
  Tcalc crabbc[3], crbccd[3], scr[3];
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  if (tcalc_is_double) {
    costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    costheta /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                      (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  return angleVerification(costheta, crabbc, crbccd, bc, scr);
}

} // namespace structure
} // namespace stormm
