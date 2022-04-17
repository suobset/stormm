// -*-c++-*-
namespace omni {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const Tcalc rotation_angle,
                     const Tcalc globalpos_scale_factor) {
  const int natom = moving_atoms.size();

  // Center and image the coordinates
  const Tcoord center_x = xcrd[atom_j];
  const Tcoord center_y = ycrd[atom_j];
  const Tcoord center_z = zcrd[atom_j];
  xcrd[atom_i] -= center_x;
  ycrd[atom_i] -= center_y;
  zcrd[atom_i] -= center_z;
  xcrd[atom_j] = static_cast<Tcoord>(0.0);
  ycrd[atom_j] = static_cast<Tcoord>(0.0);
  zcrd[atom_j] = static_cast<Tcoord>(0.0);
  for (int i = 0; i < natom; i++) {
    const int mk = moving_atoms[i];
    xcrd[mk] -= center_x;
    ycrd[mk] -= center_y;
    zcrd[mk] -= center_z;
  }
  
  // Define the vector of rotation, then the matrix.  The coordinate data type is assumed to be
  // either a real numbered type (float, double), or a signed integral type (int, long long int),
  // a fact which will be enforced by the CoordinateSeries object, the only templated coordinate
  // representation.  All other coordinate representations work in double or long long int.
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  Tcalc dx, dy, dz, invdr, cos_ra, sin_ra;
  const Tcalc value_one = 1.0;
  const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;  
  if (tcoord_is_sgnint) {
    dx = static_cast<Tcalc>(xcrd[atom_j] - xcrd[atom_i]) * inv_globalpos_scale_factor;
    dy = static_cast<Tcalc>(ycrd[atom_j] - ycrd[atom_i]) * inv_globalpos_scale_factor;
    dz = static_cast<Tcalc>(zcrd[atom_j] - zcrd[atom_i]) * inv_globalpos_scale_factor;
  }
  else {
    dx = xcrd[atom_j] - xcrd[atom_i];
    dy = ycrd[atom_j] - ycrd[atom_i];
    dz = zcrd[atom_j] - zcrd[atom_i];
  }
  if (tcalc_is_double) {
    invdr = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
    cos_ra = cos(rotation_angle);
    sin_ra = sin(rotation_angle);
  }
  else {
    invdr = 1.0f / sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    cos_ra = cosf(rotation_angle);
    sin_ra = sinf(rotation_angle);
  }
  dx *= invdr;
  dy *= invdr;
  dz *= invdr;
  const std::vector<Tcalc> rmat = { cos_ra + (dx * dx * (value_one - cos_ra)),
                                    (dy * dx * (1.0 - cos_ra)) + (dz * sin_ra),
                                    (dz * dx * (1.0 - cos_ra)) - (dy * sin_ra),
                                    (dx * dy * (1.0 - cos_ra)) - (dz * sin_ra),
                                    cos_ra + (dy * dy * (value_one - cos_ra)),
                                    (dz * dy * (1.0 - cos_ra)) + (dx * sin_ra),
                                    (dx * dz * (1.0 - cos_ra)) + (dy * sin_ra),
                                    (dy * dz * (1.0 - cos_ra)) - (dx * sin_ra),
                                    cos_ra + (dz * dz * (value_one - cos_ra)) };
  
  // Restore the original imaging and location of the bond atoms
  xcrd[atom_j] = center_x;
  ycrd[atom_j] = center_y;
  zcrd[atom_j] = center_z;
  xcrd[atom_i] += center_x;
  ycrd[atom_i] += center_y;
  zcrd[atom_i] += center_z;

  // Loop over all moving particles, rotate about the vector, and re-apply their original
  // translations relative to the origin.
  for (int i = 0; i < natom; i++) {
    const int mk = moving_atoms[i];
    const Tcalc nx = (rmat[0] * xcrd[mk]) + (rmat[3] * ycrd[mk]) + (rmat[6] * zcrd[mk]);
    const Tcalc ny = (rmat[1] * xcrd[mk]) + (rmat[4] * ycrd[mk]) + (rmat[7] * zcrd[mk]);
    const Tcalc nz = (rmat[2] * xcrd[mk]) + (rmat[5] * ycrd[mk]) + (rmat[8] * zcrd[mk]);
    if (tcoord_is_sgnint) {
      xcrd[mk] = llround(nx * globalpos_scale_factor);
      ycrd[mk] = llround(ny * globalpos_scale_factor);
      zcrd[mk] = llround(nz * globalpos_scale_factor);
    }
    else {
      xcrd[mk] = nx;
      ycrd[mk] = ny;
      zcrd[mk] = nz;
    }
    xcrd[mk] += center_x;
    ycrd[mk] += center_y;
    zcrd[mk] += center_z;    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeries<Tcoord> *cs, const int frame_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const Tcalc rotation_angle) {
  rotateAboutBond<Tcoord, Tcalc>(cs->data(), frame_index, atom_i, atom_j, moving_atoms,
                                 rotation_angle);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeriesWriter<Tcoord> csw, const int frame_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const Tcalc rotation_angle) {
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  rotateAboutBond<Tcoord, Tcalc>(&csw.xcrd[fidx_zu * natom_zu], &csw.ycrd[fidx_zu * natom_zu],
                                 &csw.zcrd[fidx_zu * natom_zu], atom_i, atom_j, moving_atoms,
                                 rotation_angle, csw.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups,
                      const Tcalc globalpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;  
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      // Unpack the appropriate inversion group
      const int root_a = inversion_groups[center_idx].root_atom;
      const int root_b = inversion_groups[center_idx].pivot_atom;
      const int ccen   = chiral_centers[center_idx];

      // Find the bisector of the root_a : chiral_center : root_b angle.  Shift the root_b atom to
      // lie along the line of the bisector, rotate the moving atoms 180 degrees about this "bond,"
      // then replace the root_b atom.
      const Tcoord orig_bx = xcrd[root_b];
      const Tcoord orig_by = ycrd[root_b];
      const Tcoord orig_bz = zcrd[root_b];
      Tcalc dax, day, daz, dbx, dby, dbz, ccenx, cceny, ccenz, invra, invrb;
      if (tcoord_is_sgnint) {
        dax = static_cast<Tcalc>(xcrd[root_a] - xcrd[ccen]) * inv_globalpos_scale_factor;
        day = static_cast<Tcalc>(ycrd[root_a] - ycrd[ccen]) * inv_globalpos_scale_factor;
        daz = static_cast<Tcalc>(zcrd[root_a] - zcrd[ccen]) * inv_globalpos_scale_factor;
        dbx = static_cast<Tcalc>(orig_bx - xcrd[ccen]) * inv_globalpos_scale_factor;
        dby = static_cast<Tcalc>(orig_by - ycrd[ccen]) * inv_globalpos_scale_factor;
        dbz = static_cast<Tcalc>(orig_bz - zcrd[ccen]) * inv_globalpos_scale_factor;
        ccenx = static_cast<Tcalc>(xcrd[ccen]) * inv_globalpos_scale_factor;
        cceny = static_cast<Tcalc>(ycrd[ccen]) * inv_globalpos_scale_factor;
        ccenz = static_cast<Tcalc>(zcrd[ccen]) * inv_globalpos_scale_factor;
      }
      else {

        // Make no assumptions about the precision of the calculation relative to the precision
        // of the coordinates.  It is possible, and even common in MD programs, to store
        // coordinates in a higher degree of precision than the calculation.  It would be faster,
        // otherwise, to set ccen(x,y,z) at the front and re-use the values in computing da or db.
        dax = xcrd[root_a] - xcrd[ccen];
        day = ycrd[root_a] - ycrd[ccen];
        daz = zcrd[root_a] - zcrd[ccen];
        dbx = orig_bx - xcrd[ccen];
        dby = orig_by - ycrd[ccen];
        dbz = orig_bz - zcrd[ccen];
        ccenx = xcrd[ccen];
        cceny = ycrd[ccen];
        ccenz = zcrd[ccen];
      }
      if (tcalc_is_double) {
        invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
        invra = 1.0 / sqrt((dax * dax) + (day * day) + (daz * daz));
      }
      else {
        invrb = 1.0f / sqrtf((dbx * dbx) + (dby * dby) + (dbz * dbz));
        invra = 1.0f / sqrtf((dax * dax) + (day * day) + (daz * daz));
      }
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      dax = ccenx + (dax * invra);
      day = cceny + (day * invra);
      daz = ccenz + (daz * invra);
      const Tcalc midpoint_x = 0.5 * (dbx + dax);
      const Tcalc midpoint_y = 0.5 * (dby + day);
      const Tcalc midpoint_z = 0.5 * (dbz + daz);
      if (tcoord_is_sgnint) {
        xcrd[root_b] = llround(midpoint_x * globalpos_scale_factor);
        ycrd[root_b] = llround(midpoint_y * globalpos_scale_factor);
        zcrd[root_b] = llround(midpoint_z * globalpos_scale_factor);
      }
      else {
        xcrd[root_b] = midpoint_x;
        ycrd[root_b] = midpoint_y;
        zcrd[root_b] = midpoint_z;
      }
      rotateAboutBond(xcrd, ycrd, zcrd, root_b, chiral_centers[center_idx],
                      inversion_groups[center_idx].rotatable_atoms,
                      static_cast<Tcalc>(symbols::pi), globalpos_scale_factor);
      xcrd[root_b] = orig_bx;
      ycrd[root_b] = orig_by;
      zcrd[root_b] = orig_bz;
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      // Find the molecule home of the present center and flip those atoms only.
      const int natom = inversion_groups[center_idx].rotatable_atoms.size();
      for (int i = 0; i < natom; i++) {
        const int atom_idx = inversion_groups[center_idx].rotatable_atoms[i];
        xcrd[atom_idx] = -xcrd[atom_idx];
      }

      // All centers have been flipped by the reflection.  Loop over all other chiral centers and
      // perform chiral inversions in order to flip the others back, where possible.  Infinite
      // recursion is limited by the fact that at most one chiral center in a molecule can be given
      // the protocol "REFLECT," and reflection only triggers subsequent rotations.
      const int nchirals = chiral_protocols.size();
      for (int i = 0; i < nchirals; i++) {
        if (chiral_protocols[i] == ChiralInversionProtocol::ROTATE) {
          flipChiralCenter(xcrd, ycrd, zcrd, i, chiral_centers, chiral_protocols,
                           inversion_groups, globalpos_scale_factor);
        }
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeries<Tcoord> *cs, const int frame_index,
                      const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter<Tcoord, Tcalc>(cs->data(), frame_index, center_idx, chiral_centers,
                                  chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeriesWriter<Tcoord> csw, const int frame_index,
                      const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  flipChiralCenter<Tcoord, Tcalc>(&csw.xcrd[fidx_zu * natom_zu], &csw.ycrd[fidx_zu * natom_zu],
                                  &csw.zcrd[fidx_zu * natom_zu], center_idx, chiral_centers,
                                  chiral_protocols, inversion_groups, csw.gpos_scale);
}
  
} // namespace structure
} // namespace omni
