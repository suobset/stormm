// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicStretch(const int i_atom, const int j_atom, const Tcalc stiffness,
                          const Tcalc equilibrium, const Tcoord* xcrd, const Tcoord* ycrd,
                          const Tcoord* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                          const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                          const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[j_atom] - xcrd[i_atom];
    dy = ycrd[j_atom] - ycrd[i_atom];
    dz = zcrd[j_atom] - zcrd[i_atom];
  }
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const Tcalc dr = (tcalc_ct == double_type_index) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                                     sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  const Tcalc dl = dr - equilibrium;

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc fmag = 2.0 * stiffness * dl / dr;
    const Tcalc fmag_dx = fmag * dx;
    const Tcalc fmag_dy = fmag * dy;
    const Tcalc fmag_dz = fmag * dz;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifmag_dx = llround(fmag_dx * force_factor);
      const Tforce ifmag_dy = llround(fmag_dy * force_factor);
      const Tforce ifmag_dz = llround(fmag_dz * force_factor);
      xfrc[i_atom] += ifmag_dx;
      yfrc[i_atom] += ifmag_dy;
      zfrc[i_atom] += ifmag_dz;
      xfrc[j_atom] -= ifmag_dx;
      yfrc[j_atom] -= ifmag_dy;
      zfrc[j_atom] -= ifmag_dz;
    }
    else {
      xfrc[i_atom] += fmag_dx;
      yfrc[i_atom] += fmag_dy;
      zfrc[i_atom] += fmag_dz;
      xfrc[j_atom] -= fmag_dx;
      yfrc[j_atom] -= fmag_dy;
      zfrc[j_atom] -= fmag_dz;
    }
  }
  return stiffness * dl * dl;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicBend(const int i_atom, const int j_atom, const int k_atom,
                       const Tcalc stiffness, const Tcalc equilibrium, const Tcoord* xcrd,
                       const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                       Tforce* yfrc, Tforce* zfrc, const EvaluateForce eval_force,
                       const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;
  
  // Compute displacements
  Tcalc ba[3], bc[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ba[0] = static_cast<Tcalc>(xcrd[i_atom] - xcrd[j_atom]) * inv_gpos_factor;
    ba[1] = static_cast<Tcalc>(ycrd[i_atom] - ycrd[j_atom]) * inv_gpos_factor;
    ba[2] = static_cast<Tcalc>(zcrd[i_atom] - zcrd[j_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
  }
  else {
    ba[0] = xcrd[i_atom] - xcrd[j_atom];
    ba[1] = ycrd[i_atom] - ycrd[j_atom];
    ba[2] = zcrd[i_atom] - zcrd[j_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
  }
  imageCoordinates(&ba[0], &ba[1], &ba[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // On to the angle force computation
  const Tcalc mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
  const Tcalc mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
  const Tcalc invbabc = (tcalc_is_double) ? 1.0 / sqrt(mgba * mgbc) :
                                            value_one / sqrtf(mgba * mgbc);
  Tcalc costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  const Tcalc theta = (tcalc_is_double) ? acos(costheta) : acosf(costheta);
  const Tcalc dtheta = theta - equilibrium;

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc dA = (tcalc_is_double) ?
                     -2.0 * stiffness * dtheta / sqrt(1.0 - (costheta * costheta)) :
                     -2.0f * stiffness * dtheta / sqrtf(value_one - (costheta * costheta));
    const Tcalc sqba = dA / mgba;
    const Tcalc sqbc = dA / mgbc;
    const Tcalc mbabc = dA * invbabc;
    if (isSignedIntegralScalarType<Tforce>()) {
      Tforce iadf[3], icdf[3];
      for (int i = 0; i < 3; i++) {
        iadf[i] = llround(((bc[i] * mbabc) - (costheta * ba[i] * sqba)) * force_factor);
        icdf[i] = llround(((ba[i] * mbabc) - (costheta * bc[i] * sqbc)) * force_factor);
      }      
      xfrc[i_atom] -= iadf[0];
      yfrc[i_atom] -= iadf[1];
      zfrc[i_atom] -= iadf[2];
      xfrc[j_atom] += iadf[0] + icdf[0];
      yfrc[j_atom] += iadf[1] + icdf[1];
      zfrc[j_atom] += iadf[2] + icdf[2];
      xfrc[k_atom] -= icdf[0];
      yfrc[k_atom] -= icdf[1];
      zfrc[k_atom] -= icdf[2];
    }
    else {
      Tcalc adf[3], cdf[3];
      for (int i = 0; i < 3; i++) {
        adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
        cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
      }
      xfrc[i_atom] -= adf[0];
      yfrc[i_atom] -= adf[1];
      zfrc[i_atom] -= adf[2];
      xfrc[j_atom] += adf[0] + cdf[0];
      yfrc[j_atom] += adf[1] + cdf[1];
      zfrc[j_atom] += adf[2] + cdf[2];
      xfrc[k_atom] -= cdf[0];
      yfrc[k_atom] -= cdf[1];
      zfrc[k_atom] -= cdf[2];
    }
  }
  return stiffness * dtheta * dtheta;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalDihedralTwist(const int i_atom, const int j_atom, const int k_atom, const int l_atom,
                        const Tcalc amplitude, const Tcalc phase_angle, const Tcalc frequency,
                        const DihedralStyle kind, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                        const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;

  // Compute displacements
  Tcalc ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ab[0] = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    ab[1] = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    ab[2] = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
    cd[0] = static_cast<Tcalc>(xcrd[l_atom] - xcrd[k_atom]) * inv_gpos_factor;
    cd[1] = static_cast<Tcalc>(ycrd[l_atom] - ycrd[k_atom]) * inv_gpos_factor;
    cd[2] = static_cast<Tcalc>(zcrd[l_atom] - zcrd[k_atom]) * inv_gpos_factor;
  }
  else {
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
  }
  imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);

  // Compute cross products and then the angle between the planes
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
  Tcalc theta;
  if (tcalc_is_double) {
    theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ? acos(costheta) : -acos(costheta);
  }
  else {
    theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0f) ?  acosf(costheta) :
                                                                  -acosf(costheta);
  }
  Tcalc sangle;
  switch (kind) {
  case DihedralStyle::COSINE:
    sangle = (frequency * theta) - phase_angle;
    break;
  case DihedralStyle::HARMONIC:
    sangle = theta - phase_angle;
    break;
  }

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {
    Tcalc fr;
    switch (kind) {
    case DihedralStyle::COSINE:
      fr = amplitude * frequency * sin(sangle);
      break;
    case DihedralStyle::HARMONIC:
      fr = -2.0 * amplitude * sangle;
      break;
    }
    Tcalc mgab, mgbc, mgcd;
    if (tcalc_is_double) {
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    else {
      mgab = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrtf(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrtf(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    const Tcalc invab = value_one / mgab;
    const Tcalc invbc = value_one / mgbc;
    const Tcalc invcd = value_one / mgcd;
    const Tcalc cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const Tcalc cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    Tcalc isinb2, isinc2;
    if (tcalc_is_double) {
      isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
      isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
    }
    else {
      isinb2 = (cosb * cosb < asymptotic_to_one_f) ?
               fr / (value_one - (cosb * cosb)) : fr * inverse_one_minus_asymptote_f;
      isinc2 = (cosc * cosc < asymptotic_to_one_f) ?
               fr / (value_one - (cosc * cosc)) : fr * inverse_one_minus_asymptote_f;
    }
    const Tcalc invabc = invab * invbc;
    const Tcalc invbcd = invbc * invcd;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
    }

    // Transform the rotational derivatives to cartesian coordinates
    const Tcalc fa = -invab * isinb2;
    const Tcalc fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
    const Tcalc fb2 = cosc * invbc * isinc2;
    const Tcalc fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
    const Tcalc fc2 = cosb * invbc * isinb2;
    const Tcalc fd = -invcd * isinc2;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifrc_ix = llround(crabbc[0] * fa * force_factor);
      const Tforce ifrc_jx = llround(((fb1 * crabbc[0]) - (fb2 * crbccd[0])) * force_factor);
      const Tforce ifrc_lx = llround(-fd * crbccd[0] * force_factor);
      xfrc[i_atom] += ifrc_ix;
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] -= ifrc_ix + ifrc_jx + ifrc_lx;
      xfrc[l_atom] += ifrc_lx;
      const Tforce ifrc_iy = llround(crabbc[1] * fa * force_factor);
      const Tforce ifrc_jy = llround(((fb1 * crabbc[1]) - (fb2 * crbccd[1])) * force_factor);
      const Tforce ifrc_ly = llround(-fd * crbccd[1] * force_factor);
      yfrc[i_atom] += ifrc_iy;
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] -= ifrc_iy + ifrc_jy + ifrc_ly;
      yfrc[l_atom] += ifrc_ly;
      const Tforce ifrc_iz = llround(crabbc[2] * fa * force_factor);
      const Tforce ifrc_jz = llround(((fb1 * crabbc[2]) - (fb2 * crbccd[2])) * force_factor);
      const Tforce ifrc_lz = llround(-fd * crbccd[2] * force_factor);
      zfrc[i_atom] += ifrc_iz;
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] -= ifrc_iz + ifrc_jz + ifrc_lz;
      zfrc[l_atom] += ifrc_lz;
    }
    else {
      const Tforce frc_ix = crabbc[0] * fa;
      const Tforce frc_jx = (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      const Tforce frc_lx = -fd * crbccd[0];
      xfrc[i_atom] += frc_ix;
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] -= frc_ix + frc_jx + frc_lx;
      xfrc[l_atom] += frc_lx;
      const Tforce frc_iy = crabbc[1] * fa;
      const Tforce frc_jy = (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      const Tforce frc_ly = -fd * crbccd[1];
      yfrc[i_atom] += frc_iy;
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] -= frc_iy + frc_jy + frc_ly;
      yfrc[l_atom] += frc_ly;
      const Tforce frc_iz = crabbc[2] * fa;
      const Tforce frc_jz = (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      const Tforce frc_lz = -fd * crbccd[2];
      zfrc[i_atom] += frc_iz;
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] -= frc_iz + frc_jz + frc_lz;
      zfrc[l_atom] += frc_lz;
    }
  }
  switch (kind) {
  case DihedralStyle::COSINE:
    return (tcalc_is_double) ? amplitude * (1.0 + cos(sangle)) :
                               amplitude * (value_one + cosf(sangle));
  case DihedralStyle::HARMONIC:
    return amplitude * sangle * sangle;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalCmap(const Tcalc* cmap_patches, const int* cmap_patch_bounds, const int surf_idx,
               const int surf_dim, const int i_atom, const int j_atom, const int k_atom,
               const int l_atom, const int m_atom, const Tcoord* xcrd, const Tcoord* ycrd,
               const Tcoord* zcrd, const double* umat, const double* invu,
               const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
               const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
               const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;

  // Compute displacements
  std::vector<Tcalc> acoef(16);
  Tcalc ab[3], bc[3], cd[3], de[3], crabbc[3], crbccd[3], crcdde[3], scr_phi[3], scr_psi[3];
  Tcalc phi_progression[4], psi_progression[4], acoef_psi[4];
  for (int i = 0; i < 4; i++) {
    acoef_psi[i] = 0.0;
  }
  if (isSignedIntegralScalarType<Tcoord>()) {
    ab[0] = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    ab[1] = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    ab[2] = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
    cd[0] = static_cast<Tcalc>(xcrd[l_atom] - xcrd[k_atom]) * inv_gpos_factor;
    cd[1] = static_cast<Tcalc>(ycrd[l_atom] - ycrd[k_atom]) * inv_gpos_factor;
    cd[2] = static_cast<Tcalc>(zcrd[l_atom] - zcrd[k_atom]) * inv_gpos_factor;
    de[0] = static_cast<Tcalc>(xcrd[m_atom] - xcrd[l_atom]) * inv_gpos_factor;
    de[1] = static_cast<Tcalc>(ycrd[m_atom] - ycrd[l_atom]) * inv_gpos_factor;
    de[2] = static_cast<Tcalc>(zcrd[m_atom] - zcrd[l_atom]) * inv_gpos_factor;
  }
  else {
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
    de[0] = xcrd[m_atom] - xcrd[l_atom];
    de[1] = ycrd[m_atom] - ycrd[l_atom];
    de[2] = zcrd[m_atom] - zcrd[l_atom];
  }
  imageCoordinates(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates(&de[0], &de[1], &de[2], umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
    
  // Compute the first dihedral
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc cos_phi = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  if (tcalc_is_double) {
    cos_phi /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                    (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    cos_phi /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr_phi);
  cos_phi = (cos_phi < -value_one) ? -value_one : (cos_phi > value_one) ? value_one : cos_phi;
  Tcalc phi;
  if (tcalc_is_double) {
    phi = (scr_phi[0]*bc[0] + scr_phi[1]*bc[1] + scr_phi[2]*bc[2] > 0.0) ?  acos(cos_phi) :
                                                                           -acos(cos_phi);
    phi += pi;
    phi = (phi < 0.0) ? phi + twopi : (phi >= twopi) ? phi - twopi : phi;
  }
  else {
    phi = (scr_phi[0]*bc[0] + scr_phi[1]*bc[1] + scr_phi[2]*bc[2] > 0.0f) ?  acosf(cos_phi) :
                                                                            -acosf(cos_phi);
    phi += pi_f;
    phi = (phi < 0.0f) ? phi + twopi_f : (phi >= twopi_f) ? phi - twopi_f : phi;
  }

  // Compute the second dihedral
  crossProduct(cd, de, crcdde);
  Tcalc cos_psi = crbccd[0]*crcdde[0] + crbccd[1]*crcdde[1] + crbccd[2]*crcdde[2];
  if (tcalc_is_double) {
    cos_psi /= sqrt((crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]) *
                    (crcdde[0]*crcdde[0] + crcdde[1]*crcdde[1] + crcdde[2]*crcdde[2]));
  }
  else {
    cos_psi /= sqrtf((crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]) *
                     (crcdde[0]*crcdde[0] + crcdde[1]*crcdde[1] + crcdde[2]*crcdde[2]));
  }
  crossProduct(crbccd, crcdde, scr_psi);
  cos_psi = (cos_psi < -value_one) ? -value_one : (cos_psi > value_one) ? value_one : cos_psi;
  Tcalc psi;
  if (tcalc_is_double) {
    psi = (scr_psi[0]*bc[0] + scr_psi[1]*bc[1] + scr_psi[2]*bc[2] > 0.0) ?  acos(cos_psi) :
                                                                           -acos(cos_psi);
    psi += pi;
    psi = (psi < 0.0) ? psi + twopi : (psi >= twopi) ? psi - twopi : psi;
  }
  else {
    psi = (scr_psi[0]*bc[0] + scr_psi[1]*bc[1] + scr_psi[2]*bc[2] > 0.0f) ?  acosf(cos_psi) :
                                                                            -acosf(cos_psi);
    psi += pi_f;
    psi = (psi < 0.0f) ? psi + twopi_f : (psi >= twopi_f) ? psi - twopi_f : psi;
  }
  
  // Compute the patch index (idx_phi, idx_psi) of the CMAP
  const Tcalc dsurf_dim = static_cast<Tcalc>(surf_dim);
  Tcalc phi_grid, psi_grid;
  if (tcalc_is_double) {
    phi_grid = phi * dsurf_dim * inverse_twopi;
    psi_grid = psi * dsurf_dim * inverse_twopi;
  }
  else {
    phi_grid = phi * dsurf_dim * inverse_twopi_f;
    psi_grid = psi * dsurf_dim * inverse_twopi_f;
  }
  const int idx_phi = phi_grid;
  const int idx_psi = psi_grid;
  const Tcalc phifrac = phi_grid - idx_phi;
  const Tcalc psifrac = psi_grid - idx_psi;

  // Draw in the matrix of spline values and derivatives
  const int patch_idx = cmap_patch_bounds[surf_idx] + (((idx_psi * surf_dim) + idx_phi) * 16);
  for (int i = 0; i < 16; i++) {
    acoef[i] = cmap_patches[patch_idx + i];
  }

  // Perform the matrix multiplications to obtain the bicubic spline coefficients
  phi_progression[0] = 1.0;
  psi_progression[0] = 1.0;
  for (int i = 1; i < 4; i++) {
    phi_progression[i] = phi_progression[i - 1] * phifrac;
    psi_progression[i] = psi_progression[i - 1] * psifrac;
  }
  matrixVectorMultiply(acoef.data(), psi_progression, acoef_psi, 4, 4, 1.0, 1.0, 0.0);
  const Tcalc contrib = (phi_progression[0] * acoef_psi[0]) + (phi_progression[1] * acoef_psi[1]) +
                        (phi_progression[2] * acoef_psi[2]) + (phi_progression[3] * acoef_psi[3]);

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {

    // The derivatives along phi and psi follow from the energy expression, evaluation of the
    // guiding matrix equation for the bicubic spline potential:
    //
    //                                                    [ a00 a01 a02 a03     [ 1
    // Energy = [ 1 + phifrac + phifrac^2 + phifrac^3 ] *   a10 a11 a12 a13  *    psifrac
    //                                                      a20 a21 a22 a23       psifrac^2
    //                                                      a30 a31 a32 a33 ]     psifrac^3 ]
    Tcalc dphi, dpsi, mgab, mgbc, mgcd, mgde;
    if (tcalc_is_double) {
      dphi  = (((3.0 * acoef[15] * phifrac) + (2.0 * acoef[14])) * phifrac) + acoef[13];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[11] * phifrac) + (2.0 * acoef[10])) * phifrac) + acoef[ 9];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[ 7] * phifrac) + (2.0 * acoef[ 6])) * phifrac) + acoef[ 5];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[ 3] * phifrac) + (2.0 * acoef[ 2])) * phifrac) + acoef[ 1];
      dpsi  = (((3.0 * acoef[15] * psifrac) + (2.0 * acoef[11])) * psifrac) + acoef[ 7];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[14] * psifrac) + (2.0 * acoef[10])) * psifrac) + acoef[ 6];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[13] * psifrac) + (2.0 * acoef[ 9])) * psifrac) + acoef[ 5];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[12] * psifrac) + (2.0 * acoef[ 8])) * psifrac) + acoef[ 4];
      dphi *= dsurf_dim / twopi;
      dpsi *= dsurf_dim / twopi;
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      mgde = sqrt(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
    }
    else {
      dphi  = (((3.0f * acoef[15] * phifrac) + (2.0f * acoef[14])) * phifrac) + acoef[13];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[11] * phifrac) + (2.0f * acoef[10])) * phifrac) + acoef[ 9];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[ 7] * phifrac) + (2.0f * acoef[ 6])) * phifrac) + acoef[ 5];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[ 3] * phifrac) + (2.0f * acoef[ 2])) * phifrac) + acoef[ 1];
      dpsi  = (((3.0f * acoef[15] * psifrac) + (2.0f * acoef[11])) * psifrac) + acoef[ 7];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[14] * psifrac) + (2.0f * acoef[10])) * psifrac) + acoef[ 6];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[13] * psifrac) + (2.0f * acoef[ 9])) * psifrac) + acoef[ 5];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[12] * psifrac) + (2.0f * acoef[ 8])) * psifrac) + acoef[ 4];
      dphi *= dsurf_dim / twopi;
      dpsi *= dsurf_dim / twopi;
      mgab = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrtf(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrtf(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      mgde = sqrtf(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
    }

    // With the derivative in hand, evaluate the transformation of coordinates for either the
    // phi or psi dihedrals.
    const Tcalc invab = value_one / mgab;
    const Tcalc invbc = value_one / mgbc;
    const Tcalc invcd = value_one / mgcd;
    const Tcalc invde = value_one / mgde;
    const Tcalc invabc = invab * invbc;
    const Tcalc invbcd = invbc * invcd;
    const Tcalc invcde = invcd * invde;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
      crcdde[i] *= invcde;
    }

    // Feed the gradient, negative of the derivative, into the functions below
    dphi *= -value_one;
    dpsi *= -value_one;

    // Phi accumulation: transform the rotational derivatives to cartesian coordinates
    const Tcalc phi_cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const Tcalc phi_cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    Tcalc phi_isinb2, phi_isinc2;
    if (tcalc_is_double) {
      phi_isinb2 = (phi_cosb * phi_cosb < asymptotic_to_one_lf) ?
                   dphi / (1.0 - (phi_cosb * phi_cosb)) : dphi * inverse_one_minus_asymptote_lf;
      phi_isinc2 = (phi_cosc * phi_cosc < asymptotic_to_one_lf) ?
                   dphi / (1.0 - (phi_cosc * phi_cosc)) : dphi * inverse_one_minus_asymptote_lf;
    }
    else {
      phi_isinb2 = (phi_cosb * phi_cosb < asymptotic_to_one_f) ?
                   dphi / (value_one - (phi_cosb * phi_cosb)) :
                   dphi * inverse_one_minus_asymptote_f;
      phi_isinc2 = (phi_cosc * phi_cosc < asymptotic_to_one_f) ?
                   dphi / (value_one - (phi_cosc * phi_cosc)) :
                   dphi * inverse_one_minus_asymptote_f;      
    }
    const Tcalc phi_fa = -invab * phi_isinb2;
    const Tcalc phi_fb1 = (mgbc - (mgab * phi_cosb)) * invabc * phi_isinb2;
    const Tcalc phi_fb2 = phi_cosc * invbc * phi_isinc2;
    const Tcalc phi_fc1 = (mgbc - (mgcd * phi_cosc)) * invbcd * phi_isinc2;
    const Tcalc phi_fc2 = phi_cosb * invbc * phi_isinb2;
    const Tcalc phi_fd = -invcd * phi_isinc2;

    // Psi accumulation: transform the rotational derivatives to cartesian coordinates
    const Tcalc psi_cosb = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const Tcalc psi_cosc = -(cd[0]*de[0] + cd[1]*de[1] + cd[2]*de[2]) * invcd * invde;
    Tcalc psi_isinb2, psi_isinc2;
    if (tcalc_is_double) {
      psi_isinb2 = (psi_cosb * psi_cosb < asymptotic_to_one_lf) ?
                   dpsi / (1.0 - (psi_cosb * psi_cosb)) : dpsi * inverse_one_minus_asymptote_lf;
      psi_isinc2 = (psi_cosc * psi_cosc < asymptotic_to_one_lf) ?
                   dpsi / (1.0 - (psi_cosc * psi_cosc)) : dpsi * inverse_one_minus_asymptote_lf;
    }
    else {
      psi_isinb2 = (psi_cosb * psi_cosb < asymptotic_to_one_f) ?
                   dpsi / (value_one - (psi_cosb * psi_cosb)) :
                   dpsi * inverse_one_minus_asymptote_f;
      psi_isinc2 = (psi_cosc * psi_cosc < asymptotic_to_one_f) ?
                   dpsi / (value_one - (psi_cosc * psi_cosc)) :
                   dpsi * inverse_one_minus_asymptote_f;
    }
    const Tcalc psi_fa = -invbc * psi_isinb2;
    const Tcalc psi_fb1 = (mgcd - (mgbc * psi_cosb)) * invbcd * psi_isinb2;
    const Tcalc psi_fb2 = psi_cosc * invcd * psi_isinc2;
    const Tcalc psi_fc1 = (mgcd - (mgde * psi_cosc)) * invcde * psi_isinc2;
    const Tcalc psi_fc2 = psi_cosb * invcd * psi_isinb2;
    const Tcalc psi_fd = -invde * psi_isinc2;

    if (isSignedIntegralScalarType<Tforce>()) {

      // Accumulate the phi dihedral forces
      Tforce ifrc_ix = llround(crabbc[0] * phi_fa * force_factor);
      Tforce ifrc_jx = llround(((phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0])) * force_factor);
      Tforce ifrc_lx = llround(-phi_fd * crbccd[0] * force_factor);
      xfrc[i_atom] += ifrc_ix;
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] -= ifrc_ix + ifrc_jx + ifrc_lx;
      xfrc[l_atom] += ifrc_lx;
      Tforce ifrc_iy = llround(crabbc[1] * phi_fa * force_factor);
      Tforce ifrc_jy = llround(((phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1])) * force_factor);
      Tforce ifrc_ly = llround(-phi_fd * crbccd[1] * force_factor);
      yfrc[i_atom] += ifrc_iy;
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] -= ifrc_iy + ifrc_jy + ifrc_ly;
      yfrc[l_atom] += ifrc_ly;
      Tforce ifrc_iz = llround(crabbc[2] * phi_fa * force_factor);
      Tforce ifrc_jz = llround(((phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2])) * force_factor);
      Tforce ifrc_lz = llround(-phi_fd * crbccd[2] * force_factor);
      zfrc[i_atom] += ifrc_iz;
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] -= ifrc_iz + ifrc_jz + ifrc_lz;
      zfrc[l_atom] += ifrc_lz;

      // Accumulate the psi dihedral forces
      ifrc_jx = llround(crbccd[0] * psi_fa * force_factor);
      Tforce ifrc_kx = llround(((psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0])) * force_factor);
      Tforce ifrc_mx = llround(-psi_fd * crcdde[0] * force_factor);
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] += ifrc_kx;
      xfrc[l_atom] -= ifrc_jx + ifrc_kx + ifrc_mx;
      xfrc[m_atom] += ifrc_mx;
      ifrc_jy = llround(crbccd[1] * psi_fa * force_factor);
      Tforce ifrc_ky = llround(((psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1])) * force_factor);
      Tforce ifrc_my = llround(-psi_fd * crcdde[1] * force_factor);
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] += ifrc_ky;
      yfrc[l_atom] -= ifrc_jy + ifrc_ky + ifrc_my;
      yfrc[m_atom] += ifrc_my;
      ifrc_jz = llround(crbccd[2] * psi_fa * force_factor);
      Tforce ifrc_kz = llround(((psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2])) * force_factor);
      Tforce ifrc_mz = llround(-psi_fd * crcdde[2] * force_factor);
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] += ifrc_kz;
      zfrc[l_atom] -= ifrc_jz + ifrc_kz + ifrc_mz;
      zfrc[m_atom] += ifrc_mz;
    }
    else {
      
      // Accumulate the phi dihedral forces
      Tforce frc_ix = crabbc[0] * phi_fa;
      Tforce frc_jx = (phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0]);
      Tforce frc_lx = -phi_fd * crbccd[0];
      xfrc[i_atom] += frc_ix;
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] -= frc_ix + frc_jx + frc_lx;
      xfrc[l_atom] += frc_lx;
      Tforce frc_iy = crabbc[1] * phi_fa;
      Tforce frc_jy = (phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1]);
      Tforce frc_ly = -phi_fd * crbccd[1];
      yfrc[i_atom] += frc_iy;
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] -= frc_iy + frc_jy + frc_ly;
      yfrc[l_atom] += frc_ly;
      Tforce frc_iz = crabbc[2] * phi_fa;
      Tforce frc_jz = (phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2]);
      Tforce frc_lz = -phi_fd * crbccd[2];
      zfrc[i_atom] += frc_iz;
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] -= frc_iz + frc_jz + frc_lz;
      zfrc[l_atom] += frc_lz;

      // Accumulate the psi dihedral forces
      frc_jx = crbccd[0] * psi_fa;
      Tforce frc_kx = (psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0]);
      Tforce frc_mx = -psi_fd * crcdde[0];
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] += frc_kx;
      xfrc[l_atom] -= frc_jx + frc_kx + frc_mx;
      xfrc[m_atom] += frc_mx;
      frc_jy = crbccd[1] * psi_fa;
      Tforce frc_ky = (psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1]);
      Tforce frc_my = -psi_fd * crcdde[1];
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] += frc_ky;
      yfrc[l_atom] -= frc_jy + frc_ky + frc_my;
      yfrc[m_atom] += frc_my;
      frc_jz = crbccd[2] * psi_fa;
      Tforce frc_kz = (psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2]);
      Tforce frc_mz = -psi_fd * crcdde[2];
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] += frc_kz;
      zfrc[l_atom] -= frc_jz + frc_kz + frc_mz;
      zfrc[m_atom] += frc_mz;
    }
  }

  // Return the double-precision energy result.  It will be converted to a fixed-precision
  // quantity later.
  return contrib;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Vec2<Tcalc> evaluateAttenuated14Pair(const int i_atom, const int l_atom, const int attn_idx,
                                     const Tcalc coulomb_constant, const Tcalc* charges,
                                     const int* lj_param_idx, const Tcalc* attn14_elec_factors,
                                     const Tcalc* attn14_vdw_factors, const Tcalc* lja_14_coeff,
                                     const Tcalc* ljb_14_coeff, const int n_lj_types,
                                     const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                     const double* umat, const double* invu,
                                     const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                     Tforce* zfrc, const EvaluateForce eval_elec_force,
                                     const EvaluateForce eval_vdw_force,
                                     const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;
  const int ilj_t = lj_param_idx[i_atom];
  const int jlj_t = lj_param_idx[l_atom];
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = static_cast<Tcalc>(xcrd[l_atom] - xcrd[i_atom]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[l_atom] - ycrd[i_atom]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[l_atom] - zcrd[i_atom]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[l_atom] - xcrd[i_atom];
    dy = ycrd[l_atom] - ycrd[i_atom];
    dz = zcrd[l_atom] - zcrd[i_atom];
  }
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE);
  const Tcalc invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
  const Tcalc invr = (tcalc_is_double) ? sqrt(invr2) : sqrtf(invr2);
  const Tcalc invr4 = invr2 * invr2;
  const Tcalc ele_scale = attn14_elec_factors[attn_idx];
  const Tcalc vdw_scale = attn14_vdw_factors[attn_idx];
  const Tcalc qiqj = (coulomb_constant * charges[i_atom] * charges[l_atom]) / ele_scale;
  const Tcalc lja = lja_14_coeff[(ilj_t * n_lj_types) + jlj_t] / vdw_scale;
  const Tcalc ljb = ljb_14_coeff[(ilj_t * n_lj_types) + jlj_t] / vdw_scale;
  const Tcalc ele_contrib = qiqj * invr;
  const Tcalc vdw_contrib = (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);

  // Evaluate the force, if requested
  if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
    Tcalc fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) : 0.0;
    if (eval_vdw_force == EvaluateForce::YES) {
      if (tcalc_is_double) {
        fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
      }
      else {
        fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
      }
    }
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifmag_dx = llround(fmag * dx * force_factor);
      const Tforce ifmag_dy = llround(fmag * dy * force_factor);
      const Tforce ifmag_dz = llround(fmag * dz * force_factor);
      xfrc[i_atom] += ifmag_dx;
      yfrc[i_atom] += ifmag_dy;
      zfrc[i_atom] += ifmag_dz;
      xfrc[l_atom] -= ifmag_dx;
      yfrc[l_atom] -= ifmag_dy;
      zfrc[l_atom] -= ifmag_dz;
    }
    else {
      const Tcalc fmag_dx = fmag * dx;
      const Tcalc fmag_dy = fmag * dy;
      const Tcalc fmag_dz = fmag * dz;
      xfrc[i_atom] += fmag_dx;
      yfrc[i_atom] += fmag_dy;
      zfrc[i_atom] += fmag_dz;
      xfrc[l_atom] -= fmag_dx;
      yfrc[l_atom] -= fmag_dy;
      zfrc[l_atom] -= fmag_dz;
    }
  } 
  return Vec2<Tcalc>(ele_contrib, vdw_contrib);
}

} // namespace energy
} // namespace omni
