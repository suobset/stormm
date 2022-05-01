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

} // namespace energy
} // namespace omni
