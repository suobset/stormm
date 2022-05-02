// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateNonbondedEnergy(const NonbondedKit<Tcalc> nbk, const StaticExclusionMaskReader ser,
                                const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                const double* umat, const double* invu,
                                const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                Tforce* zfrc, ScoreCard *ecard,
                                const EvaluateForce eval_elec_force,
                                const EvaluateForce eval_vdw_force, const int system_index,
                                const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  const Tcalc value_one = 1.0;
  
  // Initialize the energy result as two separate accumulators for each of a pair of quantities
  double ele_energy = 0.0;
  double vdw_energy = 0.0;
  llint ele_acc = 0LL;
  llint vdw_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Allocate arrays for tile coordinates and accumulated forces, akin to GPU protocols.
  std::vector<Tcalc> cachi_xcrd(tile_length), cachi_ycrd(tile_length), cachi_zcrd(tile_length);
  std::vector<Tcalc> cachj_xcrd(tile_length), cachj_ycrd(tile_length), cachj_zcrd(tile_length);
  std::vector<Tcalc> cachi_q(tile_length), cachj_q(tile_length);
  std::vector<int> cachi_ljidx(tile_length), cachj_ljidx(tile_length);
  std::vector<Tcalc> cachi_xfrc(tile_length), cachi_yfrc(tile_length), cachi_zfrc(tile_length);
  std::vector<Tcalc> cachj_xfrc(tile_length), cachj_yfrc(tile_length), cachj_zfrc(tile_length);
  
  // Perform nested loops over all supertiles and all tiles within them
  for (int sti = 0; sti < ser.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {

      // Tile dimensions and locations
      const int stni_atoms = std::min(ser.atom_count - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(ser.atom_count - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;

      // Access the supertile's map index: if zero, there are no exclusions to worry about
      const int stij_map_index = ser.supertile_map_idx[(stj * ser.supertile_stride_count) + sti];
      const int diag_supertile = (sti == stj);

      // The outer loops can proceed until the branch about exclusions
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);

          // Pre-cache the atom positions to avoid conversions in the inner loops.  Pre-cache
          // properties and local force accumulators to mimic GPU activity.
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tcoord_is_sgnint) {
              cachi_xcrd[i] = static_cast<Tcalc>(xcrd[atom_i]) * inv_gpos_factor;
              cachi_ycrd[i] = static_cast<Tcalc>(ycrd[atom_i]) * inv_gpos_factor;
              cachi_zcrd[i] = static_cast<Tcalc>(zcrd[atom_i]) * inv_gpos_factor;
            }
            else {
              cachi_xcrd[i] = xcrd[atom_i];
              cachi_ycrd[i] = ycrd[atom_i];
              cachi_zcrd[i] = zcrd[atom_i];
            }
            cachi_xfrc[i] = 0.0;
            cachi_yfrc[i] = 0.0;
            cachi_zfrc[i] = 0.0;
            cachi_q[i] = nbk.coulomb_constant * nbk.charge[atom_i];
            cachi_ljidx[i] = nbk.lj_idx[atom_i];
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tcoord_is_sgnint) {
              cachj_xcrd[j] = static_cast<Tcalc>(xcrd[atom_j]) * inv_gpos_factor;
              cachj_ycrd[j] = static_cast<Tcalc>(ycrd[atom_j]) * inv_gpos_factor;
              cachj_zcrd[j] = static_cast<Tcalc>(zcrd[atom_j]) * inv_gpos_factor;
            }
            else {
              cachj_xcrd[j] = xcrd[atom_j];
              cachj_ycrd[j] = ycrd[atom_j];
              cachj_zcrd[j] = zcrd[atom_j];
            }
            cachj_xfrc[j] = 0.0;
            cachj_yfrc[j] = 0.0;
            cachj_zfrc[j] = 0.0;
            cachj_q[j] = nbk.charge[atom_j];
            cachj_ljidx[j] = nbk.lj_idx[atom_j];              
          }

          // Branch for different types of tiles
          if (stij_map_index == 0) {

            // Exclusion-free loops
            for (int i = 0; i < ni_atoms; i++) {
              const Tcalc qi = cachi_q[i];
              const int ljt_i = cachi_ljidx[i];
              const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
              const Tcalc atomi_x = cachi_xcrd[i];
              const Tcalc atomi_y = cachi_ycrd[i];
              const Tcalc atomi_z = cachi_zcrd[i];
              for (int j = 0; j < jlim; j++) {
                const Tcalc dx = cachj_xcrd[j] - atomi_x;
                const Tcalc dy = cachj_ycrd[j] - atomi_y;
                const Tcalc dz = cachj_zcrd[j] - atomi_z;
                Tcalc invr2, invr;
                if (tcalc_is_double) {
                  invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
                  invr = sqrt(invr2);
                }
                else {
                  invr2 = value_one / ((dx * dx) + (dy * dy) + (dz * dz));
                  invr = sqrtf(invr2);
                }
                const Tcalc invr4 = invr2 * invr2;
                const Tcalc qiqj = qi * cachj_q[j];
                const Tcalc ele_contrib = qiqj * invr;
                ele_energy += ele_contrib;
                ele_acc += static_cast<llint>(llround(ele_contrib * nrg_scale_factor));
                const int ljt_j = cachj_ljidx[j];
                const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                const Tcalc vdw_contrib = ((lja * invr4 * invr4) - (ljb * invr2)) * invr4;
                vdw_energy += vdw_contrib;
                vdw_acc += static_cast<llint>(llround(vdw_contrib * nrg_scale_factor));

                // Evaluate the force, if requested
                if (eval_elec_force == EvaluateForce::YES ||
                    eval_vdw_force == EvaluateForce::YES) {
                  Tcalc fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) :
                                                                         0.0;
                  if (eval_vdw_force == EvaluateForce::YES) {
                    if (tcalc_is_double) {
                      fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
                    }
                    else {
                      fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
                    }
                  }
                  const Tcalc fmag_dx = fmag * dx;
                  const Tcalc fmag_dy = fmag * dy;
                  const Tcalc fmag_dz = fmag * dz;
                  cachi_xfrc[i] += fmag_dx;
                  cachi_yfrc[i] += fmag_dy;
                  cachi_zfrc[i] += fmag_dz;
                  cachj_xfrc[j] -= fmag_dx;
                  cachj_yfrc[j] -= fmag_dy;
                  cachj_zfrc[j] -= fmag_dz;
                }
              }
            }
          }
          else {

            // Get the tile's mask and check exclusions with each interaction
            const int tij_map_index = ser.tile_map_idx[stij_map_index +
                                                       (tj * tile_lengths_per_supertile) + ti];
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              const uint mask_i = ser.mask_data[tij_map_index + i];
              const Tcalc atomi_x = cachi_xcrd[i];
              const Tcalc atomi_y = cachi_ycrd[i];
              const Tcalc atomi_z = cachi_zcrd[i];
              const Tcalc qi = cachi_q[i];
              const int ljt_i = cachi_ljidx[i];
              const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
              for (int j = 0; j < jlim; j++) {
                if ((mask_i >> j) & 0x1) {
                  continue;
                }
                const Tcalc dx = cachj_xcrd[j] - atomi_x;
                const Tcalc dy = cachj_ycrd[j] - atomi_y;
                const Tcalc dz = cachj_zcrd[j] - atomi_z;
                Tcalc invr2, invr;
                if (tcalc_is_double) {
                  invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
                  invr = sqrt(invr2);
                }
                else {
                  invr2 = value_one / ((dx * dx) + (dy * dy) + (dz * dz));
                  invr = sqrtf(invr2);
                }
                const Tcalc invr4 = invr2 * invr2;
                const Tcalc qiqj = qi * cachj_q[j];
                const Tcalc ele_contrib = qiqj * invr;
                ele_energy += ele_contrib;
                ele_acc += static_cast<llint>(llround(ele_contrib * nrg_scale_factor));
                const int ljt_j = cachj_ljidx[j];
                const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                const Tcalc vdw_contrib = ((lja * invr4 * invr4) - (ljb * invr2)) * invr4;
                vdw_energy += vdw_contrib;
                vdw_acc += static_cast<llint>(llround(vdw_contrib * nrg_scale_factor));

                // Evaluate the force, if requested
                if (eval_elec_force == EvaluateForce::YES ||
                    eval_vdw_force == EvaluateForce::YES) {
                  Tcalc fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) :
                                                                          0.0;
                  if (eval_vdw_force == EvaluateForce::YES) {
                    if (tcalc_is_double) {
                      fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
                    }
                    else {
                      fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
                    }
                  }
                  const Tcalc fmag_dx = fmag * dx;
                  const Tcalc fmag_dy = fmag * dy;
                  const Tcalc fmag_dz = fmag * dz;
                  cachi_xfrc[i] += fmag_dx;
                  cachi_yfrc[i] += fmag_dy;
                  cachi_zfrc[i] += fmag_dz;
                  cachj_xfrc[j] -= fmag_dx;
                  cachj_yfrc[j] -= fmag_dy;
                  cachj_zfrc[j] -= fmag_dz;
                }                
              }
            }
          }

          // Contribute cached forces back to global arrays.
          if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_i] += llround(cachi_xfrc[i] * force_factor);              
                yfrc[atom_i] += llround(cachi_yfrc[i] * force_factor);
                zfrc[atom_i] += llround(cachi_zfrc[i] * force_factor);
              }
              else {
                xfrc[atom_i] += cachi_xfrc[i];
                yfrc[atom_i] += cachi_yfrc[i];
                zfrc[atom_i] += cachi_zfrc[i];
              }
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_j] += llround(cachj_xfrc[j] * force_factor);              
                yfrc[atom_j] += llround(cachj_yfrc[j] * force_factor);
                zfrc[atom_j] += llround(cachj_zfrc[j] * force_factor);
              }
              else {
                xfrc[atom_j] += cachj_xfrc[j];
                yfrc[atom_j] += cachj_yfrc[j];
                zfrc[atom_j] += cachj_zfrc[j];
              }
            }
          }
        }
      }
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::ELECTROSTATIC, ele_acc, system_index);
  ecard->contribute(StateVariable::VDW, vdw_acc, system_index);

  // Return the double-precision energy sums, if of interest
  return { ele_energy, vdw_energy };
}

//-------------------------------------------------------------------------------------------------
#if 0
double evaluateGeneralizedBornEnergy(const NonbondedKit<double> nbk,
                                     const ImplicitSolventKit<double> isk,
                                     const NeckGeneralizedBornTable &ngb_tables,
                                     const double* xcrd, const double* ycrd, const double* zcrd,
                                     const double* umat, const double* invu,
                                     const UnitCellType unit_cell, double* xfrc, double* yfrc,
                                     double* zfrc, ScoreCard *ecard,
                                     const EvaluateForce eval_force, const int system_index) {

  // Complete the implicit solvent model and initialize the energy result
  const ImplicitSolventRecipe<double> isr(isk, ngb_tables.getDoublePrecisionAbstract());
  double egb_energy = 0.0;
  llint egb_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Return zero if there is no implicit solvent model.  This is a reference calculation and
  // intended to keep clean code.
  if (isr.igb == ImplicitSolventModel::NONE) {
    ecard->contribute(StateVariable::GENERALIZED_BORN, egb_acc, system_index);
    return egb_energy;
  }

  // Allocate space for the calculation
  std::vector<double> effective_gb_radii(nbk.natom);
  std::vector<double> psi(nbk.natom, 0.0);
  std::vector<double> sumdeijda(nbk.natom);

  // Shorthand for some constants defined in Constants/generalized_born.h
  const double gta  = gb_taylor_a_lf;
  const double gtb  = gb_taylor_b_lf;
  const double gtc  = gb_taylor_c_lf;
  const double gtd  = gb_taylor_d_lf;
  const double gtdd = gb_taylor_dd_lf;
  const double gte  = gb_taylor_e_lf;
  const double gtf  = gb_taylor_f_lf;
  const double gtg  = gb_taylor_g_lf;
  const double gth  = gb_taylor_h_lf;
  const double gthh = gb_taylor_hh_lf;

  // Perform nested loops over all atoms to compute the GB radii
  for (int i = 1; i < nbk.natom; i++) {
    const double atomi_x = xcrd[i];
    const double atomi_y = ycrd[i];
    const double atomi_z = zcrd[i];
    const double atomi_radius = isr.pb_radii[i] - isr.gb_offset;
    const double atomi_inv_radius = 1.0 / atomi_radius;
    for (int j = 0; j < i; j++) {
      const double dx   = atomi_x - xcrd[j];
      const double dy   = atomi_y - ycrd[j];
      const double dz   = atomi_z - zcrd[j];
      const double r2   = (dx * dx) + (dy * dy) + (dz * dz);
      const double r    = sqrt(r2);
      const double invr = 1.0 / r;
      const double atomj_radius = isr.pb_radii[j] - isr.gb_offset;
      const double atomj_inv_radius = 1.0 / atomj_radius;

      // First computation: atom I -> atom J
      const double sj = isr.gb_screen[j] * atomj_radius;
      const double sj2 = sj * sj;
      if (r > 4.0 * sj) {
        const double invr2 = invr * invr;
        const double tmpsd = sj2 * invr2;
        const double dumbo = gta + tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd * gtdd))); 
        psi[i] -= sj * tmpsd * invr2 * dumbo;
      }
      else if (r > atomi_radius + sj) {
        psi[i] -= 0.5 * ((sj / (r2 - sj2)) + (0.5 * invr * log((r - sj) / (r + sj))));
      }
      else if (r > fabs(atomi_radius - sj)) {
        const double theta = 0.5 * atomi_inv_radius * invr *
                             (r2 + (atomi_radius * atomi_radius) - sj2);
        const double uij   = 1.0 / (r + sj);
        psi[i] -= 0.25 * (atomi_inv_radius * (2.0 - theta) -
                            uij + invr * log(atomi_radius * uij));
      }
      else if (atomi_radius < sj) {
        psi[i] -= 0.5 * ((sj / (r2 - sj2)) + (2.0 * atomi_inv_radius) +
                         (0.5 * invr * log((sj - r) / (sj + r))));
      }

      // Second computation: atom J -> atom I
      const double si = isr.gb_screen[i] * atomi_radius;
      const double si2 = si * si;
      if (r > 4.0 * si) {
        const double invr2  = invr * invr;
        const double tmpsd  = si2 * invr2;
        const double dumbo  = gta + tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd * gtdd)));
        psi[j] -= si * tmpsd * invr2 * dumbo;
      }
      else if (r > atomj_radius + si) {
        psi[j] -= 0.5 * ((si / (r2 - si2)) + (0.5 * invr * log((r - si) / (r + si))));
      }
      else if (r > fabs(atomj_radius - si)) {
        const double theta = 0.5 * atomj_inv_radius * invr *
                             (r2 + (atomj_radius * atomj_radius) - si2);
        const double uij   = 1.0 / (r + si);
        psi[j] -= 0.25 * (atomj_inv_radius * (2.0 - theta) - uij +
                          invr * log(atomj_radius * uij));
      }
      else if (atomj_radius < si) {
        psi[j] -= 0.5 * ((si / (r2 - si2)) + (2.0 * atomj_inv_radius) +
                         (0.5 * invr * log((si - r) / (si + r))));
      }

      // Neck GB contribution
      if ((isr.igb == ImplicitSolventModel::NECK_GB ||
           isr.igb == ImplicitSolventModel::NECK_GB_II) &&
          r < isr.pb_radii[i] + isr.pb_radii[j] + isr.gb_neckcut) {

        // First computation: atom I -> atom J
        const int ij_table_idx = (isr.table_size * isr.neck_gb_idx[j]) + isr.neck_gb_idx[i];
        double mdist  = r - isr.neck_max_sep[ij_table_idx];
        double mdist2 = mdist * mdist;
        double mdist6 = mdist2 * mdist2 * mdist2;
        const double ij_neck = isr.neck_max_val[ij_table_idx] / (1.0 + mdist2 + (0.3 * mdist6));
        psi[i] -= isr.gb_neckscale * ij_neck;

        // Second computation: atom J -> atom I
        const int ji_table_idx = (isr.table_size * isr.neck_gb_idx[i]) + isr.neck_gb_idx[j];
        mdist  = r - isr.neck_max_sep[ji_table_idx];
        mdist2 = mdist * mdist;
        mdist6 = mdist2 * mdist2 * mdist2;
        const double ji_neck = isr.neck_max_val[ji_table_idx] / (1.0 + mdist2 + (0.3 * mdist6));
        psi[j] -= isr.gb_neckscale * ji_neck;
      }
    }
  }

  // Make a second pass to finalize the effective GB radii
  for (int i = 0; i < nbk.natom; i++) {
    switch (isr.igb) {
    case ImplicitSolventModel::HCT_GB:
      {
        // Original (Hawkins-Craemer-Truhlar) effective radii
        const double atomi_inv_radius = 1.0 / (isr.pb_radii[i] - isr.gb_offset);
        effective_gb_radii[i] = 1.0 / (atomi_inv_radius + psi[i]);
        if (effective_gb_radii[i] < 0.0) {
          effective_gb_radii[i] = 30.0;
        }
      }
      break;
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      {
        // "GBAO" formulas
        const double atomi_radius = isr.pb_radii[i] - isr.gb_offset;
        const double atomi_inv_radius = 1.0 / atomi_radius;
        const double fipsi = psi[i] * (-atomi_radius);
        effective_gb_radii[i] = 1.0 / (atomi_inv_radius -
                                       tanh((isr.gb_alpha[i] - (isr.gb_beta[i] * fipsi) +
                                             (isr.gb_gamma[i] * fipsi * fipsi)) * fipsi) / 
                                       isr.pb_radii[i]);
      }
      break;
    case ImplicitSolventModel::NONE:
      break;
    }
  }

  // Compute inherent Generalized Born energies and initialize an array for solvent forces
  for (int i = 0; i < nbk.natom; i++) {
    const double atomi_q = nbk.charge[i];
    const double atomi_radius = effective_gb_radii[i];
    const double expmkf = exp(-gb_kscale * isr.kappa * atomi_radius) / isr.dielectric;
    const double dielfac = 1.0 - expmkf;
    const double atmq2h = 0.5 * atomi_q * atomi_q * nbk.coulomb_constant;
    const double atmqd2h = atmq2h * dielfac;
    const double contrib = -atmqd2h / atomi_radius;
    egb_energy += contrib;
    egb_acc += static_cast<llint>(contrib * nrg_scale_factor);
    if (eval_force == EvaluateForce::YES) {
      sumdeijda[i] = atmqd2h - (gb_kscale * isr.kappa * atmq2h * expmkf * atomi_radius);
    }
  }

  // Due to the lack of exclusions, the Generalized Born reference calculation is a much simpler
  // pair of nested loops over all atoms without self-interactions or double-counting.
  for (int i = 1; i < nbk.natom; i++) {
    const double atomi_x = xcrd[i];
    const double atomi_y = ycrd[i];
    const double atomi_z = zcrd[i];
    const double atomi_q = nbk.coulomb_constant * nbk.charge[i];
    const double atomi_radius = effective_gb_radii[i];
    for (int j = 0; j < i; j++) {
      const double dx = xcrd[j] - atomi_x;
      const double dy = ycrd[j] - atomi_y;
      const double dz = zcrd[j] - atomi_z;
      const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
      const double qiqj = atomi_q * nbk.charge[j];
      const double ij_born_radius = atomi_radius * effective_gb_radii[j];
      const double efac = exp(-r2 / (4.0 * ij_born_radius));
      const double fgbi = 1.0 / sqrt(r2 + ij_born_radius * efac);
      const double fgbk = -isr.kappa * gb_kscale / fgbi;
      const double expmkf = exp(fgbk) / isr.dielectric;
      const double dielfac = 1.0 - expmkf;
      const double contrib = -qiqj * dielfac * fgbi;
      egb_energy += contrib;
      egb_acc += static_cast<llint>(contrib * nrg_scale_factor);
      if (eval_force == EvaluateForce::YES) {
        const double temp4 = fgbi * fgbi * fgbi;
        const double temp6 = qiqj * temp4 * (dielfac + fgbk * expmkf);
        const double fmag = temp6 * (1.0 - (0.25 * efac));
        const double temp5 = 0.5 * efac * temp6 * (ij_born_radius + (0.25 * r2));
        sumdeijda[i] += atomi_radius * temp5;
        sumdeijda[j] += effective_gb_radii[j] * temp5;
        xfrc[i] += fmag * dx;
        yfrc[i] += fmag * dy;
        zfrc[i] += fmag * dz;
        xfrc[j] -= fmag * dx;
        yfrc[j] -= fmag * dy;
        zfrc[j] -= fmag * dz;
      }
    }
  }

  // A second pair of nested loops over all atoms is needed to fold in derivatives of the
  // effective Born radii to the forces on each atom.
  if (eval_force == EvaluateForce::YES) {
    switch (isr.igb) {
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      for (int i = 0; i < nbk.natom; i++) {
        const double atomi_radius = isr.pb_radii[i] - isr.gb_offset;
        const double fipsi = psi[i] * (-atomi_radius);
        const double thi = tanh((isr.gb_alpha[i] -
                                 (isr.gb_beta[i] - (isr.gb_gamma[i] * fipsi)) * fipsi) * fipsi);
        sumdeijda[i] *= (isr.gb_alpha[i] -
                         ((2.0 * isr.gb_beta[i]) - (3.0 * isr.gb_gamma[i] * fipsi)) * fipsi) *
                        (1.0 - thi * thi) * atomi_radius / isr.pb_radii[i];
      }
    }
    for (int i = 1; i < nbk.natom; i++) {
      const double atomi_x = xcrd[i];
      const double atomi_y = ycrd[i];
      const double atomi_z = zcrd[i];
      const double atomi_radius = isr.pb_radii[i] - isr.gb_offset;
      const double atomi_inv_radius = 1.0 / atomi_radius;
      for (int j = 0; j < i; j++) {
        const double dx = xcrd[j] - atomi_x;
        const double dy = ycrd[j] - atomi_y;
        const double dz = zcrd[j] - atomi_z;
        const double atomj_radius = isr.pb_radii[j] - isr.gb_offset;
        const double atomj_inv_radius = 1.0 / atomj_radius;
        const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
        const double invr = 1.0 / sqrt(r2);
        const double invr2 = invr * invr;
        const double r = r2 * invr;

        // First computation: atom I -> atom J
        const double sj = isr.gb_screen[j] * atomj_radius;
        const double sj2 = sj * sj;
        double datmpi, datmpj;
        if (r > 4.0 * sj) {
          const double tmpsd  = sj2 * invr2;
          const double dumbo  = gte + tmpsd*(gtf + tmpsd*(gtg + tmpsd*(gth + tmpsd*gthh)));
          datmpi = tmpsd * sj * invr2 * invr2 * dumbo;
        }
        else if (r > atomi_radius + sj) {
          const double temp1  = 1.0 / (r2 - sj2);
          datmpi = (temp1 * sj * (-0.5 * invr2 + temp1)) +
                   (0.25 * invr * invr2 * log((r - sj) / (r + sj)));
        }
        else if (r > fabs(atomi_radius - sj)) {
          const double temp1  = 1.0 / (r + sj);
          const double invr3  = invr2 * invr;
          datmpi = -0.25 * ((-0.5 * (r2 - atomi_radius * atomi_radius + sj2) *
                             invr3 * atomi_inv_radius * atomi_inv_radius) +
                            (invr * temp1 * (temp1 - invr)) - (invr3 * log(atomi_radius * temp1)));
        }
        else if (atomi_radius < sj) {
          const double temp1  = 1.0 / (r2 - sj2);
          datmpi = -0.5 * ((sj * invr2 * temp1) - (2.0 * sj * temp1 * temp1) -
                           (0.5 * invr2 * invr * log((sj - r) / (sj + r))));
        }
        else {
          datmpi = 0.0;
        }

        // Second computation: atom J -> atom I
        const double si = isr.gb_screen[i] * atomi_radius;
        const double si2 = si * si;
        if (r > 4.0 * si) {
          const double tmpsd  = si2 * invr2;
          const double dumbo  = gte + tmpsd*(gtf + tmpsd*(gtg + tmpsd*(gth + tmpsd*gthh)));
          datmpj = tmpsd * si * invr2 * invr2 * dumbo;
        }
        else if (r > atomj_radius + si) {
          const double temp1  = 1.0 / (r2 - si2);
          datmpj = (temp1 * si * (-0.5 * invr2 + temp1)) +
                   (0.25 * invr * invr2 * log((r - si) / (r + si)));
        }
        else if (r > fabs(atomj_radius - si)) {
          const double temp1 = 1.0 / (r + si);
          const double invr3 = invr2 * invr;
          datmpj = -0.25 * ((-0.5 * (r2 - atomj_radius * atomj_radius + si2) *
                             invr3 * atomj_inv_radius * atomj_inv_radius) +
                            (invr * temp1 * (temp1 - invr)) - (invr3 * log(atomj_radius * temp1)));
        }
        else if (atomj_radius < si) {
          const double temp1  = 1.0 / (r2 - si2);
          datmpj = -0.5 * ((si * invr2 * temp1) - (2.0 * si * temp1 * temp1) -
                           (0.5 * invr2 * invr * log((si - r) / (si + r))));
        }
        else {
          datmpj = 0.0;
        }

        // Neck GB contributions
        if ((isr.igb == ImplicitSolventModel::NECK_GB ||
             isr.igb == ImplicitSolventModel::NECK_GB_II) &&
            r < isr.pb_radii[i] + isr.pb_radii[j] + isr.gb_neckcut) {

          // First computation: atom I -> atom J
          const int ij_table_idx = (isr.table_size * isr.neck_gb_idx[j]) + isr.neck_gb_idx[i];
          double mdist = r - isr.neck_max_sep[ij_table_idx];
          double mdist2 = mdist * mdist;
          double mdist6 = mdist2 * mdist2 * mdist2;
          double temp1 = 1.0 + mdist2 + (0.3 * mdist6);
          temp1 = temp1 * temp1 * r;
          datmpi += (((2.0 * mdist) + (1.8 * mdist2 * mdist2 * mdist)) *
                     isr.neck_max_val[ij_table_idx] * isr.gb_neckscale) / temp1;

          // Second computation: atom J -> atom I
          const int ji_table_idx = (isr.table_size * isr.neck_gb_idx[i]) + isr.neck_gb_idx[j];
          mdist = r - isr.neck_max_sep[ji_table_idx];
          mdist2 = mdist * mdist;
          mdist6 = mdist2 * mdist2 * mdist2;
          temp1 = 1.0 + mdist2 + 0.3*mdist6;
          temp1 = temp1 * temp1 * r;
          datmpj += (((2.0 * mdist) + (1.8 * mdist2 * mdist2 * mdist)) *
                     isr.neck_max_val[ji_table_idx] * isr.gb_neckscale) / temp1;
        }

        // Contribute the derivatives to the force arrays
        const double fmag = (datmpi * sumdeijda[i]) + (datmpj * sumdeijda[j]);
        xfrc[i] -= fmag * dx;
        yfrc[i] -= fmag * dy;
        zfrc[i] -= fmag * dz;
        xfrc[j] += fmag * dx;
        yfrc[j] += fmag * dy;
        zfrc[j] += fmag * dz;
      }
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::GENERALIZED_BORN, egb_acc, system_index);

  // Return the double-precision energy sum, if of interest
  return egb_energy;
}
#endif

} // namespace energy
} // namespace omni
