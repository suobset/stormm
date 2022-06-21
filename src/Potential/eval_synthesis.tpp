// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void synthesisVwuEvaluation(const SyValenceKit<Tcalc> syvk,
                            const SyRestraintKit<Tcalc, Tcalc2, Tcalc4> syrk,
                            const Tcalc* sh_charges, const int* sh_lj_idx, llint* sh_xcrd,
                            llint* sh_ycrd, llint* sh_zcrd, llint* sh_xvel, llint* sh_yvel,
                            llint* sh_zvel, llint* sh_xfrc, llint* sh_yfrc, llint* sh_zfrc,
                            const double inv_gpos_scale, const double force_scale,
                            ScoreCard *ecard, const int vwu_idx, const EvaluateForce eval_force,
                            const VwuTask activity, const VwuGoal purpose, const int step_number) {

  // Initialize energy accumulators
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();
  llint bond_acc = 0LL;
  llint angl_acc = 0LL;
  llint dihe_acc = 0LL;
  llint impr_acc = 0LL;
  llint ubrd_acc = 0LL;
  llint cimp_acc = 0LL;
  llint cmap_acc = 0LL;
  llint qq14_acc = 0LL;
  llint lj14_acc = 0LL;
  llint rest_acc = 0LL;
  const int sysid = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                       static_cast<int>(VwuAbstractMap::SYSTEM_ID)].x;
  
  // Evaluate bonds and Urey-Bradley terms
  if (activity == VwuTask::BOND || activity == VwuTask::UBRD || activity == VwuTask::ALL_TASKS) {
    const int2 cbnd_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::CBND)];
    const int2 cbnd_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                    static_cast<int>(VwuAbstractMap::CBND_NRG)];
    for (int pos = cbnd_limits.x; pos < cbnd_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syvk.cbnd_acc[cbnd_nrg_limits.x], pos - cbnd_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.cbnd_acc[cbnd_nrg_limits.x], pos - cbnd_limits.x);
        break;
      }
      const uint2 tinsr = syvk.cbnd_insr[pos];
      const bool is_urey_bradley = ((tinsr.x >> 20) & 0x1);

      // Skip Urey-Bradley interactions if only bonds are desired, and skip bonds if only
      // Urey-Bradley interactions are of interest.
      if ((activity == VwuTask::BOND && is_urey_bradley) ||
          (activity == VwuTask::UBRD && (! is_urey_bradley))) {
        continue;
      }
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int param_idx = tinsr.y;
      const Tcalc keq = (is_urey_bradley) ? syvk.ubrd_keq[param_idx] : syvk.bond_keq[param_idx];
      const Tcalc leq = (is_urey_bradley) ? syvk.ubrd_leq[param_idx] :
                                            std::abs(syvk.bond_leq[param_idx]);
      const Tcalc du =
        evalHarmonicStretch<llint, llint, Tcalc>(i_atom, j_atom, keq, leq, sh_xcrd, sh_ycrd,
                                                 sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                                 sh_xfrc, sh_yfrc, sh_zfrc, eval_force,
                                                 inv_gpos_scale, force_scale);
      if (log_term) {
        if (is_urey_bradley) {
          ubrd_acc += llround(du * nrg_scale_factor);
        }
        else {
          bond_acc += llround(du * nrg_scale_factor);
        }
      }
    }
  }
  
  // Evaluate harmonic bond angles
  if (activity == VwuTask::ANGL || activity == VwuTask::ALL_TASKS) {
    const int2 angl_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::ANGL)];
    const int2 angl_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                    static_cast<int>(VwuAbstractMap::ANGL_NRG)];
    for (int pos = angl_limits.x; pos < angl_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syvk.angl_acc[angl_nrg_limits.x], pos - angl_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.angl_acc[angl_nrg_limits.x], pos - angl_limits.x);
        break;
      }
      const uint2 tinsr = syvk.angl_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = tinsr.y;
      const double keq = syvk.angl_keq[param_idx];
      const double theta0 = syvk.angl_theta[param_idx];
      const double du =
        evalHarmonicBend<llint, llint, Tcalc>(i_atom, j_atom, k_atom, keq, theta0, sh_xcrd,
                                              sh_ycrd, sh_zcrd, nullptr, nullptr,
                                              UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc,
                                              eval_force, inv_gpos_scale, force_scale);
      if (log_term) {
        angl_acc += llround(du * nrg_scale_factor);
      }
    }
  }

  // Evaluate cosine-based dihedrals and CHARMM improper dihedrals
  if (activity == VwuTask::DIHE || activity == VwuTask::CIMP || activity == VwuTask::INFR14 ||
      activity == VwuTask::ALL_TASKS) {
    const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
    const bool tcalc_is_double = (tcalc_ct == double_type_index);
    const Tcalc value_one = 1.0;
    const int2 cdhe_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::CDHE)];
    const int2 cdhe_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                    static_cast<int>(VwuAbstractMap::CDHE_NRG)];
    const int ljabc_offset = syvk.ljabc_offsets[sysid];
    const int nljt = syvk.n_lj_types[sysid];
    for (int pos = cdhe_limits.x; pos < cdhe_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syvk.cdhe_acc[cdhe_nrg_limits.x], pos - cdhe_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.cdhe_acc[cdhe_nrg_limits.x], pos - cdhe_limits.x);
        break;
      }
      const uint2 tinsr = syvk.cdhe_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
        const int attn_idx = ((tinsr.y >> 10) & 0x1f);
        if (attn_idx > 0) {
          const Vec2<double> uc =
            evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, syvk.coulomb, sh_charges, sh_lj_idx,
                                     syvk.attn14_elec, syvk.attn14_vdw, syvk.lja_14_coeff,
                                     syvk.ljb_14_coeff, ljabc_offset, nljt, sh_xcrd, sh_ycrd,
                                     sh_zcrd, nullptr, nullptr, UnitCellType::NONE, sh_xfrc,
                                     sh_yfrc, sh_zfrc, eval_force, eval_force, inv_gpos_scale,
                                     force_scale);
          if (log_term) {
            qq14_acc += llround(uc.x * nrg_scale_factor);
            lj14_acc += llround(uc.y * nrg_scale_factor);
          }
        }
        if (activity == VwuTask::INFR14) {
          continue;
        }
      }
      const bool is_charmm_improper = ((tinsr.x >> 30) & 0x1);

      // Skip CHARMM improper interactions if only standard dihedrals are of interest, and
      // skip standard dihedrals if only CHARMM impropers are of interest.
      if ((activity == VwuTask::DIHE && is_charmm_improper) ||
          (activity == VwuTask::CIMP && (! is_charmm_improper))) {
        continue;
      }
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = ((tinsr.y >> 16) & 0xffff);
      const bool second_term = ((tinsr.y >> 15) & 0x1);
      const TorsionKind kind = (tinsr.x >> 31) ? TorsionKind::IMPROPER : TorsionKind::PROPER;
      DihedralStyle dihe_style;
      Tcalc ampl, phi, freq;
      if (is_charmm_improper) {
        dihe_style = DihedralStyle::HARMONIC;
        ampl = syvk.cimp_keq[param_idx];
        phi  = syvk.cimp_phi[param_idx];
        freq = value_one;
      }
      else {
        dihe_style = DihedralStyle::COSINE;
        ampl = syvk.dihe_amp[param_idx];
        phi  = syvk.dihe_phi[param_idx];
        freq = syvk.dihe_freq[param_idx];
      }
      double du =
        evalDihedralTwist<llint, llint, Tcalc>(i_atom, j_atom, k_atom, l_atom, ampl, phi, freq,
                                               dihe_style, sh_xcrd, sh_ycrd, sh_zcrd, nullptr,
                                               nullptr, UnitCellType::NONE, sh_xfrc, sh_yfrc,
                                               sh_zfrc, eval_force, inv_gpos_scale, force_scale);
      if (second_term) {

        // Any second term will, by construction, be a cosine-based dihedral.  Computing it in
        // this way foregoes the efficiency that could be gained by combining both terms into a
        // single twist applied to the same four atoms.  GPU kernels will exploit the opportunity.
        const uint xinsr = syvk.cdhe_ovrt_insr[pos];
        const int param_ii_idx = (xinsr & 0xfffff);
        const Tcalc ampl_ii = syvk.dihe_amp[param_ii_idx];
        const Tcalc phi_ii  = syvk.dihe_phi[param_ii_idx];
        const Tcalc freq_ii = syvk.dihe_freq[param_ii_idx];
        du += evalDihedralTwist<llint, llint, Tcalc>(i_atom, j_atom, k_atom, l_atom, ampl_ii,
                                                     phi_ii, freq_ii, dihe_style, sh_xcrd, sh_ycrd,
                                                     sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                                     sh_xfrc, sh_yfrc, sh_zfrc, eval_force,
                                                     inv_gpos_scale, force_scale);
      }
      if (log_term) {
        if (is_charmm_improper) {
          cimp_acc += llround(du * nrg_scale_factor);
        }
        else {
          if (kind == TorsionKind::PROPER) {
            dihe_acc += llround(du * nrg_scale_factor);
          }
          else {
            impr_acc += llround(du * nrg_scale_factor);
          }
        }
      }
    }
  }

  // Evaluate CMAP interactions
  if (activity == VwuTask::CMAP || activity == VwuTask::ALL_TASKS) {
    const int2 cmap_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::CMAP)];
    const int2 cmap_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                    static_cast<int>(VwuAbstractMap::CMAP_NRG)];
    for (int pos = cmap_limits.x; pos < cmap_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syvk.cmap_acc[cmap_nrg_limits.x], pos - cmap_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.cmap_acc[cmap_nrg_limits.x], pos - cmap_limits.x);
        break;
      }
      const uint2 tinsr = syvk.cmap_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      const int m_atom = ((tinsr.y >> 10) & 0x3ff);
      const int surf_idx = (tinsr.y >> 20);
      const double du =
        evalCmap<llint, llint, Tcalc>(syvk.cmap_patches, syvk.cmap_patch_bounds, surf_idx,
                                      syvk.cmap_dim[surf_idx], i_atom, j_atom, k_atom, l_atom,
                                      m_atom, sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr,
                                      UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc, eval_force,
                                      inv_gpos_scale, force_scale);
      if (log_term) {
        cmap_acc += llround(du * nrg_scale_factor);
      }
    }
  }
  
  // Evaluate any remaining 1:4 interactions
  if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
    const int2 infr14_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                  static_cast<int>(VwuAbstractMap::INFR14)];
    const int2 infr14_nrg_limits =
      syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                         static_cast<int>(VwuAbstractMap::INFR14_NRG)];
    const int ljabc_offset = syvk.ljabc_offsets[sysid];
    const int nljt = syvk.n_lj_types[sysid];    
    for (int pos = infr14_limits.x; pos < infr14_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syvk.infr14_acc[infr14_nrg_limits.x], pos - infr14_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.infr14_acc[infr14_nrg_limits.x], pos - infr14_limits.x);
        break;
      }
      const uint tinsr = syvk.infr14_insr[pos];
      const int i_atom = (tinsr & 0x3ff);
      const int l_atom = ((tinsr >> 10) & 0x3ff);
      const int attn_idx = (tinsr >> 20);
      if (attn_idx == 0) {
        continue;
      }
      const Vec2<double> uc =
        evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, syvk.coulomb, sh_charges,
                                 sh_lj_idx, syvk.attn14_elec, syvk.attn14_vdw, syvk.lja_14_coeff,
                                 syvk.ljb_14_coeff, ljabc_offset, nljt, sh_xcrd, sh_ycrd, sh_zcrd,
                                 nullptr, nullptr, UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc,
                                 eval_force, eval_force, inv_gpos_scale, force_scale);
      if (log_term) {
        qq14_acc += llround(uc.x * nrg_scale_factor);
        lj14_acc += llround(uc.y * nrg_scale_factor);
      }
    }
  }

  // Evaluate positional restraints
  if (activity == VwuTask::RPOSN || activity == VwuTask::ALL_TASKS) {
    const int2 rposn_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                 static_cast<int>(VwuAbstractMap::RPOSN)];
    const int2 rposn_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::RPOSN_NRG)];
    for (int pos = rposn_limits.x; pos < rposn_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syrk.rposn_acc[rposn_nrg_limits.x], pos - rposn_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syrk.rposn_acc[rposn_nrg_limits.x], pos - rposn_limits.x);
        break;
      }
      const uint2 tinsr = syrk.rposn_insr[pos];
      const int p_atom = (tinsr.x & 0x3ff);
      const int kr_param_idx = ((tinsr.x >> 10) & 0x1fffff);
      const int xyz_param_idx = tinsr.y;
      const double contrib =
        evalPosnRestraint(p_atom, step_number, syrk.rposn_step_bounds[kr_param_idx].x,
                          syrk.rposn_step_bounds[kr_param_idx].y,
                          syrk.rposn_init_xy[xyz_param_idx], syrk.rposn_finl_xy[xyz_param_idx],
                          syrk.rposn_init_z[xyz_param_idx], syrk.rposn_finl_z[xyz_param_idx],
                          syrk.rposn_init_k[kr_param_idx], syrk.rposn_finl_k[kr_param_idx],
                          syrk.rposn_init_r[kr_param_idx], syrk.rposn_finl_r[kr_param_idx],
                          sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE, sh_xfrc,
                          sh_yfrc, sh_zfrc, eval_force, inv_gpos_scale, force_scale);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }

  // Evaluate distance restraints
  if (activity == VwuTask::RBOND || activity == VwuTask::ALL_TASKS) {
    const int2 rbond_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                 static_cast<int>(VwuAbstractMap::RBOND)];
    const int2 rbond_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::RBOND_NRG)];
    for (int pos = rbond_limits.x; pos < rbond_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syrk.rbond_acc[rbond_nrg_limits.x], pos - rbond_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syrk.rbond_acc[rbond_nrg_limits.x], pos - rbond_limits.x);
        break;
      }
      const uint2 tinsr = syrk.rbond_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int param_idx = tinsr.y;
      const double contrib =
        evalBondRestraint<llint,
                          llint,
                          double,
                          double2,
                          double4>(i_atom, j_atom, step_number,
                                   syrk.rbond_step_bounds[param_idx].x,
                                   syrk.rbond_step_bounds[param_idx].y,
                                   syrk.rbond_init_k[param_idx], syrk.rbond_finl_k[param_idx],
                                   syrk.rbond_init_r[param_idx], syrk.rbond_finl_r[param_idx],
                                   sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                   sh_xfrc, sh_yfrc, sh_zfrc, eval_force, inv_gpos_scale,
                                   force_scale);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }

  // Evaluate angle restraints
  if (activity == VwuTask::RANGL || activity == VwuTask::ALL_TASKS) {
    const int2 rangl_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                 static_cast<int>(VwuAbstractMap::RANGL)];
    const int2 rangl_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::RANGL_NRG)];
    for (int pos = rangl_limits.x; pos < rangl_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syrk.rangl_acc[rangl_nrg_limits.x], pos - rangl_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syrk.rangl_acc[rangl_nrg_limits.x], pos - rangl_limits.x);
        break;
      }
      const uint2 tinsr = syrk.rangl_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = tinsr.y;
      const double contrib =
        evalAnglRestraint<llint,
                          llint,
                          double,
                          double2,
                          double4>(i_atom, j_atom, k_atom, step_number,
                                   syrk.rangl_step_bounds[param_idx].x,
                                   syrk.rangl_step_bounds[param_idx].y,
                                   syrk.rangl_init_k[param_idx], syrk.rangl_finl_k[param_idx],
                                   syrk.rangl_init_r[param_idx], syrk.rangl_finl_r[param_idx],
                                   sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                   sh_xfrc, sh_yfrc, sh_zfrc, eval_force, inv_gpos_scale,
                                   force_scale);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }      
    }
  }

  // Evaluate dihedral restraints
  if (activity == VwuTask::RDIHE || activity == VwuTask::ALL_TASKS) {
    const int2 rdihe_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                 static_cast<int>(VwuAbstractMap::RDIHE)];
    const int2 rdihe_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::RDIHE_NRG)];
    for (int pos = rdihe_limits.x; pos < rdihe_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(&syrk.rdihe_acc[rdihe_nrg_limits.x], pos - rdihe_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syrk.rdihe_acc[rdihe_nrg_limits.x], pos - rdihe_limits.x);
        break;
      }
      const uint2 tinsr = syrk.rdihe_insr[pos];
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      const int param_idx = ((tinsr.y >> 10) & 0x3fffff);
      const double contrib =
        evalDiheRestraint<llint,
                          llint,
                          double,
                          double2,
                          double4>(i_atom, j_atom, k_atom, l_atom, step_number,
                                   syrk.rdihe_step_bounds[param_idx].x,
                                   syrk.rdihe_step_bounds[param_idx].y,
                                   syrk.rdihe_init_k[param_idx], syrk.rdihe_finl_k[param_idx],
                                   syrk.rdihe_init_r[param_idx], syrk.rdihe_finl_r[param_idx],
                                   sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                   sh_xfrc, sh_yfrc, sh_zfrc, eval_force, inv_gpos_scale,
                                   force_scale);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }
  
  // Contribute results as the instantaneous states
  commitVwuEnergies(bond_acc, angl_acc, dihe_acc, impr_acc, ubrd_acc, cimp_acc, cmap_acc,
                    qq14_acc, lj14_acc, rest_acc, sysid, activity, ecard);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalSyValenceEnergy(const SyValenceKit<Tcalc> syvk,
                         const SyRestraintKit<Tcalc, Tcalc2, Tcalc4> syrk, PsSynthesisWriter psyw,
                         ScoreCard *ecard, const EvaluateForce eval_force, const VwuTask activity,
                         const VwuGoal purpose, const int step_number) {

  // Loop over all valence work units within the topology synthesis.  Each holds within it the
  // system index number and a list of atom imports.
  std::vector<Tcalc> sh_charges(maximum_valence_work_unit_atoms);
  std::vector<int> sh_lj_idx(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xvel(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yvel(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zvel(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zfrc(maximum_valence_work_unit_atoms);
  for (int i = 0; i < syvk.nvwu; i++) {

    // Import atoms and clear forces for this work unit
    const int2 atom_limits = syvk.vwu_abstracts[(i * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::IMPORT)];
    for (int j = atom_limits.x; j < atom_limits.y; j++) {
      const int atom_idx = syvk.vwu_imports[j];
      const size_t j_sh = j - atom_limits.x;
      sh_charges[j_sh] = syvk.charges[atom_idx];
      sh_lj_idx[j_sh]  = syvk.lj_idx[atom_idx];
      sh_xcrd[j_sh]    = psyw.xcrd[atom_idx];
      sh_ycrd[j_sh]    = psyw.ycrd[atom_idx];
      sh_zcrd[j_sh]    = psyw.zcrd[atom_idx];
      sh_xfrc[j_sh]    = 0LL;
      sh_yfrc[j_sh]    = 0LL;
      sh_zfrc[j_sh]    = 0LL;
    }
    synthesisVwuEvaluation(syvk, syrk, sh_charges.data(), sh_lj_idx.data(), sh_xcrd.data(),
                           sh_ycrd.data(), sh_zcrd.data(), sh_xvel.data(), sh_yvel.data(),
                           sh_zvel.data(), sh_xfrc.data(), sh_yfrc.data(), sh_zfrc.data(),
                           psyw.inv_gpos_scale, psyw.frc_scale, ecard, i, eval_force, activity,
                           purpose, step_number);
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
      for (int j = atom_limits.x; j < atom_limits.y; j++) {
        const int atom_idx = syvk.vwu_imports[j];
        const size_t j_sh = j - atom_limits.x;
        psyw.xfrc[atom_idx] += sh_xfrc[j_sh];
        psyw.yfrc[atom_idx] += sh_yfrc[j_sh];
        psyw.zfrc[atom_idx] += sh_zfrc[j_sh];
      }
      break;
    case VwuGoal::MOVE_PARTICLES:
      {
        const int2 manip_limits = syvk.vwu_abstracts[(i * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::MANIPULATE)];
        for (int j = atom_limits.x; j < atom_limits.y; j++) {
          const int atom_idx = syvk.vwu_imports[j];
          const int j_sh = j - atom_limits.x;
          const int manip_idx = manip_limits.x + (j_sh / uint_bit_count_int);
          const uint2 manip_mask = syvk.vwu_manip[manip_idx];
          const int manip_bit = j_sh - ((j_sh / uint_bit_count_int) * uint_bit_count_int);
          if ((manip_mask.y >> manip_bit) & 0x1) {
            psyw.xcrd[atom_idx] = sh_xcrd[j_sh];
            psyw.ycrd[atom_idx] = sh_ycrd[j_sh];
            psyw.zcrd[atom_idx] = sh_zcrd[j_sh];
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void evalSyNonbondedTileGroups(const SyNonbondedKit<Tcalc> synbk, PsSynthesisWriter psyw,
                               ScoreCard *ecard, const EvaluateForce eval_force) {

  // Loop over all non-bonded work units within the topology synthesis.  Each holds within it a
  // list of atom imports and a series of tiles.  Within each tile, all atoms pertain to the same
  // system, but within each work unit, depending on the boundary conditions and the non-bonded
  // list type, different tiles may correspond to different systems.
  std::vector<int> nbwu_abstract(tile_groups_wu_abstract_length);
  std::vector<Tcalc> sh_charges(maximum_valence_work_unit_atoms);
  std::vector<int> sh_lj_idx(maximum_valence_work_unit_atoms);
  std::vector<int> sh_n_lj_types(small_block_max_imports);
  std::vector<int> sh_ljabc_offsets(small_block_max_imports);
  std::vector<llint> sh_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zfrc(maximum_valence_work_unit_atoms);
  for (int i = 0; i < synbk.nnbwu; i++) {

    // Import the abstract.
    for (int j = 0; j < 48; j++) {
      nbwu_abstract[j] = synbk[(tile_groups_wu_abstract_length * i) + j];
    }
    const int ntile_sides = nbwu_abstract[0];
    for (int j = 0; j < ntile_sides; j++) {
      const int system_idx = nbwu_abstract[j + 28];
      sh_n_lj_types[j]    = synbk.n_lj_types[system_idx];
      sh_ljabc_offsets[j] = synbk.ljabc_offsets[system_idx];
    }

    // Import atoms into the appropriate arrays
    for (int j = 0; j < ntile_sides; j++) {
      const int atom_start_idx = nbwu_abstract[j + 1];
      const int system_idx     = nbwu_abstract[j + 28];
      const int key_idx        = j / 4;
      const int key_pos        = j - (key_idx * 4);
      const int tside_count    = ((nbwu_abstract[21 + key_idx] >> (8 * key_pos)) & 0xff);
      for (int k = 0; k < tside_count; k++) {
        const size_t localpos = (tile_length * j) + k;
        const size_t synthpos = atom_start_idx + k;
        sh_xcrd[localpos]    = psyw.xcrd[synthpos];
        sh_ycrd[localpos]    = psyw.ycrd[synthpos];
        sh_zcrd[localpos]    = psyw.zcrd[synthpos];
        sh_xfrc[localpos]    = 0LL;
        sh_yfrc[localpos]    = 0LL;
        sh_zfrc[localpos]    = 0LL;
        sh_charges[localpos] = synbk.charge[synthpos];

        // The Lennard-Jones indices recorded here are specific to each system.  For each tile,
        // it is still critical to know the relevant system's total number of Lennard-Jones types
        // as well as its table offset.  On the GPU, these pieces of information will be imported
        // into __shared__ memory for convenient access.
        sh_lj_idx[localpos]  = synbk.lj_idx[synthpos];
      }
    }

  }
}

} // namespace energy
} // namespace omni
