// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void synthesisVwuEvaluation(const SyValenceKit<Tcalc> syvk, const Tcalc* sh_charges,
                            const int* sh_lj_idx, llint* sh_xcrd, llint* sh_ycrd, llint* sh_zcrd,
                            llint* sh_xvel, llint* sh_yvel, llint* sh_zvel, llint* sh_xfrc,
                            llint* sh_yfrc, llint* sh_zfrc, const double inv_gpos_scale,
                            const double force_scale, ScoreCard *ecard, const int vwu_idx,
                            const EvaluateForce eval_force, const VwuTask activity,
                            const VwuGoal purpose, const int step_number) {

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
      case VwuGoal::ACCUMULATE_FORCES:
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
      case VwuGoal::ACCUMULATE_FORCES:
        if (readBitFromMask(&syvk.angl_acc[angl_limits.x], pos - angl_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.angl_acc[angl_limits.x], pos - angl_limits.x);
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
    const int2 cdhe_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                static_cast<int>(VwuAbstractMap::CDHE)];
    const int2 cdhe_nrg_limits = syvk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                    static_cast<int>(VwuAbstractMap::CDHE_NRG)];
    for (int pos = cdhe_limits.x; pos < cdhe_limits.y; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE_FORCES:
        if (readBitFromMask(&syvk.cdhe_acc[cdhe_limits.x], pos - cdhe_limits.x) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(&syvk.cdhe_acc[cdhe_limits.x], pos - cdhe_limits.x);
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
                                     syvk.ljb_14_coeff, syvk.ljabc_offsets[sysid],
                                     syvk.n_lj_types[sysid], sh_xcrd, sh_ycrd, sh_zcrd, nullptr,
                                     nullptr, UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc,
                                     eval_force, eval_force, inv_gpos_scale, force_scale);
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

    }
  }
  
  // Contribute results as the instantaneous states
  commitVwuEnergies(bond_acc, angl_acc, dihe_acc, impr_acc, ubrd_acc, cimp_acc, cmap_acc,
                    qq14_acc, lj14_acc, rest_acc, sysid, activity, ecard);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void evalSyValenceEnergy(const SyValenceKit<Tcalc> syvk, PsSynthesisWriter psyw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const VwuTask activity,
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
    synthesisVwuEvaluation(syvk, sh_charges.data(), sh_lj_idx.data(), sh_xcrd.data(),
                           sh_ycrd.data(), sh_zcrd.data(), sh_xvel.data(), sh_yvel.data(),
                           sh_zvel.data(), sh_xfrc.data(), sh_yfrc.data(), sh_zfrc.data(),
                           psyw.inv_gpos_scale, psyw.frc_scale, ecard, i, eval_force, activity,
                           purpose, step_number);
  }
}

} // namespace energy
} // namespace omni
