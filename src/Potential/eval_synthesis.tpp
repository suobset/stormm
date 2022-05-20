// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void synthesisVwuEvaluation(const SyValenceKit<Tcalc> syvk, llint* sh_xcrd, llint* sh_ycrd,
                            llint* sh_zcrd, llint* sh_xvel, llint* sh_yvel, llint* sh_zvel,
                            llint* sh_xfrc, llint* sh_yfrc, llint* sh_zfrc,
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
      sh_xcrd[j - atom_limits.x] = psyw.xcrd[atom_idx];
      sh_ycrd[j - atom_limits.x] = psyw.ycrd[atom_idx];
      sh_zcrd[j - atom_limits.x] = psyw.zcrd[atom_idx];
      sh_xfrc[j - atom_limits.x] = 0LL;
      sh_yfrc[j - atom_limits.x] = 0LL;
      sh_zfrc[j - atom_limits.x] = 0LL;
    }
    synthesisVwuEvaluation(syvk, sh_xcrd.data(), sh_ycrd.data(), sh_zcrd.data(), sh_xvel.data(),
                           sh_yvel.data(), sh_zvel.data(), sh_xfrc.data(), sh_yfrc.data(),
                           sh_zfrc.data(), psyw.inv_gpos_scale, psyw.frc_scale, ecard, i,
                           eval_force, activity, purpose, step_number);
  }
}

} // namespace energy
} // namespace omni
