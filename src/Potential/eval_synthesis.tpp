// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2>
void evalSyNonbondedTileGroups(const SyNonbondedKit<Tcalc, Tcalc2> synbk,
                               const SeMaskSynthesisReader syse,
                               PsSynthesisWriter psyw, ScoreCard *ecard,
                               const EvaluateForce eval_elec_force,
                               const EvaluateForce eval_vdw_force) {

  // Critical indexing and offsets for each work unit
  std::vector<int> sh_nbwu_abstract(tile_groups_wu_abstract_length);
  std::vector<int> sh_system_indices(small_block_max_imports);
  std::vector<int> sh_n_lj_types(small_block_max_imports);
  std::vector<int> sh_ljabc_offsets(small_block_max_imports);

  // Pre-computed information for rapidly manipulating the particles of any one tile
  std::vector<Tcalc> sh_tile_xcog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_ycog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_zcog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_tpts(small_block_max_imports);

  // L1-cached coordinates of particles.  These will be accessed repeatedly as the list of tile
  // instructions gets processed.
  std::vector<llint> lc_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> lc_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> lc_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> lc_xcrd_overflow(maximum_valence_work_unit_atoms);
  std::vector<int> lc_ycrd_overflow(maximum_valence_work_unit_atoms);
  std::vector<int> lc_zcrd_overflow(maximum_valence_work_unit_atoms);
  std::vector<Tcalc> lc_charge(maximum_valence_work_unit_atoms);
  std::vector<int>   lc_lj_idx(maximum_valence_work_unit_atoms);

  // Local force accumulators, stored in __shared__ on the GPU to do atomic operations to L1.
  std::vector<llint> sh_xfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zfrc(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xfrc_overflow(maximum_valence_work_unit_atoms);
  std::vector<int> sh_yfrc_overflow(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zfrc_overflow(maximum_valence_work_unit_atoms);

  // Arrays for mocking the register content of various threads
  std::vector<Tcalc> reg_xcrd(tile_length * 2);
  std::vector<Tcalc> reg_ycrd(tile_length * 2);
  std::vector<Tcalc> reg_zcrd(tile_length * 2);
  std::vector<Tcalc> reg_xfrc(tile_length * 2);
  std::vector<Tcalc> reg_yfrc(tile_length * 2);
  std::vector<Tcalc> reg_zfrc(tile_length * 2);
  std::vector<uint>  reg_excl(tile_length * 2);
  std::vector<int>   reg_lj_idx(tile_length * 2);
  std::vector<Tcalc> reg_charge(tile_length * 2);

  // Constants needed for the calculation precision
  const Tcalc value_one = 1.0;
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc sqrt_coul = (tcalc_is_double) ? sqrt(synbk.coulomb) : sqrtf(synbk.coulomb);
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();
  const bool do_elec_force = (eval_elec_force == EvaluateForce::YES);
  const bool do_vdw_force  = (eval_vdw_force == EvaluateForce::YES);
  const bool do_either_force = (do_elec_force || do_vdw_force);
  
  // Initialize the appropriate energy terms in all systems
  for (int i = 0; i < psyw.system_count; i++) {
    ecard->initialize(StateVariable::ELECTROSTATIC, i);
    ecard->initialize(StateVariable::VDW, i);
  }
  
  // Loop over all non-bonded work units within the topology synthesis.  Each holds within it a
  // list of atom imports and a series of tiles.  Within each tile, all atoms pertain to the same
  // system, but within each work unit, depending on the boundary conditions and the non-bonded
  // list type, different tiles may correspond to different systems.
  for (int nbwu_idx = 0; nbwu_idx < synbk.nnbwu; nbwu_idx++) {

    // Import the abstract.
    for (int pos = 0; pos < 48; pos++) {
      sh_nbwu_abstract[pos] = synbk.nbwu_abstracts[(tile_groups_wu_abstract_length * nbwu_idx) +
                                                   pos];
    }
    const int ntile_sides = sh_nbwu_abstract[0];
    for (int pos = 0; pos < ntile_sides; pos++) {
      const int system_idx   = sh_nbwu_abstract[pos + 28];
      sh_system_indices[pos] = system_idx;
      sh_n_lj_types[pos]     = synbk.n_lj_types[system_idx];
      sh_ljabc_offsets[pos]  = synbk.ljabc_offsets[system_idx];
    }
    const int tile_insr_start = sh_nbwu_abstract[small_block_max_imports + 6];
    const int tile_insr_end   = sh_nbwu_abstract[small_block_max_imports + 7];
    
    // Import atoms into the appropriate arrays.  Prepare to compute the center of geometry for
    // all imported atoms.
    for (int pos = 0; pos < ntile_sides; pos++) {
      const int atom_start_idx = sh_nbwu_abstract[pos + 1];
      const int system_idx     = sh_system_indices[pos];
      const int key_idx        = pos / 4;
      const int key_pos        = pos - (key_idx * 4);
      const int tside_count    = ((sh_nbwu_abstract[small_block_max_imports + 1 + key_idx] >>
                                   (8 * key_pos)) & 0xff);
      
      // Pre-compute the centers of geometry for each batch of tile_length atoms, storing the
      // results (totals plus weights) in floating-point format.  When it comes time to do actual
      // tiles, combine the results for the abscissa and ordinate atoms, divide by the combined
      // weight, and use that number to shift the tile atoms to the optimal, precision-preserving
      // locations in the fixed-precision format prior to converting to floating point numbers.
      Tcalc x_cog = 0.0;
      Tcalc y_cog = 0.0;
      Tcalc z_cog = 0.0;
      Tcalc t_pts = 0.0;
      for (int i = 0; i < tside_count; i++) {
        const size_t localpos = (tile_length * pos) + i;
        const size_t synthpos = atom_start_idx + i;
        lc_xcrd[localpos]    = psyw.xcrd[synthpos];
        lc_ycrd[localpos]    = psyw.ycrd[synthpos];
        lc_zcrd[localpos]    = psyw.zcrd[synthpos];
        if (psyw.gpos_bits > globalpos_scale_nonoverflow_bits) {
          lc_xcrd_overflow[localpos] = psyw.xcrd_ovrf[synthpos];
          lc_ycrd_overflow[localpos] = psyw.ycrd_ovrf[synthpos];
          lc_zcrd_overflow[localpos] = psyw.zcrd_ovrf[synthpos];
        }
        sh_xfrc[localpos]    = 0LL;
        sh_yfrc[localpos]    = 0LL;
        sh_zfrc[localpos]    = 0LL;
        if (psyw.frc_bits > force_scale_nonoverflow_bits) {
          sh_xfrc_overflow[localpos] = 0;
          sh_yfrc_overflow[localpos] = 0;
          sh_zfrc_overflow[localpos] = 0;
        }
        
        // Pre-scale all charges by the square root of Coulomb's constant so that it will carry
        // through in all subsequent electrostatic calculations.  On the GPU, the latency involved
        // in this step will likely hide all of the cost of the extra multiplication.
        lc_charge[localpos]  = synbk.charge[synthpos] * sqrt_coul;

        // The Lennard-Jones indices recorded here are specific to each system.  For each tile,
        // it is still critical to know the relevant system's total number of Lennard-Jones types
        // as well as its table offset.  On the GPU, these pieces of information will be imported
        // into __shared__ memory for convenient access.
        lc_lj_idx[localpos]  = synbk.lj_idx[synthpos];

        // Center of geometry computation--coordinates are not scaled to real units at this stage
        if (psyw.gpos_bits > globalpos_scale_nonoverflow_bits) {
          x_cog += hostInt95ToDouble(lc_xcrd[localpos], lc_xcrd_overflow[localpos]);
          y_cog += hostInt95ToDouble(lc_ycrd[localpos], lc_ycrd_overflow[localpos]);
          z_cog += hostInt95ToDouble(lc_zcrd[localpos], lc_zcrd_overflow[localpos]);
        }
        else {
          x_cog += static_cast<Tcalc>(lc_xcrd[localpos]);
          y_cog += static_cast<Tcalc>(lc_ycrd[localpos]);
          z_cog += static_cast<Tcalc>(lc_zcrd[localpos]);
        }
        t_pts += value_one;
      }
      sh_tile_xcog[pos] = x_cog;
      sh_tile_ycog[pos] = y_cog;
      sh_tile_zcog[pos] = z_cog;
      sh_tile_tpts[pos] = t_pts;
    }
    
    // Loop over tile instructions
    for (int pos = tile_insr_start; pos < tile_insr_end; pos++) {
      uint2 tinsr = synbk.nbwu_insr[pos];
      for (int i = 0; i < 2 * tile_length; i++) {
        reg_excl[i] = syse.mask_data[tinsr.y + i];
      }
      const int local_absc_start = (tinsr.x & 0xffff);
      const int local_ordi_start = ((tinsr.x >> 16) & 0xffff);
      const int absc_import_idx = local_absc_start >> tile_length_bits;
      const int ordi_import_idx = local_ordi_start >> tile_length_bits;

      // The system index is needed in order to know where to accumulate the resulting energy
      const int system_idx = sh_system_indices[absc_import_idx];
      
      // On the GPU, the atomic coordinates stored in signed long long integers will be converted
      // to floating point numbers at the beginning of the tile calculation, but to preserve as
      // much information as possible that conversion will also involve centering those atoms on
      // the tile's center of geometry.
      const Tcalc inv_tile_pts = value_one /
                                 (sh_tile_tpts[absc_import_idx] + sh_tile_tpts[ordi_import_idx]);
      const Tcalc tx_cog = (sh_tile_xcog[absc_import_idx] + sh_tile_xcog[ordi_import_idx]) *
                           inv_tile_pts;
      const Tcalc ty_cog = (sh_tile_ycog[absc_import_idx] + sh_tile_ycog[ordi_import_idx]) *
                           inv_tile_pts;
      const Tcalc tz_cog = (sh_tile_zcog[absc_import_idx] + sh_tile_zcog[ordi_import_idx]) *
                           inv_tile_pts;
      if (psyw.gpos_bits > globalpos_scale_nonoverflow_bits) {
        const int95_t x_center = hostDoubleToInt95(tx_cog);
        const int95_t y_center = hostDoubleToInt95(ty_cog);
        const int95_t z_center = hostDoubleToInt95(tz_cog);

        // Based on knowledge of the proper centering, re-import the coordinates after making the
        // shift in precision-preserving integer arithmetic.
        for (int i = 0; i < tile_length; i++) {
          const size_t ilabsc = i + local_absc_start;
          const size_t ilordi = i + local_ordi_start;
          const size_t iplust = i + tile_length;
          const int95_t ix_absc = hostInt95Subtract(lc_xcrd[ilabsc], lc_xcrd_overflow[ilabsc],
                                                    x_center.x, x_center.y);
          const int95_t iy_absc = hostInt95Subtract(lc_ycrd[ilabsc], lc_ycrd_overflow[ilabsc],
                                                    y_center.x, y_center.y);
          const int95_t iz_absc = hostInt95Subtract(lc_zcrd[ilabsc], lc_zcrd_overflow[ilabsc],
                                                    z_center.x, z_center.y);
          reg_xcrd[i]   = hostInt95ToDouble(ix_absc) * psyw.inv_gpos_scale;
          reg_ycrd[i]   = hostInt95ToDouble(iy_absc) * psyw.inv_gpos_scale;
          reg_zcrd[i]   = hostInt95ToDouble(iz_absc) * psyw.inv_gpos_scale;
          reg_lj_idx[i] = lc_lj_idx[ilabsc];
          reg_charge[i] = lc_charge[ilabsc];
          const int95_t ix_ordi = hostInt95Subtract(lc_xcrd[ilordi], lc_xcrd_overflow[ilordi],
                                                    x_center.x, x_center.y);
          const int95_t iy_ordi = hostInt95Subtract(lc_ycrd[ilordi], lc_ycrd_overflow[ilordi],
                                                    y_center.x, y_center.y);
          const int95_t iz_ordi = hostInt95Subtract(lc_zcrd[ilordi], lc_zcrd_overflow[ilordi],
                                                    z_center.x, z_center.y);
          reg_xcrd[iplust] = hostInt95ToDouble(ix_ordi) * psyw.inv_gpos_scale;
          reg_ycrd[iplust] = hostInt95ToDouble(iy_ordi) * psyw.inv_gpos_scale;
          reg_zcrd[iplust] = hostInt95ToDouble(iz_ordi) * psyw.inv_gpos_scale;
          reg_lj_idx[iplust] = lc_lj_idx[ilordi];
          reg_charge[iplust] = lc_charge[ilordi];
        }
      }
      else {
        const llint x_center = static_cast<llint>(tx_cog);
        const llint y_center = static_cast<llint>(ty_cog);
        const llint z_center = static_cast<llint>(tz_cog);

        // Based on knowledge of the proper centering, re-import the coordinates after making the
        // shift in precision-preserving integer arithemtic.
        for (int i = 0; i < tile_length; i++) {
          const size_t ilabsc = i + local_absc_start;
          const size_t ilordi = i + local_ordi_start;
          const size_t iplust = i + tile_length;
          reg_xcrd[i]   = static_cast<Tcalc>(lc_xcrd[ilabsc] - x_center) * psyw.inv_gpos_scale;
          reg_ycrd[i]   = static_cast<Tcalc>(lc_ycrd[ilabsc] - y_center) * psyw.inv_gpos_scale;
          reg_zcrd[i]   = static_cast<Tcalc>(lc_zcrd[ilabsc] - z_center) * psyw.inv_gpos_scale;
          reg_lj_idx[i] = lc_lj_idx[ilabsc];
          reg_charge[i] = lc_charge[ilabsc];
          reg_xcrd[iplust] = static_cast<Tcalc>(lc_xcrd[ilordi] - x_center) * psyw.inv_gpos_scale;
          reg_ycrd[iplust] = static_cast<Tcalc>(lc_ycrd[ilordi] - y_center) * psyw.inv_gpos_scale;
          reg_zcrd[iplust] = static_cast<Tcalc>(lc_zcrd[ilordi] - z_center) * psyw.inv_gpos_scale;
          reg_lj_idx[iplust] = lc_lj_idx[ilordi];
          reg_charge[iplust] = lc_charge[ilordi];
        }
      }
      for (int i = 0; i < 2 * tile_length; i++) {
        reg_xfrc[i] = 0.0;
        reg_yfrc[i] = 0.0;
        reg_zfrc[i] = 0.0;
      }

      // Scan over the entire tile--if the tile runs past the end of the system's atoms, there
      // will be a solid mask of excluded interactions.
      const int nljt = sh_n_lj_types[absc_import_idx];
      const int lj_offset = sh_ljabc_offsets[absc_import_idx];
      Tcalc elec_nrg = 0.0;
      Tcalc vdw_nrg  = 0.0;
      for (int i = 0; i < tile_length; i++) {
        const uint i_mask = reg_excl[i];
        const Tcalc xi = reg_xcrd[i];
        const Tcalc yi = reg_ycrd[i];
        const Tcalc zi = reg_zcrd[i];
        const Tcalc qi = reg_charge[i];
        const int ilj_idx = (nljt * reg_lj_idx[i]) + lj_offset;
        for (int j = tile_length; j < 2 * tile_length; j++) {
          if ((i_mask >> (j - tile_length)) & 0x1) {
            continue;
          }
          const Tcalc dx       = reg_xcrd[j] - xi;
          const Tcalc dy       = reg_ycrd[j] - yi;
          const Tcalc dz       = reg_zcrd[j] - zi;
          const Tcalc dr       = (tcalc_is_double) ?  sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                                     sqrtf((dx * dx) + (dy * dy) + (dz * dz));
          const Tcalc invr     = value_one / dr;
          const Tcalc invr2    = invr * invr;
          const Tcalc invr4    = invr2 * invr2;
          const Tcalc qqij     = reg_charge[j] * qi;
          const int   ij_ljidx = reg_lj_idx[j] + ilj_idx;
          const Tcalc2 ljab    = synbk.ljab_coeff[ij_ljidx];

          // Log the energy.  This is obligatory on the CPU, but the GPU may or may not do it.
          elec_nrg += qqij * invr;
          vdw_nrg  += ((ljab.x * invr4 * invr2) - ljab.y) * invr4 * invr2;
          
          // Compute the forces and contribute them to accumulators.
          if (do_either_force) {
            Tcalc fmag = (do_elec_force) ? -qqij * invr * invr2 : 0.0;
            if (do_vdw_force) {
              if (tcalc_is_double) {
                fmag += ((6.0 * ljab.y) - (12.0 * ljab.x * invr2 * invr4)) * invr4 * invr4;
              }
              else {
                fmag += ((6.0f * ljab.y) - (12.0f * ljab.x * invr2 * invr4)) * invr4 * invr4;
              }
            }
            const Tcalc fmag_dx = fmag * dx;
            const Tcalc fmag_dy = fmag * dy;
            const Tcalc fmag_dz = fmag * dz;
            reg_xfrc[i] += fmag_dx;
            reg_yfrc[i] += fmag_dy;
            reg_zfrc[i] += fmag_dz;
            reg_xfrc[j] -= fmag_dx;
            reg_yfrc[j] -= fmag_dy;
            reg_zfrc[j] -= fmag_dz;
          }
        }
      }

      // There is no need to test whether each work unit is responsible for accumulating the
      // energy computed in a tile.  However, each tile will need to contribute its result to
      // the energy accumulator for a particular system, due to the fact that work units can
      // contain tiles from different systems.
      const llint elec_acc = llround(elec_nrg * nrg_scale_factor);
      ecard->add(StateVariable::ELECTROSTATIC, elec_acc, system_idx);
      const llint vdw_acc  = llround(vdw_nrg * nrg_scale_factor);
      ecard->add(StateVariable::VDW, vdw_acc, system_idx);
      if (do_either_force) {
        for (int i = 0; i < tile_length; i++) {
          const size_t ilabsc = i + local_absc_start;
          const size_t ilordi = i + local_ordi_start;
          const size_t iplust = i + tile_length;
          if (psyw.frc_bits > force_scale_nonoverflow_bits) {
            const int95_t nfx = hostInt95Sum(sh_xfrc[ilabsc], sh_xfrc_overflow[ilabsc],
                                             reg_xfrc[i] * psyw.frc_scale);
            const int95_t nfy = hostInt95Sum(sh_yfrc[ilabsc], sh_yfrc_overflow[ilabsc],
                                             reg_yfrc[i] * psyw.frc_scale);
            const int95_t nfz = hostInt95Sum(sh_zfrc[ilabsc], sh_zfrc_overflow[ilabsc],
                                             reg_zfrc[i] * psyw.frc_scale);
            sh_xfrc[ilabsc] = nfx.x;
            sh_yfrc[ilabsc] = nfy.x;
            sh_zfrc[ilabsc] = nfz.x;
            sh_xfrc_overflow[ilabsc] = nfx.y;
            sh_yfrc_overflow[ilabsc] = nfy.y;
            sh_zfrc_overflow[ilabsc] = nfz.y;
            const int95_t pfx = hostInt95Sum(sh_xfrc[ilordi], sh_xfrc_overflow[ilordi],
                                             reg_xfrc[iplust] * psyw.frc_scale);
            const int95_t pfy = hostInt95Sum(sh_yfrc[ilordi], sh_yfrc_overflow[ilordi],
                                             reg_yfrc[iplust] * psyw.frc_scale);
            const int95_t pfz = hostInt95Sum(sh_zfrc[ilordi], sh_zfrc_overflow[ilordi],
                                             reg_zfrc[iplust] * psyw.frc_scale);
            sh_xfrc[ilordi] = pfx.x;
            sh_yfrc[ilordi] = pfy.x;
            sh_zfrc[ilordi] = pfz.x;
            sh_xfrc_overflow[ilordi] = pfx.y;
            sh_yfrc_overflow[ilordi] = pfy.y;
            sh_zfrc_overflow[ilordi] = pfz.y;
          }
          else {
            sh_xfrc[ilabsc] += llround(reg_xfrc[i] * psyw.frc_scale);
            sh_yfrc[ilabsc] += llround(reg_yfrc[i] * psyw.frc_scale);
            sh_zfrc[ilabsc] += llround(reg_zfrc[i] * psyw.frc_scale);
            sh_xfrc[ilordi] += llround(reg_xfrc[iplust] * psyw.frc_scale);
            sh_yfrc[ilordi] += llround(reg_yfrc[iplust] * psyw.frc_scale);
            sh_zfrc[ilordi] += llround(reg_zfrc[iplust] * psyw.frc_scale);
          }
        }
      }
    }

    // Contribute local force accumulators back to global
    for (int pos = 0; pos < ntile_sides; pos++) {
      const int atom_start_idx = sh_nbwu_abstract[pos + 1];
      const int key_idx        = pos / 4;
      const int key_pos        = pos - (key_idx * 4);
      const int tside_count    = ((sh_nbwu_abstract[21 + key_idx] >> (8 * key_pos)) & 0xff);
      for (int i = 0; i < tside_count; i++) {
        const size_t localpos = (tile_length * pos) + i;
        const size_t synthpos = atom_start_idx + i;
        if (psyw.frc_bits > force_scale_nonoverflow_bits) {
          const int95_t nfx = hostInt95Sum(psyw.xfrc[synthpos], psyw.xfrc_ovrf[synthpos],
                                           sh_xfrc[localpos], sh_xfrc_overflow[localpos]);
          const int95_t nfy = hostInt95Sum(psyw.yfrc[synthpos], psyw.yfrc_ovrf[synthpos],
                                           sh_yfrc[localpos], sh_yfrc_overflow[localpos]);
          const int95_t nfz = hostInt95Sum(psyw.zfrc[synthpos], psyw.zfrc_ovrf[synthpos],
                                           sh_zfrc[localpos], sh_zfrc_overflow[localpos]);
          psyw.xfrc[synthpos] = nfx.x;
          psyw.yfrc[synthpos] = nfy.x;
          psyw.zfrc[synthpos] = nfz.x;
          psyw.xfrc_ovrf[synthpos] = nfx.y;
          psyw.yfrc_ovrf[synthpos] = nfy.y;
          psyw.zfrc_ovrf[synthpos] = nfz.y;
        }
        else {
          psyw.xfrc[synthpos] += sh_xfrc[localpos];
          psyw.yfrc[synthpos] += sh_yfrc[localpos];
          psyw.zfrc[synthpos] += sh_zfrc[localpos];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void evalSyNonbondedEnergy(const AtomGraphSynthesis &poly_ag,
                           const StaticExclusionMaskSynthesis &poly_se,
                           PhaseSpaceSynthesis *poly_ps, ScoreCard *ecard,
                           const EvaluateForce eval_elec_force,
                           const EvaluateForce eval_vdw_force) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  if (tcalc_is_double) {
    evalSyNonbondedTileGroups(poly_ag.getDoublePrecisionNonbondedKit(), poly_se.data(),
                              poly_ps->data(), ecard, eval_elec_force, eval_vdw_force);
  }
  else {
    evalSyNonbondedTileGroups(poly_ag.getSinglePrecisionNonbondedKit(), poly_se.data(),
                              poly_ps->data(), ecard, eval_elec_force, eval_vdw_force);
  }
}

} // namespace energy
} // namespace stormm
