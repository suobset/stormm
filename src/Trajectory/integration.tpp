// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void velocityVerletVelocityUpdate(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                  const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                  const int natom, const Tcalc* masses, Tcoord* vxalt,
                                  Tcoord* vyalt, Tcoord* vzalt,
                                  const ThermostatReader<Tcalc> &tstr,
                                  const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                                  const int* xfrc_ovrf, const int* yfrc_ovrf, const int* zfrc_ovrf,
                                  int* vxalt_ovrf, int* vyalt_ovrf, int* vzalt_ovrf,
                                  const int atom_offset, const Tcalc vel_scale_factor,
                                  const Tcalc frc_scale_factor) {
  const Tcalc zero = 0.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  Tcalc dt_factor = (tcalc_is_double) ? 0.5  * tstr.dt * kcal_to_gafs :
                                        0.5f * tstr.dt * kcal_to_gafs_f;
  if (tcoord_is_integral) {
    dt_factor *= vel_scale_factor / frc_scale_factor;
  }
  const int atom_limit = atom_offset + natom;
  switch (tstr.kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::BERENDSEN:
    for (int i = atom_offset; i < atom_limit; i++) {
      const Tcalc mss_i = masses[i];
      const Tcalc hmdt = (mss_i > constants::small) ? dt_factor / mss_i : zero;
      if (tcoord_is_integral) {
        Tcalc xpush, ypush, zpush;
        if (xfrc_ovrf == nullptr) {
          xpush = hmdt * static_cast<Tcalc>(xfrc[i]);
          ypush = hmdt * static_cast<Tcalc>(yfrc[i]);
          zpush = hmdt * static_cast<Tcalc>(zfrc[i]);
        }
        else {
          xpush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(xfrc[i], xfrc_ovrf[i]));
          ypush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(yfrc[i], yfrc_ovrf[i]));
          zpush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(zfrc[i], zfrc_ovrf[i]));
        }
        if (xvel_ovrf == nullptr) {
          vxalt[i] = xvel[i] + static_cast<Tcoord>(llround(xpush));
          vyalt[i] = yvel[i] + static_cast<Tcoord>(llround(ypush));
          vzalt[i] = zvel[i] + static_cast<Tcoord>(llround(zpush));
        }
        else {
          const int95_t vx_next = hostInt95Sum(xvel[i], xvel_ovrf[i], xpush);
          const int95_t vy_next = hostInt95Sum(yvel[i], yvel_ovrf[i], ypush);
          const int95_t vz_next = hostInt95Sum(zvel[i], zvel_ovrf[i], zpush);
          vxalt[i] = vx_next.x;
          vyalt[i] = vy_next.x;
          vzalt[i] = vz_next.x;
          vxalt_ovrf[i] = vx_next.y;
          vyalt_ovrf[i] = vy_next.y;
          vzalt_ovrf[i] = vz_next.y;
        }
      }
      else {
        vxalt[i] = xvel[i] + (hmdt * xfrc[i]);
        vyalt[i] = yvel[i] + (hmdt * yfrc[i]);
        vzalt[i] = zvel[i] + (hmdt * zfrc[i]);
      }
    }
    break;
  case ThermostatKind::LANGEVIN:
    switch (tstr.layout) {
    case ThermostatPartition::COMMON:
      {
        Tcalc sdfac;
        if (tcalc_is_double) {
          sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                       getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        else {
          sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                        getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        if (tcoord_is_integral) {
          sdfac *= frc_scale_factor;
        }
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = atom_offset; i < atom_limit; i++) {
          const Tcalc mss_i = masses[i];
          Tcalc rsd, hmdt;
          if (mss_i > constants::small) {
            rsd = (tcalc_is_double) ? sdfac * sqrt(mss_i) : sdfac * sqrtf(mss_i);
            hmdt = dt_factor / mss_i;
          }
          else {
            rsd = zero;
            hmdt = zero;
          }
          Tcalc xbump, ybump, zbump;
          switch (tstr.rng_mode) {
          case PrecisionModel::DOUBLE:
            xbump = rsd * tstr.cache[rnd_offset +                           i];
            ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          case PrecisionModel::SINGLE:
            xbump = rsd * tstr.sp_cache[rnd_offset +                           i];
            ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          }
          
          // The forces acting on each particle will not be modified with the Langevin bumps in
          // this case.  The same composite forces will be re-computed for the second half update.
          if (tcoord_is_integral) {
            Tcalc fx_comp, fy_comp, fz_comp;
            if (xfrc_ovrf == nullptr) {
              fx_comp = static_cast<Tcalc>(xfrc[i]) + xbump;
              fy_comp = static_cast<Tcalc>(yfrc[i]) + ybump;
              fz_comp = static_cast<Tcalc>(zfrc[i]) + zbump;
            }
            else {
              fx_comp = static_cast<Tcalc>(hostInt95ToDouble(xfrc[i], xfrc_ovrf[i])) + xbump;
              fy_comp = static_cast<Tcalc>(hostInt95ToDouble(yfrc[i], yfrc_ovrf[i])) + ybump;
              fz_comp = static_cast<Tcalc>(hostInt95ToDouble(zfrc[i], zfrc_ovrf[i])) + zbump;
            }
            if (xvel_ovrf == nullptr) {
              vxalt[i] = llround((static_cast<Tcalc>(xvel[i]) + (hmdt * fx_comp)) *
                                 tstr.ln_implicit);
              vyalt[i] = llround((static_cast<Tcalc>(yvel[i]) + (hmdt * fy_comp)) *
                                 tstr.ln_implicit);
              vzalt[i] = llround((static_cast<Tcalc>(zvel[i]) + (hmdt * fz_comp)) *
                                 tstr.ln_implicit);
            }
            else {
              const Tcalc dvx = static_cast<Tcalc>(hostInt95ToDouble(xvel[i], xvel_ovrf[i])) +
                                (hmdt * fx_comp);
              const Tcalc dvy = static_cast<Tcalc>(hostInt95ToDouble(yvel[i], yvel_ovrf[i])) +
                                (hmdt * fy_comp);
              const Tcalc dvz = static_cast<Tcalc>(hostInt95ToDouble(zvel[i], zvel_ovrf[i])) +
                                (hmdt * fz_comp);
              const int95_t vx_next = hostDoubleToInt95(dvx * tstr.ln_implicit);
              const int95_t vy_next = hostDoubleToInt95(dvy * tstr.ln_implicit);
              const int95_t vz_next = hostDoubleToInt95(dvz * tstr.ln_implicit);
              vxalt[i] = vx_next.x;
              vyalt[i] = vy_next.x;
              vzalt[i] = vz_next.x;
              vxalt_ovrf[i] = vx_next.y;
              vyalt_ovrf[i] = vy_next.y;
              vzalt_ovrf[i] = vz_next.y;
            }
          }
          else {
            vxalt[i] = (xvel[i] + (hmdt * (xfrc[i] + xbump))) * tstr.ln_implicit;
            vyalt[i] = (yvel[i] + (hmdt * (yfrc[i] + ybump))) * tstr.ln_implicit;
            vzalt[i] = (zvel[i] + (hmdt * (zfrc[i] + zbump))) * tstr.ln_implicit;
          }
        }
      }
      break;
    case ThermostatPartition::SYSTEMS:
    case ThermostatPartition::ATOMS:
      {
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < tstr.npart; i++) {
          const int4 prtn = tstr.partition_bounds[i];

          // Skip segments of the thermostating that are not within the current limits.
          if (prtn.y < atom_offset || prtn.x >= atom_limit) {
            continue;
          }
          Tcalc sdfac;
          if (tcalc_is_double) {
            sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                         getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          else {
            sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                          getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          if (tcoord_is_integral) {
            sdfac *= frc_scale_factor;
          }
          for (int j = prtn.x; j < prtn.y; j++) {

            // Check that each atom in the partition is within the current system before applying
            // the thermostat effect to it.
            if (j < atom_offset || j >= atom_limit) {
              continue;
            }
            const Tcalc mss_j = masses[j];
            Tcalc rsd, hmdt;
            if (mss_j > constants::small) {
              rsd = (tcalc_is_double) ? sdfac * sqrt(mss_j) : sdfac * sqrtf(mss_j);
              hmdt = dt_factor / mss_j;
            }
            else {
              rsd = zero;
              hmdt = zero;
            }
            Tcalc xbump, ybump, zbump;
            switch (tstr.rng_mode) {
            case PrecisionModel::DOUBLE:
              xbump = rsd * tstr.cache[rnd_offset +                           j];
              ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            case PrecisionModel::SINGLE:
              xbump = rsd * tstr.sp_cache[rnd_offset +                           j];
              ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            }
            if (tcoord_is_integral) {
              Tcalc fx_comp, fy_comp, fz_comp;
              if (xfrc_ovrf == nullptr) {
                fx_comp = static_cast<Tcalc>(xfrc[j]) + xbump;
                fy_comp = static_cast<Tcalc>(yfrc[j]) + ybump;
                fz_comp = static_cast<Tcalc>(zfrc[j]) + zbump;
              }
              else {
                fx_comp = static_cast<Tcalc>(hostInt95ToDouble(xfrc[j], xfrc_ovrf[j])) + xbump;
                fy_comp = static_cast<Tcalc>(hostInt95ToDouble(yfrc[j], yfrc_ovrf[j])) + ybump;
                fz_comp = static_cast<Tcalc>(hostInt95ToDouble(zfrc[j], zfrc_ovrf[j])) + zbump;
              }
              if (xvel_ovrf == nullptr) {
                vxalt[j] = llround((static_cast<Tcalc>(xvel[j]) + (hmdt * fx_comp)) *
                                   tstr.ln_implicit);
                vyalt[j] = llround((static_cast<Tcalc>(yvel[j]) + (hmdt * fy_comp)) *
                                   tstr.ln_implicit);
                vzalt[j] = llround((static_cast<Tcalc>(zvel[j]) + (hmdt * fz_comp)) *
                                   tstr.ln_implicit);
              }
              else {
                const Tcalc dvx = static_cast<Tcalc>(hostInt95ToDouble(xvel[j], xvel_ovrf[j])) +
                                  (hmdt * fx_comp);
                const Tcalc dvy = static_cast<Tcalc>(hostInt95ToDouble(yvel[j], yvel_ovrf[j])) +
                                  (hmdt * fy_comp);
                const Tcalc dvz = static_cast<Tcalc>(hostInt95ToDouble(zvel[j], zvel_ovrf[j])) +
                                  (hmdt * fz_comp);
                const int95_t vx_next = hostDoubleToInt95(dvx * tstr.ln_implicit);
                const int95_t vy_next = hostDoubleToInt95(dvy * tstr.ln_implicit);
                const int95_t vz_next = hostDoubleToInt95(dvz * tstr.ln_implicit);
                vxalt[j] = vx_next.x;
                vyalt[j] = vy_next.x;
                vzalt[j] = vz_next.x;
                vxalt_ovrf[j] = vx_next.y;
                vyalt_ovrf[j] = vy_next.y;
                vzalt_ovrf[j] = vz_next.y;
              }
            }
            else {
              vxalt[j] = (xvel[j] + (hmdt * (xfrc[j] + xbump))) * tstr.ln_implicit;
              vyalt[j] = (yvel[j] + (hmdt * (yfrc[j] + ybump))) * tstr.ln_implicit;
              vzalt[j] = (zvel[j] + (hmdt * (zfrc[j] + zbump))) * tstr.ln_implicit;
            }
          }
        }
      }
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void velocityVerletVelocityUpdate(PsSynthesisWriter *poly_psw,
                                  const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                                  const ThermostatReader<T> &tstr) {
  const int last_sys = poly_psw->system_count - 1;
  const int natom = poly_psw->atom_starts[last_sys] + poly_psw->atom_counts[last_sys];
  if (poly_psw->frc_bits <= force_scale_nonoverflow_bits) {
    if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
      velocityVerletVelocityUpdate<llint, T>(poly_psw->xvel, poly_psw->yvel, poly_psw->zvel,
                                             poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                             natom, poly_auk.masses, poly_psw->vxalt,
                                             poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                             nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                             nullptr, nullptr, 0, poly_psw->vel_scale,
                                             poly_psw->frc_scale);
    }
    else {
      velocityVerletVelocityUpdate<llint, T>(poly_psw->xvel, poly_psw->yvel, poly_psw->zvel,
                                             poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                             natom, poly_auk.masses, poly_psw->vxalt,
                                             poly_psw->vyalt, poly_psw->vzalt, tstr,
                                             poly_psw->xvel_ovrf, poly_psw->yvel_ovrf,
                                             poly_psw->zvel_ovrf, nullptr, nullptr, nullptr,
                                             poly_psw->vxalt_ovrf, poly_psw->vyalt_ovrf,
                                             poly_psw->vzalt_ovrf, 0, poly_psw->vel_scale,
                                             poly_psw->frc_scale);
    }
  }
  else {
    if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
      velocityVerletVelocityUpdate<llint, T>(poly_psw->xvel, poly_psw->yvel, poly_psw->zvel,
                                             poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                             natom, poly_auk.masses, poly_psw->vxalt,
                                             poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                             nullptr, nullptr, poly_psw->xfrc_ovrf,
                                             poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf, nullptr,
                                             nullptr, nullptr, 0, poly_psw->vel_scale,
                                             poly_psw->frc_scale);
    }
    else {
      velocityVerletVelocityUpdate<llint, T>(poly_psw->xvel, poly_psw->yvel, poly_psw->zvel,
                                             poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                             natom, poly_auk.masses, poly_psw->vxalt,
                                             poly_psw->vyalt, poly_psw->vzalt, tstr,
                                             poly_psw->xvel_ovrf, poly_psw->yvel_ovrf,
                                             poly_psw->zvel_ovrf, poly_psw->xfrc_ovrf,
                                             poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf,
                                             poly_psw->vxalt_ovrf, poly_psw->vyalt_ovrf,
                                             poly_psw->vzalt_ovrf, 0, poly_psw->vel_scale,
                                             poly_psw->frc_scale);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcoord4>
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg,
                                  const AtomGraphSynthesis &poly_ag, const Thermostat &tst,
                                  const PrecisionModel prec) {
  cg.contributeForces(poly_ps);
  velocityVerletVelocityUpdate(poly_ps, poly_ag.getSelfPointer(), tst.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcoord4>
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg_qq,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg_lj,
                                  const AtomGraphSynthesis &poly_ag, const Thermostat &tst,
                                  const PrecisionModel prec) {
  cg_qq.contributeForces(poly_ps);
  cg_lj.contributeForces(poly_ps);
  velocityVerletVelocityUpdate(poly_ps, poly_ag.getSelfPointer(), tst.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void velocityVerletCoordinateUpdate(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                    const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                    const int natom, const Tcalc* masses, Tcoord* xalt,
                                    Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt,
                                    Tcoord* vzalt, const ThermostatReader<Tcalc> &tstr,
                                    const int* xcrd_ovrf, const int* ycrd_ovrf,
                                    const int* zcrd_ovrf, const int* xfrc_ovrf,
                                    const int* yfrc_ovrf, const int* zfrc_ovrf,
                                    int* xalt_ovrf, int* yalt_ovrf, int* zalt_ovrf,
                                    int* vxalt_ovrf, int* vyalt_ovrf, int* vzalt_ovrf,
                                    const int atom_offset, const Tcalc gpos_scale_factor,
                                    const Tcalc vel_scale_factor, const Tcalc frc_scale_factor) {
  const Tcalc zero = 0.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  Tcalc dt_factor = (tcalc_is_double) ? 0.5  * tstr.dt * kcal_to_gafs :
                                        0.5f * tstr.dt * kcal_to_gafs_f;
  if (tcoord_is_integral) {
    dt_factor *= vel_scale_factor / frc_scale_factor;
  }
  const int atom_limit = atom_offset + natom;
  switch (tstr.kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::BERENDSEN:
    for (int i = atom_offset; i < atom_limit; i++) {
      const Tcalc mss_i = masses[i];
      const Tcalc hmdt = (mss_i > constants::small) ? dt_factor / mss_i : zero;
      if (tcoord_is_integral) {
        Tcalc xpush, ypush, zpush;
        if (xfrc_ovrf == nullptr) {
          xpush = hmdt * static_cast<Tcalc>(xfrc[i]);
          ypush = hmdt * static_cast<Tcalc>(yfrc[i]);
          zpush = hmdt * static_cast<Tcalc>(zfrc[i]);
        }
        else {
          xpush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(xfrc[i], xfrc_ovrf[i]));
          ypush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(yfrc[i], yfrc_ovrf[i]));
          zpush = hmdt * static_cast<Tcalc>(hostInt95ToDouble(zfrc[i], zfrc_ovrf[i]));
        }
        if (vxalt_ovrf == nullptr) {
          vxalt[i] += static_cast<Tcoord>(llround(xpush));
          vyalt[i] += static_cast<Tcoord>(llround(ypush));
          vzalt[i] += static_cast<Tcoord>(llround(zpush));
        }
        else {
          const int95_t vx_next = hostInt95Sum(vxalt[i], vxalt_ovrf[i], xpush);
          const int95_t vy_next = hostInt95Sum(vyalt[i], vyalt_ovrf[i], ypush);
          const int95_t vz_next = hostInt95Sum(vzalt[i], vzalt_ovrf[i], zpush);
          vxalt[i] = vx_next.x;
          vyalt[i] = vy_next.x;
          vzalt[i] = vz_next.x;
          vxalt_ovrf[i] = vx_next.y;
          vyalt_ovrf[i] = vy_next.y;
          vzalt_ovrf[i] = vz_next.y;
        }
      }
      else {
        vxalt[i] += hmdt * xfrc[i];
        vyalt[i] += hmdt * yfrc[i];
        vzalt[i] += hmdt * zfrc[i];
      }
    }
    break;
  case ThermostatKind::LANGEVIN:
    switch (tstr.layout) {
    case ThermostatPartition::COMMON:
      {
        Tcalc sdfac;
        if (tcalc_is_double) {
          sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                       getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        else {
          sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                        getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        if (tcoord_is_integral) {
          sdfac *= frc_scale_factor;
        }
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = atom_offset; i < atom_limit; i++) {
          const Tcalc mss_i = masses[i];
          Tcalc rsd, hmdt;
          if (mss_i > constants::small) {
            rsd = (tcalc_is_double) ? sdfac * sqrt(mss_i) : sdfac * sqrtf(mss_i);
            hmdt = dt_factor / mss_i;
          }
          else {
            rsd = zero;
            hmdt = zero;
          }
          Tcalc xbump, ybump, zbump;
          switch (tstr.rng_mode) {
          case PrecisionModel::DOUBLE:
            xbump = rsd * tstr.cache[rnd_offset +                           i];
            ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          case PrecisionModel::SINGLE:
            xbump = rsd * tstr.sp_cache[rnd_offset +                           i];
            ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          }
          if (tcoord_is_integral) {
            Tcalc fx_comp, fy_comp, fz_comp;
            if (xfrc_ovrf == nullptr) {
              fx_comp = static_cast<Tcalc>(xfrc[i]) + xbump;
              fy_comp = static_cast<Tcalc>(yfrc[i]) + ybump;
              fz_comp = static_cast<Tcalc>(zfrc[i]) + zbump;
            }
            else {
              fx_comp = static_cast<Tcalc>(hostInt95ToDouble(xfrc[i], xfrc_ovrf[i])) + xbump;
              fy_comp = static_cast<Tcalc>(hostInt95ToDouble(yfrc[i], yfrc_ovrf[i])) + ybump;
              fz_comp = static_cast<Tcalc>(hostInt95ToDouble(zfrc[i], zfrc_ovrf[i])) + zbump;
            }
            if (vxalt_ovrf == nullptr) {
              vxalt[i] = llround((static_cast<Tcalc>(vxalt[i]) * tstr.ln_explicit) +
                                 (hmdt * fx_comp));
              vyalt[i] = llround((static_cast<Tcalc>(vyalt[i]) * tstr.ln_explicit) +
                                 (hmdt * fy_comp));
              vzalt[i] = llround((static_cast<Tcalc>(vzalt[i]) * tstr.ln_explicit) +
                                 (hmdt * fz_comp));
            }
            else {
              const Tcalc vx_next = (static_cast<Tcalc>(hostInt95ToDouble(vxalt[i],
                                                                          vxalt_ovrf[i])) *
                                     tstr.ln_explicit) + (hmdt * fx_comp);
              const Tcalc vy_next = (static_cast<Tcalc>(hostInt95ToDouble(vyalt[i],
                                                                          vyalt_ovrf[i])) *
                                     tstr.ln_explicit) + (hmdt * fy_comp);
              const Tcalc vz_next = (static_cast<Tcalc>(hostInt95ToDouble(vzalt[i],
                                                                          vzalt_ovrf[i])) *
                                     tstr.ln_explicit) + (hmdt * fz_comp);
              const int95_t ivx_next = hostDoubleToInt95(vx_next);
              const int95_t ivy_next = hostDoubleToInt95(vy_next);
              const int95_t ivz_next = hostDoubleToInt95(vz_next);
              vxalt[i] = ivx_next.x;
              vyalt[i] = ivy_next.x;
              vzalt[i] = ivz_next.x;
              vxalt_ovrf[i] = ivx_next.y;
              vyalt_ovrf[i] = ivy_next.y;
              vzalt_ovrf[i] = ivz_next.y;
            }
          }
          else {
            vxalt[i] = (vxalt[i] * tstr.ln_explicit) + (hmdt * (xfrc[i] + xbump));
            vyalt[i] = (vyalt[i] * tstr.ln_explicit) + (hmdt * (yfrc[i] + ybump));
            vzalt[i] = (vzalt[i] * tstr.ln_explicit) + (hmdt * (zfrc[i] + zbump));
          }
        }
      }
      break;
    case ThermostatPartition::SYSTEMS:
    case ThermostatPartition::ATOMS:
      {
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < tstr.npart; i++) {
          const int4 prtn = tstr.partition_bounds[i];

          // Skip segments of the thermostating that are not within the current limits.
          if (prtn.y < atom_offset || prtn.x >= atom_limit) {
            continue;
          }
          Tcalc sdfac;
          if (tcalc_is_double) {
            sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                         getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          else {
            sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                          getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          if (tcoord_is_integral) {
            sdfac *= frc_scale_factor;
          }
          for (int j = prtn.x; j < prtn.y; j++) {

            // Check that each atom in the partition is within the current system before applying
            // the thermostat effect to it.
            if (j < atom_offset || j >= atom_limit) {
              continue;
            }
            const Tcalc mss_j = masses[j];
            Tcalc rsd, hmdt;
            if (mss_j > constants::small) {
              rsd = (tcalc_is_double) ? sdfac * sqrt(mss_j) : sdfac * sqrtf(mss_j);
              hmdt = dt_factor / mss_j;
            }
            else {
              rsd = zero;
              hmdt = zero;
            }
            Tcalc xbump, ybump, zbump;
            switch (tstr.rng_mode) {
            case PrecisionModel::DOUBLE:
              xbump = rsd * tstr.cache[rnd_offset +                           j];
              ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            case PrecisionModel::SINGLE:
              xbump = rsd * tstr.sp_cache[rnd_offset +                           j];
              ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            }
            if (tcoord_is_integral) {
              Tcalc fx_comp, fy_comp, fz_comp;
              if (xfrc_ovrf == nullptr) {
                fx_comp = static_cast<Tcalc>(xfrc[j]) + xbump;
                fy_comp = static_cast<Tcalc>(yfrc[j]) + ybump;
                fz_comp = static_cast<Tcalc>(zfrc[j]) + zbump;
              }
              else {
                fx_comp = static_cast<Tcalc>(hostInt95ToDouble(xfrc[j], xfrc_ovrf[j])) + xbump;
                fy_comp = static_cast<Tcalc>(hostInt95ToDouble(yfrc[j], yfrc_ovrf[j])) + ybump;
                fz_comp = static_cast<Tcalc>(hostInt95ToDouble(zfrc[j], zfrc_ovrf[j])) + zbump;
              }
              if (vxalt_ovrf == nullptr) {
                vxalt[j] = llround((static_cast<Tcalc>(vxalt[j]) * tstr.ln_explicit) +
                                   (hmdt * fx_comp));
                vyalt[j] = llround((static_cast<Tcalc>(vyalt[j]) * tstr.ln_explicit) +
                                   (hmdt * fy_comp));
                vzalt[j] = llround((static_cast<Tcalc>(vzalt[j]) * tstr.ln_explicit) +
                                   (hmdt * fz_comp));
              }
              else {
                const Tcalc vx_next = (static_cast<Tcalc>(hostInt95ToDouble(vxalt[j],
                                                                            vxalt_ovrf[j])) *
                                       tstr.ln_explicit) + (hmdt * fx_comp);
                const Tcalc vy_next = (static_cast<Tcalc>(hostInt95ToDouble(vyalt[j],
                                                                            vyalt_ovrf[j])) *
                                       tstr.ln_explicit) + (hmdt * fy_comp);
                const Tcalc vz_next = (static_cast<Tcalc>(hostInt95ToDouble(vzalt[j],
                                                                            vzalt_ovrf[j])) *
                                       tstr.ln_explicit) + (hmdt * fz_comp);
                const int95_t ivx_next = hostDoubleToInt95(vx_next);
                const int95_t ivy_next = hostDoubleToInt95(vy_next);
                const int95_t ivz_next = hostDoubleToInt95(vz_next);
                vxalt[j] = ivx_next.x;
                vyalt[j] = ivy_next.x;
                vzalt[j] = ivz_next.x;
                vxalt_ovrf[j] = ivx_next.y;
                vyalt_ovrf[j] = ivy_next.y;
                vzalt_ovrf[j] = ivz_next.y;
              }
            }
            else {
              vxalt[j] = (vxalt[j] * tstr.ln_explicit) + (hmdt * (xfrc[j] + xbump));
              vyalt[j] = (vyalt[j] * tstr.ln_explicit) + (hmdt * (yfrc[j] + ybump));
              vzalt[j] = (vzalt[j] * tstr.ln_explicit) + (hmdt * (zfrc[j] + zbump));
            }
          }
        }
      }
      break;
    }
    break;
  }

  // With the complete velocity update, move the particles forward.
  if (tcoord_is_integral) {
    const Tcalc vcdt_factor = tstr.dt * gpos_scale_factor / vel_scale_factor;
    for (int i = atom_offset; i < atom_limit; i++) {
      Tcalc x_updt, y_updt, z_updt;
      if (vxalt_ovrf == nullptr) {
        x_updt = static_cast<Tcalc>(vxalt[i]) * vcdt_factor;
        y_updt = static_cast<Tcalc>(vyalt[i]) * vcdt_factor;
        z_updt = static_cast<Tcalc>(vzalt[i]) * vcdt_factor;
      }
      else {
        x_updt = static_cast<Tcalc>(hostInt95ToDouble(vxalt[i], vxalt_ovrf[i])) * vcdt_factor;
        y_updt = static_cast<Tcalc>(hostInt95ToDouble(vyalt[i], vyalt_ovrf[i])) * vcdt_factor;
        z_updt = static_cast<Tcalc>(hostInt95ToDouble(vzalt[i], vzalt_ovrf[i])) * vcdt_factor;
      }
      if (xcrd_ovrf == nullptr) {
        xalt[i] = xcrd[i] + static_cast<Tcoord>(llround(x_updt));
        yalt[i] = ycrd[i] + static_cast<Tcoord>(llround(y_updt));
        zalt[i] = zcrd[i] + static_cast<Tcoord>(llround(z_updt));
      }
      else {
        const int95_t ix_next = hostInt95Sum(xcrd[i], xcrd_ovrf[i], x_updt);
        const int95_t iy_next = hostInt95Sum(ycrd[i], ycrd_ovrf[i], y_updt);
        const int95_t iz_next = hostInt95Sum(zcrd[i], zcrd_ovrf[i], z_updt);
        xalt[i] = ix_next.x;
        yalt[i] = iy_next.x;
        zalt[i] = iz_next.x;
        xalt_ovrf[i] = ix_next.y;
        yalt_ovrf[i] = iy_next.y;
        zalt_ovrf[i] = iz_next.y;
      }
    }
  }
  else {
    for (int i = atom_offset; i < atom_limit; i++) {
      xalt[i] = xcrd[i] + (vxalt[i] * tstr.dt);
      yalt[i] = ycrd[i] + (vyalt[i] * tstr.dt);
      zalt[i] = zcrd[i] + (vzalt[i] * tstr.dt);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void velocityVerletCoordinateUpdate(PsSynthesisWriter *poly_psw,
                                    const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                                    const ThermostatReader<T> &tstr) {
  const int last_sys = poly_psw->system_count - 1;
  const int natom = poly_psw->atom_starts[last_sys] + poly_psw->atom_counts[last_sys];
  if (poly_psw->frc_bits <= force_scale_nonoverflow_bits) {
    if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
      if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                                 nullptr, 0, poly_psw->gpos_scale,
                                                 poly_psw->vel_scale, poly_psw->frc_scale);
      }
      else {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses,
                                                 poly_psw->xalt, poly_psw->yalt, poly_psw->zalt,
                                                 poly_psw->vxalt, poly_psw->vyalt, poly_psw->vzalt,
                                                 tstr, poly_psw->xcrd_ovrf, poly_psw->ycrd_ovrf,
                                                 poly_psw->zcrd_ovrf, nullptr, nullptr, nullptr,
                                                 poly_psw->xalt_ovrf, poly_psw->yalt_ovrf,
                                                 poly_psw->zalt_ovrf, nullptr, nullptr, nullptr,
                                                 0, poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
    }
    else {
      if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr,
                                                 nullptr, nullptr, nullptr, poly_psw->vxalt_ovrf,
                                                 poly_psw->vyalt_ovrf, poly_psw->vzalt_ovrf, 0,
                                                 poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
      else {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr,
                                                 poly_psw->xcrd_ovrf, poly_psw->ycrd_ovrf,
                                                 poly_psw->zcrd_ovrf, nullptr, nullptr, nullptr,
                                                 poly_psw->xalt_ovrf, poly_psw->yalt_ovrf,
                                                 poly_psw->zalt_ovrf, poly_psw->vxalt_ovrf,
                                                 poly_psw->vyalt_ovrf, poly_psw->vzalt_ovrf,
                                                 0, poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
    }
  }
  else {
    if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
      if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                                 nullptr, nullptr, poly_psw->xfrc_ovrf,
                                                 poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf, nullptr,
                                                 nullptr, nullptr, nullptr, nullptr, nullptr, 0,
                                                 poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
      else {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr,
                                                 poly_psw->xcrd_ovrf, poly_psw->ycrd_ovrf,
                                                 poly_psw->zcrd_ovrf, poly_psw->xfrc_ovrf,
                                                 poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf,
                                                 poly_psw->xalt_ovrf, poly_psw->yalt_ovrf,
                                                 poly_psw->zalt_ovrf, nullptr, nullptr, nullptr, 0,
                                                 poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
    }
    else {
      if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr, nullptr,
                                                 nullptr, nullptr, poly_psw->xfrc_ovrf,
                                                 poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf,
                                                 poly_psw->xalt_ovrf, poly_psw->yalt_ovrf,
                                                 poly_psw->zalt_ovrf, poly_psw->vxalt_ovrf,
                                                 poly_psw->vyalt_ovrf, poly_psw->vzalt_ovrf, 0,
                                                 poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
      else {
        velocityVerletCoordinateUpdate<llint, T>(poly_psw->xcrd, poly_psw->ycrd, poly_psw->zcrd,
                                                 poly_psw->xfrc, poly_psw->yfrc, poly_psw->zfrc,
                                                 natom, poly_auk.masses, poly_psw->xalt,
                                                 poly_psw->yalt, poly_psw->zalt, poly_psw->vxalt,
                                                 poly_psw->vyalt, poly_psw->vzalt, tstr,
                                                 poly_psw->xcrd_ovrf, poly_psw->ycrd_ovrf,
                                                 poly_psw->zcrd_ovrf, poly_psw->xfrc_ovrf,
                                                 poly_psw->yfrc_ovrf, poly_psw->zfrc_ovrf,
                                                 poly_psw->xalt_ovrf, poly_psw->yalt_ovrf,
                                                 poly_psw->zalt_ovrf, poly_psw->vxalt_ovrf,
                                                 poly_psw->vyalt_ovrf, poly_psw->vzalt_ovrf, 0,
                                                 poly_psw->gpos_scale, poly_psw->vel_scale,
                                                 poly_psw->frc_scale);
      }
    }
  }
}


} // namespace structure
} // namespace stormm
