// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                        const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                        const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                        const Tcalc inv_vel_scale) {
  llint result = 0LL;
  const Tcalc my_gafs_to_kcal = gafs_to_kcal;
  if (isSignedIntegralScalarType<Tcoord>()) {
    if (xvel_ovrf == nullptr || yvel_ovrf == nullptr || zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        const Tcalc vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        const Tcalc vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
      const Tcalc contrib_i = mss_i * my_gafs_to_kcal *
                              ((xvel[i] * xvel[i]) + (yvel[i] * yvel[i]) + (zvel[i] * zvel[i]));
      result += llround(contrib_i * nrg_scale_factor);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                       ScoreCard *sc, const ChemicalDetailsKit &cdk, const int system_index,
                       const Tcalc inv_vel_scale) {
  const Tcalc nrg_scale_factor = sc->getEnergyScalingFactor<Tcalc>();
  const llint acc = evalKineticEnergy<Tcoord, double, Tcalc>(xvel, yvel, zvel, nullptr, nullptr,
                                                             nullptr, cdk.masses, cdk.natom,
                                                             nrg_scale_factor, inv_vel_scale);
  sc->contribute(StateVariable::KINETIC, acc, system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tmass, typename Tmass2, typename Tmass4, typename Tcalc>
void evalKineticEnergy(const PsSynthesisReader &poly_psr, ScoreCard *sc,
                       const SyAtomUpdateKit<Tmass, Tmass2, Tmass4> &poly_auk) {
  const Tcalc nrg_scale_factor = sc->getEnergyScalingFactor<Tcalc>();
  for (int i = 0; i < poly_psr.system_count; i++) {
    const size_t iao = poly_psr.atom_starts[i];
    const llint tke = evalKineticEnergy<llint,
                                        Tmass,
                                        Tcalc>(&poly_psr.vxalt[iao], &poly_psr.vyalt[iao],
                                               &poly_psr.vzalt[iao], &poly_psr.vxalt_ovrf[iao],
                                               &poly_psr.vyalt_ovrf[iao],
                                               &poly_psr.vzalt_ovrf[iao],
                                               &poly_auk.masses[iao], poly_psr.atom_counts[i],
                                               nrg_scale_factor, poly_psr.inv_vel_scale);
    sc->contribute(StateVariable::KINETIC, tke, i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalRescaledKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                                const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                                const Tcalc inv_vel_scale) {
  llint result = 0LL;
  const Tcalc my_gafs_to_kcal = gafs_to_kcal;
  if (isSignedIntegralScalarType<Tcoord>()){
    if (xvel_ovrf == nullptr && yvel_ovrf == nullptr && zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++){
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        const Tcalc vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        const Tcalc vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
    }   
    }
    else { 
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
      const Tcalc contrib_i = mss_i * my_gafs_to_kcal *
                              ((xvel[i] * xvel[i]) + (yvel[i] * yvel[i]) + (zvel[i] * zvel[i]));
      result += llround(contrib_i * nrg_scale_factor);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalMomenta(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                  const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                  const Tmass* masses, const int natom, const Tcalc inv_vel_scale) {
  llint result = 0LL;
  if (isSignedIntegralScalarType<Tcoord>()){
    if (xvel_ovrf == nullptr && yvel_ovrf == nullptr && zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++){
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        const Tcalc vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        const Tcalc vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
        result += mss_i * ((vx) + (vy) + (vz));
      }   
    }
    else { 
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
        result += mss_i * ((vx) + (vy) + (vz));
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
      const Tcalc contrib_i = mss_i * ((xvel[i]) + (yvel[i]) + (zvel[i]));
      result += llround(contrib_i);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
Tcalc computeTemperature(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                         const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                         const Tmass* masses, const int natom, const int ndof,
                         const Tcalc nrg_scale_factor, const Tcalc inv_vel_scale) {

  // Begin by computing the kinetic energy
  const llint ke = evalKineticEnergy<Tcoord, Tmass, Tcalc>(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf,
                                                           zvel_ovrf, masses, natom,
                                                           nrg_scale_factor, inv_vel_scale);
  const Tcalc tpre = static_cast<Tcalc>(ke) / (nrg_scale_factor * static_cast<Tcalc>(ndof));
  const Tcalc result = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) ?
                       2.0 * tpre / boltzmann_constant : 2.0f * tpre / boltzmann_constant_f;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace *ps, const AtomGraph *ag,
                         const ApplyConstraints use_cnst) {
  PhaseSpaceReader psr = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  int ndof;
  switch (use_cnst) {
  case ApplyConstraints::YES:
    ndof = cdk.cnst_dof;
    break;
  case ApplyConstraints::NO:
    ndof = cdk.free_dof;
    break;
  }
  return computeTemperature<double, double, Tcalc>(psr.xvel, psr.yvel, psr.zvel, nullptr,
                                                   nullptr, nullptr, cdk.masses, cdk.natom, ndof,
                                                   pow(2.0, 32), 1.0);
}

//------------------------------------------------------------------------------------------------
template <typename Tmass, typename Tmass2, typename Tmass4, typename Tcalc>
void computeTemperature(const PsSynthesisReader &poly_psr, ScoreCard *sc,
                        const SyAtomUpdateKit<Tmass, Tmass2, Tmass4> &poly_auk,
                        const ThermostatReader<Tcalc> &tstr, const bool ke_computed) {
  const Tcalc nrg_scale_factor = sc->getEnergyScalingFactor<Tcalc>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  for (int i = 0; i < poly_psr.system_count; i++) {
    const size_t iao = poly_psr.atom_starts[i];
    const int ndof = (tstr.cnst_geom) ? poly_auk.cnst_dof[i] : poly_auk.free_dof[i];
    Tcalc ttemp;
    if (ke_computed) {
      const llint ke = sc->reportInstantaneousStates(StateVariable::KINETIC, i);
      const Tcalc tpre = static_cast<Tcalc>(ke) / (nrg_scale_factor * static_cast<Tcalc>(ndof));
      ttemp = (tcalc_is_double) ? 2.0 * tpre / boltzmann_constant :
                                  2.0f * tpre / boltzmann_constant_f;
    }
    else {
      ttemp = computeTemperature<llint,
                                 Tmass, Tcalc>(&poly_psr.vxalt[iao], &poly_psr.vyalt[iao],
                                               &poly_psr.vzalt[iao], &poly_psr.vxalt_ovrf[iao],
                                               &poly_psr.vyalt_ovrf[iao],
                                               &poly_psr.vzalt_ovrf[iao], &poly_auk.masses[iao],
                                               poly_psr.atom_counts[i], ndof, nrg_scale_factor,
                                               poly_psr.inv_vel_scale);
    }
    const llint ittemp = llround(ttemp * nrg_scale_factor);
    sc->contribute(StateVariable::TEMPERATURE_ALL, ittemp, i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace &ps, const AtomGraph &ag,
                         const ApplyConstraints use_cnst) {
  return computeTemperature<Tcalc>(ps.getSelfPointer(), ag.getSelfPointer(), use_cnst);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
std::vector<llint> initiateMomentaRescale(const PhaseSpaceSynthesis *PsSynthesis,
                                          const AtomGraphSynthesis *AgSynthesis,
                                          const ScoreCard *sc, std::vector<int> swap_indices,
                                          std::vector<double> temps, bool flag) {
  std::vector<llint> rescaledMomenta;
  const PsSynthesisReader pssr = PsSynthesis->data();
  const std::vector<AtomGraph*> topologies = AgSynthesis->getSystemTopologyPointer();
  const int system_count = pssr.system_count;
  const int* atom_starts = pssr.atom_starts;
  for (int i = 0; i < system_count; i++) {
    llint* xvel = const_cast<llint*>(&pssr.xvel[atom_starts[i]]);
    llint* yvel = const_cast<llint*>(&pssr.yvel[atom_starts[i]]);
    llint* zvel = const_cast<llint*>(&pssr.zvel[atom_starts[i]]);
    const int* xvel_ovrf = nullptr;
    const int* yvel_ovrf = nullptr;
    const int* zvel_ovrf = nullptr;
    if (pssr.vel_bits > velocity_scale_nonoverflow_bits) {
      xvel_ovrf = &pssr.xvel_ovrf[atom_starts[i]];
      yvel_ovrf = &pssr.yvel_ovrf[atom_starts[i]];
      zvel_ovrf = &pssr.zvel_ovrf[atom_starts[i]];
    }
    const AtomGraph* top = topologies[pssr.unique_ag_idx[atom_starts[i]]];
    const ChemicalDetailsKit cdk = top->getChemicalDetailsKit();
    const int natom = pssr.atom_counts[atom_starts[i]];
    const int inv_vel_scale = pssr.inv_vel_scale;
    if (flag) {
      if (i % 2 == 0 && i < system_count - 1) {

        // Even i: rescale xvel, yvel, zvel by sqrt(temps[i + 1] / temps[i])
        double scale_factor = sqrt(temps[i + 1] / temps[i]);
        for (int j = 0; j < natom; j++) {
          xvel[j] *= scale_factor;
          yvel[j] *= scale_factor;
          zvel[j] *= scale_factor;
        }
      }
      else if (i % 2 != 0 && i > 0) {

        // Odd i: rescale xvel, yvel, zvel by sqrt(temps[i] / temps[i - 1])
        double scale_factor = sqrt(temps[i] / temps[i - 1]);
        for (int j = 0; j < natom; j++) {
          xvel[j] *= scale_factor;
          yvel[j] *= scale_factor;
          zvel[j] *= scale_factor;
        }
      }
    }
    else {
      if (i % 2 != 0 && i < system_count - 1) {

        // Odd i: rescale xvel, yvel, zvel by sqrt(temps[i + 1] / temps[i])
        double scale_factor = sqrt(temps[i + 1] / temps[i]);
        for (int j = 0; j < natom; j++) {
          xvel[j] *= scale_factor;
          yvel[j] *= scale_factor;
          zvel[j] *= scale_factor;
        }
      }
      else if (i % 2 == 0 && i > 0) {

        // Even i: rescale xvel, yvel, zvel by sqrt(temps[i] / temps[i - 1])
        double scale_factor = sqrt(temps[i] / temps[i - 1]);
        for (int j = 0; j < natom; j++) {
          xvel[j] *= scale_factor;
          yvel[j] *= scale_factor;
          zvel[j] *= scale_factor;
        }
      }
    }

    // Calculate momenta after rescaling
    rescaledMomenta.push_back(evalMomenta(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf,
                                          zvel_ovrf, cdk.masses, natom, inv_vel_scale));
  }

  // Perform swaps of momenta as per swap_indices and temperature scaling
  for (int i = 0; i < static_cast<int>(swap_indices.size()); i += 2) {
    const llint p1 = rescaledMomenta[swap_indices[i]];
    const llint p2 = rescaledMomenta[swap_indices[i + 1]];
    const double t1 = temps[swap_indices[i]];
    const double t2 = temps[swap_indices[i + 1]];
    rescaledMomenta[swap_indices[i]] = static_cast<llint>(sqrt(t2 / t1) * p1);
    rescaledMomenta[swap_indices[i + 1]] = static_cast<llint>(sqrt(t1 / t2) * p2);
  }
  return rescaledMomenta;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
std::vector<llint> initiateKineticEnergyRescale(const PhaseSpaceSynthesis *PsSynthesis,
                                                const AtomGraphSynthesis *AgSynthesis,
                                                const ScoreCard *sc, std::vector<double> temps,
                                                bool flag) {
  std::vector<llint> results;
  const PsSynthesisReader pssr = PsSynthesis->data();
  const std::vector<AtomGraph*> topologies = AgSynthesis->getSystemTopologyPointer();
  const int system_count = pssr.system_count;
  const int* atom_starts = pssr.atom_starts;

  for (int i = 0; i < system_count; i++) {
    const llint* xvel = &pssr.xvel[atom_starts[i]];
    const llint* yvel = &pssr.yvel[atom_starts[i]];
    const llint* zvel = &pssr.zvel[atom_starts[i]];
    const int* xvel_ovrf = nullptr;
    const int* yvel_ovrf = nullptr;
    const int* zvel_ovrf = nullptr;

    if (pssr.vel_bits >= velocity_scale_nonoverflow_bits) {
      xvel_ovrf = &pssr.xvel_ovrf[atom_starts[i]];
      yvel_ovrf = &pssr.yvel_ovrf[atom_starts[i]];
      zvel_ovrf = &pssr.zvel_ovrf[atom_starts[i]];
    }

    const AtomGraph* top = topologies[pssr.unique_ag_idx[atom_starts[i]]];
    const ChemicalDetailsKit cdk = top->getChemicalDetailsKit();
    const int natom = pssr.atom_counts[atom_starts[i]];
    const double inv_vel_scale = pssr.inv_vel_scale;

    double scale_factor = 1.0;

    if (flag) {
      // Even i: rescale by sqrt(temps[i + 1] / temps[i])
      if (i % 2 == 0 && i < system_count - 1) {
        scale_factor = sqrt(temps[i + 1] / temps[i]);
      }
      // Odd i: rescale by sqrt(temps[i] / temps[i - 1])
      else if (i % 2 != 0 && i > 0) {
        scale_factor = sqrt(temps[i] / temps[i - 1]);
      }
    } else {
        // Odd i: rescale by sqrt(temps[i + 1] / temps[i])
        if (i % 2 != 0 && i < system_count - 1) {
          scale_factor = sqrt(temps[i + 1] / temps[i]);
        }
        // Even i: rescale by sqrt(temps[i] / temps[i - 1])
        else if (i % 2 == 0 && i > 0) {
          scale_factor = sqrt(temps[i] / temps[i - 1]);
        }
    }

    // Compute the kinetic energy after rescaling
    llint result = evalRescaledKineticEnergy(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf, zvel_ovrf,
                                             cdk.masses, natom, pow(2.0, 32), 
                                             inv_vel_scale * scale_factor);
    results.push_back(result);
    }
  return results;
}

                      
} // namespace mm
} // namespace stormm
