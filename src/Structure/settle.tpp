// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void settleGroupVelocity(const int oxy_idx, const int hd1_idx, const int hd2_idx, const Tcalc m_o,
                         const Tcalc m_h, const Tcalc m_oh, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, Tcoord* xvel, Tcoord* yvel, Tcoord* zvel,
                         const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                         int* xvel_ovrf, int* yvel_ovrf, int* zvel_ovrf, const Tcalc gpos_scale,
                         const Tcalc vel_scale) {
  const Tcalc value_one  = 1.0;
  const Tcalc value_two  = 2.0;
  const Tcalc value_half = 0.5;
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  Tcalc doh_x, doh_y, doh_z, dhh_x, dhh_y, dhh_z, dho_x, dho_y, dho_z;
  if (tcoord_is_integral) {
    if (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr && zcrd_ovrf != nullptr) {
      const int95_t idoh_x = hostInt95Subtract(xcrd[hd1_idx], xcrd_ovrf[hd1_idx], xcrd[oxy_idx],
                                               xcrd_ovrf[oxy_idx]);
      const int95_t idoh_y = hostInt95Subtract(ycrd[hd1_idx], ycrd_ovrf[hd1_idx], ycrd[oxy_idx],
                                               ycrd_ovrf[oxy_idx]);
      const int95_t idoh_z = hostInt95Subtract(zcrd[hd1_idx], zcrd_ovrf[hd1_idx], zcrd[oxy_idx],
                                               zcrd_ovrf[oxy_idx]);
      doh_x = hostInt95ToDouble(idoh_x) / gpos_scale;
      doh_y = hostInt95ToDouble(idoh_y) / gpos_scale;
      doh_z = hostInt95ToDouble(idoh_z) / gpos_scale;
      const int95_t idhh_x = hostInt95Subtract(xcrd[hd2_idx], xcrd_ovrf[hd2_idx], xcrd[hd1_idx],
                                               xcrd_ovrf[hd1_idx]);
      const int95_t idhh_y = hostInt95Subtract(ycrd[hd2_idx], ycrd_ovrf[hd2_idx], ycrd[hd1_idx],
                                               ycrd_ovrf[hd1_idx]);
      const int95_t idhh_z = hostInt95Subtract(zcrd[hd2_idx], zcrd_ovrf[hd2_idx], zcrd[hd1_idx],
                                               zcrd_ovrf[hd1_idx]);
      dhh_x = hostInt95ToDouble(idhh_x) / gpos_scale;
      dhh_y = hostInt95ToDouble(idhh_y) / gpos_scale;
      dhh_z = hostInt95ToDouble(idhh_z) / gpos_scale;
      const int95_t idho_x = hostInt95Subtract(xcrd[oxy_idx], xcrd_ovrf[oxy_idx], xcrd[hd2_idx],
                                               xcrd_ovrf[hd2_idx]);
      const int95_t idho_y = hostInt95Subtract(ycrd[oxy_idx], ycrd_ovrf[oxy_idx], ycrd[hd2_idx],
                                               ycrd_ovrf[hd2_idx]);
      const int95_t idho_z = hostInt95Subtract(zcrd[oxy_idx], zcrd_ovrf[oxy_idx], zcrd[hd2_idx],
                                               zcrd_ovrf[hd2_idx]);
      dho_x = hostInt95ToDouble(idho_x) / gpos_scale;
      dho_y = hostInt95ToDouble(idho_y) / gpos_scale;
      dho_z = hostInt95ToDouble(idho_z) / gpos_scale;
    }
    else {
      doh_x = static_cast<Tcalc>(xcrd[hd1_idx] - xcrd[oxy_idx]) / gpos_scale;
      doh_y = static_cast<Tcalc>(ycrd[hd1_idx] - ycrd[oxy_idx]) / gpos_scale;
      doh_z = static_cast<Tcalc>(zcrd[hd1_idx] - zcrd[oxy_idx]) / gpos_scale;
      dhh_x = static_cast<Tcalc>(xcrd[hd2_idx] - xcrd[hd1_idx]) / gpos_scale;
      dhh_y = static_cast<Tcalc>(ycrd[hd2_idx] - ycrd[hd1_idx]) / gpos_scale;
      dhh_z = static_cast<Tcalc>(zcrd[hd2_idx] - zcrd[hd1_idx]) / gpos_scale;
      dho_x = static_cast<Tcalc>(xcrd[oxy_idx] - xcrd[hd2_idx]) / gpos_scale;
      dho_y = static_cast<Tcalc>(ycrd[oxy_idx] - ycrd[hd2_idx]) / gpos_scale;
      dho_z = static_cast<Tcalc>(zcrd[oxy_idx] - zcrd[hd2_idx]) / gpos_scale;
    }
  }
  else {
    doh_x = xcrd[hd1_idx] - xcrd[oxy_idx];
    doh_y = ycrd[hd1_idx] - ycrd[oxy_idx];
    doh_z = zcrd[hd1_idx] - zcrd[oxy_idx];
    dhh_x = xcrd[hd2_idx] - xcrd[hd1_idx];
    dhh_y = ycrd[hd2_idx] - ycrd[hd1_idx];
    dhh_z = zcrd[hd2_idx] - zcrd[hd1_idx];
    dho_x = xcrd[oxy_idx] - xcrd[hd2_idx];
    dho_y = ycrd[oxy_idx] - ycrd[hd2_idx];
    dho_z = zcrd[oxy_idx] - zcrd[hd2_idx];
  }
  Tcalc inv_mag_doh, inv_mag_dhh, inv_mag_dho;
  if (tcalc_is_double) {
    inv_mag_doh = value_one / sqrt((doh_x * doh_x) + (doh_y * doh_y) + (doh_z * doh_z));
    inv_mag_dhh = value_one / sqrt((dhh_x * dhh_x) + (dhh_y * dhh_y) + (dhh_z * dhh_z));
    inv_mag_dho = value_one / sqrt((dho_x * dho_x) + (dho_y * dho_y) + (dho_z * dho_z));
  }
  else {
    inv_mag_doh = value_one / sqrtf((doh_x * doh_x) + (doh_y * doh_y) + (doh_z * doh_z));
    inv_mag_dhh = value_one / sqrtf((dhh_x * dhh_x) + (dhh_y * dhh_y) + (dhh_z * dhh_z));
    inv_mag_dho = value_one / sqrtf((dho_x * dho_x) + (dho_y * dho_y) + (dho_z * dho_z));
  }
  doh_x *= inv_mag_doh;
  doh_y *= inv_mag_doh;
  doh_z *= inv_mag_doh;
  dhh_x *= inv_mag_dhh;
  dhh_y *= inv_mag_dhh;
  dhh_z *= inv_mag_dhh;
  dho_x *= inv_mag_dho;
  dho_y *= inv_mag_dho;
  dho_z *= inv_mag_dho;
  const Tcalc cos_a = -((doh_x * dho_x) + (doh_y * dho_y) + (doh_z * dho_z));
  const Tcalc cos_b = -((dhh_x * doh_x) + (dhh_y * doh_y) + (dhh_z * doh_z));
  const Tcalc cos_c = -((dho_x * dhh_x) + (dho_y * dhh_y) + (dho_z * dhh_z));
  
  // Get the velocity projections onto each displacement
  Tcalc dvoh_x, dvoh_y, dvoh_z, dvhh_x, dvhh_y, dvhh_z, dvho_x, dvho_y, dvho_z;
  if (tcoord_is_integral) {
    if (xvel_ovrf != nullptr && yvel_ovrf != nullptr && zvel_ovrf != nullptr) {
      const int95_t idvoh_x = hostInt95Subtract(xvel[hd1_idx], xvel_ovrf[hd1_idx], xvel[oxy_idx],
                                                xvel_ovrf[oxy_idx]);
      const int95_t idvoh_y = hostInt95Subtract(yvel[hd1_idx], yvel_ovrf[hd1_idx], yvel[oxy_idx],
                                                yvel_ovrf[oxy_idx]);
      const int95_t idvoh_z = hostInt95Subtract(zvel[hd1_idx], zvel_ovrf[hd1_idx], zvel[oxy_idx],
                                                zvel_ovrf[oxy_idx]);
      dvoh_x = hostInt95ToDouble(idvoh_x) / vel_scale;
      dvoh_y = hostInt95ToDouble(idvoh_y) / vel_scale;
      dvoh_z = hostInt95ToDouble(idvoh_z) / vel_scale;
      const int95_t idvhh_x = hostInt95Subtract(xvel[hd2_idx], xvel_ovrf[hd2_idx], xvel[hd1_idx],
                                                xvel_ovrf[hd1_idx]);
      const int95_t idvhh_y = hostInt95Subtract(yvel[hd2_idx], yvel_ovrf[hd2_idx], yvel[hd1_idx],
                                                yvel_ovrf[hd1_idx]);
      const int95_t idvhh_z = hostInt95Subtract(zvel[hd2_idx], zvel_ovrf[hd2_idx], zvel[hd1_idx],
                                                zvel_ovrf[hd1_idx]);
      dvhh_x = hostInt95ToDouble(idvhh_x) / vel_scale;
      dvhh_y = hostInt95ToDouble(idvhh_y) / vel_scale;
      dvhh_z = hostInt95ToDouble(idvhh_z) / vel_scale;
      const int95_t idvho_x = hostInt95Subtract(xvel[oxy_idx], xvel_ovrf[oxy_idx], xvel[hd2_idx],
                                                xvel_ovrf[hd2_idx]);
      const int95_t idvho_y = hostInt95Subtract(yvel[oxy_idx], yvel_ovrf[oxy_idx], yvel[hd2_idx],
                                                yvel_ovrf[hd2_idx]);
      const int95_t idvho_z = hostInt95Subtract(zvel[oxy_idx], zvel_ovrf[oxy_idx], zvel[hd2_idx],
                                                zvel_ovrf[hd2_idx]);
      dvho_x = hostInt95ToDouble(idvho_x) / vel_scale;
      dvho_y = hostInt95ToDouble(idvho_y) / vel_scale;
      dvho_z = hostInt95ToDouble(idvho_z) / vel_scale;
    }
    else {
      dvoh_x = static_cast<Tcalc>(xvel[hd1_idx] - xvel[oxy_idx]) / vel_scale;
      dvoh_y = static_cast<Tcalc>(yvel[hd1_idx] - yvel[oxy_idx]) / vel_scale;
      dvoh_z = static_cast<Tcalc>(zvel[hd1_idx] - zvel[oxy_idx]) / vel_scale;
      dvhh_x = static_cast<Tcalc>(xvel[hd2_idx] - xvel[hd1_idx]) / vel_scale;
      dvhh_y = static_cast<Tcalc>(yvel[hd2_idx] - yvel[hd1_idx]) / vel_scale;
      dvhh_z = static_cast<Tcalc>(zvel[hd2_idx] - zvel[hd1_idx]) / vel_scale;
      dvho_x = static_cast<Tcalc>(xvel[oxy_idx] - xvel[hd2_idx]) / vel_scale;
      dvho_y = static_cast<Tcalc>(yvel[oxy_idx] - yvel[hd2_idx]) / vel_scale;
      dvho_z = static_cast<Tcalc>(zvel[oxy_idx] - zvel[hd2_idx]) / vel_scale;
    }
  }
  else {
    dvoh_x = xvel[hd1_idx] - xvel[oxy_idx];
    dvoh_y = yvel[hd1_idx] - yvel[oxy_idx];
    dvoh_z = zvel[hd1_idx] - zvel[oxy_idx];
    dvhh_x = xvel[hd2_idx] - xvel[hd1_idx];
    dvhh_y = yvel[hd2_idx] - yvel[hd1_idx];
    dvhh_z = zvel[hd2_idx] - zvel[hd1_idx];
    dvho_x = xvel[oxy_idx] - xvel[hd2_idx];
    dvho_y = yvel[oxy_idx] - yvel[hd2_idx];
    dvho_z = zvel[oxy_idx] - zvel[hd2_idx];
  }
  const Tcalc voh = (doh_x * dvoh_x) + (doh_y * dvoh_y) + (doh_z * dvoh_z);
  const Tcalc vhh = (dhh_x * dvhh_x) + (dhh_y * dvhh_y) + (dhh_z * dvhh_z);
  const Tcalc vho = (dho_x * dvho_x) + (dho_y * dvho_y) + (dho_z * dvho_z);

  // Act on the results to adjust the velocities
  const Tcalc d = m_h / ((m_oh * m_oh) +
                         (m_h * cos_a * ((m_o * cos_b * cos_c) - (m_h * cos_a))) -
                         (value_half * m_o * m_oh * ((cos_b * cos_b) + (cos_c * cos_c))));
  const Tcalc t_oh = ((voh * ((value_two * m_oh) - (m_o * cos_c * cos_c))) +
                      (vhh * ((m_h * cos_c * cos_a) - (m_oh * cos_b))) +
                      (vho * ((m_o * cos_b * cos_c) - (value_two * m_h * cos_a)))) * m_o * d;
  const Tcalc t_hh = ((vhh * ((m_oh * m_oh) - (m_h * m_h * cos_a * cos_a))) +
                      (vho * m_o * ((m_h * cos_a * cos_b) - (m_oh * cos_c))) +
                      (voh * m_o * ((m_h * cos_c * cos_a) - (m_oh * cos_b)))) * d;
  const Tcalc t_ho = ((vho * ((value_two * m_oh) - (m_o * cos_b * cos_b))) +
                      (voh * ((m_o * cos_b * cos_c) - (value_two * m_h * cos_a))) +
                      (vhh * ((m_h * cos_a * cos_b) - (m_oh * cos_c)))) * m_o * d;
  
  // Return the results to the original arrays
  const Tcalc vox_bump = value_half * ((t_oh * doh_x) - (t_ho * dho_x)) / m_o;
  const Tcalc voy_bump = value_half * ((t_oh * doh_y) - (t_ho * dho_y)) / m_o;
  const Tcalc voz_bump = value_half * ((t_oh * doh_z) - (t_ho * dho_z)) / m_o;
  const Tcalc vh1x_bump = value_half * ((t_hh * dhh_x) - (t_oh * doh_x)) / m_h;
  const Tcalc vh1y_bump = value_half * ((t_hh * dhh_y) - (t_oh * doh_y)) / m_h;
  const Tcalc vh1z_bump = value_half * ((t_hh * dhh_z) - (t_oh * doh_z)) / m_h;
  const Tcalc vh2x_bump = value_half * ((t_ho * dho_x) - (t_hh * dhh_x)) / m_h;
  const Tcalc vh2y_bump = value_half * ((t_ho * dho_y) - (t_hh * dhh_y)) / m_h;
  const Tcalc vh2z_bump = value_half * ((t_ho * dho_z) - (t_hh * dhh_z)) / m_h;
  if (tcoord_is_integral) {
    if (xvel_ovrf != nullptr && yvel_ovrf != nullptr && zvel_ovrf != nullptr) {
      const int95_t ivox = hostInt95Sum(xvel[oxy_idx], xvel_ovrf[oxy_idx],
                                        vox_bump * vel_scale);
      const int95_t ivoy = hostInt95Sum(yvel[oxy_idx], yvel_ovrf[oxy_idx],
                                        voy_bump * vel_scale);
      const int95_t ivoz = hostInt95Sum(zvel[oxy_idx], zvel_ovrf[oxy_idx],
                                        voz_bump * vel_scale);
      const int95_t ivh1x = hostInt95Sum(xvel[hd1_idx], xvel_ovrf[hd1_idx],
                                         vh1x_bump * vel_scale);
      const int95_t ivh1y = hostInt95Sum(yvel[hd1_idx], yvel_ovrf[hd1_idx],
                                         vh1y_bump * vel_scale);
      const int95_t ivh1z = hostInt95Sum(zvel[hd1_idx], zvel_ovrf[hd1_idx],
                                         vh1z_bump * vel_scale);
      const int95_t ivh2x = hostInt95Sum(xvel[hd2_idx], xvel_ovrf[hd2_idx],
                                         vh2x_bump * vel_scale);
      const int95_t ivh2y = hostInt95Sum(yvel[hd2_idx], yvel_ovrf[hd2_idx],
                                         vh2y_bump * vel_scale);
      const int95_t ivh2z = hostInt95Sum(zvel[hd2_idx], zvel_ovrf[hd2_idx],
                                         vh2z_bump * vel_scale);
      xvel[oxy_idx] = ivox.x;
      yvel[oxy_idx] = ivoy.x;
      zvel[oxy_idx] = ivoz.x;
      xvel_ovrf[oxy_idx] = ivox.y;
      yvel_ovrf[oxy_idx] = ivoy.y;
      zvel_ovrf[oxy_idx] = ivoz.y;
      xvel[hd1_idx] = ivh1x.x;
      yvel[hd1_idx] = ivh1y.x;
      zvel[hd1_idx] = ivh1z.x;
      xvel_ovrf[hd1_idx] = ivh1x.y;
      yvel_ovrf[hd1_idx] = ivh1y.y;
      zvel_ovrf[hd1_idx] = ivh1z.y;
      xvel[hd2_idx] = ivh2x.x;
      yvel[hd2_idx] = ivh2y.x;
      zvel[hd2_idx] = ivh2z.x;
      xvel_ovrf[hd2_idx] = ivh2x.y;
      yvel_ovrf[hd2_idx] = ivh2y.y;
      zvel_ovrf[hd2_idx] = ivh2z.y;
    }
    else {
      xvel[oxy_idx] += static_cast<llint>(llround(vox_bump * vel_scale));
      yvel[oxy_idx] += static_cast<llint>(llround(voy_bump * vel_scale));
      zvel[oxy_idx] += static_cast<llint>(llround(voz_bump * vel_scale));
      xvel[hd1_idx] += static_cast<llint>(llround(vh1x_bump * vel_scale));
      yvel[hd1_idx] += static_cast<llint>(llround(vh1y_bump * vel_scale));
      zvel[hd1_idx] += static_cast<llint>(llround(vh1z_bump * vel_scale));
      xvel[hd2_idx] += static_cast<llint>(llround(vh2x_bump * vel_scale));
      yvel[hd2_idx] += static_cast<llint>(llround(vh2y_bump * vel_scale));
      zvel[hd2_idx] += static_cast<llint>(llround(vh2z_bump * vel_scale));
    }
  }
  else {
    xvel[oxy_idx] += vox_bump;
    yvel[oxy_idx] += voy_bump;
    zvel[oxy_idx] += voz_bump;
    xvel[hd1_idx] += vh1x_bump;
    yvel[hd1_idx] += vh1y_bump;
    zvel[hd1_idx] += vh1z_bump;
    xvel[hd2_idx] += vh2x_bump;
    yvel[hd2_idx] += vh2y_bump;
    zvel[hd2_idx] += vh2z_bump;
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void settleVelocities(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, Tcoord* xvel,
                      Tcoord* yvel, Tcoord* zvel, const ConstraintKit<Tcalc> &cnst,
                      const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                      int* xvel_ovrf, int* yvel_ovrf, int* zvel_ovrf, const Tcalc gpos_scale,
                      const Tcalc vel_scale) {

  // The SETTLE constraint groups indicate relevant parameter sets, and there is the provision for
  // more than one SETTLE constraint geometry.
  for (int i = 0; i < cnst.nsettle; i++) {
    const int parm_idx = cnst.settle_param_idx[i];
    settleGroupVelocity<Tcoord, Tcalc>(cnst.settle_ox_atoms[i], cnst.settle_h1_atoms[i],
                                       cnst.settle_h2_atoms[i], cnst.settle_mo[parm_idx],
                                       cnst.settle_mh[parm_idx], cnst.settle_moh[parm_idx], xcrd,
                                       ycrd, zcrd, xvel, yvel, zvel, xcrd_ovrf, ycrd_ovrf,
                                       zcrd_ovrf, xvel_ovrf, yvel_ovrf, zvel_ovrf, gpos_scale,
                                       vel_scale);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void settleVelocities(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnst) {
  settleVelocities<double, Tcalc>(psw->xcrd, psw->ycrd, psw->zcrd, psw->vxalt, psw->vyalt,
                                  psw->vzalt, cnst);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void settleVelocities(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                      const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk) {
  std::vector<llint> local_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_zcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_xvel(maximum_valence_work_unit_atoms);
  std::vector<llint> local_yvel(maximum_valence_work_unit_atoms);
  std::vector<llint> local_zvel(maximum_valence_work_unit_atoms);
  std::vector<int> local_xcrd_ovrf, local_ycrd_ovrf, local_zcrd_ovrf;
  if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
    local_xcrd_ovrf.resize(maximum_valence_work_unit_atoms);
    local_ycrd_ovrf.resize(maximum_valence_work_unit_atoms);
    local_zcrd_ovrf.resize(maximum_valence_work_unit_atoms);
  }
  std::vector<int> local_xvel_ovrf, local_yvel_ovrf, local_zvel_ovrf;
  if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
    local_xvel_ovrf.resize(maximum_valence_work_unit_atoms);
    local_yvel_ovrf.resize(maximum_valence_work_unit_atoms);
    local_zvel_ovrf.resize(maximum_valence_work_unit_atoms);
  }
  std::vector<int2> vwu_notes(vwu_abstract_length);
  for (int i = 0; i < poly_vk.nvwu; i++) {

    // Obtain the work unit's abstract
    const int abst_llim = vwu_abstract_length * i;
    const int abst_hlim = abst_llim + vwu_abstract_length;
    for (int j = abst_llim; j < abst_hlim; j++) {
      vwu_notes[j - abst_llim] = poly_vk.vwu_abstracts[j];
    }
    
    // Import atoms and velocities
    const int2 import_lims = vwu_notes[static_cast<size_t>(VwuAbstractMap::IMPORT)];
    for (int j = import_lims.x; j < import_lims.y; j++) {
      const size_t jatom = poly_vk.vwu_imports[j];
      const size_t jlocal_idx = j - import_lims.x;
      local_xcrd[jlocal_idx] = poly_psw->xcrd[jatom];
      local_ycrd[jlocal_idx] = poly_psw->ycrd[jatom];
      local_zcrd[jlocal_idx] = poly_psw->zcrd[jatom];
      if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        local_xcrd_ovrf[jlocal_idx] = poly_psw->xcrd_ovrf[jatom];
        local_ycrd_ovrf[jlocal_idx] = poly_psw->ycrd_ovrf[jatom];
        local_zcrd_ovrf[jlocal_idx] = poly_psw->zcrd_ovrf[jatom];
      }
      local_xvel[jlocal_idx] = poly_psw->vxalt[jatom];
      local_yvel[jlocal_idx] = poly_psw->vyalt[jatom];
      local_zvel[jlocal_idx] = poly_psw->vzalt[jatom];
      if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
        local_xvel_ovrf[jlocal_idx] = poly_psw->vxalt_ovrf[jatom];
        local_yvel_ovrf[jlocal_idx] = poly_psw->vyalt_ovrf[jatom];
        local_zvel_ovrf[jlocal_idx] = poly_psw->vzalt_ovrf[jatom];
      }
    }

    // Loop over all instructions
    const int2 insr_lims = vwu_notes[static_cast<size_t>(VwuAbstractMap::SETTLE)];
    int *xcrd_ovrf_ptr, *ycrd_ovrf_ptr, *zcrd_ovrf_ptr;
    if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
      xcrd_ovrf_ptr = nullptr;
      ycrd_ovrf_ptr = nullptr;
      zcrd_ovrf_ptr = nullptr;
    }
    else {
      xcrd_ovrf_ptr = local_xcrd_ovrf.data();
      ycrd_ovrf_ptr = local_ycrd_ovrf.data();
      zcrd_ovrf_ptr = local_zcrd_ovrf.data();
    }
    int *xvel_ovrf_ptr, *yvel_ovrf_ptr, *zvel_ovrf_ptr;
    if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
      xvel_ovrf_ptr = nullptr;
      yvel_ovrf_ptr = nullptr;
      zvel_ovrf_ptr = nullptr;
    }
    else {
      xvel_ovrf_ptr = local_xvel_ovrf.data();
      yvel_ovrf_ptr = local_yvel_ovrf.data();
      zvel_ovrf_ptr = local_zvel_ovrf.data();
    }
    for (int j = insr_lims.x; j < insr_lims.y; j++) {
      const uint2 insr = poly_auk.sett_insr[j];
      const int oxy_idx = (insr.x & 0x3ff);
      const int hd1_idx = ((insr.x >> 10) & 0x3ff);
      const int hd2_idx = ((insr.x >> 20) & 0x3ff);
      const Tcalc4 sgeom = poly_auk.settle_geom[insr.y];
      const Tcalc4 smass = poly_auk.settle_mass[insr.y];
      settleGroupVelocity<llint, Tcalc>(oxy_idx, hd1_idx, hd2_idx, smass.x, smass.y, sgeom.w,
                                        local_xcrd.data(), local_ycrd.data(), local_zcrd.data(),
                                        local_xvel.data(), local_yvel.data(), local_zvel.data(),
                                        xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr, xvel_ovrf_ptr,
                                        yvel_ovrf_ptr, zvel_ovrf_ptr, poly_psw->gpos_scale,
                                        poly_psw->vel_scale);
    }
    
    // Return atom velocity information which the work unit is responsible for updating
    const int manip_llim = vwu_notes[(size_t)(VwuAbstractMap::MANIPULATE)].x;
    for (int j = import_lims.x; j < import_lims.y; j++) {
      const size_t jatom = poly_vk.vwu_imports[j];
      const size_t jlocal_idx = j - import_lims.x;

      // Determine the atom ownership.  Only the update mask matters, as opposed to the movement
      // mask, as no other work will be done while these work units are temporarily laid out.
      const int mask_segment = manip_llim + (jlocal_idx >> 5);
      const int segment_bit  = jlocal_idx - (mask_segment << 5);
      if (poly_auk.vwu_manip[mask_segment].y & (0x1 << segment_bit)) {
        poly_psw->vxalt[jatom] = local_xvel[jlocal_idx];
        poly_psw->vyalt[jatom] = local_yvel[jlocal_idx];
        poly_psw->vzalt[jatom] = local_zvel[jlocal_idx];
        if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
          poly_psw->vxalt_ovrf[jatom] = local_xvel_ovrf[jlocal_idx];
          poly_psw->vyalt_ovrf[jatom] = local_yvel_ovrf[jlocal_idx];
          poly_psw->vzalt_ovrf[jatom] = local_zvel_ovrf[jlocal_idx];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void updateSettlePVC(const int atm_idx, const Tcalc p_to_v_factor,
                     const Tcalc r_adj, const Tcoord r_ref, const int r_ref_ovrf, Tcoord* r_dev, 
                     int* r_dev_ovrf, Tcoord* v_dev, int* v_dev_ovrf) {
  if (r_dev_ovrf == nullptr) {
    const Tcoord ir_next = r_ref + static_cast<Tcoord>(llround(r_adj));
    const Tcalc d_chng = ir_next - r_dev[atm_idx];
    if (v_dev_ovrf == nullptr) {
      const Tcoord dv = llround(d_chng * p_to_v_factor);
      v_dev[atm_idx] += dv;
    }
    else {
      const int95_t dv = hostDoubleToInt95(d_chng * p_to_v_factor);
      const int95_t iv_next = hostSplitFPSum(dv, v_dev[atm_idx], v_dev_ovrf[atm_idx]);
      v_dev[atm_idx] = iv_next.x;
      v_dev_ovrf[atm_idx] = iv_next.y;
    }
    r_dev[atm_idx] = ir_next;
  }
  else {
    const int95_t ir_next = hostInt95Sum(r_ref, r_ref_ovrf, r_adj);
    const int95_t i_chng = hostSplitFPSubtract(ir_next, r_dev[atm_idx], r_dev_ovrf[atm_idx]);
    const Tcalc d_chng = hostInt95ToDouble(i_chng);
    if (v_dev_ovrf == nullptr) {
      const Tcoord dv = llround(d_chng * p_to_v_factor);
      v_dev[atm_idx] += dv;
    }
    else {
      const int95_t dv = hostDoubleToInt95(d_chng * p_to_v_factor);
      const int95_t iv_next = hostSplitFPSum(dv, v_dev[atm_idx], v_dev_ovrf[atm_idx]);
      v_dev[atm_idx] = iv_next.x;
      v_dev_ovrf[atm_idx] = iv_next.y;
    }
    r_dev[atm_idx] = ir_next.x;
    r_dev_ovrf[atm_idx] = ir_next.y;
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc, typename Tcalc3>
void settleGroupPosition(const int oxy_idx, const int hd1_idx, const int hd2_idx,
                         const Tcalc frac_heavy, const Tcalc frac_light, const Tcalc ra,
                         const Tcalc rb, const Tcalc rc, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, Tcoord* xalt, Tcoord* yalt, Tcoord* zalt,
                         Tcoord* vxalt, Tcoord* vyalt, Tcoord* vzalt, const Tcalc dt,
                         const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                         int* xalt_ovrf, int* yalt_ovrf, int* zalt_ovrf, int* vxalt_ovrf,
                         int* vyalt_ovrf, int* vzalt_ovrf, const Tcalc gpos_scale,
                         const Tcalc vel_scale) {
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc value_zero = 0.0;
  const Tcalc value_one  = 1.0;
  const Tcalc value_two  = 2.0;
  
  // The above terms that follow unite under the following diagram, with the system's center of
  // mass placed at the origin for numerical convenience.  On the GPU, the best policy will be to
  // consider placing the oxygen at the origin using the fixed-precision coordinates, then moving
  // the hydrogens along with it.  The implementation is inspired by that found in NAMD, from the
  // Beckman institute at the University of Illinois at Urbana-Champaign.
  //
  //                 |               
  //                 |                       "Primed" coordinate axes:
  //              a0 x -----
  //                 |   |                       ^
  //                 |   ra                      |
  //                 |   |                       j'
  //    -----------Orig.-------------            |
  //         |       |        |                  |
  //         rb      |---rc---|                  +---- i' ---->
  //         |       |        |                k'
  //  b0 x-----------|        c0 x
  //
  
  // Take the reference coordinates of the oxygen atom as an impromptu center of coordinates, to
  // get the numbers into small scales without excessive register usage in a GPU implementation.
  Tcalc3 oxy_local_crd, hd1_local_crd, hd2_local_crd;
  if (tcoord_is_integral) {
    if (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr && zcrd_ovrf != nullptr &&
        xalt_ovrf != nullptr && yalt_ovrf != nullptr && zalt_ovrf != nullptr) {
      const int95_t ioxy_dx = hostInt95Subtract(xalt[oxy_idx], xalt_ovrf[oxy_idx],
                                                xcrd[oxy_idx], xcrd_ovrf[oxy_idx]);
      const int95_t ioxy_dy = hostInt95Subtract(yalt[oxy_idx], yalt_ovrf[oxy_idx],
                                                ycrd[oxy_idx], ycrd_ovrf[oxy_idx]);
      const int95_t ioxy_dz = hostInt95Subtract(zalt[oxy_idx], zalt_ovrf[oxy_idx],
                                                zcrd[oxy_idx], zcrd_ovrf[oxy_idx]);
      oxy_local_crd.x = hostInt95ToDouble(ioxy_dx) / gpos_scale;
      oxy_local_crd.y = hostInt95ToDouble(ioxy_dy) / gpos_scale;
      oxy_local_crd.z = hostInt95ToDouble(ioxy_dz) / gpos_scale;
      const int95_t ihd1_dx = hostInt95Subtract(xalt[hd1_idx], xalt_ovrf[hd1_idx],
                                                xcrd[oxy_idx], xcrd_ovrf[oxy_idx]);
      const int95_t ihd1_dy = hostInt95Subtract(yalt[hd1_idx], yalt_ovrf[hd1_idx],
                                                ycrd[oxy_idx], ycrd_ovrf[oxy_idx]);
      const int95_t ihd1_dz = hostInt95Subtract(zalt[hd1_idx], zalt_ovrf[hd1_idx],
                                                zcrd[oxy_idx], zcrd_ovrf[oxy_idx]);
      hd1_local_crd.x = hostInt95ToDouble(ihd1_dx) / gpos_scale;
      hd1_local_crd.y = hostInt95ToDouble(ihd1_dy) / gpos_scale;
      hd1_local_crd.z = hostInt95ToDouble(ihd1_dz) / gpos_scale;
      const int95_t ihd2_dx = hostInt95Subtract(xalt[hd2_idx], xalt_ovrf[hd2_idx],
                                                xcrd[oxy_idx], xcrd_ovrf[oxy_idx]);
      const int95_t ihd2_dy = hostInt95Subtract(yalt[hd2_idx], yalt_ovrf[hd2_idx],
                                                ycrd[oxy_idx], ycrd_ovrf[oxy_idx]);
      const int95_t ihd2_dz = hostInt95Subtract(zalt[hd2_idx], zalt_ovrf[hd2_idx],
                                                zcrd[oxy_idx], zcrd_ovrf[oxy_idx]);
      hd2_local_crd.x = hostInt95ToDouble(ihd2_dx) / gpos_scale;
      hd2_local_crd.y = hostInt95ToDouble(ihd2_dy) / gpos_scale;
      hd2_local_crd.z = hostInt95ToDouble(ihd2_dz) / gpos_scale;
    }
    else {
      oxy_local_crd.x = static_cast<Tcalc>(xalt[oxy_idx] - xcrd[oxy_idx]) / gpos_scale;
      oxy_local_crd.y = static_cast<Tcalc>(yalt[oxy_idx] - ycrd[oxy_idx]) / gpos_scale;
      oxy_local_crd.z = static_cast<Tcalc>(zalt[oxy_idx] - zcrd[oxy_idx]) / gpos_scale;
      hd1_local_crd.x = static_cast<Tcalc>(xalt[hd1_idx] - xcrd[oxy_idx]) / gpos_scale;
      hd1_local_crd.y = static_cast<Tcalc>(yalt[hd1_idx] - ycrd[oxy_idx]) / gpos_scale;
      hd1_local_crd.z = static_cast<Tcalc>(zalt[hd1_idx] - zcrd[oxy_idx]) / gpos_scale;
      hd2_local_crd.x = static_cast<Tcalc>(xalt[hd2_idx] - xcrd[oxy_idx]) / gpos_scale;
      hd2_local_crd.y = static_cast<Tcalc>(yalt[hd2_idx] - ycrd[oxy_idx]) / gpos_scale;
      hd2_local_crd.z = static_cast<Tcalc>(zalt[hd2_idx] - zcrd[oxy_idx]) / gpos_scale;
    }
  }
  else {
    oxy_local_crd.x = xalt[oxy_idx] - xcrd[oxy_idx];
    oxy_local_crd.y = yalt[oxy_idx] - ycrd[oxy_idx];
    oxy_local_crd.z = zalt[oxy_idx] - zcrd[oxy_idx];
    hd1_local_crd.x = xalt[hd1_idx] - xcrd[oxy_idx];
    hd1_local_crd.y = yalt[hd1_idx] - ycrd[oxy_idx];
    hd1_local_crd.z = zalt[hd1_idx] - zcrd[oxy_idx];
    hd2_local_crd.x = xalt[hd2_idx] - xcrd[oxy_idx];
    hd2_local_crd.y = yalt[hd2_idx] - ycrd[oxy_idx];
    hd2_local_crd.z = zalt[hd2_idx] - zcrd[oxy_idx];
  }
  const Tcalc3 curr_com = { (frac_heavy * oxy_local_crd.x) +
                            (frac_light * (hd1_local_crd.x + hd2_local_crd.x)),
                            (frac_heavy * oxy_local_crd.y) + 
                            (frac_light * (hd1_local_crd.y + hd2_local_crd.y)),
                            (frac_heavy * oxy_local_crd.z) +
                            (frac_light * (hd1_local_crd.z + hd2_local_crd.z)) };
  
  // Shift the local coordinates to place their center of mass at the origin.
  oxy_local_crd.x -= curr_com.x;
  oxy_local_crd.y -= curr_com.y;
  oxy_local_crd.z -= curr_com.z;
  hd1_local_crd.x -= curr_com.x;
  hd1_local_crd.y -= curr_com.y;
  hd1_local_crd.z -= curr_com.z;
  hd2_local_crd.x -= curr_com.x;
  hd2_local_crd.y -= curr_com.y;
  hd2_local_crd.z -= curr_com.z;
  
  // Compute the difference vectors for the reference coordinates
  Tcalc3 ref_oh1, ref_oh2;
  if (tcoord_is_integral) {
    if (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr && zcrd_ovrf != nullptr) {
      const int95_t ioh1_dx = hostInt95Subtract(xcrd[hd1_idx], xcrd_ovrf[hd1_idx],
                                                xcrd[oxy_idx], xcrd_ovrf[oxy_idx]);
      const int95_t ioh1_dy = hostInt95Subtract(ycrd[hd1_idx], ycrd_ovrf[hd1_idx],
                                                ycrd[oxy_idx], ycrd_ovrf[oxy_idx]);
      const int95_t ioh1_dz = hostInt95Subtract(zcrd[hd1_idx], zcrd_ovrf[hd1_idx],
                                                zcrd[oxy_idx], zcrd_ovrf[oxy_idx]);
      ref_oh1.x = hostInt95ToDouble(ioh1_dx) / gpos_scale;
      ref_oh1.y = hostInt95ToDouble(ioh1_dy) / gpos_scale;
      ref_oh1.z = hostInt95ToDouble(ioh1_dz) / gpos_scale;
      const int95_t ioh2_dx = hostInt95Subtract(xcrd[hd2_idx], xcrd_ovrf[hd2_idx],
                                                xcrd[oxy_idx], xcrd_ovrf[oxy_idx]);
      const int95_t ioh2_dy = hostInt95Subtract(ycrd[hd2_idx], ycrd_ovrf[hd2_idx],
                                                ycrd[oxy_idx], ycrd_ovrf[oxy_idx]);
      const int95_t ioh2_dz = hostInt95Subtract(zcrd[hd2_idx], zcrd_ovrf[hd2_idx],
                                                zcrd[oxy_idx], zcrd_ovrf[oxy_idx]);
      ref_oh2.x = hostInt95ToDouble(ioh2_dx) / gpos_scale;
      ref_oh2.y = hostInt95ToDouble(ioh2_dy) / gpos_scale;
      ref_oh2.z = hostInt95ToDouble(ioh2_dz) / gpos_scale;
    }
    else {
      ref_oh1.x = static_cast<Tcalc>(xcrd[hd1_idx] - xcrd[oxy_idx]) / gpos_scale;
      ref_oh1.y = static_cast<Tcalc>(ycrd[hd1_idx] - ycrd[oxy_idx]) / gpos_scale;
      ref_oh1.z = static_cast<Tcalc>(zcrd[hd1_idx] - zcrd[oxy_idx]) / gpos_scale;
      ref_oh2.x = static_cast<Tcalc>(xcrd[hd2_idx] - xcrd[oxy_idx]) / gpos_scale;
      ref_oh2.y = static_cast<Tcalc>(ycrd[hd2_idx] - ycrd[oxy_idx]) / gpos_scale;
      ref_oh2.z = static_cast<Tcalc>(zcrd[hd2_idx] - zcrd[oxy_idx]) / gpos_scale;
    }
  }
  else {
    ref_oh1.x = xcrd[hd1_idx] - xcrd[oxy_idx];
    ref_oh1.y = ycrd[hd1_idx] - ycrd[oxy_idx];
    ref_oh1.z = zcrd[hd1_idx] - zcrd[oxy_idx];
    ref_oh2.x = xcrd[hd2_idx] - xcrd[oxy_idx];
    ref_oh2.y = ycrd[hd2_idx] - ycrd[oxy_idx];
    ref_oh2.z = zcrd[hd2_idx] - zcrd[oxy_idx];
  }
  
  // Compute a series of cross products taking the normal of the plane of the heavy (oxygen)
  // atom's displacement to the normal of the plane of the molecule.  The resulting vector will
  // lie in the plane of the molecule running parallel to the h1 -> h2 vector.
  Tcalc3 prime_k = crossProduct(ref_oh1, ref_oh2);
  Tcalc3 prime_i = crossProduct(oxy_local_crd, prime_k);
  Tcalc3 prime_j = crossProduct(prime_k, prime_i);
  normalize<Tcalc3, Tcalc>(&prime_i);
  normalize<Tcalc3, Tcalc>(&prime_j);
  normalize<Tcalc3, Tcalc>(&prime_k);
  
  // Transform various vectors into the primed coordinate system
  rotateCoordinates<Tcalc3, Tcalc>(&ref_oh1, prime_i, prime_j, prime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&ref_oh2, prime_i, prime_j, prime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&oxy_local_crd, prime_i, prime_j, prime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&hd1_local_crd, prime_i, prime_j, prime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&hd2_local_crd, prime_i, prime_j, prime_k);

  // Compute the positions of the canonical SETTLE group, oriented in the plane of the current
  // particle set.
  const Tcalc sinphi = oxy_local_crd.z / ra;
  const Tcalc cosphi = (tcalc_is_double) ? sqrt(value_one - (sinphi * sinphi)) :
                                           sqrtf(value_one - (sinphi * sinphi));
  const Tcalc sinpsi = (hd1_local_crd.z - hd2_local_crd.z) / (value_two * rc * cosphi);
  const Tcalc cospsi = (tcalc_is_double) ? sqrt(value_one - (sinpsi * sinpsi)) :
                                           sqrtf(value_one - (sinpsi * sinpsi));
  const Tcalc rc_ss = rc * sinpsi * sinphi;
  const Tcalc rc_sc = rc * sinpsi * cosphi;
  const Tcalc3 canonical_oxy = { value_zero, ra * cosphi, ra * sinphi };
  const Tcalc3 canonical_hd1 = { -rc * cospsi, -(rc_ss + (rb * cosphi)), rc_sc - (rb * sinphi) };
  const Tcalc3 canonical_hd2 = { rc * cosphi, rc_ss - (rb * cosphi), -(rc_sc + (rb * sinphi)) };
  
  // Taking the displacements of light atoms relative to the heavy atom in the original geometry
  // facilitates derivation of the rotation matrix.
  const Tcalc alpha = (canonical_hd1.x * (ref_oh1.x - ref_oh2.x)) +
                      (canonical_hd1.y * ref_oh1.y) + (canonical_hd2.y * ref_oh2.y);
  const Tcalc beta  = (canonical_hd1.x * (ref_oh2.y - ref_oh1.y)) +
                      (canonical_hd1.y * ref_oh1.x) + (canonical_hd2.y * ref_oh2.x);
  const Tcalc gamma = (ref_oh1.x * hd1_local_crd.y) - (hd1_local_crd.x * ref_oh1.y) +
                      (ref_oh2.x * hd2_local_crd.y) - (hd2_local_crd.x * ref_oh2.y);
  const Tcalc a2b2 = (alpha * alpha) + (beta * beta);
  Tcalc sintheta, costheta;
  if (tcalc_is_double) {
    sintheta = ((alpha * gamma) - (beta * sqrt(a2b2 - (gamma * gamma)))) / a2b2;
    costheta = sqrt(value_one - (sintheta * sintheta));
  }
  else {
    sintheta = ((alpha * gamma) - (beta * sqrtf(a2b2 - (gamma * gamma)))) / a2b2;
    costheta = sqrtf(value_one - (sintheta * sintheta));
  }
  
  // Rotate the canonical water about the local coordinate system's Z axis and impart the
  // Cartesian Z elevations of the particles' raw coordinates as they have been translated to
  // the local frame.
  Tcalc3 oxy_adj_crd = { -canonical_oxy.y * sintheta, canonical_oxy.y * costheta,
                         oxy_local_crd.z };
  Tcalc3 hd1_adj_crd = { (canonical_hd1.x * costheta) - (canonical_hd1.y * sintheta),
                         (canonical_hd1.x * sintheta) + (canonical_hd1.y * costheta),
                         hd1_local_crd.z };
  Tcalc3 hd2_adj_crd = { -(canonical_hd1.x * costheta) - (canonical_hd2.y * sintheta),
                         -(canonical_hd1.x * sintheta) + (canonical_hd2.y * costheta),
                         hd2_local_crd.z };
  
  // Undo the coordinate transformation by generating new normal vectors with a three-dimensional
  // transpose.
  const Tcalc3 unprime_i = { prime_i.x, prime_j.x, prime_k.x };
  const Tcalc3 unprime_j = { prime_i.y, prime_j.y, prime_k.y };
  const Tcalc3 unprime_k = { prime_i.z, prime_j.z, prime_k.z };
  rotateCoordinates<Tcalc3, Tcalc>(&oxy_adj_crd, unprime_i, unprime_j, unprime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&hd1_adj_crd, unprime_i, unprime_j, unprime_k);
  rotateCoordinates<Tcalc3, Tcalc>(&hd2_adj_crd, unprime_i, unprime_j, unprime_k);

  // Replace the developing coordinates with the constrained system after returning its center
  // of mass to the original location and then adjusting the entire coordinate system based on
  // the position of the reference oxygen.
  oxy_adj_crd.x += curr_com.x;
  oxy_adj_crd.y += curr_com.y;
  oxy_adj_crd.z += curr_com.z;
  hd1_adj_crd.x += curr_com.x;
  hd1_adj_crd.y += curr_com.y;
  hd1_adj_crd.z += curr_com.z;
  hd2_adj_crd.x += curr_com.x;
  hd2_adj_crd.y += curr_com.y;
  hd2_adj_crd.z += curr_com.z;
  if (tcoord_is_integral) {

    // As with SHAKE on the positions, it is necessary to apply a post-hoc correction after
    // SETTLE to keep particle velocities consistent.  Adjust the velocities with the changes in
    // each particle's position.  This prefactor will reduce clutter in the fixed-precision code.
    const Tcalc p_to_v_factor = vel_scale / (dt * gpos_scale);
    const Tcoord ref_oxy_x = xcrd[oxy_idx];
    const Tcoord ref_oxy_y = ycrd[oxy_idx];
    const Tcoord ref_oxy_z = zcrd[oxy_idx];
    const Tcoord ref_oxy_ovrf_x = (xcrd_ovrf != nullptr) ? xcrd_ovrf[oxy_idx] : 0;
    const Tcoord ref_oxy_ovrf_y = (ycrd_ovrf != nullptr) ? ycrd_ovrf[oxy_idx] : 0;
    const Tcoord ref_oxy_ovrf_z = (zcrd_ovrf != nullptr) ? zcrd_ovrf[oxy_idx] : 0;
    updateSettlePVC(oxy_idx, p_to_v_factor, oxy_adj_crd.x * gpos_scale, ref_oxy_x, ref_oxy_ovrf_x,
                    xalt, xalt_ovrf, vxalt, vxalt_ovrf);
    updateSettlePVC(oxy_idx, p_to_v_factor, oxy_adj_crd.y * gpos_scale, ref_oxy_y, ref_oxy_ovrf_y,
                    yalt, yalt_ovrf, vyalt, vyalt_ovrf);
    updateSettlePVC(oxy_idx, p_to_v_factor, oxy_adj_crd.z * gpos_scale, ref_oxy_z, ref_oxy_ovrf_z,
                    zalt, zalt_ovrf, vzalt, vzalt_ovrf);
    updateSettlePVC(hd1_idx, p_to_v_factor, hd1_adj_crd.x * gpos_scale, ref_oxy_x, ref_oxy_ovrf_x,
                    xalt, xalt_ovrf, vxalt, vxalt_ovrf);
    updateSettlePVC(hd1_idx, p_to_v_factor, hd1_adj_crd.y * gpos_scale, ref_oxy_y, ref_oxy_ovrf_y,
                    yalt, yalt_ovrf, vyalt, vyalt_ovrf);
    updateSettlePVC(hd1_idx, p_to_v_factor, hd1_adj_crd.z * gpos_scale, ref_oxy_z, ref_oxy_ovrf_z,
                    zalt, zalt_ovrf, vzalt, vzalt_ovrf);
    updateSettlePVC(hd2_idx, p_to_v_factor, hd2_adj_crd.x * gpos_scale, ref_oxy_x, ref_oxy_ovrf_x,
                    xalt, xalt_ovrf, vxalt, vxalt_ovrf);
    updateSettlePVC(hd2_idx, p_to_v_factor, hd2_adj_crd.y * gpos_scale, ref_oxy_y, ref_oxy_ovrf_y,
                    yalt, yalt_ovrf, vyalt, vyalt_ovrf);
    updateSettlePVC(hd2_idx, p_to_v_factor, hd2_adj_crd.z * gpos_scale, ref_oxy_z, ref_oxy_ovrf_z,
                    zalt, zalt_ovrf, vzalt, vzalt_ovrf);
  }
  else {
    const Tcoord nx_oxy = xcrd[oxy_idx] + oxy_adj_crd.x;
    const Tcoord ny_oxy = ycrd[oxy_idx] + oxy_adj_crd.y;
    const Tcoord nz_oxy = zcrd[oxy_idx] + oxy_adj_crd.z;
    const Tcoord nx_hd1 = xcrd[oxy_idx] + hd1_adj_crd.x;
    const Tcoord ny_hd1 = ycrd[oxy_idx] + hd1_adj_crd.y;
    const Tcoord nz_hd1 = zcrd[oxy_idx] + hd1_adj_crd.z;
    const Tcoord nx_hd2 = xcrd[oxy_idx] + hd2_adj_crd.x;
    const Tcoord ny_hd2 = ycrd[oxy_idx] + hd2_adj_crd.y;
    const Tcoord nz_hd2 = zcrd[oxy_idx] + hd2_adj_crd.z;
    const Tcoord dvx_oxy = static_cast<Tcalc>(nx_oxy - xalt[oxy_idx]) / dt;
    const Tcoord dvy_oxy = static_cast<Tcalc>(ny_oxy - yalt[oxy_idx]) / dt;
    const Tcoord dvz_oxy = static_cast<Tcalc>(nz_oxy - zalt[oxy_idx]) / dt;
    const Tcoord dvx_hd1 = static_cast<Tcalc>(nx_hd1 - xalt[hd1_idx]) / dt;
    const Tcoord dvy_hd1 = static_cast<Tcalc>(ny_hd1 - yalt[hd1_idx]) / dt;
    const Tcoord dvz_hd1 = static_cast<Tcalc>(nz_hd1 - zalt[hd1_idx]) / dt;
    const Tcoord dvx_hd2 = static_cast<Tcalc>(nx_hd2 - xalt[hd2_idx]) / dt;
    const Tcoord dvy_hd2 = static_cast<Tcalc>(ny_hd2 - yalt[hd2_idx]) / dt;
    const Tcoord dvz_hd2 = static_cast<Tcalc>(nz_hd2 - zalt[hd2_idx]) / dt;
    vxalt[oxy_idx] += dvx_oxy;
    vyalt[oxy_idx] += dvy_oxy;
    vzalt[oxy_idx] += dvz_oxy;
    vxalt[hd1_idx] += dvx_hd1;
    vyalt[hd1_idx] += dvy_hd1;
    vzalt[hd1_idx] += dvz_hd1;
    vxalt[hd2_idx] += dvx_hd2;
    vyalt[hd2_idx] += dvy_hd2;
    vzalt[hd2_idx] += dvz_hd2;
    xalt[oxy_idx] = nx_oxy;
    yalt[oxy_idx] = ny_oxy;
    zalt[oxy_idx] = nz_oxy;
    xalt[hd1_idx] = nx_hd1;
    yalt[hd1_idx] = ny_hd1;
    zalt[hd1_idx] = nz_hd1;
    xalt[hd2_idx] = nx_hd2;
    yalt[hd2_idx] = ny_hd2;
    zalt[hd2_idx] = nz_hd2;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc, typename Tcalc3>
void settlePositions(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, Tcoord* xalt,
                     Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt, Tcoord* vzalt,
                     const ConstraintKit<Tcalc> &cnst, const Tcalc dt, const int* xcrd_ovrf,
                     const int* ycrd_ovrf, const int* zcrd_ovrf, int* xalt_ovrf, int* yalt_ovrf,
                     int* zalt_ovrf, int* vxalt_ovrf, int* vyalt_ovrf, int* vzalt_ovrf,
                     const Tcalc gpos_scale, const Tcalc vel_scale) {

  // The SETTLE constraint groups indicate relevant parameter sets, and there is the provision for
  // more than one SETTLE constraint geometry.
  for (int i = 0; i < cnst.nsettle; i++) {

    // Obtain the SETTLE parameter set.  On the GPU, the parameters (other than the atom indices)
    // will have to linger in L1, as taking them into registers would be unaffordable.
    const int parm_idx = cnst.settle_param_idx[i];
    const Tcalc frac_heavy = cnst.settle_mormt[parm_idx];
    const Tcalc frac_light = cnst.settle_mhrmt[parm_idx];
    settleGroupPosition<Tcoord, Tcalc, Tcalc3>(cnst.settle_ox_atoms[i], cnst.settle_h1_atoms[i],
                                               cnst.settle_h2_atoms[i], frac_heavy, frac_light,
                                               cnst.settle_ra[parm_idx], cnst.settle_rb[parm_idx],
                                               cnst.settle_rc[parm_idx], xcrd, ycrd, zcrd, xalt,
                                               yalt, zalt, vxalt, vyalt, vzalt, dt, xcrd_ovrf,
                                               ycrd_ovrf, zcrd_ovrf, xalt_ovrf, yalt_ovrf,
                                               zalt_ovrf, vxalt_ovrf, vyalt_ovrf, vzalt_ovrf,
                                               gpos_scale, vel_scale);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
void settlePositions(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnst, const Tcalc dt) {
  settlePositions<double, Tcalc, Tcalc3>(psw->xcrd, psw->ycrd, psw->zcrd, psw->xalt, psw->yalt,
                                         psw->zalt, psw->vxalt, psw->vyalt, psw->vzalt, cnst, dt);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc3, typename Tcalc4>
void settlePositions(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                     const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk, const Tcalc dt) {
  std::vector<llint> local_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_zcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> local_xalt(maximum_valence_work_unit_atoms);
  std::vector<llint> local_yalt(maximum_valence_work_unit_atoms);
  std::vector<llint> local_zalt(maximum_valence_work_unit_atoms);
  std::vector<llint> local_vxalt(maximum_valence_work_unit_atoms);
  std::vector<llint> local_vyalt(maximum_valence_work_unit_atoms);
  std::vector<llint> local_vzalt(maximum_valence_work_unit_atoms);
  std::vector<int> local_xcrd_ovrf, local_xalt_ovrf, local_vxalt_ovrf;
  std::vector<int> local_ycrd_ovrf, local_yalt_ovrf, local_vyalt_ovrf;
  std::vector<int> local_zcrd_ovrf, local_zalt_ovrf, local_vzalt_ovrf;
  if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
    local_xcrd_ovrf.resize(maximum_valence_work_unit_atoms);
    local_ycrd_ovrf.resize(maximum_valence_work_unit_atoms);
    local_zcrd_ovrf.resize(maximum_valence_work_unit_atoms);
    local_xalt_ovrf.resize(maximum_valence_work_unit_atoms);
    local_yalt_ovrf.resize(maximum_valence_work_unit_atoms);
    local_zalt_ovrf.resize(maximum_valence_work_unit_atoms);
  }
  if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
    local_vxalt_ovrf.resize(maximum_valence_work_unit_atoms);
    local_vyalt_ovrf.resize(maximum_valence_work_unit_atoms);
    local_vzalt_ovrf.resize(maximum_valence_work_unit_atoms);
  }
  std::vector<int2> vwu_notes(vwu_abstract_length);
  for (int i = 0; i < poly_vk.nvwu; i++) {

    // Obtain the work unit's abstract
    const int abst_llim = vwu_abstract_length * i;
    const int abst_hlim = abst_llim + vwu_abstract_length;
    for (int j = abst_llim; j < abst_hlim; j++) {
      vwu_notes[j - abst_llim] = poly_vk.vwu_abstracts[j];
    }
    
    // Import atoms and velocities
    const int2 import_lims = vwu_notes[static_cast<size_t>(VwuAbstractMap::IMPORT)];
    for (int j = import_lims.x; j < import_lims.y; j++) {
      const size_t jatom = poly_vk.vwu_imports[j];
      const size_t jlocal_idx = j - import_lims.x;
      local_xcrd[jlocal_idx] = poly_psw->xcrd[jatom];
      local_ycrd[jlocal_idx] = poly_psw->ycrd[jatom];
      local_zcrd[jlocal_idx] = poly_psw->zcrd[jatom];
      local_xalt[jlocal_idx] = poly_psw->xalt[jatom];
      local_yalt[jlocal_idx] = poly_psw->yalt[jatom];
      local_zalt[jlocal_idx] = poly_psw->zalt[jatom];
      local_vxalt[jlocal_idx] = poly_psw->vxalt[jatom];
      local_vyalt[jlocal_idx] = poly_psw->vyalt[jatom];
      local_vzalt[jlocal_idx] = poly_psw->vzalt[jatom];
      if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        local_xcrd_ovrf[jlocal_idx] = poly_psw->xcrd_ovrf[jatom];
        local_ycrd_ovrf[jlocal_idx] = poly_psw->ycrd_ovrf[jatom];
        local_zcrd_ovrf[jlocal_idx] = poly_psw->zcrd_ovrf[jatom];
        local_xalt_ovrf[jlocal_idx] = poly_psw->xalt_ovrf[jatom];
        local_yalt_ovrf[jlocal_idx] = poly_psw->yalt_ovrf[jatom];
        local_zalt_ovrf[jlocal_idx] = poly_psw->zalt_ovrf[jatom];
      }
      if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
        local_vxalt_ovrf[jlocal_idx] = poly_psw->vxalt_ovrf[jatom];
        local_vyalt_ovrf[jlocal_idx] = poly_psw->vyalt_ovrf[jatom];
        local_vzalt_ovrf[jlocal_idx] = poly_psw->vzalt_ovrf[jatom];
      }
    }

    // Loop over all instructions
    const int2 insr_lims = vwu_notes[static_cast<size_t>(VwuAbstractMap::SETTLE)];
    for (int j = insr_lims.x; j < insr_lims.y; j++) {
      const uint2 insr = poly_auk.sett_insr[j];
      const int oxy_idx = (insr.x & 0x3ff);
      const int hd1_idx = ((insr.x >> 10) & 0x3ff);
      const int hd2_idx = ((insr.x >> 20) & 0x3ff);
      const Tcalc4 sgeom = poly_auk.settle_geom[insr.y];
      const Tcalc4 smass = poly_auk.settle_mass[insr.y];
      settleGroupPosition<llint, Tcalc, Tcalc3>(oxy_idx, hd1_idx, hd2_idx, smass.z, smass.w,
                                                sgeom.x, sgeom.y, sgeom.z, local_xcrd.data(),
                                                local_ycrd.data(), local_zcrd.data(),
                                                local_xalt.data(), local_yalt.data(),
                                                local_zalt.data(), local_vxalt.data(),
                                                local_vyalt.data(), local_vzalt.data(), dt,
                                                local_xcrd_ovrf.data(), local_ycrd_ovrf.data(),
                                                local_zcrd_ovrf.data(), local_xalt_ovrf.data(),
                                                local_yalt_ovrf.data(), local_zalt_ovrf.data(),
                                                local_vxalt_ovrf.data(), local_vyalt_ovrf.data(),
                                                local_vzalt_ovrf.data(), poly_psw->gpos_scale,
                                                poly_psw->vel_scale);
    }
    
    // Return atom position information which the work unit is responsible for updating
    const int manip_llim = vwu_notes[(size_t)(VwuAbstractMap::MANIPULATE)].x;
    for (int j = import_lims.x; j < import_lims.y; j++) {
      const size_t jatom = poly_vk.vwu_imports[j];
      const size_t jlocal_idx = j - import_lims.x;

      // Determine the atom ownership.  Only the update mask matters, as opposed to the movement
      // mask, as no other work will be done while these work units are temporarily laid out.
      const int mask_segment = manip_llim + (jlocal_idx >> 5);
      const int segment_bit  = jlocal_idx - (mask_segment << 5);
      if (poly_auk.vwu_manip[mask_segment].y & (0x1 << segment_bit)) {
        poly_psw->xalt[jatom] = local_xalt[jlocal_idx];
        poly_psw->yalt[jatom] = local_yalt[jlocal_idx];
        poly_psw->zalt[jatom] = local_zalt[jlocal_idx];
        poly_psw->vxalt[jatom] = local_vxalt[jlocal_idx];
        poly_psw->vyalt[jatom] = local_vyalt[jlocal_idx];
        poly_psw->vzalt[jatom] = local_vzalt[jlocal_idx];
        if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
          poly_psw->xalt_ovrf[jatom] = local_xalt_ovrf[jlocal_idx];
          poly_psw->yalt_ovrf[jatom] = local_yalt_ovrf[jlocal_idx];
          poly_psw->zalt_ovrf[jatom] = local_zalt_ovrf[jlocal_idx];
        }
        if (poly_psw->vel_bits > velocity_scale_nonoverflow_bits) {
          poly_psw->vxalt_ovrf[jatom] = local_vxalt_ovrf[jlocal_idx];
          poly_psw->vyalt_ovrf[jatom] = local_vyalt_ovrf[jlocal_idx];
          poly_psw->vzalt_ovrf[jatom] = local_vzalt_ovrf[jlocal_idx];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void idealizeSettleReference(PhaseSpaceWriter *psw, const ConstraintKit<T> &cnst) {
  const T value_half = 0.5;
  const T value_two = 2.0;
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  for (int i = 0; i < cnst.nsettle; i++) {

    // Locate the various atoms
    const int oxy_idx = cnst.settle_ox_atoms[i];
    const int hd1_idx = cnst.settle_h1_atoms[i];
    const int hd2_idx = cnst.settle_h2_atoms[i];

    // Find the ideal distances between SETTLE oxygen and hydrogen atoms.
    const int parm_idx = cnst.settle_param_idx[i];
    const T hh = cnst.settle_rc[parm_idx] * value_two;
    const T t1 = value_half * cnst.settle_mo[parm_idx] / cnst.settle_mh[parm_idx];
    const T tri_oh = cnst.settle_ra[parm_idx] * (1.0 + t1);
    const T tri_oh2 = tri_oh * tri_oh;
    const T oh = sqrt(tri_oh2 + (cnst.settle_rc[parm_idx] * cnst.settle_rc[parm_idx]));
    const T oxy_x = psw->xcrd[oxy_idx];
    const T oxy_y = psw->ycrd[oxy_idx];
    const T oxy_z = psw->zcrd[oxy_idx];

    // Compute the H-O-H angle, set up a system of local axes, and derive the ideal positions of
    // the hydrogens.
    const T ang_hoh = (tcalc_is_double) ? value_two * asin(value_half * hh / oh) :
                                          value_two * asinf(value_half * hh / oh);
    const double3 oh1 = { psw->xcrd[hd1_idx] - psw->xcrd[oxy_idx],
                          psw->ycrd[hd1_idx] - psw->ycrd[oxy_idx],
                          psw->zcrd[hd1_idx] - psw->zcrd[oxy_idx] };
    const double3 oh2 = { psw->xcrd[hd2_idx] - psw->xcrd[oxy_idx],
                          psw->ycrd[hd2_idx] - psw->ycrd[oxy_idx],
                          psw->zcrd[hd2_idx] - psw->zcrd[oxy_idx] };
    double3 plane_normal = crossProduct(oh1, oh2);
    double3 water_xaxis  = oh1;
    normalize<double3, double>(&plane_normal);
    normalize<double3, double>(&water_xaxis);
    const double3 water_yaxis  = crossProduct(plane_normal, water_xaxis);
    const double hd2_waterx = oh * cos(ang_hoh);
    const double hd2_watery = oh * sin(ang_hoh);
    psw->xcrd[hd1_idx] = psw->xcrd[oxy_idx] + (oh * water_xaxis.x);
    psw->ycrd[hd1_idx] = psw->ycrd[oxy_idx] + (oh * water_xaxis.y);
    psw->zcrd[hd1_idx] = psw->zcrd[oxy_idx] + (oh * water_xaxis.z);
    psw->xcrd[hd2_idx] = psw->xcrd[oxy_idx] +
                         (hd2_waterx * water_xaxis.x) + (hd2_watery * water_yaxis.x);
    psw->ycrd[hd2_idx] = psw->ycrd[oxy_idx] +
                         (hd2_waterx * water_xaxis.y) + (hd2_watery * water_yaxis.y);
    psw->zcrd[hd2_idx] = psw->zcrd[oxy_idx] +
                         (hd2_waterx * water_xaxis.z) + (hd2_watery * water_yaxis.z);
  }
}

} // namespace structure
} // namespace stormm
