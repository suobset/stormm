#include "Math/summation.h"
#include "valence_workunit.h"

namespace omni {
namespace synthesis {

using math::prefixSumInPlace;
using math::PrefixSumType;
using topology::ConstraintKit;
using topology::ValenceKit;
using topology::VirtualSiteKind;
using topology::VirtualSiteKit;
  
//-------------------------------------------------------------------------------------------------
ValenceDelegator::ValenceDelegator(const AtomGraph &ag, const RestraintApparatus &ra) :
    atom_count{ag.getAtomCount()},
    bond_i_presence{2 * ag.getBondTermCount(), -1},
    bond_j_presence{2 * ag.getBondTermCount(), -1},
    angl_i_presence{2 * ag.getAngleTermCount(), -1},
    angl_j_presence{2 * ag.getAngleTermCount(), -1},
    angl_k_presence{2 * ag.getAngleTermCount(), -1},
    dihe_i_presence{2 * ag.getDihedralTermCount(), -1},
    dihe_j_presence{2 * ag.getDihedralTermCount(), -1},
    dihe_k_presence{2 * ag.getDihedralTermCount(), -1},
    dihe_l_presence{2 * ag.getDihedralTermCount(), -1},
    ubrd_i_presence{2 * ag.getUreyBradleyTermCount(), -1},
    ubrd_k_presence{2 * ag.getUreyBradleyTermCount(), -1},
    cimp_i_presence{2 * ag.getCharmmImprTermCount(), -1},
    cimp_j_presence{2 * ag.getCharmmImprTermCount(), -1},
    cimp_k_presence{2 * ag.getCharmmImprTermCount(), -1},
    cimp_l_presence{2 * ag.getCharmmImprTermCount(), -1},
    cmap_i_presence{2 * ag.getCmapTermCount(), -1},
    cmap_j_presence{2 * ag.getCmapTermCount(), -1},
    cmap_k_presence{2 * ag.getCmapTermCount(), -1},
    cmap_l_presence{2 * ag.getCmapTermCount(), -1},
    cmap_m_presence{2 * ag.getCmapTermCount(), -1},
    vs_presence{2 * ag.getVirtualSiteCount(), -1},
    vsf1_presence{2 * ag.getVirtualSiteCount(), -1},
    vsf2_presence{2 * ag.getVirtualSiteCount(), -1},
    vsf3_presence{2 * ag.getVirtualSiteCount(), -1},
    vsf4_presence{2 * ag.getVirtualSiteCount(), -1},
    cnst_n_presence{2 * ag.getConstraintGroupTotalSize(), -1},
    sett_ox_presence{2 * ag.getRigidWaterCount(), -1},
    sett_h1_presence{2 * ag.getRigidWaterCount(), -1},
    sett_h2_presence{2 * ag.getRigidWaterCount(), -1},
    rposn_i_presence{2 * ra.getPositionalRestraintCount(), -1},
    rbond_i_presence{2 * ra.getDistanceRestraintCount(), -1},
    rbond_j_presence{2 * ra.getDistanceRestraintCount(), -1},
    rangl_i_presence{2 * ra.getAngleRestraintCount(), -1},
    rangl_j_presence{2 * ra.getAngleRestraintCount(), -1},
    rangl_k_presence{2 * ra.getAngleRestraintCount(), -1},
    rdihe_i_presence{2 * ra.getDihedralRestraintCount(), -1},
    rdihe_j_presence{2 * ra.getDihedralRestraintCount(), -1},
    rdihe_k_presence{2 * ra.getDihedralRestraintCount(), -1},
    rdihe_l_presence{2 * ra.getDihedralRestraintCount(), -1},
    bond_affector_list{ag.getBondTermCount() * 2},
    bond_affector_bounds{atom_count + 1, 0},
    angl_affector_list{ag.getAngleTermCount() * 3},
    angl_affector_bounds{atom_count + 1, 0},
    dihe_affector_list{ag.getDihedralTermCount() * 4},
    dihe_affector_bounds{atom_count + 1, 0},
    ubrd_affector_list{ag.getUreyBradleyTermCount() * 2},
    ubrd_affector_bounds{atom_count + 1, 0},
    cimp_affector_list{ag.getCharmmImprTermCount() * 4},
    cimp_affector_bounds{atom_count + 1, 0},
    cmap_affector_list{ag.getCmapTermCount() * 5},
    cmap_affector_bounds{atom_count + 1, 0},
    vste_affector_list{ag.getVirtualSiteCount() * 4},
    vste_affector_bounds{atom_count + 1, 0},
    cnst_affector_list{ag.getConstraintGroupTotalSize()},
    cnst_affector_bounds{atom_count + 1, 0},
    sett_affector_list{ag.getRigidWaterCount() * 3},
    sett_affector_bounds{atom_count + 1, 0},
    work_unit_assignments{atom_count, 0},
    work_unit_presence{atom_count * 4, 0}
{
  // Pass through the topology, filling out the valence term affector arrays and the virtual site
  // frame atom arrays.
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const ConstraintKit cnk = ag.getConstraintKit();
  const RestraintApparatusDpReader rar = ra.dpData();
  for (int pos = 0; pos < vk.nbond; pos++) {
    bond_affector_bounds[vk.bond_i_atoms[pos]] += 1;
    bond_affector_bounds[vk.bond_j_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vk.nangl; pos++) {
    angl_affector_bounds[vk.angl_i_atoms[pos]] += 1;
    angl_affector_bounds[vk.angl_j_atoms[pos]] += 1;
    angl_affector_bounds[vk.angl_k_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vk.ndihe; pos++) {
    dihe_affector_bounds[vk.dihe_i_atoms[pos]] += 1;
    dihe_affector_bounds[vk.dihe_j_atoms[pos]] += 1;
    dihe_affector_bounds[vk.dihe_k_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vk.nubrd; pos++) {
    ubrd_affector_bounds[vk.ubrd_i_atoms[pos]] += 1;
    ubrd_affector_bounds[vk.ubrd_k_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vk.ncimp; pos++) {
    cimp_affector_bounds[vk.cimp_i_atoms[pos]] += 1;
    cimp_affector_bounds[vk.cimp_j_atoms[pos]] += 1;
    cimp_affector_bounds[vk.cimp_k_atoms[pos]] += 1;
    cimp_affector_bounds[vk.cimp_l_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vk.ncmap; pos++) {
    cmap_affector_bounds[vk.cmap_i_atoms[pos]] += 1;
    cmap_affector_bounds[vk.cmap_j_atoms[pos]] += 1;
    cmap_affector_bounds[vk.cmap_k_atoms[pos]] += 1;
    cmap_affector_bounds[vk.cmap_l_atoms[pos]] += 1;
    cmap_affector_bounds[vk.cmap_m_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < vsk.nsite; pos++) {
    vste_affector_bounds[vsk.frame1_idx[pos]] += 1;
    vste_affector_bounds[vsk.frame2_idx[pos]] += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      break;
    case VirtualSiteKind::FIXED_4:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      vste_affector_bounds[vsk.frame4_idx[pos]] += 1;
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      cnst_affector_bounds[cnk.group_list[j]] += 1;
    }
  }
  for (int pos = 0; pos < cnk.nsettle; pos++) {
    sett_affector_bounds[cnk.settle_ox_atoms[pos]] += 1;
    sett_affector_bounds[cnk.settle_h1_atoms[pos]] += 1;
    sett_affector_bounds[cnk.settle_h2_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < rar.nposn; pos++) {
    rposn_affector_bounds[rar.rposn_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < rar.nbond; pos++) {
    rbond_affector_bounds[rar.rbond_i_atoms[pos]] += 1;
    rbond_affector_bounds[rar.rbond_j_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < rar.nangl; pos++) {
    rangl_affector_bounds[rar.rangl_i_atoms[pos]] += 1;
    rangl_affector_bounds[rar.rangl_j_atoms[pos]] += 1;
    rangl_affector_bounds[rar.rangl_k_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < rar.ndihe; pos++) {
    rdihe_affector_bounds[rar.rdihe_i_atoms[pos]] += 1;
    rdihe_affector_bounds[rar.rdihe_j_atoms[pos]] += 1;
    rdihe_affector_bounds[rar.rdihe_k_atoms[pos]] += 1;
    rdihe_affector_bounds[rar.rdihe_l_atoms[pos]] += 1;
  }
  prefixSumInPlace<int>(&bond_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&angl_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&dihe_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&ubrd_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cimp_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cmap_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&vste_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cnst_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&sett_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&rposn_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&rbond_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&rangl_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&rdihe_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int i_atom = vk.bond_i_atoms[pos];
    const int j_atom = vk.bond_j_atoms[pos];
    const int list_i_idx = bond_affector_bounds[i_atom];
    const int list_j_idx = bond_affector_bounds[j_atom];
    bond_affector_list[list_i_idx] = pos;
    bond_affector_list[list_j_idx] = pos;
    bond_affector_bounds[i_atom] += 1;
    bond_affector_bounds[j_atom] += 1;
  }
  for (int pos = 0; pos < vk.nangl; pos++) {
    const int i_atom = vk.angl_i_atoms[pos];
    const int j_atom = vk.angl_j_atoms[pos];
    const int k_atom = vk.angl_k_atoms[pos];
    const int list_i_idx = angl_affector_bounds[i_atom];
    const int list_j_idx = angl_affector_bounds[j_atom];
    const int list_k_idx = angl_affector_bounds[k_atom];
    angl_affector_list[list_i_idx] = pos;
    angl_affector_list[list_j_idx] = pos;
    angl_affector_list[list_k_idx] = pos;
    angl_affector_bounds[i_atom] += 1;
    angl_affector_bounds[j_atom] += 1;
    angl_affector_bounds[k_atom] += 1;
  }
  for (int pos = 0; pos < vk.ndihe; pos++) {
    const int i_atom = vk.dihe_i_atoms[pos];
    const int j_atom = vk.dihe_j_atoms[pos];
    const int k_atom = vk.dihe_k_atoms[pos];
    const int l_atom = vk.dihe_l_atoms[pos];
    const int list_i_idx = dihe_affector_bounds[i_atom];
    const int list_j_idx = dihe_affector_bounds[j_atom];
    const int list_k_idx = dihe_affector_bounds[k_atom];
    const int list_l_idx = dihe_affector_bounds[l_atom];
    dihe_affector_list[list_i_idx] = pos;
    dihe_affector_list[list_j_idx] = pos;
    dihe_affector_list[list_k_idx] = pos;
    dihe_affector_list[list_l_idx] = pos;
    dihe_affector_bounds[i_atom] += 1;
    dihe_affector_bounds[j_atom] += 1;
    dihe_affector_bounds[k_atom] += 1;
    dihe_affector_bounds[l_atom] += 1;
  }
  for (int pos = 0; pos < vk.nubrd; pos++) {
    const int i_atom = vk.ubrd_i_atoms[pos];
    const int k_atom = vk.ubrd_k_atoms[pos];
    const int list_i_idx = ubrd_affector_bounds[i_atom];
    const int list_k_idx = ubrd_affector_bounds[k_atom];
    ubrd_affector_list[list_i_idx] = pos;
    ubrd_affector_list[list_k_idx] = pos;
    ubrd_affector_bounds[i_atom] += 1;
    ubrd_affector_bounds[k_atom] += 1;
  }
  for (int pos = 0; pos < vk.ncimp; pos++) {
    const int i_atom = vk.cimp_i_atoms[pos];
    const int j_atom = vk.cimp_j_atoms[pos];
    const int k_atom = vk.cimp_k_atoms[pos];
    const int l_atom = vk.cimp_l_atoms[pos];
    const int list_i_idx = cimp_affector_bounds[i_atom];
    const int list_j_idx = cimp_affector_bounds[j_atom];
    const int list_k_idx = cimp_affector_bounds[k_atom];
    const int list_l_idx = cimp_affector_bounds[l_atom];
    cimp_affector_list[list_i_idx] = pos;
    cimp_affector_list[list_j_idx] = pos;
    cimp_affector_list[list_k_idx] = pos;
    cimp_affector_list[list_l_idx] = pos;
    cimp_affector_bounds[i_atom] += 1;
    cimp_affector_bounds[j_atom] += 1;
    cimp_affector_bounds[k_atom] += 1;
    cimp_affector_bounds[l_atom] += 1;
  }
  for (int pos = 0; pos < vk.ncmap; pos++) {
    const int i_atom = vk.cmap_i_atoms[pos];
    const int j_atom = vk.cmap_j_atoms[pos];
    const int k_atom = vk.cmap_k_atoms[pos];
    const int l_atom = vk.cmap_l_atoms[pos];
    const int m_atom = vk.cmap_m_atoms[pos];
    const int list_i_idx = cmap_affector_bounds[i_atom];
    const int list_j_idx = cmap_affector_bounds[j_atom];
    const int list_k_idx = cmap_affector_bounds[k_atom];
    const int list_l_idx = cmap_affector_bounds[l_atom];
    const int list_m_idx = cmap_affector_bounds[m_atom];
    cmap_affector_list[list_i_idx] = pos;
    cmap_affector_list[list_j_idx] = pos;
    cmap_affector_list[list_k_idx] = pos;
    cmap_affector_list[list_l_idx] = pos;
    cmap_affector_list[list_m_idx] = pos;
    cmap_affector_bounds[i_atom] += 1;
    cmap_affector_bounds[j_atom] += 1;
    cmap_affector_bounds[k_atom] += 1;
    cmap_affector_bounds[l_atom] += 1;
    cmap_affector_bounds[m_atom] += 1;
  }
  for (int pos = 0; pos < vsk.nsite; pos++) {
    const int parent_atom = vsk.frame1_idx[pos];
    const int frame2_atom = vsk.frame2_idx[pos];
    const int list_p_idx  = vste_affector_bounds[parent_atom];
    const int list_f2_idx = vste_affector_bounds[frame2_atom];
    vste_affector_list[list_p_idx ] = pos;
    vste_affector_list[list_f2_idx] = pos;
    vste_affector_bounds[parent_atom] += 1;
    vste_affector_bounds[frame2_atom] += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int frame4_atom = vsk.frame4_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        const int list_f4_idx = vste_affector_bounds[frame4_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_list[list_f4_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
        vste_affector_bounds[frame4_atom] += 1;
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      const int c_atom = cnk.group_list[j];
      const int list_c_idx = cnst_affector_bounds[c_atom];
      cnst_affector_list[list_c_idx] = pos;
      cnst_affector_bounds[c_atom] += 1;
    }
  }
  for (int pos = 0; pos < cnk.nsettle; pos++) {
    const int ox_atom = cnk.settle_ox_atoms[pos];
    const int h1_atom = cnk.settle_h1_atoms[pos];
    const int h2_atom = cnk.settle_h2_atoms[pos];
    const int list_ox_idx = sett_affector_bounds[ox_atom];
    const int list_h1_idx = sett_affector_bounds[h1_atom];
    const int list_h2_idx = sett_affector_bounds[h2_atom];
    sett_affector_list[list_ox_idx] = pos;
    sett_affector_list[list_h1_idx] = pos;
    sett_affector_list[list_h2_idx] = pos;
    sett_affector_bounds[ox_atom] += 1;
    sett_affector_bounds[h1_atom] += 1;
    sett_affector_bounds[h2_atom] += 1;
  }
  for (int pos = 0; pos < rar.nposn; pos++) {
    const int pr_atom = rar.rposn_atoms[pos];
    const int list_pr_idx = rposn_affector_bounds[pr_atom];
    rposn_affector_list[list_pr_idx] = pos;
    rposn_affector_bounds[pr_atom] += 1;
  }
  for (int pos = 0; pos < rar.nbond; pos++) {
    const int i_atom = rar.rbond_i_atoms[pos];
    const int j_atom = rar.rbond_j_atoms[pos];
    const int list_i_idx = rbond_affector_bounds[i_atom];
    const int list_j_idx = rbond_affector_bounds[j_atom];
    rbond_affector_list[list_i_idx] = pos;
    rbond_affector_list[list_j_idx] = pos;
    rbond_affector_bounds[i_atom] += 1;
    rbond_affector_bounds[j_atom] += 1;
  }
  for (int pos = 0; pos < rar.nangl; pos++) {
    const int i_atom = rar.rangl_i_atoms[pos];
    const int j_atom = rar.rangl_j_atoms[pos];
    const int k_atom = rar.rangl_k_atoms[pos];
    const int list_i_idx = rangl_affector_bounds[i_atom];
    const int list_j_idx = rangl_affector_bounds[j_atom];
    const int list_k_idx = rangl_affector_bounds[k_atom];
    rangl_affector_list[list_i_idx] = pos;
    rangl_affector_list[list_j_idx] = pos;
    rangl_affector_list[list_k_idx] = pos;
    rangl_affector_bounds[i_atom] += 1;
    rangl_affector_bounds[j_atom] += 1;
    rangl_affector_bounds[k_atom] += 1;
  }
  for (int pos = 0; pos < rar.ndihe; pos++) {
    const int i_atom = rar.rdihe_i_atoms[pos];
    const int j_atom = rar.rdihe_j_atoms[pos];
    const int k_atom = rar.rdihe_k_atoms[pos];
    const int l_atom = rar.rdihe_l_atoms[pos];
    const int list_i_idx = rdihe_affector_bounds[i_atom];
    const int list_j_idx = rdihe_affector_bounds[j_atom];
    const int list_k_idx = rdihe_affector_bounds[k_atom];
    const int list_l_idx = rdihe_affector_bounds[l_atom];
    rdihe_affector_list[list_i_idx] = pos;
    rdihe_affector_list[list_j_idx] = pos;
    rdihe_affector_list[list_k_idx] = pos;
    rdihe_affector_list[list_l_idx] = pos;
    rdihe_affector_bounds[i_atom] += 1;
    rdihe_affector_bounds[j_atom] += 1;
    rdihe_affector_bounds[k_atom] += 1;
    rdihe_affector_bounds[l_atom] += 1;
  }

  // Cast the atom count to a constant, just to avoid the *this pointer in subsequent loops.
  // Rewind the prefix sums after using them to populate the various affector lists.
  const int natom = atom_count;
  for (int i = natom; i > 0; i--) {
    bond_affector_bounds[i]  = bond_affector_bounds[i - 1];
    angl_affector_bounds[i]  = angl_affector_bounds[i - 1];
    dihe_affector_bounds[i]  = dihe_affector_bounds[i - 1];
    ubrd_affector_bounds[i]  = ubrd_affector_bounds[i - 1];
    cimp_affector_bounds[i]  = cimp_affector_bounds[i - 1];
    cmap_affector_bounds[i]  = cmap_affector_bounds[i - 1];
    vste_affector_bounds[i]  = vste_affector_bounds[i - 1];
    cnst_affector_bounds[i]  = cnst_affector_bounds[i - 1];
    sett_affector_bounds[i]  = sett_affector_bounds[i - 1];
    rposn_affector_bounds[i] = rposn_affector_bounds[i - 1];
    rbond_affector_bounds[i] = rbond_affector_bounds[i - 1];
    rangl_affector_bounds[i] = rangl_affector_bounds[i - 1];
    rdihe_affector_bounds[i] = rdihe_affector_bounds[i - 1];
  }
  bond_affector_bounds[0]  = 0;
  angl_affector_bounds[0]  = 0;
  dihe_affector_bounds[0]  = 0;
  ubrd_affector_bounds[0]  = 0;
  cimp_affector_bounds[0]  = 0;
  cmap_affector_bounds[0]  = 0;
  vste_affector_bounds[0]  = 0;
  cnst_affector_bounds[0]  = 0;
  sett_affector_bounds[0]  = 0;
  rposn_affector_bounds[0] = 0;
  rbond_affector_bounds[0] = 0;
  rangl_affector_bounds[0] = 0;
  rdihe_affector_bounds[0] = 0;
}


//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::ValenceWorkUnit(const AtomGraph &ag, ValenceDelegator *vdel, const int seed_atom,
                                 const int max_atoms) :
    atom_count{ag.getAtomCount()},
    atom_import_list{}
{}

//-------------------------------------------------------------------------------------------------
#if 0
ValenceWorkUnit::addNewAtom(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                            const ConstraintKit<double> &cnsk, const RestraintKit<double> &rstk) {

}
#endif

} // namespace topology
} // namespace omni
