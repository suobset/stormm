#include "Math/summation.h"
#include "valence_workunit.h"

namespace omni {
namespace topology {

using math::prefixSumInPlace;
using math::PrefixSumType;
  
//-------------------------------------------------------------------------------------------------
ValenceDelegator::ValenceDelegator(const AtomGraph &ag) :
    atom_count{ag.getAtomCount()},
    bond_i_atoms{ag.getBondTermCount(), -1},
    bond_j_atoms{ag.getBondTermCount(), -1},
    angl_i_atoms{ag.getAngleTermCount(), -1},
    angl_j_atoms{ag.getAngleTermCount(), -1},
    angl_k_atoms{ag.getAngleTermCount(), -1},
    dihe_i_atoms{ag.getDihedralTermCount(), -1},
    dihe_j_atoms{ag.getDihedralTermCount(), -1},
    dihe_k_atoms{ag.getDihedralTermCount(), -1},
    dihe_l_atoms{ag.getDihedralTermCount(), -1},
    ubrd_i_atoms{ag.getUreyBradleyTermCount(), -1},
    ubrd_k_atoms{ag.getUreyBradleyTermCount(), -1},
    cimp_i_atoms{ag.getCharmmImprTermCount(), -1},
    cimp_j_atoms{ag.getCharmmImprTermCount(), -1},
    cimp_k_atoms{ag.getCharmmImprTermCount(), -1},
    cimp_l_atoms{ag.getCharmmImprTermCount(), -1},
    cmap_i_atoms{ag.getCmapTermCount(), -1},
    cmap_j_atoms{ag.getCmapTermCount(), -1},
    cmap_k_atoms{ag.getCmapTermCount(), -1},
    cmap_l_atoms{ag.getCmapTermCount(), -1},
    cmap_m_atoms{ag.getCmapTermCount(), -1},
    virtual_site_placement{ag.getVirtualSiteCount(), -1},
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
    work_unit_assignments{atom_count, 0},
    work_unit_presence{atom_count * 4}
{
  // Pass through the topology, filling out the valence term affector arrays and the virtual site
  // frame atom arrays.
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
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
    if (vsk.frame3_idx[pos] >= 0) {
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
    }
    if (vsk.frame4_idx[pos] >= 0) {
      vste_affector_bounds[vsk.frame4_idx[pos]] += 1;
    }
  }
  prefixSumInPlace<int>(&bond_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&angl_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&dihe_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&ubrd_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cimp_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cmap_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&vste_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
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
  
  // Cast the atom count to a constant, just to avoid the *this pointer in subsequent loops.
  // Rewind the prefix sums after using them to populate the various affector lists.
  const int natom = atom_count;
  for (int i = natom; i > 0; i--) {
    bond_affector_bounds[i] = bond_affector_bounds[i - 1];
    angl_affector_bounds[i] = angl_affector_bounds[i - 1];
    dihe_affector_bounds[i] = dihe_affector_bounds[i - 1];
    ubrd_affector_bounds[i] = ubrd_affector_bounds[i - 1];
    cimp_affector_bounds[i] = cimp_affector_bounds[i - 1];
    cmap_affector_bounds[i] = cmap_affector_bounds[i - 1];
    vste_affector_bounds[i] = vste_affector_bounds[i - 1];
  }
  bond_affector_bounds[0] = 0;
  angl_affector_bounds[0] = 0;
  dihe_affector_bounds[0] = 0;
  ubrd_affector_bounds[0] = 0;
  cimp_affector_bounds[0] = 0;
  cmap_affector_bounds[0] = 0;
  vste_affector_bounds[0] = 0;
}

//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::addNewAtom(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                            const ConstraintKit<double> &cnsk, const RestraintKit<double> &rstk) {

}

//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::ValenceWorkUnit(const AtomGraph &ag, ValenceDelegator *vdel, const int seed_atom,
                                 const int max_atoms) :
    atom_count{ag.getAtomCount()},
    atom_import_list{}
{}
  
} // namespace topology
} // namespace omni
