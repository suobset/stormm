#include "valence_workunit.h"

namespace omni {
namespace topology {

//-------------------------------------------------------------------------------------------------
ValenceDelegation::ValenceDelegation(const AtomGraph &ag) :
    atom_count{ag.getAtomCount()},
    bond_i_atoms{ag.getBondTermCount()},
    bond_j_atoms{ag.getBondTermCount()},
    angl_i_atoms{ag.getAngleTermCount()},
    angl_j_atoms{ag.getAngleTermCount()},
    angl_k_atoms{ag.getAngleTermCount()},
    dihe_i_atoms{ag.getDihedralTermCount()},
    dihe_j_atoms{ag.getDihedralTermCount()},
    dihe_k_atoms{ag.getDihedralTermCount()},
    dihe_l_atoms{ag.getDihedralTermCount()},
    ubrd_i_atoms{ag.getUreyBradleyTermCount()},
    ubrd_k_atoms{ag.getUreyBradleyTermCount()},
    cimp_i_atoms{ag.getCharmmImprTermCount()},
    cimp_j_atoms{ag.getCharmmImprTermCount()},
    cimp_k_atoms{ag.getCharmmImprTermCount()},
    cimp_l_atoms{ag.getCharmmImprTermCount()},
    cmap_i_atoms{ag.getCmapTermCount()},
    cmap_j_atoms{ag.getCmapTermCount()},
    cmap_k_atoms{ag.getCmapTermCount()},
    cmap_l_atoms{ag.getCmapTermCount()},
    cmap_m_atoms{ag.getCmapTermCount()},
    virtual_site_placement{ag.getVirtualSiteCount()},
    bond_affector_list{ag.getBondTermCount() * 2},
    bond_affector_bounds{atom_count + 1},
    angl_affector_list{ag.getAngleTermCount() * 3},
    angl_affector_bounds{atom_count + 1},
    dihe_affector_list{ag.getDihedralTermCount() * 4},
    dihe_affector_bounds{atom_count + 1},
    ubrd_affector_list{ag.getUreyBradleyTermCount() * 2},
    ubrd_affector_bounds{atom_count + 1},
    cimp_affector_list{ag.getCharmmImprTermCount() * 4},
    cimp_affector_bounds{atom_count + 1},
    cmap_affector_list{ag.getCmapTermCount() * 5},
    cmap_affector_bounds{atom_count + 1},
    vste_affector_list{ag.getVirtualSiteCount() * 4},
    vste_affector_bounds{atom_count + 1}
{
  // Pass through the topology, filling out the valence term affector arrays and the virtual site
  // frame atom arrays.
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  for (int i = 0; i < vk.nbond; i++) {

  }
}

//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::ValenceWorkUnit(const AtomGraph &ag, ValenceDelegation *vdel, const int seed_atom,
                                 const int max_atoms) :
    atom_count{ag.getAtomCount()},
    atom_import_list{}
{}
  
} // namespace topology
} // namespace omni
