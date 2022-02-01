#include "atomgraph_synthesis.h"

namespace omni {
namespace topology {

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<int> topology_indices_in) :
  topology_count{static_cast<int>(topologies_in.size())},
  system_count{static_cast<int>(topology_indices_in.size())},
  total_atoms{0}, total_virtual_sites{0}, total_bond_terms{0}, total_angl_terms,
  total_dihe_terms{0}, total_ubrd_terms{0}, total_cimp_terms{0}, total_cmap_terms{0},
  total_atom_types{0}, total_charge_types{0}, total_bond_params{0}, total_angl_params{0},
  total_dihe_params{0}, total_ubrd_params{0}, total_cimp_params{0}, total_cmap_surfaces{0},
  periodic_box_class{UnitCellType::NONE}, gb_style{ImplicitSolventModel::NONE},
  dielectric_constant{80.0}, 
{

  // Loop over all systems and compute the sizes of various arrays

}

} // namespace topology
} // namespace omni
