// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
SyPermutorKit<T>::SyPermutorKit(int system_count_in, int perm_map_count_in,
                                const int* perm_map_idx_in, const int2* perm_elements_in,
                                const int* perm_element_bounds_in, const int* system_settings_in,
                                const int* rot_grp_atoms_in, const int* rot_grp_bounds_in,
                                const int* prm_rot_grp_bounds_in, const int* ctx_grp_atoms_in,
                                const int* ctx_grp_bounds_in, const int* prm_ctx_grp_bounds_in,
                                const int* inv_grp_atoms_in, const int* inv_grp_bounds_in,
                                const int* prm_inv_grp_bounds_in, const int4* rot_bond_markers_in,
                                const int4* ctx_bond_markers_in, const int4* chiral_markers_in,
                                const T* rot_bond_settings_in, const T* ctx_bond_settings_in,
                                const int* rot_bond_settings_bounds_in,
                                const int* ctx_bond_settings_bounds_in,
                                const int* chiral_settings_in,
                                const int* chiral_settings_bounds_in) :
    system_count{system_count_in}, perm_map_count{perm_map_count_in},
    perm_map_idx{perm_map_idx_in}, perm_elements{perm_elements_in},
    perm_element_bounds{perm_element_bounds_in}, system_settings{system_settings_in},
    rot_grp_atoms{rot_grp_atoms_in}, rot_grp_bounds{rot_grp_bounds_in},
    prm_rot_grp_bounds{prm_rot_grp_bounds_in}, ctx_grp_atoms{ctx_grp_atoms_in},
    ctx_grp_bounds{ctx_grp_bounds_in}, prm_ctx_grp_bounds{prm_ctx_grp_bounds_in},
    inv_grp_atoms{inv_grp_atoms_in}, inv_grp_bounds{inv_grp_bounds_in},
    prm_inv_grp_bounds{prm_inv_grp_bounds_in}, rot_bond_markers{rot_bond_markers_in},
    ctx_bond_markers{ctx_bond_markers_in}, chiral_markers{chiral_markers_in},
    rot_bond_settings{rot_bond_settings_in}, ctx_bond_settings{ctx_bond_settings_in},
    rot_bond_settings_bounds{rot_bond_settings_bounds_in},
    ctx_bond_settings_bounds{ctx_bond_settings_bounds_in},
    chiral_settings_bounds{chiral_settings_bounds_in}, chiral_settings{chiral_settings_in}
{}

} // namespace synthesis
} // namespace stormm
