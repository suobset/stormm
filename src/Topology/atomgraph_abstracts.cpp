#include "DataTypes/omni_vector_types.h"
#include "atomgraph_abstracts.h"

namespace omni {
namespace topology {

//-------------------------------------------------------------------------------------------------
ChemicalDetailsKit::ChemicalDetailsKit(const int natom_in, const int nres_in, const int nmol_in,
                                       const char4* atom_names_in, const char4* res_names_in,
                                       const char4* atom_types_in, const int* z_numbers_in,
                                       const int* res_limits_in, const int* atom_numbers_in,
                                       const int* res_numbers_in, const int* mol_home_in,
                                       const int* mol_contents_in, const int* mol_limits_in,
                                       const double* masses_in, const float* sp_masses_in) :
    natom{natom_in}, nres{nres_in}, nmol{nmol_in}, atom_names{atom_names_in},
    res_names{res_names_in}, atom_types{atom_types_in}, z_numbers{z_numbers_in},
    res_limits{res_limits_in}, atom_numbers{atom_numbers_in}, res_numbers{res_numbers_in},
    mol_home{mol_home_in}, mol_contents{mol_contents_in}, mol_limits{mol_limits_in},
    masses{masses_in}, sp_masses{sp_masses_in}
{}

//-------------------------------------------------------------------------------------------------
ConstraintKit::ConstraintKit(const int ngroup_in, const int nsettle_in,
                             const int* group_list_in, const int* group_bounds_in,
                             const double* group_targets_in, const double* group_inv_masses_in,
                             const int* settle_ox_atoms_in, const int* settle_h1_atoms_in,
                             const int* settle_h2_atoms_in) :
    ngroup{ngroup_in}, nsettle{nsettle_in}, group_list{group_list_in},
    group_bounds{group_bounds_in}, group_targets{group_targets_in},
    group_inv_masses{group_inv_masses_in}, settle_ox_atoms{settle_ox_atoms_in},
    settle_h1_atoms{settle_h1_atoms_in}, settle_h2_atoms{settle_h2_atoms_in}
{}

} // namespace topology
} // namespace omni
