#include "Math/rounding.h"
#include "atomgraph_synthesis.h"
#include "atomgraph_abstracts.h"

namespace omni {
namespace topology {

using math::roundUp;
  
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
    dielectric_constant{1.0}, salt_concentration{0.0}, coulomb_constant{accepted_coulomb_constant},
    pb_radii_set{""}, use_bond_constraints{ShakeSetting::OFF}, use_settle{SettleSetting::OFF},
    water_residue_name{' ', ' ', ' ', ' '},
    topologies{topologies_in},
    topology_indices{HybridKind::POINTER, "tpsyn_indices"},
    atom_counts{HybridKind::POINTER, "typsyn_atom_counts"},
    residue_counts{HybridKind::POINTER, "typsyn_res_counts"},
    molecule_counts{HybridKind::POINTER, "typsyn_mol_counts"},
    largest_residue_sizes{HybridKind::POINTER, "typsyn_max_res"},
    last_solute_residues{HybridKind::POINTER, "typsyn_last_sol_res"},
    last_solute_atoms{HybridKind::POINTER, "typsyn_last_sol_atm"},
    first_solvent_molecules{HybridKind::POINTER, "typsyn_1st_solv_mol"},
    urey_bradley_term_counts{HybridKind::POINTER, "typsyn_ubrd_counts"},
    charmm_impr_term_counts{HybridKind::POINTER, "typsyn_cimp_counts"},
    cmap_term_counts{HybridKind::POINTER, "typsyn_cmap_counts"},
    bond_term_counts{HybridKind::POINTER, "typsyn_bond_counts"},
    angl_term_counts{HybridKind::POINTER, "typsyn_angl_counts"},
    dihe_term_counts{HybridKind::POINTER, "typsyn_dihe_counts"},
    virtual_site_counts{HybridKind::POINTER, "typsyn_vsite_counts"},
    atom_type_counts{HybridKind::POINTER, "typsyn_atype_counts"},
    total_exclusion_counts{HybridKind::POINTER, "typsyn_excl_counts"},
    rigid_water_counts{HybridKind::POINTER, "typsyn_rwat_counts"},
    bond_constraint_counts{HybridKind::POINTER, "typsyn_bcnst_counts"},
    degrees_of_freedom{HybridKind::POINTER, "typsyn_deg_freedom"},
    nonrigid_particle_counts{HybridKind::POINTER, "typsyn_n_nonrigid"},
    atom_offsets{HybridKind::POINTER, "typsyn_atom_offsets"},
    urey_bradley_term_offsets{HybridKind::POINTER, "typsyn_ubrd_offset"},
    charmm_impr_term_offsets{HybridKind::POINTER, "typsyn_cimp_offset"},
    cmap_term_offsets{HybridKind::POINTER, "typsyn_cmap_offset"},
    bond_term_offsets{HybridKind::POINTER, "typsyn_bond_offset"},
    angl_term_offsets{HybridKind::POINTER, "typsyn_angl_offset"},
    dihe_term_offsets{HybridKind::POINTER, "typsyn_dihe_offset"},
    virtual_site_offsets{HybridKind::POINTER, "typsyn_vsite_offset"},
    nb_exclusion_offsets{HybridKind::POINTER, "typsyn_nbexcl_offset"},
    int_system_data{HybridKind::ARRAY, "tpsyn_int_data"},
    residue_limits{HybridKind::ARRAY, "tpsyn_res_lims"},
    atom_struc_numbers{HybridKind::ARRAY, "tpsyn_atom_struc_nums"},
    residue_numbers{HybridKind::ARRAY, "tpsyn_res_numbers"},
    molecule_limits{HybridKind::ARRAY, "tpsyn_mol_limits"},
    atomic_numbers{HybridKind::ARRAY, "tpsyn_znum"},
    mobile_atoms{HybridKind::ARRAY, "tpsyn_belly"},
    molecule_membership{HybridKind::ARRAY, "tpsyn_molnum"},
    molecule_contents{HybridKind::ARRAY, "tpsyn_mol_contents"},
    atomic_charges{HybridKind::ARRAY, "tpsyn_atomq"},
    atomic_masses{HybridKind::ARRAY, "tpsyn_mass"},
    inverse_atomic_masses{HybridKind::ARRAY, "tpsyn_invmass"},
    sp_atomic_charges{HybridKind::ARRAY, "tpsyn_atomq_sp"},
    sp_atomic_masses{HybridKind::ARRAY, "tpsyn_mass_sp"},
    sp_inverse_atomic_masses{HybridKind::ARRAY, "tpsyn_invmass_sp"},
    atom_names{HybridKind::ARRAY, "tpsyn_atom_names"},
    atom_types{HybridKind::ARRAY, "tpsyn_atom_types"},
    res_names{HybridKind::ARRAY, "tpsyn_res_names"},
    urey_bradley_parameter_indicesHybridKind::ARRAY, "tpsyn_urey_parm",
    charmm_impr_parameter_indicesHybridKind::ARRAY, "tpsyn_cimp_parm",
    cmap_surface_indicesHybridKind::ARRAY, "tpsyn_cmap_parm",
    cmap_surface_dimensions{HybridKind::ARRAY, "tpsyn_"},
    urey_bradley_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    urey_bradley_equilibria{HybridKind::ARRAY, "tpsyn_"},
    charmm_impr_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    charmm_impr_phase_angles{HybridKind::ARRAY, "tpsyn_"},
    cmap_surfaces{HybridKind::ARRAY, "tpsyn_"},
    sp_urey_bradley_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    sp_urey_bradley_equilibria{HybridKind::ARRAY, "tpsyn_"},
    sp_charmm_impr_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    sp_charmm_impr_phase_angles{HybridKind::ARRAY, "tpsyn_"},
    sp_cmap_surfaces{HybridKind::ARRAY, "tpsyn_"},
    bond_parameter_indices{HybridKind::ARRAY, "tpsyn_"},
    angl_parameter_indices{HybridKind::ARRAY, "tpsyn_"},
    dihe_parameter_indices{HybridKind::ARRAY, "tpsyn_"},
    bond_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    bond_equilibria{HybridKind::ARRAY, "tpsyn_"},
    angl_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    angl_equilibria{HybridKind::ARRAY, "tpsyn_"},
    dihe_amplitudes{HybridKind::ARRAY, "tpsyn_"},
    dihe_periodicities{HybridKind::ARRAY, "tpsyn_"},
    dihe_phase_angles{HybridKind::ARRAY, "tpsyn_"},
    sp_bond_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    sp_bond_equilibria{HybridKind::ARRAY, "tpsyn_"},
    sp_angl_stiffnesses{HybridKind::ARRAY, "tpsyn_"},
    sp_angl_equilibria{HybridKind::ARRAY, "tpsyn_"},
    sp_dihe_amplitudes{HybridKind::ARRAY, "tpsyn_"},
    sp_dihe_periodicities{HybridKind::ARRAY, "tpsyn_"},
    sp_dihe_phase_angles{HybridKind::ARRAY, "tpsyn_"},
    nmr_initial_steps{HybridKind::ARRAY, "tpsyn_"},
    nmr_final_steps{HybridKind::ARRAY, "tpsyn_"},
    nmr_increments{HybridKind::ARRAY, "tpsyn_"},
    nmr_k_initial_values{HybridKind::ARRAY, "tpsyn_"},
    nmr_r_initial_values{HybridKind::ARRAY, "tpsyn_"},
    nmr_k_final_values{HybridKind::ARRAY, "tpsyn_"},
    nmr_r_final_values{HybridKind::ARRAY, "tpsyn_"},
    sp_nmr_k_initial_values{HybridKind::ARRAY, "tpsyn_"},
    sp_nmr_r_initial_values{HybridKind::ARRAY, "tpsyn_"},
    sp_nmr_r_final_values{HybridKind::ARRAY, "tpsyn_"},
    sp_nmr_r_final_values{HybridKind::ARRAY, "tpsyn_"},
    atom_imports{HybridKind::ARRAY, "tpsyn_"},
    vwu_instruction_sets{HybridKind::ARRAY, "tpsyn_"},
    bond_instructions{HybridKind::ARRAY, "tpsyn_"},
    angl_instructions{HybridKind::ARRAY, "tpsyn_"},
    dihe_instructions{HybridKind::ARRAY, "tpsyn_"},
    cmap_instructions{HybridKind::ARRAY, "tpsyn_"},
    nmr2_instructions{HybridKind::ARRAY, "tpsyn_"},
    nmr3_instructions{HybridKind::ARRAY, "tpsyn_"},
    nmr4_instructions{HybridKind::ARRAY, "tpsyn_"}
{
  // Allocate memory and set POINTER-kind arrays for the small packets of data
  const int padded_system_count = roundUp(system_count, warp_size_int);
  int_data.resize(32 * padded_system_count);
  int pivot = 0;
  topology_indices.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  atom_counts.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  residue_counts.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  molecule_counts.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  largest_residue_sizes.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_residues.setPointer(int_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_atoms.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  first_solvent_molecules.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  urey_bradley_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  charmm_impr_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  cmap_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  bond_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  angl_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  dihe_term_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  virtual_site_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  atom_type_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  total_exclusion_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  rigid_water_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  bond_constraint_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  degrees_of_freedom.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  nonrigid_particle_counts.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  atom_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  residue_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  molecule_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  urey_bradley_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  charmm_impr_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  cmap_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  bond_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  angl_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  dihe_term_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  virtual_site_offsets.setPointer(int_data, 0, system_count);
  pivot += padded_system_count;
  nb_exclusion_offsets.setPointer(int_data, 0, system_count);
  
  // Loop over all systems, fill in the above details, and compute the sizes of various arrays
  int atom_offset = 0;
  int resi_offset = 0;
  int mole_offset = 0;
  int ubrd_offset = 0;
  int cimp_offset = 0;
  int cmap_offset = 0;
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  int vste_offset = 0;
  int excl_offset = 0;
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* ag_ptr = topologies[topology_indices_in[i]];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const NonbondedKit<float> nbk_sp = ag_ptr->getSinglePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const ValenceKit<float> vk_sp    = ag_ptr->getSinglePrecisionValenceKit();
    total_atoms += cdk.natom;
    atom_counts.putHost(cdk.natom, i);
    residue_counts.putHost(cdk.nres, i);
    molecule_counts.putHost(cdk.nmol, i);
    largest_residue_sizes.putHost(ag_ptr->getLargestResidueSize(), i);
    last_solute_residues.putHost(ag_ptr->getLastSoluteResidue(), i);
    last_solute_atoms.putHost(ag_ptr->getLastSoluteAtom(), i);
    first_solvent_molecules.putHost(ag_ptr->getFirstSolventMolecule(), i);
    urey_bradley_term_counts.putHost(vk.nubrd, i);
    charmm_impr_term_counts.putHost(vk.ncimp, i);
    cmap_term_counts.putHost(vk.ncmap, i);
    bond_term_counts.putHost(vk.nbond, i);
    angl_term_counts.putHost(vk.nangl, i);
    dihe_term_counts.putHost(vk.ndihe, i);
    virtual_site_counts.putHost(ag_ptr->getVirtualSiteCount(), i);
    atom_type_counts.putHost(nbk.n_lj_types, i);
    total_exclusion_counts.putHost(ag_ptr->getTotalExclusions(), i);
    rigid_water_counts.putHost(ag_ptr->getRigidWaterCount(), i);
    bond_constraint_counts.putHost(ag_ptr->getBondConstraintCount(), i);
    degrees_of_freedom.putHost(ag_ptr->getDegreesOfFreedom(), i);
    nonrigid_particle_counts.putHost(ag_ptr->getNonrigidParticleCount(), i);

    // Record various offsets, and increment the counters
    atom_offsets.putHost(atom_offset, i);
    atom_offset += roundUp(cdk.natom, warp_size_int);
    residue_offsets.putHost(resi_offset, i);
    resi_offset += roundUp(cdk.nres + 1, warp_size_int);
    molecule_offsets.putHost(mole_offset, i);
    mole_offset += roundUp(cdk.nmol + 1, warp_size_int);
    urey_bradley_term_offsets.putHost(ubrd_offset, i);
    ubrd_offset += roundUp(vk.nubrd, warp_size_int);
    charmm_impr_term_offsets.putHost(cimp_offset, i);
    cimp_offset += roundUp(vk.ncimp, warp_size_int);
    cmap_term_offsets.putHost(cmap_offset, i);
    cmap_offset += roundUp(vk.ncmap, warp_size_int);
    bond_term_offsets.putHost(bond_offset, i);
    bond_offset += roundUp(vk.nbond, warp_size_int);
    angl_term_offsets.putHost(angl_offset, i);
    angl_offset += roundUp(vk.nangl, warp_size_int);
    dihe_term_offsets.putHost(dihe_offset, i);
    dihe_offset += roundUp(vk.ndihe, warp_size_int);
    virtual_site_offsets.putHost(vste_offset, i);
    vste_offset += roundUp(ag_ptr->getVirtaulSiteCount(), warp_size_int);
    nb_exclusion_offsets.putHost(excl_offset, i);
    excl_offset += roundUp(ag_ptr->getTotalExclusions());
  }

  // Allocate detailed arrays for each descriptor, collating all topologies
  residue_limits.resize(resi_offset);
  atom_struc_numbers.resize(atom_offset);
  residue_numbers.resize(resi_offset);
  molecule_limits.resize(mole_offset);
  atomic_numbers.resize(atom_offset);
  mobile_atoms.resize(atom_offset);
  molecule_membership.resize(atom_offset);
  molecule_contents.resize(atom_offset);
  atomic_charges.resize(atom_offset);
  atomic_masses.resize(atom_offset);
  inverse_atomic_masses.resize(atom_offset);
  sp_atomic_charges.resize(atom_offset);
  sp_atomic_masses.resize(atom_offset);
  sp_inverse_atomic_masses.resize(atom_offset);
  atom_names.resize(atom_offset);
  atom_types.resize(atom_offset);
  residue_names.resize(resi_offset);

  // Compute the numbers of unique parameters
  std::vector<bool> covered(atom_offset, false);
  urey_bradley_parameter_indices.resize(ubrd_offset);
  charmm_impr_parameter_indices.resize(cimp_offset);
  cmap_surface_indices.resize(cmap_offset);
  cmap_surface_dimensions.resize(cmap_offset);
  urey_bradley_parameter_stiffnesses.resize(ubrd_offset);
  urey_bradley_parameter_stiffnesses.resize(ubrd_offset);
  
}

} // namespace topology
} // namespace omni
