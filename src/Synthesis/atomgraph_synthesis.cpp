#include "Math/rounding.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/atomgraph_refinement.h"
#include "atomgraph_synthesis.h"

namespace omni {
namespace synthesis {

using card::HybridKind;
using math::roundUp;
using math::maxAbsoluteDifference;
using math::maxValue;
using math::minValue;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::sum;
using parse::realToString;
using restraints::RestraintStage;
using topology::accepted_coulomb_constant;
using topology::CmapAccessories;
using topology::ComputeCmapDerivatives;
using testing::Approx;
  
//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<RestraintApparatus*> &restraints_in,
                                       const std::vector<int> &topology_indices_in,
                                       const std::vector<int> &restraint_indices_in,
                                       const ExceptionResponse policy_in) :
    policy{policy_in},
    topology_count{static_cast<int>(topologies_in.size())},
    restraint_network_count{static_cast<int>(restraints_in.size())},
    system_count{static_cast<int>(topology_indices_in.size())},
    total_atoms{0}, total_virtual_sites{0}, total_bond_terms{0}, total_angl_terms{0},
    total_dihe_terms{0}, total_ubrd_terms{0}, total_cimp_terms{0}, total_cmap_terms{0},
    total_atom_types{0}, total_charge_types{0}, total_bond_params{0}, total_angl_params{0},
    total_dihe_params{0}, total_ubrd_params{0}, total_cimp_params{0}, total_cmap_surfaces{0},
    total_position_restraints{0}, total_distance_restraints{0}, total_angle_restraints{0},
    total_dihedral_restraints{0}, periodic_box_class{UnitCellType::NONE},
    gb_style{ImplicitSolventModel::NONE}, dielectric_constant{1.0}, salt_concentration{0.0},
    coulomb_constant{accepted_coulomb_constant}, use_bond_constraints{ShakeSetting::OFF},
    use_settle{SettleSetting::OFF}, water_residue_name{' ', ' ', ' ', ' '}, pb_radii_sets{},
    topologies{topologies_in},
    restraint_networks{restraints_in},
    topology_indices{HybridKind::POINTER, "tpsyn_top_indices"},
    restraint_indices{HybridKind::POINTER, "tpsyn_rst_indices"},
    atom_counts{HybridKind::POINTER, "tpsyn_atom_counts"},
    residue_counts{HybridKind::POINTER, "tpsyn_res_counts"},
    molecule_counts{HybridKind::POINTER, "tpsyn_mol_counts"},
    largest_residue_sizes{HybridKind::POINTER, "tpsyn_max_res"},
    last_solute_residues{HybridKind::POINTER, "tpsyn_last_sol_res"},
    last_solute_atoms{HybridKind::POINTER, "tpsyn_last_sol_atm"},
    first_solvent_molecules{HybridKind::POINTER, "tpsyn_1st_solv_mol"},
    ubrd_term_counts{HybridKind::POINTER, "tpsyn_ubrd_counts"},
    cimp_term_counts{HybridKind::POINTER, "tpsyn_cimp_counts"},
    cmap_term_counts{HybridKind::POINTER, "tpsyn_cmap_counts"},
    bond_term_counts{HybridKind::POINTER, "tpsyn_bond_counts"},
    angl_term_counts{HybridKind::POINTER, "tpsyn_angl_counts"},
    dihe_term_counts{HybridKind::POINTER, "tpsyn_dihe_counts"},
    virtual_site_counts{HybridKind::POINTER, "tpsyn_vsite_counts"},
    posn_restraint_counts{HybridKind::POINTER, "tpsyn_rposn_counts"},
    bond_restraint_counts{HybridKind::POINTER, "tpsyn_rbond_counts"},
    angl_restraint_counts{HybridKind::POINTER, "tpsyn_rangl_counts"},
    dihe_restraint_counts{HybridKind::POINTER, "tpsyn_rdihe_counts"},
    atom_type_counts{HybridKind::POINTER, "tpsyn_atype_counts"},
    total_exclusion_counts{HybridKind::POINTER, "tpsyn_excl_counts"},
    rigid_water_counts{HybridKind::POINTER, "tpsyn_rwat_counts"},
    bond_constraint_counts{HybridKind::POINTER, "tpsyn_bcnst_counts"},
    degrees_of_freedom{HybridKind::POINTER, "tpsyn_deg_freedom"},
    nonrigid_particle_counts{HybridKind::POINTER, "tpsyn_n_nonrigid"},
    atom_offsets{HybridKind::POINTER, "tpsyn_atom_offsets"},
    atom_bit_offsets{HybridKind::POINTER, "tpsyn_abit_offsets"},
    residue_offsets{HybridKind::POINTER, "tpsyn_res_offsets"},
    molecule_offsets{HybridKind::POINTER, "tpsyn_res_offsets"},
    ubrd_term_offsets{HybridKind::POINTER, "tpsyn_ubrd_offset"},
    cimp_term_offsets{HybridKind::POINTER, "tpsyn_cimp_offset"},
    cmap_term_offsets{HybridKind::POINTER, "tpsyn_cmap_offset"},
    bond_term_offsets{HybridKind::POINTER, "tpsyn_bond_offset"},
    angl_term_offsets{HybridKind::POINTER, "tpsyn_angl_offset"},
    dihe_term_offsets{HybridKind::POINTER, "tpsyn_dihe_offset"},
    virtual_site_offsets{HybridKind::POINTER, "tpsyn_vsite_offset"},
    posn_restraint_offsets{HybridKind::POINTER, "tpsyn_rposn_offset"},
    bond_restraint_offsets{HybridKind::POINTER, "tpsyn_rbond_offset"},
    angl_restraint_offsets{HybridKind::POINTER, "tpsyn_rangl_offset"},
    dihe_restraint_offsets{HybridKind::POINTER, "tpsyn_rdihe_offset"},
    nb_exclusion_offsets{HybridKind::POINTER, "tpsyn_nbexcl_offset"},
    lennard_jones_abc_offsets{HybridKind::POINTER, "tpsyn_ljtable_offset"},
    int_system_data{HybridKind::ARRAY, "tpsyn_int_data"},
    residue_limits{HybridKind::POINTER, "tpsyn_res_lims"},
    atom_struc_numbers{HybridKind::POINTER, "tpsyn_atom_struc_nums"},
    residue_numbers{HybridKind::POINTER, "tpsyn_res_numbers"},
    molecule_limits{HybridKind::POINTER, "tpsyn_mol_limits"},
    atomic_numbers{HybridKind::POINTER, "tpsyn_znum"},
    mobile_atoms{HybridKind::POINTER, "tpsyn_belly"},
    molecule_membership{HybridKind::POINTER, "tpsyn_molnum"},
    molecule_contents{HybridKind::POINTER, "tpsyn_mol_contents"},
    atomic_charges{HybridKind::POINTER, "tpsyn_atomq"},
    atomic_masses{HybridKind::POINTER, "tpsyn_mass"},
    inverse_atomic_masses{HybridKind::POINTER, "tpsyn_invmass"},
    sp_atomic_charges{HybridKind::POINTER, "tpsyn_atomq_sp"},
    sp_atomic_masses{HybridKind::POINTER, "tpsyn_mass_sp"},
    sp_inverse_atomic_masses{HybridKind::POINTER, "tpsyn_invmass_sp"},
    atom_names{HybridKind::POINTER, "tpsyn_atom_names"},
    atom_types{HybridKind::POINTER, "tpsyn_atom_types"},
    residue_names{HybridKind::POINTER, "tpsyn_res_names"},
    chem_int_data{HybridKind::ARRAY, "tpsyn_chem_ints"},
    chem_double_data{HybridKind::ARRAY, "tpsyn_chem_doubles"},
    chem_float_data{HybridKind::ARRAY, "tpsyn_chem_floats"},
    chem_char4_data{HybridKind::ARRAY, "tpsyn_chem_char4s"},
    ubrd_stiffnesses{HybridKind::POINTER, "tpsyn_ub_stiff"},
    ubrd_equilibria{HybridKind::POINTER, "tpsyn_ub_equil"},
    cimp_stiffnesses{HybridKind::POINTER, "tpsyn_cimp_stiff"},
    cimp_phase_angles{HybridKind::POINTER, "tpsyn_cimp_equil"},
    cmap_surface_dimensions{HybridKind::POINTER, "tpsyn_cmap_dims"},
    cmap_surface_bounds{HybridKind::POINTER, "tpsyn_cmap_bounds"},
    cmap_patch_bounds{HybridKind::POINTER, "tpsyn_cmpatch_bounds"},
    cmap_surfaces{HybridKind::POINTER, "tpsyn_cmap_surf"},
    cmap_patches{HybridKind::POINTER, "tpsyn_cmap_patch"},
    sp_ubrd_stiffnesses{HybridKind::POINTER, "tpsyn_ub_stiff_sp"},
    sp_ubrd_equilibria{HybridKind::POINTER, "tpsyn_ub_equil_sp"},
    sp_cimp_stiffnesses{HybridKind::POINTER, "tpsyn_cimp_stiff_sp"},
    sp_cimp_phase_angles{HybridKind::POINTER, "tpsyn_cimp_equil_sp"},
    sp_cmap_surfaces{HybridKind::POINTER, "tpsyn_cmap_surf_sp"},
    sp_cmap_patches{HybridKind::POINTER, "tpsyn_cmap_patch_sp"},
    bond_stiffnesses{HybridKind::POINTER, "tpsyn_bondk"},
    bond_equilibria{HybridKind::POINTER, "tpsyn_bondl0"},
    angl_stiffnesses{HybridKind::POINTER, "tpsyn_anglk"},
    angl_equilibria{HybridKind::POINTER, "tpsyn_anglt0"},
    dihe_amplitudes{HybridKind::POINTER, "tpsyn_dihek"},
    dihe_periodicities{HybridKind::POINTER, "tpsyn_dihen"},
    dihe_phase_angles{HybridKind::POINTER, "tpsyn_dihepsi"},
    sp_bond_stiffnesses{HybridKind::POINTER, "tpsyn_bondk_sp"},
    sp_bond_equilibria{HybridKind::POINTER, "tpsyn_bondl0_sp"},
    sp_angl_stiffnesses{HybridKind::POINTER, "tpsyn_anglk_sp"},
    sp_angl_equilibria{HybridKind::POINTER, "tpsyn_anglt0_sp"},
    sp_dihe_amplitudes{HybridKind::POINTER, "tpsyn_dihek_sp"},
    sp_dihe_periodicities{HybridKind::POINTER, "tpsyn_dihen_sp"},
    sp_dihe_phase_angles{HybridKind::POINTER, "tpsyn_dihepsi_sp"},
    valparam_double_data{HybridKind::ARRAY, "tpsyn_vparm_dbl"},
    valparam_float_data{HybridKind::ARRAY, "tpsyn_vparm_flt"},
    valparam_int_data{HybridKind::ARRAY, "tpsyn_vparm_int"},
    ubrd_i_atoms{HybridKind::POINTER, "tpsyn_ubrd_i"},
    ubrd_k_atoms{HybridKind::POINTER, "tpsyn_ubrd_k"},
    ubrd_param_idx{HybridKind::POINTER, "tpsyn_ubrd_idx"},
    cimp_i_atoms{HybridKind::POINTER, "tpsyn_cimp_i"},
    cimp_j_atoms{HybridKind::POINTER, "tpsyn_cimp_j"},
    cimp_k_atoms{HybridKind::POINTER, "tpsyn_cimp_k"},
    cimp_l_atoms{HybridKind::POINTER, "tpsyn_cimp_l"},
    cimp_param_idx{HybridKind::POINTER, "tpsyn_cimp_idx"},
    cmap_i_atoms{HybridKind::POINTER, "tpsyn_cmap_i"},
    cmap_j_atoms{HybridKind::POINTER, "tpsyn_cmap_j"},
    cmap_k_atoms{HybridKind::POINTER, "tpsyn_cmap_k"},
    cmap_l_atoms{HybridKind::POINTER, "tpsyn_cmap_l"},
    cmap_m_atoms{HybridKind::POINTER, "tpsyn_cmap_m"},
    cmap_param_idx{HybridKind::POINTER, "tpsyn_cmap_idx"},
    bond_i_atoms{HybridKind::POINTER, "tpsyn_bond_i"},
    bond_j_atoms{HybridKind::POINTER, "tpsyn_bond_j"},
    bond_param_idx{HybridKind::POINTER, "tpsyn_bond_idx"},
    angl_i_atoms{HybridKind::POINTER, "tpsyn_angl_i"},
    angl_j_atoms{HybridKind::POINTER, "tpsyn_angl_j"},
    angl_k_atoms{HybridKind::POINTER, "tpsyn_angl_k"},
    angl_param_idx{HybridKind::POINTER, "tpsyn_angl_idx"},
    dihe_i_atoms{HybridKind::POINTER, "tpsyn_dihe_i"},
    dihe_j_atoms{HybridKind::POINTER, "tpsyn_dihe_j"},
    dihe_k_atoms{HybridKind::POINTER, "tpsyn_dihe_k"},
    dihe_l_atoms{HybridKind::POINTER, "tpsyn_dihe_l"},
    dihe_param_idx{HybridKind::POINTER, "tpsyn_dihe_idx"},
    valence_int_data{HybridKind::ARRAY, "tpsyn_val_ints"},
    charge_indices{HybridKind::POINTER, "tpsyn_q_idx"},
    lennard_jones_indices{HybridKind::POINTER, "tpsyn_lj_idx"},
    lennard_jones_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_ab"},
    lennard_jones_c_coeff{HybridKind::ARRAY, "tpsyn_lj_c"},
    lennard_jones_14_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_14_ab"},
    lennard_jones_14_c_coeff{HybridKind::ARRAY, "tpsyn_lj_14_c"},
    sp_lennard_jones_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_ab_sp"},
    sp_lennard_jones_c_coeff{HybridKind::ARRAY, "tpsyn_lj_c_sp"},
    sp_lennard_jones_14_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_14_ab_sp"},
    sp_lennard_jones_14_c_coeff{HybridKind::ARRAY, "tpsyn_lj_14_c_sp"},
    rposn_step_bounds{HybridKind::POINTER, "tpsyn_rposn_steps"},
    rbond_step_bounds{HybridKind::POINTER, "tpsyn_rbond_steps"},
    rangl_step_bounds{HybridKind::POINTER, "tpsyn_rangl_steps"},
    rdihe_step_bounds{HybridKind::POINTER, "tpsyn_rdihe_steps"},
    rposn_init_k{HybridKind::POINTER, "tpsyn_rposn_init_k"},
    rposn_final_k{HybridKind::POINTER, "tpsyn_rposn_finl_k"},
    rposn_init_r{HybridKind::POINTER, "tpsyn_rposn_init_k"},
    rposn_final_r{HybridKind::POINTER, "tpsyn_rposn_finl_k"},
    rposn_init_xy{HybridKind::POINTER, "tpsyn_rposn_init_xy"},
    rposn_init_z{HybridKind::POINTER, "tpsyn_rposn_init_z"},
    rposn_final_xy{HybridKind::POINTER, "tpsyn_rposn_final_xy"},
    rposn_final_z{HybridKind::POINTER, "tpsyn_rposn_final_z"},
    rbond_init_k{HybridKind::POINTER, "tpsyn_rbond_init_k"},
    rbond_final_k{HybridKind::POINTER, "tpsyn_rbond_finl_k"},
    rbond_init_r{HybridKind::POINTER, "tpsyn_rbond_init_k"},
    rbond_final_r{HybridKind::POINTER, "tpsyn_rbond_finl_k"},
    rangl_init_k{HybridKind::POINTER, "tpsyn_rangl_init_k"},
    rangl_final_k{HybridKind::POINTER, "tpsyn_rangl_finl_k"},
    rangl_init_r{HybridKind::POINTER, "tpsyn_rangl_init_k"},
    rangl_final_r{HybridKind::POINTER, "tpsyn_rangl_finl_k"},
    rdihe_init_k{HybridKind::POINTER, "tpsyn_rdihe_init_k"},
    rdihe_final_k{HybridKind::POINTER, "tpsyn_rdihe_finl_k"},
    rdihe_init_r{HybridKind::POINTER, "tpsyn_rdihe_init_k"},
    rdihe_final_r{HybridKind::POINTER, "tpsyn_rdihe_finl_k"},
    sp_rposn_init_k{HybridKind::POINTER, "tpsynf_rposn_init_k"},
    sp_rposn_final_k{HybridKind::POINTER, "tpsynf_rposn_finl_k"},
    sp_rposn_init_r{HybridKind::POINTER, "tpsynf_rposn_init_k"},
    sp_rposn_final_r{HybridKind::POINTER, "tpsynf_rposn_finl_k"},
    sp_rposn_init_xy{HybridKind::POINTER, "tpsynf_rposn_init_xy"},
    sp_rposn_init_z{HybridKind::POINTER, "tpsynf_rposn_init_z"},
    sp_rposn_final_xy{HybridKind::POINTER, "tpsynf_rposn_final_xy"},
    sp_rposn_final_z{HybridKind::POINTER, "tpsynf_rposn_final_z"},
    sp_rbond_init_k{HybridKind::POINTER, "tpsynf_rbond_init_k"},
    sp_rbond_final_k{HybridKind::POINTER, "tpsynf_rbond_finl_k"},
    sp_rbond_init_r{HybridKind::POINTER, "tpsynf_rbond_init_k"},
    sp_rbond_final_r{HybridKind::POINTER, "tpsynf_rbond_finl_k"},
    sp_rangl_init_k{HybridKind::POINTER, "tpsynf_rangl_init_k"},
    sp_rangl_final_k{HybridKind::POINTER, "tpsynf_rangl_finl_k"},
    sp_rangl_init_r{HybridKind::POINTER, "tpsynf_rangl_init_k"},
    sp_rangl_final_r{HybridKind::POINTER, "tpsynf_rangl_finl_k"},
    sp_rdihe_init_k{HybridKind::POINTER, "tpsynf_rdihe_init_k"},
    sp_rdihe_final_k{HybridKind::POINTER, "tpsynf_rdihe_finl_k"},
    sp_rdihe_init_r{HybridKind::POINTER, "tpsynf_rdihe_init_k"},
    sp_rdihe_final_r{HybridKind::POINTER, "tpsynf_rdihe_finl_k"},
    nmr_int2_data{HybridKind::ARRAY, "tpsyn_nmr_int2_data"},
    nmr_double_data{HybridKind::ARRAY, "tpsyn_nmr_dbl_data"},
    nmr_double2_data{HybridKind::ARRAY, "tpsyn_nmr_dbl2_data"},
    nmr_double4_data{HybridKind::ARRAY, "tpsyn_nmr_dbl4_data"},
    nmr_float_data{HybridKind::ARRAY, "tpsyn_nmr_flt_data"},
    nmr_float2_data{HybridKind::ARRAY, "tpsyn_nmr_flt2_data"},
    nmr_float4_data{HybridKind::ARRAY, "tpsyn_nmr_flt4_data"},
    rposn_atoms{HybridKind::POINTER, "tpsyn_rposn_at"},
    rposn_kr_param_idx{HybridKind::POINTER, "tpsyn_rposn_kr_idx"},
    rposn_xyz_param_idx{HybridKind::POINTER, "tpsyn_rposn_xyz_idx"},
    rbond_i_atoms{HybridKind::POINTER, "tpsyn_rbond_iat"},
    rbond_j_atoms{HybridKind::POINTER, "tpsyn_rbond_jat"},
    rbond_param_idx{HybridKind::POINTER, "tpsyn_rbond_param"},
    rangl_i_atoms{HybridKind::POINTER, "tpsyn_rangl_iat"},
    rangl_j_atoms{HybridKind::POINTER, "tpsyn_rangl_jat"},
    rangl_k_atoms{HybridKind::POINTER, "tpsyn_rangl_kat"},
    rangl_param_idx{HybridKind::POINTER, "tpsyn_rangl_param"},
    rdihe_i_atoms{HybridKind::POINTER, "tpsyn_rdihe_iat"},
    rdihe_j_atoms{HybridKind::POINTER, "tpsyn_rdihe_jat"},
    rdihe_k_atoms{HybridKind::POINTER, "tpsyn_rdihe_kat"},
    rdihe_l_atoms{HybridKind::POINTER, "tpsyn_rdihe_lat"},
    rdihe_param_idx{HybridKind::POINTER, "tpsyn_rdihe_param"},
    nmr_int_data{HybridKind::ARRAY, "tpsyn_nmr_ints"},
    virtual_site_parameters{HybridKind::ARRAY, "tpsyn_vs_params"},
    sp_virtual_site_parameters{HybridKind::ARRAY, "tpsynf_vs_params"},
    virtual_site_atoms{HybridKind::POINTER, "tpsyn_vs_atoms"},
    virtual_site_frame1_atoms{HybridKind::POINTER, "tpsyn_vs_frame1"},
    virtual_site_frame2_atoms{HybridKind::POINTER, "tpsyn_vs_frame2"},
    virtual_site_frame3_atoms{HybridKind::POINTER, "tpsyn_vs_frame3"},
    virtual_site_frame4_atoms{HybridKind::POINTER, "tpsyn_vs_frame4"},
    virtual_site_parameter_indices{HybridKind::POINTER, "tpsyn_vs_param_idx"},
    vsite_int_data{HybridKind::ARRAY, "tpsyn_vsite_ints"},
    atom_import_list{HybridKind::ARRAY, "tpsyn_atom_imports"},
    atom_manipulation_mask{HybridKind::POINTER, "tpsyn_atom_manip"},
    vwu_instruction_sets{HybridKind::ARRAY, "tpsyn_vwu_insr_sets"},
    bond_instructions{HybridKind::POINTER, "tpsyn_bond_insr"},
    angl_instructions{HybridKind::POINTER, "tpsyn_angl_insr"},
    cdhe_instructions{HybridKind::POINTER, "tpsyn_cdhe_insr"},
    cdhe_overtones{HybridKind::POINTER, "tpsyn_ovrt_insr"},
    cmap_instructions{HybridKind::POINTER, "tpsyn_cmap_insr"},
    rposn_instructions{HybridKind::POINTER, "tpsyn_nmr1_insr"},
    rbond_instructions{HybridKind::POINTER, "tpsyn_nmr2_insr"},
    rangl_instructions{HybridKind::POINTER, "tpsyn_nmr3_insr"},
    rdihe_instructions{HybridKind::POINTER, "tpsyn_nmr4_insr"},
    insr_uint_data{HybridKind::ARRAY, "tpsyn_insr_data1"},
    insr_uint2_data{HybridKind::ARRAY, "tpsyn_insr_data2"}
{
  // Setup and memory layout
  const std::vector<int> topology_index_rebase = checkTopologyList(topology_indices_in);
  const std::vector<int> restraint_index_rebase = checkRestraintList(restraint_indices_in,
                                                                     topology_index_rebase);
  checkCommonSettings();
  buildAtomAndTermArrays(topology_indices_in, topology_index_rebase, restraint_indices_in,
                         restraint_index_rebase);

  // Condense valence and charge parameters into compact tables of unique values
  condenseParameterTables();
  
  // The unique Lennard-Jones parameters are hard to map out.  To the degree that there are
  // unique parameters in each system, the Lennard-Jones tables might as well be block matrices.
  // Keep a list of the A, B, and possibly C coefficients of each Lennard-Jones type interacting
  // with all other types.
  extendLJMatrices();

  // Restraint data can now be incorporated.
  condenseRestraintNetworks();
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<int> &topology_indices_in,
                                       const ExceptionResponse policy_in) :
    AtomGraphSynthesis(topologies_in, std::vector<RestraintApparatus*>(1, nullptr),
                       topology_indices_in, std::vector<int>(topology_indices_in.size(), 0),
                       policy_in)
{}

//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::checkTopologyList(const std::vector<int> &topology_indices_in) {

  // Check that no system indexes a topology outside of the given list
  if (maxValue(topology_indices_in) >= topology_count || minValue(topology_indices_in) < 0) {
    rtErr("One or more systems references a topology index outside of the list supplied.",
          "AtomGraphSynthesis", "checkTopologyList");
  }
  
  // Check that each topology in the supplied list is being used.  Roll that result into an
  // analysis of the uniqueness of each topology.
  std::vector<bool> topology_unique(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    topology_unique[topology_indices_in[i]] = true;
  }

  // Check that all topologies are, in fact, unique.  Compact the list if necessary and update
  // the system topology indexing.
  for (int i = 0; i < topology_count; i++) {
    if (topology_unique[i]) {
      const int topi_natom = topologies[i]->getAtomCount();
      for (int j = i + 1; j < topology_count; j++) {
        if (topologies[j] == topologies[i] || topologies[j]->getAtomCount() == topi_natom &&
            topologies[i]->getFileName() == topologies[j]->getFileName()) {
          topology_unique[j] = false;
        }
      }
    }
  }

  // Warn if there are unused topologies
  int n_unused_topology = 0;
  for (int i = 0; i < topology_count; i++) {
    n_unused_topology += (topology_unique[i] == false);
  }
  if (n_unused_topology > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Out of " + std::to_string(topology_count) + " topologies, " +
             std::to_string(n_unused_topology) + " are not referenced by any systems in this "
             "synthesis.", "AtomGraphSynthesis", "checkTopologyList");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Rebase the list of topology indices
  std::vector<int> topology_index_rebase(topology_count, -1);
  int n_unique_top = 0;
  for (int i = 0; i < topology_count; i++) {
    topology_index_rebase[i] = n_unique_top;
    topologies[n_unique_top] = topologies[i];
    if (topology_unique[i]) {
      n_unique_top++;
    }
  }
  topology_count = n_unique_top;
  if (topology_count == 0) {
    rtErr("No topologies were detected to describe " + std::to_string(system_count) + " systems.",
          "AtomGraphSynthesis", "checkTopologyList");
  }
  if (system_count == 0) {
    rtErr("No systems making use of any of the " + std::to_string(topology_count) +
          " topologies were detected.", "AtomGraphSynthesis", "checkTopologyList");
  }
  return topology_index_rebase;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::checkRestraintList(const std::vector<int> &restraint_indices_in,
                                       const std::vector<int> &topology_indices) {
  
  // For each system in the synthesis, there is a restraint group and a topology.  Each restraint
  // group references a topology, to ensure consistency in future applications.  Check that the
  // restraint apparatus for each system references the same topology as the system itself.
  // A reference to a restraint apparatus index less than zero implies that the system uses no
  // restraints.
  if (static_cast<int>(restraint_indices_in.size()) != system_count) {
    rtErr("A restraint network index must be provided for each system that the synthesis is "
          "designed to describe, even if the restraint apparatus pointer it indicates is a "
          "null pointer.  The number of systems in the synthesis is determined by the length of "
          "the topology index array, so there must be as many indices into a list of restraint "
          "networks as there are indices into distinct system topologies.  " +
          std::to_string(static_cast<int>(restraint_indices_in.size())) + " restraint network "
          "pointer indices were provided for a synthesis that should contain " +
          std::to_string(system_count) + " systems.", "AtomGraphSynthesis", "checkRestraintList");
  }
  if (maxValue(restraint_indices_in) >= restraint_network_count) {
    rtErr("One or more systems references a restraint apparatus index outside of the list "
          "supplied.", "AtomGraphSynthesis", "checkRestraintList");
  }
  for (int i = 0; i < system_count; i++) {
    if (restraint_indices_in[i] < 0 || restraint_networks[restraint_indices_in[i]] == nullptr) {
      continue;
    }
    const AtomGraph* topref = topologies[topology_indices[i]];    
    const AtomGraph* rstref = restraint_networks[restraint_indices_in[i]]->getTopologyPointer();
    if (rstref != topref) {
      rtErr("Mismatch in topologies referenced by the restraint apparatus for system " +
            std::to_string(i) + " and the system itself.  Atom counts of the topologies are " +
            std::to_string(rstref->getAtomCount()) + " and " +
            std::to_string(topref->getAtomCount()) + ".", "AtomGraphSynthesis",
            "checkRestraintList");
    }
  }

  // Create a list of unique restraint apparatuses in use by the AtomGraphSynthesis, analogous
  // to the list created for topologies.
  std::vector<bool> network_unique(restraint_network_count, false);
  for (int i = 0; i < system_count; i++) {
    network_unique[restraint_indices_in[i]] = true;
  }
  for (int i = 0; i < restraint_network_count; i++) {
    if (network_unique[i]) {
      for (int j = i + 1; j < restraint_network_count; j++) {
        if (restraint_networks[j] == restraint_networks[i]) {
          network_unique[j] = false;
          continue;
        }
      }
    }
  }
  
  // Warn if there are unused restraint apparatuses
  int n_unused_network = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    n_unused_network += (network_unique[i] == false);
  }
  if (n_unused_network > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Out of " + std::to_string(restraint_network_count) + " networks, " +
             std::to_string(n_unused_network) + " are not referenced by any systems in this "
             "synthesis.", "AtomGraphSynthesis", "checkRestraintList");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Rebase the list of restraint network indices
  std::vector<int> restraint_index_rebase(restraint_network_count, -1);
  int n_unique_network = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    restraint_index_rebase[i] = n_unique_network;
    restraint_networks[n_unique_network] = restraint_networks[i];
    if (network_unique[i]) {
      n_unique_network++;
    }
  }
  restraint_network_count = n_unique_network;
  if (restraint_network_count == 0) {
    rtErr("No restraint apparatuses were detected to describe " + std::to_string(system_count) +
          " systems.", "AtomGraphSynthesis", "checkRestraintList");
  }
  return restraint_index_rebase;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::checkCommonSettings() {

  // Check that all topologies contain compatible boundary conditions and solvent models
  periodic_box_class = topologies[0]->getUnitCellType();
  gb_style = topologies[0]->getImplicitSolventModel();
  dielectric_constant = topologies[0]->getDielectricConstant();
  salt_concentration = topologies[0]->getSaltConcentration();
  coulomb_constant = topologies[0]->getCoulombConstant();
  const Approx appr_dielcon(dielectric_constant, constants::tiny);
  const Approx appr_saltcon(salt_concentration, constants::tiny);
  const Approx appr_coulomb(coulomb_constant, constants::tiny);
  bool ism_problem  = false;
  bool diel_problem = false;
  bool salt_problem = false;
  bool coul_problem = false;
  for (int i = 1; i < topology_count; i++) {
    ism_problem =  (ism_problem  || (topologies[i]->getImplicitSolventModel() != gb_style));
    diel_problem = (diel_problem ||
                    appr_dielcon.test(topologies[i]->getDielectricConstant()) == false);
    salt_problem = (salt_problem ||
                    appr_saltcon.test(topologies[i]->getSaltConcentration()) == false);
    coul_problem = (coul_problem ||
                    appr_coulomb.test(topologies[i]->getCoulombConstant()) == false);
  }
  if (ism_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must have a consistent implicit solvent model setting.  The first, "
            "set to " + getImplicitSolventModelName(gb_style) +
            ", does not agree with subsequent topologies.", "AtomGraphSynthesis",
            "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must have a consistent implicit solvent model setting.  The first, "
             "set to " + getImplicitSolventModelName(gb_style) +
             ", will be applied to all systems.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (diel_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must have a consistent dielectric constant, but values conflicting "
            "with the primary system's " + realToString(dielectric_constant) + " were found.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must have a consistent dielectric constant.  A value of " +
             realToString(dielectric_constant) + ", taken from the first topology, will be "
             "applied throughout.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (salt_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All systems must use the same salt concentration, but values conflicting with the "
            "primary system's " + realToString(salt_concentration) + " were found.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All systems must use the same salt concentration.  A value of " +
             realToString(salt_concentration) + ", found in the first topology, will be applied "
             "throughout.", "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (coul_problem) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("All topologies must make use of the same definitions of physical constants.  "
            "Coulomb's constant is defined as " + realToString(coulomb_constant, 6) +
            " in the first topology but differently in subsequent topologies.",
            "AtomGraphSynthesis", "checkCommonSettings");
    case ExceptionResponse::WARN:
      rtWarn("All topologies must make use of the same definitions of physical constants.  "
             "Coulomb's constant is defined as " + realToString(coulomb_constant, 6) +
             " in the first topology and this definition will be applied throughout.",
             "AtomGraphSynthesis", "checkCommonSettings");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::buildAtomAndTermArrays(const std::vector<int> &topology_indices_in,
                                                const std::vector<int> &topology_index_rebase,
                                                const std::vector<int> &restraint_indices_in,
                                                const std::vector<int> &restraint_index_rebase) {

  // Allocate memory and set POINTER-kind arrays for the small packets of data
  const int padded_system_count = roundUp(system_count, warp_size_int);
  int_system_data.resize(43 * padded_system_count);
  int pivot = 0;
  topology_indices.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  restraint_indices.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  residue_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  molecule_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  largest_residue_sizes.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_residues.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  last_solute_atoms.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  first_solvent_molecules.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  ubrd_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cimp_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cmap_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_term_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  virtual_site_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  posn_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_restraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_type_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  total_exclusion_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  rigid_water_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_constraint_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  degrees_of_freedom.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  nonrigid_particle_counts.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  atom_bit_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  residue_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  molecule_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  ubrd_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cimp_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  cmap_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_term_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  virtual_site_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  posn_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  bond_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  angl_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  dihe_restraint_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  nb_exclusion_offsets.setPointer(&int_system_data, pivot, system_count);
  pivot += padded_system_count;
  lennard_jones_abc_offsets.setPointer(&int_system_data, pivot, system_count);
  
  // Load the topology indexing first
  for (int i = 0; i < system_count; i++) {
    topology_indices.putHost(topology_index_rebase[topology_indices_in[i]], i);
    restraint_indices.putHost(restraint_index_rebase[restraint_indices_in[i]], i);
  }

  // Loop over all systems, fill in the above details, and compute the sizes of various arrays
  // for topological valence terms, virtual sites, exclusions, and restraints
  int atom_offset = 0;
  int abit_offset = 0;
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
  int rposn_offset = 0;
  int rbond_offset = 0;
  int rangl_offset = 0;
  int rdihe_offset = 0;
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(i)];
    const RestraintApparatus* ra_ptr     = restraint_networks[restraint_indices.readHost(i)];
    const ChemicalDetailsKit cdk         = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk       = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk          = ag_ptr->getDoublePrecisionValenceKit();
    const int nvsite = ag_ptr->getVirtualSiteCount();
    atom_counts.putHost(cdk.natom, i);
    total_atoms += cdk.natom;
    virtual_site_counts.putHost(nvsite, i);
    total_virtual_sites += nvsite;
    residue_counts.putHost(cdk.nres, i);
    molecule_counts.putHost(cdk.nmol, i);
    largest_residue_sizes.putHost(ag_ptr->getLargestResidueSize(), i);
    last_solute_residues.putHost(ag_ptr->getLastSoluteResidue(), i);
    last_solute_atoms.putHost(ag_ptr->getLastSoluteAtom(), i);
    first_solvent_molecules.putHost(ag_ptr->getFirstSolventMolecule(), i);
    ubrd_term_counts.putHost(vk.nubrd, i);
    cimp_term_counts.putHost(vk.ncimp, i);
    cmap_term_counts.putHost(vk.ncmap, i);
    bond_term_counts.putHost(vk.nbond, i);
    angl_term_counts.putHost(vk.nangl, i);
    dihe_term_counts.putHost(vk.ndihe, i);
    total_ubrd_terms += vk.nubrd;
    total_cimp_terms += vk.ncimp;
    total_cmap_terms += vk.ncmap;
    total_bond_terms += vk.nbond;
    total_angl_terms += vk.nangl;
    total_dihe_terms += vk.ndihe;
    atom_type_counts.putHost(nbk.n_lj_types, i);
    total_exclusion_counts.putHost(ag_ptr->getTotalExclusions(), i);
    rigid_water_counts.putHost(ag_ptr->getRigidWaterCount(), i);
    bond_constraint_counts.putHost(ag_ptr->getBondConstraintCount(), i);
    degrees_of_freedom.putHost(ag_ptr->getDegreesOfFreedom(), i);
    nonrigid_particle_counts.putHost(ag_ptr->getNonrigidParticleCount(), i);

    // Record various offsets, and increment the counters
    atom_offsets.putHost(atom_offset, i);
    atom_offset += roundUp(cdk.natom, warp_size_int);
    const int uint_bits = sizeof(uint) * 8;
    atom_bit_offsets.putHost(abit_offset, i);
    abit_offset += roundUp((cdk.natom + uint_bits - 1) / uint_bits, warp_size_int);
    virtual_site_offsets.putHost(vste_offset, i);
    vste_offset += roundUp(nvsite, warp_size_int);
    residue_offsets.putHost(resi_offset, i);
    resi_offset += roundUp(cdk.nres + 1, warp_size_int);
    molecule_offsets.putHost(mole_offset, i);
    mole_offset += roundUp(cdk.nmol + 1, warp_size_int);
    ubrd_term_offsets.putHost(ubrd_offset, i);
    ubrd_offset += roundUp(vk.nubrd, warp_size_int);
    cimp_term_offsets.putHost(cimp_offset, i);
    cimp_offset += roundUp(vk.ncimp, warp_size_int);
    cmap_term_offsets.putHost(cmap_offset, i);
    cmap_offset += roundUp(vk.ncmap, warp_size_int);
    bond_term_offsets.putHost(bond_offset, i);
    bond_offset += roundUp(vk.nbond, warp_size_int);
    angl_term_offsets.putHost(angl_offset, i);
    angl_offset += roundUp(vk.nangl, warp_size_int);
    dihe_term_offsets.putHost(dihe_offset, i);
    dihe_offset += roundUp(vk.ndihe, warp_size_int);
    nb_exclusion_offsets.putHost(excl_offset, i);
    excl_offset += roundUp(ag_ptr->getTotalExclusions(), warp_size_int);
    if (ra_ptr != nullptr) {
      const RestraintApparatusDpReader rar = ra_ptr->dpData();
      posn_restraint_counts.putHost(rar.nposn, i);
      bond_restraint_counts.putHost(rar.nbond, i);
      angl_restraint_counts.putHost(rar.nangl, i);
      dihe_restraint_counts.putHost(rar.ndihe, i);
      total_position_restraints += rar.nposn;
      total_distance_restraints += rar.nbond;
      total_angle_restraints    += rar.nangl;
      total_dihedral_restraints += rar.ndihe;
      rposn_offset += roundUp(rar.nposn, warp_size_int);
      posn_restraint_offsets.putHost(rposn_offset, i);
      rbond_offset += roundUp(rar.nbond, warp_size_int);
      bond_restraint_offsets.putHost(rbond_offset, i);
      rangl_offset += roundUp(rar.nangl, warp_size_int);
      angl_restraint_offsets.putHost(rangl_offset, i);
      rdihe_offset += roundUp(rar.ndihe, warp_size_int);
      dihe_restraint_offsets.putHost(rdihe_offset, i);
    }
  }

  // Allocate detailed arrays for each descriptor, then collate all topologies.  This
  // fills the "atom and residue details" arrays listed in atomgraph_synthesis.h.
  chem_int_data.resize(resi_offset + mole_offset + (8 * atom_offset));
  chem_double_data.resize(3 * atom_offset);
  chem_float_data.resize(3 * atom_offset);
  chem_char4_data.resize(resi_offset + (2 * atom_offset));
  pivot = 0;
  residue_limits.setPointer(&chem_int_data, pivot, resi_offset);
  pivot += resi_offset;
  atom_struc_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  residue_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  molecule_limits.setPointer(&chem_int_data, pivot, mole_offset);
  pivot += mole_offset;
  atomic_numbers.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  mobile_atoms.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  molecule_membership.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  molecule_contents.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  charge_indices.setPointer(&chem_int_data, pivot, atom_offset);
  pivot += atom_offset;
  lennard_jones_indices.setPointer(&chem_int_data, pivot, atom_offset);
  pivot = 0;
  atomic_charges.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  atomic_masses.setPointer(&chem_double_data, pivot, atom_offset);
  pivot += atom_offset;
  inverse_atomic_masses.setPointer(&chem_double_data, pivot, atom_offset);
  pivot = 0;
  sp_atomic_charges.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_atomic_masses.setPointer(&chem_float_data, pivot, atom_offset);
  pivot += atom_offset;
  sp_inverse_atomic_masses.setPointer(&chem_float_data, pivot, atom_offset);
  pivot = 0;
  atom_names.setPointer(&chem_char4_data, pivot, atom_offset);
  pivot += atom_offset;
  atom_types.setPointer(&chem_char4_data, pivot, atom_offset);
  pivot += atom_offset;
  residue_names.setPointer(&chem_char4_data, pivot, resi_offset);

  // Fill the above atom and residue descriptor arrays
  int* residue_limits_ptr      = residue_limits.data();
  int* atom_struc_numbers_ptr  = atom_struc_numbers.data();
  int* residue_numbers_ptr     = residue_numbers.data();
  int* molecule_limits_ptr     = molecule_limits.data();
  int* atomic_numbers_ptr      = atomic_numbers.data();
  int* molecule_membership_ptr = molecule_membership.data();
  int* molecule_contents_ptr   = molecule_contents.data();
  double* atomic_charges_ptr   = atomic_charges.data();
  float* sp_atomic_charges_ptr = sp_atomic_charges.data();
  char4* atom_names_ptr        = atom_names.data();
  char4* atom_types_ptr        = atom_types.data();
  char4* residue_names_ptr     = residue_names.data();
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(i)];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const NonbondedKit<float> nbk_sp = ag_ptr->getSinglePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const int synth_bit_base = atom_bit_offsets.readHost(i);
    const int synth_atom_base = atom_offsets.readHost(i);
    const int synth_residue_base = residue_offsets.readHost(i);
    const int synth_molecule_base = molecule_offsets.readHost(i);
    for (int j = 0; j < cdk.nres + 1; j++) {
      residue_limits_ptr[synth_residue_base + j] = cdk.res_limits[j] + synth_atom_base;
    }
    for (int j = 0; j < cdk.natom; j++) {
      atom_struc_numbers_ptr[synth_atom_base + j] = cdk.atom_numbers[j];
      residue_numbers_ptr[synth_residue_base + j] = cdk.res_numbers[j];
      atomic_numbers_ptr[synth_atom_base + j] = cdk.z_numbers[j];
      molecule_membership_ptr[synth_atom_base + j] = cdk.mol_home[j] + synth_molecule_base;
      molecule_contents_ptr[synth_atom_base + j] = cdk.mol_contents[j] + synth_atom_base;
      atomic_charges_ptr[synth_atom_base + j] = nbk.charge[j];
      sp_atomic_charges_ptr[synth_atom_base + j] = nbk_sp.charge[j];
      atom_names_ptr[synth_atom_base + j].x = cdk.atom_names[j].x;
      atom_names_ptr[synth_atom_base + j].y = cdk.atom_names[j].y;
      atom_names_ptr[synth_atom_base + j].z = cdk.atom_names[j].z;
      atom_names_ptr[synth_atom_base + j].w = cdk.atom_names[j].w;
    }
    for (int j = 0; j < cdk.nres; j++) {
      residue_names_ptr[synth_residue_base + j].x = cdk.res_names[j].x;
      residue_names_ptr[synth_residue_base + j].y = cdk.res_names[j].y;
      residue_names_ptr[synth_residue_base + j].z = cdk.res_names[j].z;
      residue_names_ptr[synth_residue_base + j].w = cdk.res_names[j].w;
    }
    for (int j = 0; j < cdk.nmol + 1; j++) {
      molecule_limits_ptr[synth_molecule_base + j] = cdk.mol_limits[j] + synth_atom_base;
    }
    const std::vector<bool> atom_mobility = ag_ptr->getAtomMobility();
    const std::vector<int> iatom_mobility(atom_mobility.begin(), atom_mobility.end());
    mobile_atoms.putHost(iatom_mobility, synth_atom_base, cdk.natom);
    atomic_masses.putHost(ag_ptr->getAtomicMass<double>(MassForm::ORDINARY), synth_atom_base,
                          cdk.natom);
    inverse_atomic_masses.putHost(ag_ptr->getAtomicMass<double>(MassForm::INVERSE),
                                  synth_atom_base, cdk.natom);
    sp_atomic_masses.putHost(ag_ptr->getAtomicMass<float>(MassForm::ORDINARY), synth_atom_base,
                             cdk.natom);
    sp_inverse_atomic_masses.putHost(ag_ptr->getAtomicMass<float>(MassForm::INVERSE),
                                     synth_atom_base, cdk.natom);
    atom_types.putHost(ag_ptr->getAtomType(), synth_atom_base, cdk.natom);
  }

  // Fill in the valence term atom indexing arrays.
  valence_int_data.resize((3 * ubrd_offset) + (5 * cimp_offset) + (6 * cmap_offset) +
                          (3 * bond_offset) + (4 * angl_offset) + (5 * dihe_offset));
  pivot = 0;
  ubrd_i_atoms.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  ubrd_k_atoms.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  ubrd_param_idx.setPointer(&valence_int_data, pivot, ubrd_offset);
  pivot += ubrd_offset;
  cimp_i_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_j_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_k_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_l_atoms.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cimp_param_idx.setPointer(&valence_int_data, pivot, cimp_offset);
  pivot += cimp_offset;
  cmap_i_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_j_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_k_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_l_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_m_atoms.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  cmap_param_idx.setPointer(&valence_int_data, pivot, cmap_offset);
  pivot += cmap_offset;
  bond_i_atoms.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  bond_j_atoms.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  bond_param_idx.setPointer(&valence_int_data, pivot, bond_offset);
  pivot += bond_offset;
  angl_i_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_j_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_k_atoms.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  angl_param_idx.setPointer(&valence_int_data, pivot, angl_offset);
  pivot += angl_offset;
  dihe_i_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_j_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_k_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_l_atoms.setPointer(&valence_int_data, pivot, dihe_offset);
  pivot += dihe_offset;
  dihe_param_idx.setPointer(&valence_int_data, pivot, dihe_offset);
  for (int sysid = 0; sysid < system_count; sysid++) {
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(sysid)];
    const ValenceKit<double> vk = ag_ptr->getDoublePrecisionValenceKit();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_ubrd_offset = ubrd_term_offsets.readHost(sysid);
    const int synth_cimp_offset = cimp_term_offsets.readHost(sysid);
    const int synth_cmap_offset = cmap_term_offsets.readHost(sysid);
    const int synth_bond_offset = bond_term_offsets.readHost(sysid);
    const int synth_angl_offset = angl_term_offsets.readHost(sysid);
    const int synth_dihe_offset = dihe_term_offsets.readHost(sysid);
    for (int pos = 0; pos < vk.nubrd; pos++) {
      ubrd_i_atoms.putHost(vk.ubrd_i_atoms[pos] + synth_atom_base, synth_ubrd_offset + pos);
      ubrd_k_atoms.putHost(vk.ubrd_k_atoms[pos] + synth_atom_base, synth_ubrd_offset + pos);
    }
    for (int pos = 0; pos < vk.ncimp; pos++) {
      cimp_i_atoms.putHost(vk.cimp_i_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_j_atoms.putHost(vk.cimp_j_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_k_atoms.putHost(vk.cimp_k_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
      cimp_l_atoms.putHost(vk.cimp_l_atoms[pos] + synth_atom_base, synth_cimp_offset + pos);
    }
    for (int pos = 0; pos < vk.ncmap; pos++) {
      cmap_i_atoms.putHost(vk.cmap_i_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_j_atoms.putHost(vk.cmap_j_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_k_atoms.putHost(vk.cmap_k_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_l_atoms.putHost(vk.cmap_l_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
      cmap_m_atoms.putHost(vk.cmap_m_atoms[pos] + synth_atom_base, synth_cmap_offset + pos);
    }
    for (int pos = 0; pos < vk.nbond; pos++) {
      bond_i_atoms.putHost(vk.bond_i_atoms[pos] + synth_atom_base, synth_bond_offset + pos);
      bond_j_atoms.putHost(vk.bond_j_atoms[pos] + synth_atom_base, synth_bond_offset + pos);
    }
    for (int pos = 0; pos < vk.nangl; pos++) {
      angl_i_atoms.putHost(vk.angl_i_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
      angl_j_atoms.putHost(vk.angl_j_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
      angl_k_atoms.putHost(vk.angl_k_atoms[pos] + synth_atom_base, synth_angl_offset + pos);
    }
    for (int pos = 0; pos < vk.ndihe; pos++) {
      dihe_i_atoms.putHost(vk.dihe_i_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_j_atoms.putHost(vk.dihe_j_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_k_atoms.putHost(vk.dihe_k_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
      dihe_l_atoms.putHost(vk.dihe_l_atoms[pos] + synth_atom_base, synth_dihe_offset + pos);
    }
  }

  // Fill in the restraint term indexing arrays
  nmr_int_data.resize((3 * rposn_offset) + (3 * rbond_offset) + (4 * rangl_offset) +
                      (5 * rdihe_offset));
  pivot = 0;
  rposn_atoms.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rposn_kr_param_idx.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rposn_xyz_param_idx.setPointer(&nmr_int_data, pivot, rposn_offset);
  pivot += rposn_offset;
  rbond_i_atoms.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rbond_j_atoms.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rbond_param_idx.setPointer(&nmr_int_data, pivot, rbond_offset);
  pivot += rbond_offset;
  rangl_i_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_j_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_k_atoms.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rangl_param_idx.setPointer(&nmr_int_data, pivot, rangl_offset);
  pivot += rangl_offset;
  rdihe_i_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_j_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_k_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_l_atoms.setPointer(&nmr_int_data, pivot, rdihe_offset);
  pivot += rdihe_offset;
  rdihe_param_idx.setPointer(&nmr_int_data, pivot, rdihe_offset);
  for (int sysid = 0; sysid < system_count; sysid++) {
    const RestraintApparatus* ra_ptr = restraint_networks[restraint_indices.readHost(sysid)];
    if (ra_ptr == nullptr) {
      continue;
    }
    const RestraintApparatusDpReader rar = ra_ptr->dpData();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_rposn_offset = posn_restraint_offsets.readHost(sysid);
    const int synth_rbond_offset = bond_restraint_offsets.readHost(sysid);
    const int synth_rangl_offset = angl_restraint_offsets.readHost(sysid);
    const int synth_rdihe_offset = dihe_restraint_offsets.readHost(sysid);
    for (int pos = 0; pos < rar.nposn; pos++) {
      rposn_atoms.putHost(rar.rposn_atoms[pos] + synth_atom_base, synth_rposn_offset + pos);
    }
    for (int pos = 0; pos < rar.nbond; pos++) {
      rbond_i_atoms.putHost(rar.rbond_i_atoms[pos] + synth_atom_base, synth_rbond_offset + pos);
      rbond_j_atoms.putHost(rar.rbond_j_atoms[pos] + synth_atom_base, synth_rbond_offset + pos);
    }
    for (int pos = 0; pos < rar.nangl; pos++) {
      rangl_i_atoms.putHost(rar.rangl_i_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
      rangl_j_atoms.putHost(rar.rangl_j_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
      rangl_k_atoms.putHost(rar.rangl_k_atoms[pos] + synth_atom_base, synth_rangl_offset + pos);
    }
    for (int pos = 0; pos < rar.ndihe; pos++) {
      rdihe_i_atoms.putHost(rar.rdihe_i_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_j_atoms.putHost(rar.rdihe_j_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_k_atoms.putHost(rar.rdihe_k_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
      rdihe_l_atoms.putHost(rar.rdihe_l_atoms[pos] + synth_atom_base, synth_rdihe_offset + pos);
    }
  }

  // Fill in the virtual site indexing arrays
  vsite_int_data.resize(6 * vste_offset);
  pivot = 0;
  virtual_site_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame1_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame2_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame3_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_frame4_atoms.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  virtual_site_parameter_indices.setPointer(&vsite_int_data, pivot, vste_offset);
  pivot += vste_offset;
  for (int sysid = 0; sysid < system_count; sysid++) {
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(sysid)];
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    const int synth_atom_base = atom_offsets.readHost(sysid);
    const int synth_vste_offset = virtual_site_offsets.readHost(sysid);
    for (int pos = 0; pos < vsk.nsite; pos++) {
      virtual_site_atoms.putHost(vsk.vs_atoms[pos], synth_vste_offset + pos);
      virtual_site_frame1_atoms.putHost(vsk.frame1_idx[pos], synth_vste_offset + pos);
      virtual_site_frame2_atoms.putHost(vsk.frame2_idx[pos], synth_vste_offset + pos);
      virtual_site_frame3_atoms.putHost(vsk.frame3_idx[pos], synth_vste_offset + pos);
      virtual_site_frame4_atoms.putHost(vsk.frame4_idx[pos], synth_vste_offset + pos);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::condenseParameterTables() {

  // Compute the numbers of unique parameters.  Take the opportunity to compute offsets (starting
  // bounds) for various sets of terms.  While similar to arrays like bond_term_offsets already
  // stored in the AtomGraphSynthesis, these arrays pertain to term offsets for a list of unique
  // topologies, not a list of all systems (which could include multiple coordinate sets
  // referencing a single topology).
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  int ubrd_offset = 0;
  int cimp_offset = 0;
  int cmap_offset = 0;
  int chrg_offset = 0;
  int vste_offset = 0;
  std::vector<int> topology_bond_table_offsets(topology_count);
  std::vector<int> topology_angl_table_offsets(topology_count);
  std::vector<int> topology_dihe_table_offsets(topology_count);
  std::vector<int> topology_ubrd_table_offsets(topology_count);
  std::vector<int> topology_cimp_table_offsets(topology_count);
  std::vector<int> topology_cmap_table_offsets(topology_count);
  std::vector<int> topology_chrg_table_offsets(topology_count);
  std::vector<int> topology_vste_table_offsets(topology_count);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* ag_ptr = topologies[i];
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    topology_bond_table_offsets[i] = bond_offset;
    topology_angl_table_offsets[i] = angl_offset;
    topology_dihe_table_offsets[i] = dihe_offset;
    topology_ubrd_table_offsets[i] = ubrd_offset;
    topology_cimp_table_offsets[i] = cimp_offset;
    topology_cmap_table_offsets[i] = cmap_offset;
    topology_chrg_table_offsets[i] = chrg_offset;
    topology_vste_table_offsets[i] = vste_offset;
    bond_offset += vk.nbond;
    angl_offset += vk.nangl;
    dihe_offset += vk.ndihe;
    ubrd_offset += vk.nubrd;
    cimp_offset += vk.ncimp;
    cmap_offset += vk.ncmap;
    chrg_offset += nbk.natom;
    vste_offset += vsk.nsite;
  }

  // Pre-compute some quantities relating to CMAPs that will help distinguish these gargantuan
  // "parameters" in the inner loops that follow.
  std::vector<double> cmap_parameter_sums(cmap_offset);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    const ValenceKit<double> i_vk = iag_ptr->getDoublePrecisionValenceKit();
    for (int j = 0; j < i_vk.ncmap_surf; j++) {
      cmap_parameter_sums[topology_cmap_table_offsets[i] + j] =
        sum<double>(&i_vk.cmap_surf[i_vk.cmap_surf_bounds[j]],
                    i_vk.cmap_dim[j] * i_vk.cmap_dim[j]);
    }
  }
  
  // Create lists of unique parameters for the valence and non-bonded calculations.
  std::vector<int> bond_synthesis_index(bond_offset, -1);
  std::vector<int> angl_synthesis_index(angl_offset, -1);
  std::vector<int> dihe_synthesis_index(dihe_offset, -1);
  std::vector<int> ubrd_synthesis_index(ubrd_offset, -1);
  std::vector<int> cimp_synthesis_index(cimp_offset, -1);
  std::vector<int> cmap_synthesis_index(cmap_offset, -1);
  std::vector<int> chrg_synthesis_index(chrg_offset, -1);
  std::vector<int> vste_synthesis_index(vste_offset, -1);
  std::vector<double> filtered_bond_keq;
  std::vector<double> filtered_bond_leq;
  std::vector<float> sp_filtered_bond_keq;
  std::vector<float> sp_filtered_bond_leq;
  std::vector<double> filtered_angl_keq;
  std::vector<double> filtered_angl_theta;
  std::vector<float> sp_filtered_angl_keq;
  std::vector<float> sp_filtered_angl_theta;
  std::vector<double> filtered_dihe_amp;
  std::vector<double> filtered_dihe_freq;
  std::vector<double> filtered_dihe_phi;
  std::vector<float> sp_filtered_dihe_amp;
  std::vector<float> sp_filtered_dihe_freq;
  std::vector<float> sp_filtered_dihe_phi;
  std::vector<double> filtered_ubrd_keq;
  std::vector<double> filtered_ubrd_leq;
  std::vector<float> sp_filtered_ubrd_keq;
  std::vector<float> sp_filtered_ubrd_leq;
  std::vector<double> filtered_cimp_keq;
  std::vector<double> filtered_cimp_phi;
  std::vector<float> sp_filtered_cimp_keq;
  std::vector<float> sp_filtered_cimp_phi;
  std::vector<int> filtered_cmap_dim;
  std::vector<int> filtered_cmap_surf_bounds(1, 0);
  std::vector<double> filtered_cmap_surf;
  std::vector<float> sp_filtered_cmap_surf;
  std::vector<double> filtered_chrg;
  std::vector<float> sp_filtered_chrg;
  std::vector<double4> filtered_vste_params;
  std::vector<float4> sp_filtered_vste_params;
  int n_unique_bond = 0;
  int n_unique_angl = 0;
  int n_unique_dihe = 0;
  int n_unique_ubrd = 0;
  int n_unique_cimp = 0;
  int n_unique_cmap = 0;
  int n_unique_chrg = 0;
  int n_unique_vste = 0;
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    const ChemicalDetailsKit i_cdk       = iag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> i_nbk     = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk        = iag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> i_vsk   = iag_ptr->getDoublePrecisionVirtualSiteKit();
    const NonbondedKit<float> i_nbk_sp   = iag_ptr->getSinglePrecisionNonbondedKit();
    const ValenceKit<float> i_vk_sp      = iag_ptr->getSinglePrecisionValenceKit();
    const VirtualSiteKit<float> i_vsk_sp = iag_ptr->getSinglePrecisionVirtualSiteKit();

    // Seek out unique bond parameters
    for (int j = 0; j < i_vk.nbond_param; j++) {
      if (bond_synthesis_index[topology_bond_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_bond_keq(i_vk.bond_keq[j], constants::verytiny);
      const Approx ij_bond_leq(i_vk.bond_leq[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nbond_param; m++) {
          if (bond_synthesis_index[topology_bond_table_offsets[k] + m] < 0 &&
              ij_bond_keq.test(k_vk.bond_keq[m]) && ij_bond_leq.test(k_vk.bond_leq[m])) {
            bond_synthesis_index[topology_bond_table_offsets[k] + m] = n_unique_bond;
          }
        }
      }

      // Catalog this unique bond and increment the counter
      filtered_bond_keq.push_back(i_vk.bond_keq[j]);
      filtered_bond_leq.push_back(i_vk.bond_leq[j]);
      sp_filtered_bond_keq.push_back(i_vk_sp.bond_keq[j]);
      sp_filtered_bond_leq.push_back(i_vk_sp.bond_leq[j]);
      n_unique_bond++;
    }

    // Seek out unique angle parameters
    for (int j = 0; j < i_vk.nangl_param; j++) {
      if (angl_synthesis_index[topology_angl_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_angl_keq(i_vk.angl_keq[j], constants::verytiny);
      const Approx ij_angl_theta(i_vk.angl_theta[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nangl_param; m++) {
          if (angl_synthesis_index[topology_angl_table_offsets[k] + m] < 0 &&
              ij_angl_keq.test(k_vk.angl_keq[m]) && ij_angl_theta.test(k_vk.angl_theta[m])) {
            angl_synthesis_index[topology_angl_table_offsets[k] + m] = n_unique_angl;
          }
        }
      }

      // Catalog this unique angle and increment the counter
      filtered_angl_keq.push_back(i_vk.angl_keq[j]);
      filtered_angl_theta.push_back(i_vk.angl_theta[j]);
      sp_filtered_angl_keq.push_back(i_vk_sp.angl_keq[j]);
      sp_filtered_angl_theta.push_back(i_vk_sp.angl_theta[j]);
      n_unique_angl++;
    }

    // Seek out unique dihedral parameters
    for (int j = 0; j < i_vk.ndihe_param; j++) {
      if (dihe_synthesis_index[topology_dihe_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_dihe_amp(i_vk.dihe_amp[j], constants::verytiny);
      const Approx ij_dihe_freq(i_vk.dihe_freq[j], constants::verytiny);
      const Approx ij_dihe_phi(i_vk.dihe_phi[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ndihe_param; m++) {
          if (dihe_synthesis_index[topology_dihe_table_offsets[k] + m] < 0 &&
              ij_dihe_amp.test(k_vk.dihe_amp[m]) && ij_dihe_freq.test(k_vk.dihe_freq[m]) &&
              ij_dihe_phi.test(k_vk.dihe_phi[m])) {
            dihe_synthesis_index[topology_dihe_table_offsets[k] + m] = n_unique_dihe;
          }
        }
      }

      // Catalog this unique dihedral and increment the counter
      filtered_dihe_amp.push_back(i_vk.dihe_amp[j]);
      filtered_dihe_freq.push_back(i_vk.dihe_freq[j]);
      filtered_dihe_phi.push_back(i_vk.dihe_phi[j]);
      sp_filtered_dihe_amp.push_back(i_vk_sp.dihe_amp[j]);
      sp_filtered_dihe_freq.push_back(i_vk_sp.dihe_freq[j]);
      sp_filtered_dihe_phi.push_back(i_vk_sp.dihe_phi[j]);
      n_unique_dihe++;
    }
    
    // Seek out unique Urey-Bradley parameters
    for (int j = 0; j < i_vk.nubrd_param; j++) {
      if (ubrd_synthesis_index[topology_ubrd_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_ubrd_keq(i_vk.ubrd_keq[j], constants::verytiny);
      const Approx ij_ubrd_leq(i_vk.ubrd_leq[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.nubrd_param; m++) {
          if (ubrd_synthesis_index[topology_ubrd_table_offsets[k] + m] < 0 &&
              ij_ubrd_keq.test(k_vk.ubrd_keq[m]) && ij_ubrd_leq.test(k_vk.ubrd_leq[m])) {
            ubrd_synthesis_index[topology_ubrd_table_offsets[k] + m] = n_unique_ubrd;
          }
        }
      }

      // Catalog this unique Urey-Bradley parameter set and increment the counter
      filtered_ubrd_keq.push_back(i_vk.ubrd_keq[j]);
      filtered_ubrd_leq.push_back(i_vk.ubrd_leq[j]);
      sp_filtered_ubrd_keq.push_back(i_vk_sp.ubrd_keq[j]);
      sp_filtered_ubrd_leq.push_back(i_vk_sp.ubrd_leq[j]);
      n_unique_ubrd++;
    }

    // Seek out unique CHARMM improper dihedral parameters
    for (int j = 0; j < i_vk.ncimp_param; j++) {
      if (cimp_synthesis_index[topology_cimp_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_cimp_keq(i_vk.cimp_keq[j], constants::verytiny);
      const Approx ij_cimp_phi(i_vk.cimp_phi[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk   = kag_ptr->getDoublePrecisionValenceKit();
        const ValenceKit<float> k_vk_sp = kag_ptr->getSinglePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ncimp_param; m++) {
          if (cimp_synthesis_index[topology_cimp_table_offsets[k] + m] < 0 &&
              ij_cimp_keq.test(k_vk.cimp_keq[m]) && ij_cimp_phi.test(k_vk.cimp_phi[m])) {
            cimp_synthesis_index[topology_cimp_table_offsets[k] + m] = n_unique_cimp;
          }
        }
      }

      // Catalog this unique CHARMM improper dihedral and increment the counter
      filtered_cimp_keq.push_back(i_vk.cimp_keq[j]);
      filtered_cimp_phi.push_back(i_vk.cimp_phi[j]);
      sp_filtered_cimp_keq.push_back(i_vk_sp.cimp_keq[j]);
      sp_filtered_cimp_phi.push_back(i_vk_sp.cimp_phi[j]);
      n_unique_cimp++;
    }

    // Seek out unique CMAP surfaces
    for (int j = 0; j < i_vk.ncmap_surf; j++) {
      if (cmap_synthesis_index[topology_cmap_table_offsets[i] + j] >= 0) {
        continue;
      }

      // Comparing CMAP surfaces is a much more involved process than comparing other parameters.
      // Use a pre-arranged table of surface sums to expedite the process.
      const Approx ij_cmap_sum(cmap_parameter_sums[topology_cmap_table_offsets[i] + j],
                               constants::verytiny);
      const int ij_cmap_dim = i_vk.cmap_dim[j];
      const double* ij_surf_ptr = &i_vk.cmap_surf[i_vk.cmap_surf_bounds[j]];
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const ValenceKit<double> k_vk = kag_ptr->getDoublePrecisionValenceKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vk.ncmap_surf; m++) {
          if (cmap_synthesis_index[topology_cmap_table_offsets[k] + m] < 0 &&
              ij_cmap_dim != k_vk.cmap_dim[m] &&
              ij_cmap_sum.test(cmap_parameter_sums[topology_cmap_table_offsets[k] + m]) &&
              maxAbsoluteDifference(ij_surf_ptr, &k_vk.cmap_surf[k_vk.cmap_surf_bounds[m]],
                                    ij_cmap_dim * ij_cmap_dim) < constants::verytiny) {
            cmap_synthesis_index[topology_cmap_table_offsets[k] + m] = n_unique_cmap;
          }
        }
      }

      // Catalog this unique CMAP and increment the counter
      for (int k = 0; k < ij_cmap_dim * ij_cmap_dim; k++) {
        filtered_cmap_surf.push_back(i_vk.cmap_surf[k]);
        sp_filtered_cmap_surf.push_back(i_vk_sp.cmap_surf[k]);
      }
      filtered_cmap_dim.push_back(ij_cmap_dim);
      filtered_cmap_surf_bounds.push_back(filtered_cmap_surf.size());
      n_unique_cmap++;
    }

    // Seek out unique charges
    for (int j = 0; j < i_nbk.n_q_types; j++) {
      if (chrg_synthesis_index[topology_chrg_table_offsets[i] + j] >= 0) {
        continue;
      }
      const Approx ij_chrg(i_nbk.q_parameter[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const NonbondedKit<double> k_nbk   = kag_ptr->getDoublePrecisionNonbondedKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_nbk.n_q_types; m++) {
          if (chrg_synthesis_index[topology_chrg_table_offsets[k] + m] < 0 &&
              ij_chrg.test(k_nbk.q_parameter[m])) {
            chrg_synthesis_index[topology_chrg_table_offsets[k] + m] = n_unique_chrg;
          }
        }
      }

      // Catalog this unique charge and increment the counter
      filtered_chrg.push_back(i_nbk.q_parameter[j]);
      sp_filtered_chrg.push_back(i_nbk_sp.q_parameter[j]);
      n_unique_chrg++;
    }

    // Seek out unique virtual site frames
    for (int j = 0; j < i_vsk.nframe_set; j++) {
      if (vste_synthesis_index[topology_vste_table_offsets[i] + j] >= 0) {
        continue;
      }
      const int ij_frame_type = i_vsk.vs_types[j];
      const Approx ij_dim1(i_vsk.dim1[j], constants::verytiny);
      const Approx ij_dim2(i_vsk.dim2[j], constants::verytiny);
      const Approx ij_dim3(i_vsk.dim3[j], constants::verytiny);
      for (int k = i; k < topology_count; k++) {
        const AtomGraph* kag_ptr = topologies[k];
        const VirtualSiteKit<double> k_vsk    = kag_ptr->getDoublePrecisionVirtualSiteKit();
        const int mstart = (k == i) ? j : 0;
        for (int m = mstart; m < k_vsk.nframe_set; m++) {
          if (vste_synthesis_index[topology_vste_table_offsets[k] + m] < 0 &&
              ij_frame_type == k_vsk.vs_types[m] && ij_dim1.test(k_vsk.dim1[m]) &&
              ij_dim2.test(k_vsk.dim2[m]) && ij_dim3.test(k_vsk.dim3[m])) {
            vste_synthesis_index[topology_vste_table_offsets[k] + m] = n_unique_vste;
          }
        }
      }

      // Catalog this unique virtual site frame and increment the counter
      filtered_vste_params.push_back({ i_vsk.dim1[j], i_vsk.dim2[j],
                                       i_vsk.dim3[j], static_cast<double>(ij_frame_type) });
      sp_filtered_vste_params.push_back({ i_vsk_sp.dim1[j], i_vsk_sp.dim2[j],
                                          i_vsk_sp.dim3[j], static_cast<float>(ij_frame_type) });
      n_unique_vste++;
    }
  }

  // Post-process the condensed array of CMAP surfaces into patches
  const CmapAccessories cmap_digest = ComputeCmapDerivatives(n_unique_cmap, filtered_cmap_dim,
                                                             filtered_cmap_surf_bounds,
                                                             filtered_cmap_surf);
  
  // Record synthesis parameters found thus far.  Pad the relevant ARRAY-kind Hybrid objects in
  // case some developer tries to have all lanes of a warp access an element.
  total_bond_params   = n_unique_bond;
  total_angl_params   = n_unique_angl;
  total_dihe_params   = n_unique_dihe;
  total_ubrd_params   = n_unique_ubrd;
  total_cimp_params   = n_unique_cimp;
  total_cmap_surfaces = n_unique_cmap;
  total_charge_types  = n_unique_chrg;
  total_vste_params   = n_unique_vste;
  const int rn_space = (2 * roundUp(total_bond_params, warp_size_int)) +
                       (2 * roundUp(total_angl_params, warp_size_int)) +
                       (3 * roundUp(total_dihe_params, warp_size_int)) +
                       (2 * roundUp(total_ubrd_params, warp_size_int)) +
                       (2 * roundUp(total_cimp_params, warp_size_int)) +
                       roundUp(static_cast<int>(filtered_cmap_surf.size()), warp_size_int) +
                       roundUp(static_cast<int>(cmap_digest.patch_matrix_form.size()),
                               warp_size_int);
  valparam_double_data.resize(rn_space);
  valparam_float_data.resize(rn_space);
  valparam_int_data.resize(roundUp(total_cmap_surfaces, warp_size_int) +
                           (2 * roundUp(total_cmap_surfaces + 1, warp_size_int)));
  size_t ic = 0LLU;
  ic = bond_stiffnesses.putHost(&valparam_double_data, filtered_bond_keq, ic, warp_size_zu);
  ic = bond_equilibria.putHost(&valparam_double_data, filtered_bond_leq, ic, warp_size_zu);
  ic = angl_stiffnesses.putHost(&valparam_double_data, filtered_angl_keq, ic, warp_size_zu);
  ic = angl_equilibria.putHost(&valparam_double_data, filtered_angl_theta, ic, warp_size_zu);
  ic = dihe_amplitudes.putHost(&valparam_double_data, filtered_dihe_amp, ic, warp_size_zu);
  ic = dihe_periodicities.putHost(&valparam_double_data, filtered_dihe_freq, ic, warp_size_zu);
  ic = dihe_phase_angles.putHost(&valparam_double_data, filtered_dihe_phi, ic, warp_size_zu);
  ic = ubrd_stiffnesses.putHost(&valparam_double_data, filtered_ubrd_keq, ic, warp_size_zu);
  ic = ubrd_equilibria.putHost(&valparam_double_data, filtered_ubrd_leq, ic, warp_size_zu);
  ic = cimp_stiffnesses.putHost(&valparam_double_data, filtered_cimp_keq, ic, warp_size_zu);
  ic = cimp_phase_angles.putHost(&valparam_double_data, filtered_cimp_phi, ic, warp_size_zu);
  ic = cmap_surfaces.putHost(&valparam_double_data, filtered_cmap_surf, ic, warp_size_zu);
  ic = cmap_patches.putHost(&valparam_double_data, cmap_digest.patch_matrix_form, ic,
                            warp_size_zu);
  ic = 0LLU;
  ic = sp_bond_stiffnesses.putHost(&valparam_float_data, sp_filtered_bond_keq, ic, warp_size_zu);
  ic = sp_bond_equilibria.putHost(&valparam_float_data, sp_filtered_bond_leq, ic, warp_size_zu);
  ic = sp_angl_stiffnesses.putHost(&valparam_float_data, sp_filtered_angl_keq, ic, warp_size_zu);
  ic = sp_angl_equilibria.putHost(&valparam_float_data, sp_filtered_angl_theta, ic, warp_size_zu);
  ic = sp_dihe_amplitudes.putHost(&valparam_float_data, sp_filtered_dihe_amp, ic, warp_size_zu);
  ic = sp_dihe_periodicities.putHost(&valparam_float_data, sp_filtered_dihe_freq, ic,
                                     warp_size_zu);
  ic = sp_dihe_phase_angles.putHost(&valparam_float_data, sp_filtered_dihe_phi, ic, warp_size_zu);
  ic = sp_ubrd_stiffnesses.putHost(&valparam_float_data, sp_filtered_ubrd_keq, ic, warp_size_zu);
  ic = sp_ubrd_equilibria.putHost(&valparam_float_data, sp_filtered_ubrd_leq, ic, warp_size_zu);
  ic = sp_cimp_stiffnesses.putHost(&valparam_float_data, sp_filtered_cimp_keq, ic, warp_size_zu);
  ic = sp_cimp_phase_angles.putHost(&valparam_float_data, sp_filtered_cimp_phi, ic, warp_size_zu);
  ic = sp_cmap_surfaces.putHost(&valparam_float_data, sp_filtered_cmap_surf, ic, warp_size_zu);
  ic = sp_cmap_patches.putHost(&valparam_float_data,
                               std::vector<float>(cmap_digest.patch_matrix_form.begin(),
                                                  cmap_digest.patch_matrix_form.end()), ic,
                               warp_size_zu);
  ic = 0LLU;
  ic = cmap_surface_dimensions.putHost(&valparam_int_data, filtered_cmap_dim, ic, warp_size_zu);
  ic = cmap_surface_bounds.putHost(&valparam_int_data, filtered_cmap_surf_bounds, ic,
                                   warp_size_zu);
  ic = cmap_patch_bounds.putHost(&valparam_int_data, cmap_digest.patch_matrix_bounds, ic,
                                 warp_size_zu);

  // Loop back over all systems and copy the known mapping of individual topologies to the
  // synthesis as a whole
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[topology_indices.readHost(i)];
    const NonbondedKit<double> i_nbk   = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk      = iag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> i_vsk = iag_ptr->getDoublePrecisionVirtualSiteKit();
    const int tp_index = topology_indices.readHost(i);
    const int ag_bond_table_offset = topology_bond_table_offsets[tp_index];
    const int ag_angl_table_offset = topology_angl_table_offsets[tp_index];
    const int ag_dihe_table_offset = topology_dihe_table_offsets[tp_index];
    const int ag_ubrd_table_offset = topology_ubrd_table_offsets[tp_index];
    const int ag_cimp_table_offset = topology_cimp_table_offsets[tp_index];
    const int ag_cmap_table_offset = topology_cmap_table_offsets[tp_index];
    const int ag_vste_table_offset = topology_vste_table_offsets[tp_index];
    const int ag_chrg_table_offset = topology_chrg_table_offsets[tp_index];
    const int synth_bond_table_offset = bond_term_offsets.readHost(i);
    const int synth_angl_table_offset = angl_term_offsets.readHost(i);
    const int synth_dihe_table_offset = dihe_term_offsets.readHost(i);
    const int synth_ubrd_table_offset = ubrd_term_offsets.readHost(i);
    const int synth_cimp_table_offset = cimp_term_offsets.readHost(i);
    const int synth_cmap_table_offset = cmap_term_offsets.readHost(i);
    const int synth_vste_table_offset = virtual_site_offsets.readHost(i);
    const int synth_atom_offset = atom_offsets.readHost(i);
    for (int j = 0; j < i_vk.nbond; j++) {
      bond_param_idx.putHost(bond_synthesis_index[ag_bond_table_offset + j],
                             synth_bond_table_offset + j);
    }
    for (int j = 0; j < i_vk.nangl; j++) {
      angl_param_idx.putHost(angl_synthesis_index[ag_angl_table_offset + j],
                             synth_angl_table_offset + j);
    }
    for (int j = 0; j < i_vk.ndihe; j++) {
      dihe_param_idx.putHost(dihe_synthesis_index[ag_dihe_table_offset + j],
                             synth_dihe_table_offset + j);
    }
    for (int j = 0; j < i_vk.nubrd; j++) {
      ubrd_param_idx.putHost(ubrd_synthesis_index[ag_ubrd_table_offset + j],
                             synth_ubrd_table_offset + j);
    }
    for (int j = 0; j < i_vk.ncimp; j++) {
      cimp_param_idx.putHost(cimp_synthesis_index[ag_cimp_table_offset + j],
                             synth_cimp_table_offset + j);
    }
    for (int j = 0; j < i_vk.ncmap; j++) {
      cmap_param_idx.putHost(cmap_synthesis_index[ag_cmap_table_offset + j],
                             synth_cmap_table_offset + j);
    }
    for (int j = 0; j < i_nbk.natom; j++) {
      charge_indices.putHost(chrg_synthesis_index[ag_chrg_table_offset + j],
                             synth_atom_offset + j);
    }
    for (int j = 0; j < i_vsk.nsite; j++) {
      virtual_site_parameter_indices.putHost(vste_synthesis_index[ag_vste_table_offset + j],
                                             synth_vste_table_offset + j);
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::extendLJMatrices() {

  // Identify systems with identical Lennard-Jones tables and have their offsets point to the
  // same table in the synthesis.  This may conserve some L1 cache if everything else is done
  // right.  Track the necessary array size to hold all tables.
  int table_acc = 0;
  std::vector<NonbondedKit<double>> nbkvec;
  nbkvec.reserve(topology_count);
  std::vector<int> table_idx(topology_count, -1);
  for (int i = 0; i < topology_count; i++) {
    nbkvec.push_back(topologies[i]->getDoublePrecisionNonbondedKit());
  }
  int n_unique_tables = 0;
  for (int i = 0; i < topology_count; i++) {
    const int ni_lj_types = nbkvec[i].n_lj_types;
    for (int j = i + 1; j < topology_count; j++) {
      if (ni_lj_types != nbkvec[j].n_lj_types) {
        continue;
      }
      bool diags_identical = true;
      for (int k = 0; k < ni_lj_types; k++) {
        const int dgidx = k * (ni_lj_types + 1);
        diags_identical = (diags_identical &&
                           nbkvec[i].lja_coeff[dgidx] == nbkvec[j].lja_coeff[dgidx] &&
                           nbkvec[i].ljb_coeff[dgidx] == nbkvec[j].ljb_coeff[dgidx] &&
                           nbkvec[i].ljc_coeff[dgidx] == nbkvec[j].ljc_coeff[dgidx] &&
                           nbkvec[i].lja_14_coeff[dgidx] == nbkvec[j].lja_14_coeff[dgidx] &&
                           nbkvec[i].ljb_14_coeff[dgidx] == nbkvec[j].ljb_14_coeff[dgidx] &&
                           nbkvec[i].ljc_14_coeff[dgidx] == nbkvec[j].ljc_14_coeff[dgidx]);
      }
      if (diags_identical == false) {
        continue;
      }
      bool tables_identical = true;
      for (int k = 0; k < ni_lj_types * ni_lj_types; k++) {
        tables_identical = (tables_identical &&
                            nbkvec[i].lja_coeff[k] == nbkvec[j].lja_coeff[k] &&
                            nbkvec[i].ljb_coeff[k] == nbkvec[j].ljb_coeff[k] &&
                            nbkvec[i].ljc_coeff[k] == nbkvec[j].ljc_coeff[k] &&
                            nbkvec[i].lja_14_coeff[k] == nbkvec[j].lja_14_coeff[k] &&
                            nbkvec[i].ljb_14_coeff[k] == nbkvec[j].ljb_14_coeff[k] &&
                            nbkvec[i].ljc_14_coeff[k] == nbkvec[j].ljc_14_coeff[k]);
      }
      if (tables_identical) {
        table_idx[j] = n_unique_tables;      
      }
    }
    table_idx[i] = n_unique_tables;
    n_unique_tables++;
  }

  // Loop back over systems to obtain the necessary allocation base on the unique tables
  int seek_idx = 0;
  int alloc_size = 0;
  std::vector<int> tmp_table_offsets(topology_count, 0);
  for (int i = 0; i < topology_count; i++) {
    for (int j = i; j < topology_count; j++) {
      if (table_idx[i] == seek_idx) {
        tmp_table_offsets[j] = alloc_size;
      }
    }
    if (table_idx[i] == seek_idx) {
      alloc_size += roundUp(nbkvec[i].n_lj_types * nbkvec[i].n_lj_types, warp_size_int);
      seek_idx++;
    }
  }

  // Allocate and fill the Lennard-Jones tables
  lennard_jones_ab_coeff.resize(alloc_size);
  lennard_jones_c_coeff.resize(alloc_size);
  lennard_jones_14_ab_coeff.resize(alloc_size);
  lennard_jones_14_c_coeff.resize(alloc_size);
  sp_lennard_jones_ab_coeff.resize(alloc_size);
  sp_lennard_jones_c_coeff.resize(alloc_size);
  sp_lennard_jones_14_ab_coeff.resize(alloc_size);
  sp_lennard_jones_14_c_coeff.resize(alloc_size);
  double2* ab_ptr = lennard_jones_ab_coeff.data();
  double* c_ptr = lennard_jones_c_coeff.data();
  double2* ab_14_ptr = lennard_jones_14_ab_coeff.data();
  double* c_14_ptr = lennard_jones_14_c_coeff.data();
  float2* sp_ab_ptr = sp_lennard_jones_ab_coeff.data();
  float* sp_c_ptr = sp_lennard_jones_c_coeff.data();
  float2* sp_ab_14_ptr = sp_lennard_jones_14_ab_coeff.data();
  float* sp_c_14_ptr = sp_lennard_jones_14_c_coeff.data();
  seek_idx = 0;
  for (int i = 0; i < topology_count; i++) {
    if (table_idx[i] == seek_idx) {
      const int offset = tmp_table_offsets[i];
      for (int j = 0; j < nbkvec[i].n_lj_types * nbkvec[i].n_lj_types; j++) {
        const int joffset = offset + j;
        ab_ptr[joffset].x = nbkvec[i].lja_coeff[j];
        ab_ptr[joffset].y = nbkvec[i].ljb_coeff[j];
        c_ptr[joffset] = nbkvec[i].ljc_coeff[j];
        ab_14_ptr[joffset].x = nbkvec[i].lja_14_coeff[j];
        ab_14_ptr[joffset].y = nbkvec[i].ljb_14_coeff[j];
        c_14_ptr[joffset] = nbkvec[i].ljc_14_coeff[j];
        sp_ab_ptr[joffset].x = nbkvec[i].lja_coeff[j];
        sp_ab_ptr[joffset].y = nbkvec[i].ljb_coeff[j];
        sp_c_ptr[joffset] = nbkvec[i].ljc_coeff[j];
        sp_ab_14_ptr[joffset].x = nbkvec[i].lja_14_coeff[j];
        sp_ab_14_ptr[joffset].y = nbkvec[i].ljb_14_coeff[j];
        sp_c_14_ptr[joffset] = nbkvec[i].ljc_14_coeff[j];
      }
      seek_idx++;
    }
  }

  // Expand the topology-oriented offsets to become system-oriented offsets
  std::vector<int> tmp_lj_system_offsets(system_count, 0);
  for (int i = 0; i < system_count; i++) {
    tmp_lj_system_offsets[i] = tmp_table_offsets[topology_indices.readHost(i)];
  }
  lennard_jones_abc_offsets.putHost(tmp_table_offsets);
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::mapUniqueRestraintKRSeries(const int order,
                                                   const std::vector<int> &network_table_offsets,
                                                   std::vector<int> *synthesis_index,
                                                   std::vector<int2> *filtered_step_bounds,
                                                   std::vector<double2> *filtered_init_keq,
                                                   std::vector<double2> *filtered_finl_keq,
                                                   std::vector<double4> *filtered_init_r,
                                                   std::vector<double4> *filtered_finl_r) {
  int n_unique_term = 0;
  int* synthesis_index_ptr = synthesis_index->data();
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ira_ptr = restraint_networks[i];
    if (ira_ptr == nullptr) {
      continue;
    }
    const RestraintStage istage = RestraintStage::INITIAL;
    const RestraintStage fstage = RestraintStage::FINAL;
    const int* irstr_init_step    = ira_ptr->getApplicationStepPointer(order, istage);
    const int* irstr_finl_step    = ira_ptr->getApplicationStepPointer(order, fstage);
    const double2* irstr_init_keq = ira_ptr->getHarmonicStiffnessPointer(order, istage);
    const double2* irstr_finl_keq = ira_ptr->getHarmonicStiffnessPointer(order, fstage);
    const double4* irstr_init_r   = ira_ptr->getDisplacementPointer(order, istage);
    const double4* irstr_finl_r   = ira_ptr->getDisplacementPointer(order, fstage);
    int jmax;
    switch (order) {
    case 1:
      jmax = ira_ptr->getPositionalRestraintCount();
      break;
    case 2:
      jmax = ira_ptr->getDistanceRestraintCount();
      break;
    case 3:
      jmax = ira_ptr->getAngleRestraintCount();
      break;
    case 4:
      jmax = ira_ptr->getDihedralRestraintCount();
      break;
    default:
      break;
    }
    
    // Seek out unique positional restraints.  This process differs from filtering unique valence
    // parameters because with the restraints, there is just a list of parameters and atoms to
    // which they apply, whereas with the valence parameters there were lists of atoms to which
    // each term applied, a corresponding index into a parameter table, and then the separate
    // parameter tables.  This process, in contrast, must just loop over all restraints, rather
    // than all entries in the tables, to determine which have unique parameters.  In many cases,
    // all restraints will be unique, but the instructions for doing the restraints will have atoms
    // and a restraint parameter index regardless.
    for (int j = 0; j < jmax; j++) {

      // Time dependence, and the exact steps at which that time dependence starts, can signify
      // that both k / r series as well as x / y / z targets are unique.
      const int ij_init_step = irstr_init_step[j];
      const int ij_finl_step = irstr_finl_step[j];
      const bool ij_time_dep = (ij_finl_step == 0);

      // Skip restraints that have already been determined to have a known k / r series
      if (synthesis_index_ptr[network_table_offsets[i] + j] < 0) {
        const Approx ij_posn_init_k2(irstr_init_keq[j].x, constants::verytiny);
        const Approx ij_posn_init_k3(irstr_init_keq[j].y, constants::verytiny);
        const Approx ij_posn_finl_k2(irstr_finl_keq[j].x, constants::verytiny);
        const Approx ij_posn_finl_k3(irstr_finl_keq[j].y, constants::verytiny);
        const Approx ij_posn_init_r1(irstr_init_r[j].x, constants::verytiny);
        const Approx ij_posn_init_r2(irstr_init_r[j].y, constants::verytiny);
        const Approx ij_posn_init_r3(irstr_init_r[j].z, constants::verytiny);
        const Approx ij_posn_init_r4(irstr_init_r[j].w, constants::verytiny);
        const Approx ij_posn_finl_r1(irstr_finl_r[j].x, constants::verytiny);
        const Approx ij_posn_finl_r2(irstr_finl_r[j].y, constants::verytiny);
        const Approx ij_posn_finl_r3(irstr_finl_r[j].z, constants::verytiny);
        const Approx ij_posn_finl_r4(irstr_finl_r[j].w, constants::verytiny);
        for (int k = i; k < restraint_network_count; k++) {
          const RestraintApparatus* kra_ptr = restraint_networks[k];
          if (kra_ptr == nullptr) {
            continue;
          }

          // Another switch to define constants and pointers based on the order of restraint
          // being handled, this time for the abstract of the kth restraint apparatus
          const int* krstr_init_step    = kra_ptr->getApplicationStepPointer(order, istage);
          const int* krstr_finl_step    = kra_ptr->getApplicationStepPointer(order, fstage);
          const double2* krstr_init_keq = kra_ptr->getHarmonicStiffnessPointer(order, istage);
          const double2* krstr_finl_keq = kra_ptr->getHarmonicStiffnessPointer(order, fstage);
          const double4* krstr_init_r   = kra_ptr->getDisplacementPointer(order, istage);
          const double4* krstr_finl_r   = kra_ptr->getDisplacementPointer(order, fstage);          
          int mmax;
          switch (order) {
          case 1:
            mmax = kra_ptr->getPositionalRestraintCount();
            break;
          case 2:
            mmax = kra_ptr->getDistanceRestraintCount();
            break;
          case 3:
            mmax = kra_ptr->getAngleRestraintCount();
            break;
          case 4:
            mmax = kra_ptr->getDihedralRestraintCount();
            break;
          default:
            break;
          }          
          const int mstart = (k == i) ? j : 0;
          for (int m = mstart; m < mmax; m++) {
            if (synthesis_index_ptr[network_table_offsets[k] + m] >= 0) {
              continue;
            }
            const bool km_time_dep = (krstr_finl_step[m] == 0);
            if (ij_time_dep != km_time_dep) {
              continue;
            }
            if (ij_time_dep) {
              if (krstr_init_step[m] == ij_init_step &&
                  krstr_finl_step[m] == ij_finl_step &&
                  ij_posn_init_k2.test(krstr_init_keq[m].x) &&
                  ij_posn_init_k3.test(krstr_init_keq[m].y) &&
                  ij_posn_finl_k2.test(krstr_finl_keq[m].x) &&
                  ij_posn_finl_k3.test(krstr_finl_keq[m].y) &&
                  ij_posn_init_r1.test(krstr_init_r[m].x) &&
                  ij_posn_init_r2.test(krstr_init_r[m].y) &&
                  ij_posn_init_r3.test(krstr_init_r[m].z) &&
                  ij_posn_init_r4.test(krstr_init_r[m].w) &&
                  ij_posn_finl_r1.test(krstr_finl_r[m].x) &&
                  ij_posn_finl_r2.test(krstr_finl_r[m].y) &&
                  ij_posn_finl_r3.test(krstr_finl_r[m].z) &&
                  ij_posn_finl_r4.test(krstr_finl_r[m].w)) {
                synthesis_index_ptr[network_table_offsets[k] + m] = n_unique_term;
              }
            }
            else {
              if (ij_posn_init_k2.test(krstr_init_keq[m].x) &&
                  ij_posn_init_k3.test(krstr_init_keq[m].y) &&
                  ij_posn_init_r1.test(krstr_init_r[m].x) &&
                  ij_posn_init_r2.test(krstr_init_r[m].y) &&
                  ij_posn_init_r3.test(krstr_init_r[m].z) &&
                  ij_posn_init_r4.test(krstr_init_r[m].w)) {
                synthesis_index_ptr[network_table_offsets[k] + m] = n_unique_term;
              }
            }
          }
        }

        // Catalog this unique positional restraint parameter set (the x, y, and z coordinates
        // are a separate parameter set and will be considered next);
        filtered_step_bounds->push_back({ij_init_step, ij_finl_step});
        filtered_init_keq->push_back(irstr_init_keq[j]);
        filtered_finl_keq->push_back(irstr_finl_keq[j]);
        filtered_init_r->push_back(irstr_init_r[j]);
        filtered_finl_r->push_back(irstr_finl_r[j]);
        n_unique_term++;
      }
    }
  }
  return n_unique_term;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::condenseRestraintNetworks() {

  // As with topology parameters, compute the numbers of unique positional, distance, angle, and
  // dihedral angle restraints.  Take the opportunity to compute offsets (starting bounds) for
  // various sets of terms.  One difference in this implementation is that the restraint apparatus
  // pointer for any particular network may be the null pointer, indicating that there is no
  // restraint system in place for that case (or perhaps any cases).
  int posn_offset = 0;
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  std::vector<int> network_rposn_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rbond_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rangl_table_offsets(restraint_network_count, 0);
  std::vector<int> network_rdihe_table_offsets(restraint_network_count, 0);
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ra_ptr = restraint_networks[i];
    network_rposn_table_offsets[i] = posn_offset;
    network_rbond_table_offsets[i] = bond_offset;
    network_rangl_table_offsets[i] = angl_offset;
    network_rdihe_table_offsets[i] = dihe_offset;
    if (ra_ptr != nullptr) {
      posn_offset += ra_ptr->getPositionalRestraintCount();
      bond_offset += ra_ptr->getDistanceRestraintCount();
      angl_offset += ra_ptr->getAngleRestraintCount();
      dihe_offset += ra_ptr->getDihedralRestraintCount();
    }
  }

  // Create lists of unique parameters for the various types of restraints.
  std::vector<int> rposn_synthesis_kr_index(posn_offset, -1);  
  std::vector<int> rposn_synthesis_xyz_index(posn_offset, -1);  
  std::vector<int> rbond_synthesis_index(bond_offset, -1);
  std::vector<int> rangl_synthesis_index(angl_offset, -1);
  std::vector<int> rdihe_synthesis_index(dihe_offset, -1);
  std::vector<int2> filtered_rposn_step_bounds, filtered_rbond_step_bounds;
  std::vector<int2> filtered_rangl_step_bounds, filtered_rdihe_step_bounds;
  std::vector<double2> filtered_rposn_init_keq, filtered_rposn_finl_keq;
  std::vector<double4> filtered_rposn_init_r, filtered_rposn_finl_r;
  std::vector<double2> filtered_rposn_init_xy, filtered_rposn_finl_xy;
  std::vector<double> filtered_rposn_init_z, filtered_rposn_finl_z;
  std::vector<double2> filtered_rbond_init_keq, filtered_rbond_finl_keq;
  std::vector<double4> filtered_rbond_init_r, filtered_rbond_finl_r;
  std::vector<double2> filtered_rangl_init_keq, filtered_rangl_finl_keq;
  std::vector<double4> filtered_rangl_init_r, filtered_rangl_finl_r;
  std::vector<double2> filtered_rdihe_init_keq, filtered_rdihe_finl_keq;
  std::vector<double4> filtered_rdihe_init_r, filtered_rdihe_finl_r;
  const int n_unique_posn_kr = mapUniqueRestraintKRSeries(1, network_rposn_table_offsets,
                                                          &rposn_synthesis_kr_index,
                                                          &filtered_rposn_step_bounds,
                                                          &filtered_rposn_init_keq,
                                                          &filtered_rposn_finl_keq,
                                                          &filtered_rposn_init_r,
                                                          &filtered_rposn_finl_r);
  const int n_unique_bond    = mapUniqueRestraintKRSeries(1, network_rbond_table_offsets,
                                                          &rbond_synthesis_index,
                                                          &filtered_rbond_step_bounds,
                                                          &filtered_rbond_init_keq,
                                                          &filtered_rbond_finl_keq,
                                                          &filtered_rbond_init_r,
                                                          &filtered_rbond_finl_r);
  const int n_unique_angl    = mapUniqueRestraintKRSeries(1, network_rangl_table_offsets,
                                                          &rangl_synthesis_index,
                                                          &filtered_rangl_step_bounds,
                                                          &filtered_rangl_init_keq,
                                                          &filtered_rangl_finl_keq,
                                                          &filtered_rangl_init_r,
                                                          &filtered_rangl_finl_r);
  const int n_unique_dihe    = mapUniqueRestraintKRSeries(1, network_rdihe_table_offsets,
                                                          &rdihe_synthesis_index,
                                                          &filtered_rdihe_step_bounds,
                                                          &filtered_rdihe_init_keq,
                                                          &filtered_rdihe_finl_keq,
                                                          &filtered_rdihe_init_r,
                                                          &filtered_rdihe_finl_r);

  // Seek out unique positional restraint x, y, and z targets.  This process differs from
  // filtering unique valence parameters because with the restraints, there is just a list of
  // parameters and atoms to which they apply, whereas with the valence parameters there were
  // lists of atoms to which each term applied, a corresponding index into a parameter table,
  // and then the separate parameter tables.  This process, in contrast, must just loop over all
  // restraints, rather than all entries in the tables, to determine which have unique
  // parameters.  In many cases, all restraints will be unique, but the instructions for doing
  // the restraints will have atoms and a restraint parameter index regardless.
  int n_unique_posn_xyz = 0;
  for (int i = 0; i < restraint_network_count; i++) {
    const RestraintApparatus* ira_ptr = restraint_networks[i];
    if (ira_ptr == nullptr) {
      continue;
    }
    const RestraintApparatusDpReader irar_dp = ira_ptr->dpData();
    for (int j = 0; j < irar_dp.nposn; j++) {

      // Skip restraints that have been determined to have unique x / y / z targets.  This
      if (rposn_synthesis_xyz_index[network_rposn_table_offsets[i] + j] < 0) {
        const int ij_init_step = irar_dp.rposn_init_step[j];
        const int ij_finl_step = irar_dp.rposn_finl_step[j];
        const bool ij_time_dep = (ij_finl_step == 0);
        const Approx ij_posn_init_x(irar_dp.rposn_init_xy[j].x, constants::verytiny);
        const Approx ij_posn_init_y(irar_dp.rposn_init_xy[j].y, constants::verytiny);
        const Approx ij_posn_init_z(irar_dp.rposn_init_z[j],    constants::verytiny);
        const Approx ij_posn_finl_x(irar_dp.rposn_finl_xy[j].x, constants::verytiny);
        const Approx ij_posn_finl_y(irar_dp.rposn_finl_xy[j].y, constants::verytiny);
        const Approx ij_posn_finl_z(irar_dp.rposn_finl_z[j],    constants::verytiny);
        for (int k = i; k < restraint_network_count; k++) {
          const RestraintApparatus* kra_ptr = restraint_networks[k];
          if (kra_ptr == nullptr) {
            continue;
          }
          const RestraintApparatusDpReader krar_dp = kra_ptr->dpData();
          const int mstart = (k == i) ? j : 0;
          for (int m = mstart; m < krar_dp.nposn; m++) {
            if (rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] >= 0) {
              continue;
            }
            const bool km_time_dep = (krar_dp.rposn_finl_step[m] == 0);
            if (ij_time_dep != km_time_dep) {
              continue;
            }
            if (ij_time_dep) {
              if (ij_posn_init_x.test(krar_dp.rposn_init_xy[m].x) &&
                  ij_posn_init_y.test(krar_dp.rposn_init_xy[m].y) &&
                  ij_posn_init_z.test(krar_dp.rposn_init_z[m]) &&
                  ij_posn_finl_x.test(krar_dp.rposn_finl_xy[m].x) &&
                  ij_posn_finl_y.test(krar_dp.rposn_finl_xy[m].y) &&
                  ij_posn_finl_z.test(krar_dp.rposn_finl_z[m])) {
                rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] = n_unique_posn_xyz;
              }
            }
            else {
              if (ij_posn_init_x.test(krar_dp.rposn_init_xy[m].x) &&
                  ij_posn_init_y.test(krar_dp.rposn_init_xy[m].y) &&
                  ij_posn_init_z.test(krar_dp.rposn_init_z[m])) {
                rposn_synthesis_xyz_index[network_rposn_table_offsets[k] + m] = n_unique_posn_xyz;
              }
            }
          }
        }

        // Catalog this unique positional restraint target coordinate set
        filtered_rposn_init_xy.push_back({irar_dp.rposn_init_xy[j].x, irar_dp.rposn_init_xy[j].y});
        filtered_rposn_init_z.push_back(irar_dp.rposn_init_z[j]);
        filtered_rposn_finl_xy.push_back({irar_dp.rposn_finl_xy[j].x, irar_dp.rposn_finl_xy[j].y});
        filtered_rposn_finl_z.push_back(irar_dp.rposn_finl_z[j]);
        n_unique_posn_xyz++;
      }
    }
  }

  // Resize the Hybrid arrays, set pointers, and fill Hybrid objects in the AtomGraphSynthesis
  const size_t i2c = roundUp(n_unique_posn_kr, warp_size_int) +
                     roundUp(n_unique_bond, warp_size_int) +
                     roundUp(n_unique_angl, warp_size_int) +
                     roundUp(n_unique_dihe, warp_size_int);
  nmr_int2_data.resize(i2c);
  nmr_double_data.resize(2 * roundUp(n_unique_posn_xyz, warp_size_int));
  nmr_double2_data.resize(2 * (i2c + n_unique_posn_xyz));
  nmr_double4_data.resize(2 * i2c);
  nmr_float_data.resize(2 * roundUp(n_unique_posn_xyz, warp_size_int));
  nmr_float2_data.resize(2 * (i2c + n_unique_posn_xyz));
  nmr_float4_data.resize(2 * i2c);
  size_t ic = 0LLU;
  ic = rposn_step_bounds.putHost(&nmr_int2_data, filtered_rposn_step_bounds, ic, warp_size_zu);
  ic = rbond_step_bounds.putHost(&nmr_int2_data, filtered_rbond_step_bounds, ic, warp_size_zu);
  ic = rangl_step_bounds.putHost(&nmr_int2_data, filtered_rangl_step_bounds, ic, warp_size_zu);
  ic = rdihe_step_bounds.putHost(&nmr_int2_data, filtered_rdihe_step_bounds, ic, warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_z.putHost(&nmr_double_data, filtered_rposn_init_z, ic, warp_size_zu);
  ic = rposn_final_z.putHost(&nmr_double_data, filtered_rposn_finl_z, ic, warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_k.putHost(&nmr_double2_data, filtered_rposn_init_keq, ic, warp_size_zu);
  ic = rposn_final_k.putHost(&nmr_double2_data, filtered_rposn_finl_keq, ic, warp_size_zu);
  ic = rposn_init_xy.putHost(&nmr_double2_data, filtered_rposn_init_xy, ic, warp_size_zu);
  ic = rposn_final_xy.putHost(&nmr_double2_data, filtered_rposn_finl_xy, ic, warp_size_zu);
  ic = rbond_init_k.putHost(&nmr_double2_data, filtered_rbond_init_keq, ic, warp_size_zu);
  ic = rbond_final_k.putHost(&nmr_double2_data, filtered_rbond_finl_keq, ic, warp_size_zu);
  ic = rangl_init_k.putHost(&nmr_double2_data, filtered_rangl_init_keq, ic, warp_size_zu);
  ic = rangl_final_k.putHost(&nmr_double2_data, filtered_rangl_finl_keq, ic, warp_size_zu);
  ic = rdihe_init_k.putHost(&nmr_double2_data, filtered_rdihe_init_keq, ic, warp_size_zu);
  ic = rdihe_final_k.putHost(&nmr_double2_data, filtered_rdihe_finl_keq, ic, warp_size_zu);
  ic = 0LLU;
  ic = rposn_init_r.putHost(&nmr_double4_data, filtered_rposn_init_r, ic, warp_size_zu);
  ic = rposn_final_r.putHost(&nmr_double4_data, filtered_rposn_finl_r, ic, warp_size_zu);  
  ic = rbond_init_r.putHost(&nmr_double4_data, filtered_rbond_init_r, ic, warp_size_zu);
  ic = rbond_final_r.putHost(&nmr_double4_data, filtered_rbond_finl_r, ic, warp_size_zu);
  ic = rangl_init_r.putHost(&nmr_double4_data, filtered_rangl_init_r, ic, warp_size_zu);
  ic = rangl_final_r.putHost(&nmr_double4_data, filtered_rangl_finl_r, ic, warp_size_zu);
  ic = rdihe_init_r.putHost(&nmr_double4_data, filtered_rdihe_init_r, ic, warp_size_zu);
  ic = rdihe_final_r.putHost(&nmr_double4_data, filtered_rdihe_finl_r, ic, warp_size_zu);

  // Create temporary arrays for floating-point data.  This is done more in line with what happens
  // in an AtomGraph, with each POINTER-kind Hybrid object targeting a handful of larger arrays.
  // The likelihood that most simulations will have few to no restraints makes it advantageous to
  // condense these arrays as much as possible.
  const std::vector<float> spfil_rposn_init_z(filtered_rposn_init_z.begin(),
                                              filtered_rposn_init_z.end());
  const std::vector<float> spfil_rposn_finl_z(filtered_rposn_finl_z.begin(),
                                              filtered_rposn_finl_z.end());
  const std::vector<float2> spfil_rposn_init_xy = vtConv2f(filtered_rposn_init_xy);
  const std::vector<float2> spfil_rposn_finl_xy = vtConv2f(filtered_rposn_finl_xy);
  const std::vector<float2> spfil_rposn_init_keq = vtConv2f(filtered_rposn_init_keq);
  const std::vector<float2> spfil_rposn_finl_keq = vtConv2f(filtered_rposn_finl_keq);
  const std::vector<float2> spfil_rbond_init_keq = vtConv2f(filtered_rbond_init_keq);
  const std::vector<float2> spfil_rbond_finl_keq = vtConv2f(filtered_rbond_finl_keq);
  const std::vector<float2> spfil_rangl_init_keq = vtConv2f(filtered_rangl_init_keq);
  const std::vector<float2> spfil_rangl_finl_keq = vtConv2f(filtered_rangl_finl_keq);
  const std::vector<float2> spfil_rdihe_init_keq = vtConv2f(filtered_rdihe_init_keq);
  const std::vector<float2> spfil_rdihe_finl_keq = vtConv2f(filtered_rdihe_finl_keq);
  const std::vector<float4> spfil_rposn_init_r = vtConv4f(filtered_rposn_init_r);
  const std::vector<float4> spfil_rposn_finl_r = vtConv4f(filtered_rposn_finl_r);
  const std::vector<float4> spfil_rbond_init_r = vtConv4f(filtered_rbond_init_r);
  const std::vector<float4> spfil_rbond_finl_r = vtConv4f(filtered_rbond_finl_r);
  const std::vector<float4> spfil_rangl_init_r = vtConv4f(filtered_rangl_init_r);
  const std::vector<float4> spfil_rangl_finl_r = vtConv4f(filtered_rangl_finl_r);
  const std::vector<float4> spfil_rdihe_init_r = vtConv4f(filtered_rdihe_init_r);
  const std::vector<float4> spfil_rdihe_finl_r = vtConv4f(filtered_rdihe_finl_r);
  ic = 0LLU;
  ic = sp_rposn_init_z.putHost(&nmr_float_data, spfil_rposn_init_z, ic, warp_size_zu);
  ic = sp_rposn_final_z.putHost(&nmr_float_data, spfil_rposn_finl_z, ic, warp_size_zu);
  ic = 0LLU;
  ic = sp_rposn_init_k.putHost(&nmr_float2_data, spfil_rposn_init_keq, ic, warp_size_zu);
  ic = sp_rposn_final_k.putHost(&nmr_float2_data, spfil_rposn_finl_keq, ic, warp_size_zu);
  ic = sp_rposn_init_xy.putHost(&nmr_float2_data, spfil_rposn_init_xy, ic, warp_size_zu);
  ic = sp_rposn_final_xy.putHost(&nmr_float2_data, spfil_rposn_finl_xy, ic, warp_size_zu);
  ic = sp_rbond_init_k.putHost(&nmr_float2_data, spfil_rbond_init_keq, ic, warp_size_zu);
  ic = sp_rbond_final_k.putHost(&nmr_float2_data, spfil_rbond_finl_keq, ic, warp_size_zu);
  ic = sp_rangl_init_k.putHost(&nmr_float2_data, spfil_rangl_init_keq, ic, warp_size_zu);
  ic = sp_rangl_final_k.putHost(&nmr_float2_data, spfil_rangl_finl_keq, ic, warp_size_zu);
  ic = sp_rdihe_init_k.putHost(&nmr_float2_data, spfil_rdihe_init_keq, ic, warp_size_zu);
  ic = sp_rdihe_final_k.putHost(&nmr_float2_data, spfil_rdihe_finl_keq, ic, warp_size_zu);
  ic = 0LLU;
  ic = sp_rposn_init_r.putHost(&nmr_float4_data, spfil_rposn_init_r, ic, warp_size_zu);
  ic = sp_rposn_final_r.putHost(&nmr_float4_data, spfil_rposn_finl_r, ic, warp_size_zu);  
  ic = sp_rbond_init_r.putHost(&nmr_float4_data, spfil_rbond_init_r, ic, warp_size_zu);
  ic = sp_rbond_final_r.putHost(&nmr_float4_data, spfil_rbond_finl_r, ic, warp_size_zu);
  ic = sp_rangl_init_r.putHost(&nmr_float4_data, spfil_rangl_init_r, ic, warp_size_zu);
  ic = sp_rangl_final_r.putHost(&nmr_float4_data, spfil_rangl_finl_r, ic, warp_size_zu);
  ic = sp_rdihe_init_r.putHost(&nmr_float4_data, spfil_rdihe_init_r, ic, warp_size_zu);
  ic = sp_rdihe_final_r.putHost(&nmr_float4_data, spfil_rdihe_finl_r, ic, warp_size_zu);

  // With the restraint parameter tables assembled, mark the restraints for each system in
  // terms of the synthesis tables.
  for (int sysid = 0; sysid < system_count; sysid++) {
    const int ra_index = restraint_indices.readHost(sysid);
    const RestraintApparatus* ra_ptr = restraint_networks[ra_index];
    if (ra_ptr == nullptr) {
      continue;
    }
    const RestraintApparatusDpReader rar = ra_ptr->dpData();
    const int ra_posn_table_offset = network_rposn_table_offsets[ra_index];
    const int ra_bond_table_offset = network_rbond_table_offsets[ra_index];
    const int ra_angl_table_offset = network_rangl_table_offsets[ra_index];
    const int ra_dihe_table_offset = network_rdihe_table_offsets[ra_index];
    const int synth_posn_table_offset = posn_restraint_offsets.readHost(sysid);
    const int synth_bond_table_offset = bond_restraint_offsets.readHost(sysid);
    const int synth_angl_table_offset = angl_restraint_offsets.readHost(sysid);
    const int synth_dihe_table_offset = dihe_restraint_offsets.readHost(sysid);
    for (int j = 0; j < rar.nposn; j++) {
      rposn_kr_param_idx.putHost(rposn_synthesis_kr_index[ra_posn_table_offset + j],
                                 synth_posn_table_offset + j);
      rposn_xyz_param_idx.putHost(rposn_synthesis_xyz_index[ra_posn_table_offset + j],
                                  synth_posn_table_offset + j);
    }
    for (int j = 0; j < rar.nbond; j++) {
      rbond_param_idx.putHost(rbond_synthesis_index[ra_bond_table_offset + j],
                              synth_bond_table_offset + j);
    }
    for (int j = 0; j < rar.nangl; j++) {
      rangl_param_idx.putHost(rangl_synthesis_index[ra_angl_table_offset + j],
                              synth_angl_table_offset + j);
    }
    for (int j = 0; j < rar.ndihe; j++) {
      rdihe_param_idx.putHost(rdihe_synthesis_index[ra_dihe_table_offset + j],
                              synth_dihe_table_offset + j);
    }
  }
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getTopologyCount() const {
  return topology_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAtomCount() const {
  return total_atoms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getVirtualSiteCount() const {
  return total_virtual_sites;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getBondTermCount() const {
  return total_bond_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAngleTermCount() const {
  return total_angl_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getDihedralTermCount() const {
  return total_dihe_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getUreyBradleyTermCount() const {
  return total_ubrd_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCharmmImproperTermCount() const {
  return total_cimp_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCmapTermCount() const {
  return total_cmap_terms;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAtomTypeCount() const {
  return total_atom_types;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getChargeTypeCount() const {
  return total_charge_types;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getBondParameterCount() const {
  return total_bond_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getAngleParameterCount() const {
  return total_angl_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getDihedralParameterCount() const {
  return total_dihe_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getUreyBradleyParameterCount() const {
  return total_ubrd_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCharmmImproperParameterCount() const {
  return total_cimp_params;
}

//-------------------------------------------------------------------------------------------------
int AtomGraphSynthesis::getCmapSurfaceCount() const {
  return total_cmap_surfaces;
}

//-------------------------------------------------------------------------------------------------
UnitCellType AtomGraphSynthesis::getUnitCellType() const {
  return periodic_box_class;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel AtomGraphSynthesis::getImplicitSolventModel() const {
  return gb_style;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getDielectricConstant() const {
  return dielectric_constant;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getSaltConcentration() const {
  return salt_concentration;
}

//-------------------------------------------------------------------------------------------------
double AtomGraphSynthesis::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> AtomGraphSynthesis::getPBRadiiSet() const {
  return pb_radii_sets;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> AtomGraphSynthesis::getPBRadiiSet(const int low_limit,
                                                           const int high_limit) const {
  if (high_limit < low_limit || high_limit < 0 || low_limit < 0 || high_limit > system_count ||
      low_limit > system_count) {
    rtErr("The requested range " + std::to_string(low_limit) + " - " + std::to_string(high_limit) +
          " is invalid for a collection of " + std::to_string(system_count) + " systems.",
          "AtomGraphSynthesis", "getPBRadiiSet");
  }
  std::vector<std::string> result(high_limit - low_limit);
  for (int i = low_limit; i < high_limit; i++) {
    result[i - low_limit] = pb_radii_sets[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string AtomGraphSynthesis::getPBRadiiSet(const int index) const {
  if (index < 0 || index > system_count) {
    rtErr("The index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "AtomGraphSynthesis", "getPBRadiiSet");
  }
  return pb_radii_sets[index];
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::upload() {
  int_system_data.upload();
  chem_int_data.upload();
  chem_double_data.upload();
  chem_float_data.upload();
  chem_char4_data.upload();
  valparam_double_data.upload();
  valparam_float_data.upload();
  valparam_int_data.upload();
  valence_int_data.upload();
  lennard_jones_ab_coeff.upload();
  lennard_jones_c_coeff.upload();
  lennard_jones_14_ab_coeff.upload();
  lennard_jones_14_c_coeff.upload();
  sp_lennard_jones_ab_coeff.upload();
  sp_lennard_jones_c_coeff.upload();
  sp_lennard_jones_14_ab_coeff.upload();
  sp_lennard_jones_14_c_coeff.upload();
  nmr_int2_data.upload();
  nmr_double_data.upload();
  nmr_double2_data.upload();
  nmr_double4_data.upload();
  nmr_float_data.upload();
  nmr_float2_data.upload();
  nmr_float4_data.upload();
  nmr_int_data.upload();
  virtual_site_parameters.upload();
  sp_virtual_site_parameters.upload();
  vsite_int_data.upload();
  atom_imports.upload();
  vwu_instruction_sets.upload();
  insr_uint_data.upload();
  insr_uint2_data.upload();
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::download() {
  int_system_data.download();
  chem_int_data.download();
  chem_double_data.download();
  chem_float_data.download();
  chem_char4_data.download();
  valparam_double_data.download();
  valparam_float_data.download();
  valparam_int_data.download();
  valence_int_data.download();
  lennard_jones_ab_coeff.download();
  lennard_jones_c_coeff.download();
  lennard_jones_14_ab_coeff.download();
  lennard_jones_14_c_coeff.download();
  sp_lennard_jones_ab_coeff.download();
  sp_lennard_jones_c_coeff.download();
  sp_lennard_jones_14_ab_coeff.download();
  sp_lennard_jones_14_c_coeff.download();
  nmr_int2_data.download();
  nmr_double_data.download();
  nmr_double2_data.download();
  nmr_double4_data.download();
  nmr_float_data.download();
  nmr_float2_data.download();
  nmr_float4_data.download();
  nmr_int_data.download();
  virtual_site_parameters.download();
  sp_virtual_site_parameters.download();
  vsite_int_data.download();
  atom_imports.download();
  vwu_instruction_sets.download();
  insr_uint_data.download();
  insr_uint2_data.download();
}
#endif

} // namespace synthesis
} // namespace omni
