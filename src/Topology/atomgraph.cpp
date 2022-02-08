#include <cmath>
#include <cstdio>
#include <climits>
#include "Constants/generalized_born.h"
#include "Constants/scaling.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/unit_test.h"
#include "amber_prmtop_util.h"
#include "atomgraph.h"

namespace omni {
namespace topology {

using constants::PrecisionModel;
using card::HybridTargetLevel;
using card::HybridKind;
using math::addScalarToVector;
using math::roundUp;
using parse::CaseSensitivity;
using parse::char4ToString;
using parse::extractFormattedNumber;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::strncmpCased;
using parse::TextFile;
using parse::TextOrigin;
using parse::verifyNumberFormat;
using parse::operator==;
using testing::Approx;
using namespace generalized_born_defaults;

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph() :
    version_stamp{""},
    date{0, 0, 0, 0, 0, 0, 0, 0, 0},
    title{""},
    source{""},
    force_fields{},

    // Counts of atoms, residues, and other parts of the system
    atom_count{0}, residue_count{0}, molecule_count{0}, largest_residue_size{0},
    last_solute_residue{0}, last_solute_atom{0}, first_solvent_molecule{0},
    last_atom_before_cap{0}, implicit_copy_count{0}, largest_molecule_size{0},
    descriptors{HybridKind::POINTER, "tp_desc"},
    residue_limits{HybridKind::POINTER, "tp_res_limits"},
    atom_struc_numbers{HybridKind::POINTER, "tp_atom_struc_nums"},
    residue_numbers{HybridKind::POINTER, "tp_res_struc_nums"},
    molecule_limits{HybridKind::POINTER, "tp_mol_limits"},

    // Atom and residue details
    atomic_numbers{HybridKind::POINTER, "tp_znum"},
    mobile_atoms{HybridKind::POINTER, "tp_belly"},
    molecule_membership{HybridKind::POINTER, "tp_molnum"},
    molecule_contents{HybridKind::POINTER, "tp_mol_contents"},
    atomic_charges{HybridKind::POINTER, "tp_atomq"},
    atomic_masses{HybridKind::POINTER, "tp_mass"},
    inverse_atomic_masses{HybridKind::POINTER, "tp_invmass"},
    sp_atomic_charges{HybridKind::POINTER, "tp_atomq_sp"},
    sp_atomic_masses{HybridKind::POINTER, "tp_mass_sp"},
    sp_inverse_atomic_masses{HybridKind::POINTER, "tp_invmass_sp"},
    atom_names{HybridKind::POINTER, "tp_atom_names"},
    atom_types{HybridKind::POINTER, "tp_atom_types"},
    residue_names{HybridKind::POINTER, "tp_res_names"},

    // CHARMM force field family parameters
    urey_bradley_term_count{0}, charmm_impr_term_count{0}, cmap_term_count{0},
    urey_bradley_parameter_count{0}, charmm_impr_parameter_count{0}, cmap_surface_count{0},
    urey_bradley_pert_term_count{0}, charmm_impr_pert_term_count{0}, cmap_pert_term_count{0},
    urey_bradleys_in_perturbed_group{0}, charmm_imprs_in_perturbed_group{0},
    cmaps_in_perturbed_group{0},
    urey_bradley_i_atoms{HybridKind::POINTER, "tp_urey_i"},
    urey_bradley_k_atoms{HybridKind::POINTER, "tp_urey_k"},
    urey_bradley_parameter_indices{HybridKind::POINTER, "tp_urey_parm"},
    urey_bradley_assigned_atoms{HybridKind::POINTER, "tp_urey_asatoms"},
    urey_bradley_assigned_index{HybridKind::POINTER, "tp_urey_asindex"},
    urey_bradley_assigned_terms{HybridKind::POINTER, "tp_urey_asterms"},
    urey_bradley_assigned_bounds{HybridKind::POINTER, "tp_urey_asbounds"},
    charmm_impr_i_atoms{HybridKind::POINTER, "tp_charmm_impr_i"},
    charmm_impr_j_atoms{HybridKind::POINTER, "tp_charmm_impr_j"},
    charmm_impr_k_atoms{HybridKind::POINTER, "tp_charmm_impr_k"},
    charmm_impr_l_atoms{HybridKind::POINTER, "tp_charmm_impr_l"},
    charmm_impr_parameter_indices{HybridKind::POINTER, "tp_charmm_impr_parm"},
    charmm_impr_assigned_atoms{HybridKind::POINTER, "tp_cimpr_asatoms"},
    charmm_impr_assigned_index{HybridKind::POINTER, "tp_cimpr_asindex"},
    charmm_impr_assigned_terms{HybridKind::POINTER, "tp_cimpr_asterms"},
    charmm_impr_assigned_bounds{HybridKind::POINTER, "tp_cimpr_asbounds"},
    cmap_i_atoms{HybridKind::POINTER, "tp_cmap_i"},
    cmap_j_atoms{HybridKind::POINTER, "tp_cmap_j"},
    cmap_k_atoms{HybridKind::POINTER, "tp_cmap_k"},
    cmap_l_atoms{HybridKind::POINTER, "tp_cmap_l"},
    cmap_m_atoms{HybridKind::POINTER, "tp_cmap_m"},
    cmap_surface_dimensions{HybridKind::POINTER, "tp_cmap_dims"},
    cmap_surface_bounds{HybridKind::POINTER, "tp_cmap_bounds"},
    cmap_patch_bounds{HybridKind::POINTER, "tp_cmap_patch_bounds"},
    cmap_surface_indices{HybridKind::POINTER, "tp_cmap_parm"},
    cmap_assigned_atoms{HybridKind::POINTER, "tp_cmap_asatoms"},
    cmap_assigned_index{HybridKind::POINTER, "tp_cmap_asindex"},
    cmap_assigned_terms{HybridKind::POINTER, "tp_cmap_asterms"},
    cmap_assigned_bounds{HybridKind::POINTER, "tp_cmap_asbounds"},
    urey_bradley_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness"},
    urey_bradley_equilibria{HybridKind::POINTER, "tp_ub_equilibria"},
    charmm_impr_stiffnesses{HybridKind::POINTER, "tp_cimp_stiffness"},
    charmm_impr_phase_angles{HybridKind::POINTER, "tp_cimp_equilibria"},
    cmap_surfaces{HybridKind::POINTER, "tp_cmap_surface"},
    cmap_phi_derivatives{HybridKind::POINTER, "tp_cmap_dphi"},
    cmap_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi"},
    cmap_phi_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi_dpsi"},
    cmap_patches{HybridKind::POINTER, "tp_cmap_patches"},
    sp_urey_bradley_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness_sp"},
    sp_urey_bradley_equilibria{HybridKind::POINTER, "tp_ub_equilibria_sp"},
    sp_charmm_impr_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness_sp"},
    sp_charmm_impr_phase_angles{HybridKind::POINTER, "tp_ub_equilibria_sp"},
    sp_cmap_surfaces{HybridKind::POINTER, "tp_cmap_surface_sp"},
    sp_cmap_phi_derivatives{HybridKind::POINTER, "tp_cmap_dphi"},
    sp_cmap_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi"},
    sp_cmap_phi_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi_dpsi"},
    sp_cmap_patches{HybridKind::POINTER, "tp_cmap_patches"},

    // Relevant information for the bonded calculation
    bond_term_with_hydrogen{0}, angl_term_with_hydrogen{0}, dihe_term_with_hydrogen{0},
    bond_term_without_hydrogen{0}, angl_term_without_hydrogen{0}, dihe_term_without_hydrogen{0},
    bond_term_count{0}, angl_term_count{0}, dihe_term_count{0}, bond_parameter_count{0},
    angl_parameter_count{0}, dihe_parameter_count{0}, bond_perturbation_term_count{0},
    angl_perturbation_term_count{0}, dihe_perturbation_term_count{0}, bonds_in_perturbed_group{0},
    angls_in_perturbed_group{0}, dihes_in_perturbed_group{0}, bonded_group_count{0},
    bond_stiffnesses{HybridKind::POINTER, "tp_bondk"},
    bond_equilibria{HybridKind::POINTER, "tp_bondl0"},
    angl_stiffnesses{HybridKind::POINTER, "tp_anglk"},
    angl_equilibria{HybridKind::POINTER, "tp_anglt0"},
    dihe_amplitudes{HybridKind::POINTER, "tp_dihek"},
    dihe_periodicities{HybridKind::POINTER, "tp_dihen"},
    dihe_phase_angles{HybridKind::POINTER, "tp_dihepsi"},
    sp_bond_stiffnesses{HybridKind::POINTER, "tp_bondk_sp"},
    sp_bond_equilibria{HybridKind::POINTER, "tp_bondl0_sp"},
    sp_angl_stiffnesses{HybridKind::POINTER, "tp_anglk_sp"},
    sp_angl_equilibria{HybridKind::POINTER, "tp_anglt0_sp"},
    sp_dihe_amplitudes{HybridKind::POINTER, "tp_dihek_sp"},
    sp_dihe_periodicities{HybridKind::POINTER, "tp_dihen_sp"},
    sp_dihe_phase_angles{HybridKind::POINTER, "tp_dihepsi_sp"},
    bond_i_atoms{HybridKind::POINTER, "tp_bond_i"},
    bond_j_atoms{HybridKind::POINTER, "tp_bond_j"},
    bond_parameter_indices{HybridKind::POINTER, "tp_bond_parm"},
    bond_assigned_atoms{HybridKind::POINTER, "tp_bond_asatoms"},
    bond_assigned_index{HybridKind::POINTER, "tp_bond_asindex"},
    bond_assigned_terms{HybridKind::POINTER, "tp_bond_asterms"},
    bond_assigned_bounds{HybridKind::POINTER, "tp_bond_asbounds"},
    angl_i_atoms{HybridKind::POINTER, "tp_angl_i"},
    angl_j_atoms{HybridKind::POINTER, "tp_angl_j"},
    angl_k_atoms{HybridKind::POINTER, "tp_angl_k"},
    angl_parameter_indices{HybridKind::POINTER, "tp_angl_parm"},
    angl_assigned_atoms{HybridKind::POINTER, "tp_angl_asatoms"},
    angl_assigned_index{HybridKind::POINTER, "tp_angl_asindex"},
    angl_assigned_terms{HybridKind::POINTER, "tp_angl_asterms"},
    angl_assigned_bounds{HybridKind::POINTER, "tp_angl_asbounds"},
    dihe_i_atoms{HybridKind::POINTER, "tp_dihe_i"},
    dihe_j_atoms{HybridKind::POINTER, "tp_dihe_j"},
    dihe_k_atoms{HybridKind::POINTER, "tp_dihe_k"},
    dihe_l_atoms{HybridKind::POINTER, "tp_dihe_l"},
    dihe_parameter_indices{HybridKind::POINTER, "tp_dihe_parm"},
    dihe14_parameter_indices{HybridKind::POINTER, "tp_dihe_attn14_parm"},
    dihe_assigned_atoms{HybridKind::POINTER, "tp_dihe_asatoms"},
    dihe_assigned_index{HybridKind::POINTER, "tp_dihe_asindex"},
    dihe_assigned_terms{HybridKind::POINTER, "tp_dihe_asterms"},
    dihe_assigned_bounds{HybridKind::POINTER, "tp_dihe_asbounds"},
    bond_modifiers{HybridKind::POINTER, "tp_bond_mods"},
    angl_modifiers{HybridKind::POINTER, "tp_angl_mods"},
    dihe_modifiers{HybridKind::POINTER, "tp_dihe_mods"},
    bond_assigned_mods{HybridKind::POINTER, "tp_bond_asmods"},
    angl_assigned_mods{HybridKind::POINTER, "tp_angl_asmods"},
    dihe_assigned_mods{HybridKind::POINTER, "tp_dihe_asmods"},

    // Information relevant to virtual site placement
    virtual_site_count{0},
    virtual_site_atoms{HybridKind::POINTER, "tp_vsidx"},
    virtual_site_frame_types{HybridKind::POINTER, "tp_vsfrm"},
    virtual_site_frame1_atoms{HybridKind::POINTER, "tp_vsfr1"},
    virtual_site_frame2_atoms{HybridKind::POINTER, "tp_vsfr2"},
    virtual_site_frame3_atoms{HybridKind::POINTER, "tp_vsfr3"},
    virtual_site_frame4_atoms{HybridKind::POINTER, "tp_vsfr4"},
    virtual_site_frame_dim1{HybridKind::POINTER, "tp_vsdim1"},
    virtual_site_frame_dim2{HybridKind::POINTER, "tp_vsdim2"},
    virtual_site_frame_dim3{HybridKind::POINTER, "tp_vsdim3"},
    sp_virtual_site_frame_dim1{HybridKind::POINTER, "tp_vsdim1_sp"},
    sp_virtual_site_frame_dim2{HybridKind::POINTER, "tp_vsdim2_sp"},
    sp_virtual_site_frame_dim3{HybridKind::POINTER, "tp_vsdim3_sp"},

    // Relevant information for the non-bonded calculation
    charge_type_count{0}, atom_type_count{0}, total_exclusions{0}, attenuated_14_type_count{0},
    inferred_14_attenuations{0}, periodic_box_class{UnitCellType::NONE},
    gb_style{ImplicitSolventModel::NONE}, dielectric_constant{1.0}, salt_concentration{0.0},
    coulomb_constant{accepted_coulomb_constant}, pb_radii_set{""},
    charge_indices{HybridKind::POINTER, "tp_qidx"},
    lennard_jones_indices{HybridKind::POINTER, "tp_ljidx"},
    atom_exclusion_bounds{HybridKind::POINTER, "tp_nexcl"},
    atom_exclusion_list{HybridKind::POINTER, "tp_excllist"},
    nb11_exclusion_bounds{HybridKind::POINTER, "tp_nnb11"},
    nb11_exclusion_list{HybridKind::POINTER, "tp_nb11_excl"},
    nb12_exclusion_bounds{HybridKind::POINTER, "tp_nnb12"},
    nb12_exclusion_list{HybridKind::POINTER, "tp_nb12_excl"},
    nb13_exclusion_bounds{HybridKind::POINTER, "tp_nnb13"},
    nb13_exclusion_list{HybridKind::POINTER, "tp_nb13_excl"},
    nb14_exclusion_bounds{HybridKind::POINTER, "tp_nnb14"},
    nb14_exclusion_list{HybridKind::POINTER, "tp_nb14_excl"},
    infr14_i_atoms{HybridKind::POINTER, "tp_inferred14_i"},
    infr14_j_atoms{HybridKind::POINTER, "tp_inferred14_j"},
    infr14_parameter_indices{HybridKind::POINTER, "tp_inferred14_param"},
    neck_gb_indices{HybridKind::POINTER, "tp_gbneck_idx"},
    charge_parameters{HybridKind::POINTER, "tp_qparam"},
    lj_a_values{HybridKind::POINTER, "tp_lja"},
    lj_b_values{HybridKind::POINTER, "tp_ljb"},
    lj_c_values{HybridKind::POINTER, "tp_ljc"},
    lj_14_a_values{HybridKind::POINTER, "tp_lja"},
    lj_14_b_values{HybridKind::POINTER, "tp_ljb"},
    lj_14_c_values{HybridKind::POINTER, "tp_ljc"},
    lj_type_corrections{HybridKind::POINTER, "tp_lj_long"},
    attn14_elec_factors{HybridKind::POINTER, "tp_scee_param"},
    attn14_vdw_factors{HybridKind::POINTER, "tp_scnb_param"},
    atomic_pb_radii{HybridKind::POINTER, "tp_pbradii"},
    gb_screening_factors{HybridKind::POINTER, "tp_screen"},
    gb_alpha_parameters{HybridKind::POINTER, "tp_gb_alpha_sp"},
    gb_beta_parameters{HybridKind::POINTER, "tp_gb_beta_sp"},
    gb_gamma_parameters{HybridKind::POINTER, "tp_gb_gamma_sp"},
    sp_charge_parameters{HybridKind::POINTER, "tp_qparam_sp"},
    sp_lj_a_values{HybridKind::POINTER, "tp_lja_sp"},
    sp_lj_b_values{HybridKind::POINTER, "tp_ljb_sp"},
    sp_lj_c_values{HybridKind::POINTER, "tp_ljc_sp"},
    sp_lj_14_a_values{HybridKind::POINTER, "tp_lja_sp"},
    sp_lj_14_b_values{HybridKind::POINTER, "tp_ljb_sp"},
    sp_lj_14_c_values{HybridKind::POINTER, "tp_ljc_sp"},
    sp_lj_type_corrections{HybridKind::POINTER, "tp_lj_long_sp"},
    sp_attn14_elec_factors{HybridKind::POINTER, "tp_scee_param_sp"},
    sp_attn14_vdw_factors{HybridKind::POINTER, "tp_scnb_param_sp"},
    sp_atomic_pb_radii{HybridKind::POINTER, "tp_pbradii_sp"},
    sp_gb_screening_factors{HybridKind::POINTER, "tp_screen_sp"},
    sp_gb_alpha_parameters{HybridKind::POINTER, "tp_gb_alpha_sp"},
    sp_gb_beta_parameters{HybridKind::POINTER, "tp_gb_beta_sp"},
    sp_gb_gamma_parameters{HybridKind::POINTER, "tp_gb_gamma_sp"},

    // MD propagation algorithm directives
    use_bond_constraints{ShakeSetting::OFF}, use_settle{SettleSetting::OFF},
    use_perturbation_info{PerturbationSetting::OFF}, use_solvent_cap_option{SolventCapSetting::ON},
    use_polarization{PolarizationSetting::ON}, water_residue_name{' ', ' ', ' ', ' '},
    bond_constraint_mask{""}, bond_constraint_omit_mask{""}, rigid_water_count{0},
    bond_constraint_count{0}, degrees_of_freedom{0}, nonrigid_particle_count{0},

    // Overflow name keys
    atom_overflow_names{HybridKind::POINTER, "atom_name_xtkey"},
    atom_overflow_types{HybridKind::POINTER, "atom_type_xtkey"},
    residue_overflow_names{HybridKind::POINTER, "residue_name_xtkey"},

    // Information currently unused
    unused_nhparm{0}, unused_nparm{0}, unused_natyp{0}, hbond_10_12_parameter_count{0},
    heavy_bonds_plus_constraints{0}, heavy_angls_plus_constraints{0},
    heavy_dihes_plus_constraints{0},
    tree_joining_info{HybridKind::POINTER, "tp_join"},
    last_rotator_info{HybridKind::POINTER, "tp_irotat"},
    solty_info{HybridKind::POINTER, "tp_solty"},
    hbond_a_values{HybridKind::POINTER, "tp_sola"},
    hbond_b_values{HybridKind::POINTER, "tp_solb"},
    hbond_cutoffs{HybridKind::POINTER, "tp_hbcut"},
    tree_symbols{HybridKind::POINTER, "tp_symbl"},

    // Hybrid data structures (actual arrays)
    int_data{HybridKind::ARRAY, "tp_int"},
    double_data{HybridKind::ARRAY, "tp_double"},
    float_data{HybridKind::ARRAY, "tp_float"},
    char4_data{HybridKind::ARRAY, "tp_char4"}
{}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const std::string &file_name, const ExceptionResponse policy,
                     const TopologyKind engine_format) :
    AtomGraph()
{
  switch (engine_format) {
  case TopologyKind::AMBER:
    buildFromPrmtop(file_name, policy);
    break;
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    rtErr("Construction from non-Amber format files is not yet implemented.", "AtomGraph");
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const std::string &file_name, const ExceptionResponse policy,
                     const TopologyKind engine_format, const double coulomb_constant_in,
                     const double default_elec14_screening, const double default_vdw14_screening,
                     const double charge_rounding_tol, const double charge_discretization) :
    AtomGraph()
{
  switch (engine_format) {
  case TopologyKind::AMBER:
    buildFromPrmtop(file_name, policy, coulomb_constant_in, default_elec14_screening,
                    default_vdw14_screening, charge_rounding_tol, charge_discretization);
    break;
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    rtErr("Construction from non-Amber format files is not yet implemented.", "AtomGraph");
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::buildFromPrmtop(const std::string &file_name, const ExceptionResponse policy,
                                const double coulomb_constant_in,
                                const double default_elec14_screening,
                                const double default_vdw14_screening,
                                const double charge_rounding_tol,
                                const double charge_discretization) {

  // Log the file that was the source
  source = file_name;

  // Set Coulomb's constant (to the Amber-specific value if the default setting for this function
  // remains in place)
  coulomb_constant = coulomb_constant_in;

  // Get the input file as a big text vector
  const TextFile fmem(file_name, TextOrigin::DISK, std::string(""),
                      "Prmtop-based AtomGraph constructor");

  // Begin parsing, starting with the version stamp and date
  const TextFile::Reader tfr = fmem.data();
  if (tfr.line_count > 0 && strncmp(tfr.text, "%VERSION", 8) == 0) {
    if (tfr.line_limits[1] >= 34) {
      for (int i = 26; i < 35; i++) {
        version_stamp[i - 26] = tfr.text[i];
      }
      version_stamp[9] = '\0';
    }
    else {
      sprintf(version_stamp, "UNKNOWN");
    }
    if (tfr.line_limits[1] >= 61 && strncmp(&tfr.text[37], "DATE", 4) == 0) {
      char td[4];
      td[0] = tfr.text[44];
      td[1] = tfr.text[45];
      td[2] = '\0';
      date.tm_mon = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) - 1 : 0;
      td[0] = tfr.text[47];
      td[1] = tfr.text[48];
      td[2] = '\0';
      date.tm_mday = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[50];
      td[1] = tfr.text[51];
      td[2] = '\0';
      date.tm_year = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) + 100 : 0;
      td[0] = tfr.text[54];
      td[1] = tfr.text[55];
      td[2] = '\0';
      date.tm_hour = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[57];
      td[1] = tfr.text[58];
      td[2] = '\0';
      date.tm_min = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[60];
      td[1] = tfr.text[61];
      td[2] = '\0';
      date.tm_sec = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
    }
    else {
      time_t rawtime;
      tm *timeinfo;
      time(&rawtime);
      timeinfo = localtime(&rawtime);
      date = *timeinfo;
    }
  }

  // Vector to hold the detected format from section to section (this is VERY volatile)
  std::vector<int4> dfmt;

  // Read title
  int lstart = scanToFlag(fmem, "TITLE", &dfmt, TopologyRequirement::OPTIONAL);
  if (lstart == -1) {
    title = "";
  }
  else {
    const int title_char_count = tfr.line_limits[lstart + 1] - tfr.line_limits[lstart];
    title.resize(title_char_count);
    for (int i = 0; i < title_char_count; i++) {
      title[i] = tfr.text[tfr.line_limits[lstart] + i];
    }
  }

  // Get all of the descriptors.  The current Amber topology has one optional descriptor.
  const ulint max_descriptors = static_cast<ulint>(TopologyDescriptor::N_VALUES);
  lstart = scanToFlag(fmem, "POINTERS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_desc = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               max_descriptors - 1, max_descriptors);

  // Assign descriptors to counts of atoms, residues, and other parts of the system
  atom_count = tmp_desc[0];
  residue_count = tmp_desc[11];
  largest_residue_size = tmp_desc[28];
  implicit_copy_count = (tmp_desc.size() > 31) ? tmp_desc[31] : 0;

  // Assign descriptors relevant to the bonded calculation
  bond_term_with_hydrogen = tmp_desc[2];
  angl_term_with_hydrogen = tmp_desc[4];
  dihe_term_with_hydrogen = tmp_desc[6];
  bond_term_without_hydrogen = tmp_desc[3];
  angl_term_without_hydrogen = tmp_desc[5];
  dihe_term_without_hydrogen = tmp_desc[7];
  bond_term_count = bond_term_with_hydrogen + bond_term_without_hydrogen;
  angl_term_count = angl_term_with_hydrogen + angl_term_without_hydrogen;
  dihe_term_count = dihe_term_with_hydrogen + dihe_term_without_hydrogen;
  bond_parameter_count = tmp_desc[15];
  angl_parameter_count = tmp_desc[16];
  dihe_parameter_count = tmp_desc[17];
  bond_perturbation_term_count = tmp_desc[21];
  angl_perturbation_term_count = tmp_desc[22];
  dihe_perturbation_term_count = tmp_desc[23];
  bonds_in_perturbed_group = tmp_desc[24];
  angls_in_perturbed_group = tmp_desc[25];
  dihes_in_perturbed_group = tmp_desc[26];

  // Assign descriptors relevant to virtual site placement
  virtual_site_count = tmp_desc[30];

  // Assign descriptors relevant to the non-bonded calculation
  atom_type_count = tmp_desc[1];
  total_exclusions = tmp_desc[10];
  switch (tmp_desc[27]) {
  case 0:
    periodic_box_class = UnitCellType::NONE;
    break;
  case 1:
    periodic_box_class = UnitCellType::ORTHORHOMBIC;
    break;
  case 2:
    periodic_box_class = UnitCellType::TRICLINIC;
    break;
  default:
    rtErr("The IFBOX field has an invalid value of " + std::to_string(tmp_desc[27]), "AtomGraph");
  }

  // Assign descriptors relevant to the MD propagation algorithm
  switch (tmp_desc[20]) {
  case 0:
    use_perturbation_info = PerturbationSetting::OFF;
    break;
  case 1:
    use_perturbation_info = PerturbationSetting::ON;
    break;
  default:
    rtErr("The IFPERT field has an invalid value of " + std::to_string(tmp_desc[20]), "AtomGraph");
  }
  switch (tmp_desc[29]) {
  case 0:
    use_solvent_cap_option = SolventCapSetting::OFF;
    break;
  case 1:
    use_solvent_cap_option = SolventCapSetting::ON;
    break;
  default:
    rtErr("The IFCAP field has an invalid value of " + std::to_string(tmp_desc[29]), "AtomGraph");
  }

  // Assign descriptors pertaining to deprecated or unused information
  unused_nhparm = tmp_desc[8];
  unused_nparm = tmp_desc[9];
  unused_natyp = tmp_desc[18];
  hbond_10_12_parameter_count = tmp_desc[19];
  heavy_bonds_plus_constraints = tmp_desc[12];
  heavy_angls_plus_constraints = tmp_desc[13];
  heavy_dihes_plus_constraints = tmp_desc[14];

  // Read the force field citations
  lstart = scanToFlag(fmem, "FORCE_FIELD_TYPE", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    force_fields = readForceFieldReferences(fmem, lstart);
  }

  // Read atom names
  lstart = scanToFlag(fmem, "ATOM_NAME", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_atom_names = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);

  // Read atomic numbers
  lstart = scanToFlag(fmem, "ATOMIC_NUMBER", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_atomic_numbers;
  if (lstart >= 0) {
    tmp_atomic_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read charges and convert to internal units (atomic units for charges)
  lstart = scanToFlag(fmem, "CHARGE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_charges = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                     atom_count);
  for (int i = 0; i < atom_count; i++) {
    tmp_charges[i] *= inv_amber_charge_scaling;
  }
  std::vector<double> tmp_charge_parameters(atom_count);
  std::vector<int> tmp_charge_indices(atom_count);
  smoothCharges(&tmp_charges, &tmp_charge_parameters, &tmp_charge_indices,
                &charge_type_count, charge_rounding_tol, charge_discretization, file_name);

  // Read masses
  lstart = scanToFlag(fmem, "MASS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_masses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                    atom_count);

  // Read atom type (Lennard-Jones type) indices
  lstart = scanToFlag(fmem, "ATOM_TYPE_INDEX", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_lennard_jones_indices = iAmberPrmtopData(fmem, lstart, dfmt[0].x,
                                                                dfmt[0].z, atom_count);
  addScalarToVector(&tmp_lennard_jones_indices, -1);

  // Read excluded atom counts and immediately convert to a prefix sum, then cap it to create
  // a bounds array indexing into the exclusion list.
  lstart = scanToFlag(fmem, "NUMBER_EXCLUDED_ATOMS", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_atom_exclusion_counts = iAmberPrmtopData(fmem, lstart, dfmt[0].x,
                                                                dfmt[0].z, atom_count);

  // Read non-bonded parameter indices
  lstart = scanToFlag(fmem, "NONBONDED_PARM_INDEX", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_nonbonded_parameter_index =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_type_count * atom_type_count);
  addScalarToVector(&tmp_nonbonded_parameter_index, -1);

  // Read residue labels
  lstart = scanToFlag(fmem, "RESIDUE_LABEL", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_residue_names = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, residue_count);

  // Read residue limits and cap the result with the total number of atoms
  lstart = scanToFlag(fmem, "RESIDUE_POINTER", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_residue_limits = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         residue_count);
  addScalarToVector(&tmp_residue_limits, -1);
  tmp_residue_limits.push_back(atom_count);

  // Assemble the array of atom structural numbers
  lstart = scanToFlag(fmem, "ATOMIC_STRUCTURE_NUMBERS", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<int> tmp_atom_struc_numbers;
  if (lstart >= 0) {
    tmp_atom_struc_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }
  else {
    tmp_atom_struc_numbers.resize(atom_count);
    for (int i = 0; i < atom_count; i++) {
      tmp_atom_struc_numbers[i] = i + 1;
    }
  }
  
  // Assemble the array of residue numbers for each atom
  lstart = scanToFlag(fmem, "RESIDUE_NUMBERS", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_residue_numbers(atom_count);
  if (lstart >= 0) {
    std::vector<int> stated_residue_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                               residue_count);
    for (int i = 0; i < residue_count; i++) {
      for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
        tmp_residue_numbers[j] = stated_residue_numbers[i];
      }
    }
  }
  else {
    for (int i = 0; i < residue_count; i++) {
      for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
        tmp_residue_numbers[j] = i + 1;
      }
    }
  }
  
  // Read bond parameters
  lstart = scanToFlag(fmem, "BOND_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_bond_stiffnesses =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_parameter_count);
  lstart = scanToFlag(fmem, "BOND_EQUIL_VALUE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_bond_equilibria =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_parameter_count);

  // Read angle parameters
  lstart = scanToFlag(fmem, "ANGLE_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_angl_stiffnesses =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_parameter_count);
  lstart = scanToFlag(fmem, "ANGLE_EQUIL_VALUE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_angl_equilibria =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_parameter_count);

  // Read CHARMM Urey-Bradley angle terms and parameters
  std::vector<int> tmp_ub_atoms;
  std::vector<double> tmp_ub_stiffnesses;
  std::vector<double> tmp_ub_equilibria;
  lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_COUNT", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    std::vector<int> tmp_ub_counts = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 2);
    urey_bradley_term_count = tmp_ub_counts[0];
    urey_bradley_parameter_count = tmp_ub_counts[1];
    if (urey_bradley_term_count > 0) {
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_ub_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                      urey_bradley_term_count * 3);
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_FORCE_CONSTANT", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_ub_stiffnesses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                            urey_bradley_parameter_count);
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_EQUIL_VALUE", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_ub_equilibria = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                           urey_bradley_parameter_count);
    }
  }

  // Read dihedral parameters
  lstart = scanToFlag(fmem, "DIHEDRAL_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_dihe_amplitudes =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "DIHEDRAL_PERIODICITY", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_dihe_periodicities =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "DIHEDRAL_PHASE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_dihe_phase_angles =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "SCEE_SCALE_FACTOR", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_dihe_elec_screenings;
  if (lstart >= 0) {
    tmp_dihe_elec_screenings = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                dihe_parameter_count);
  }
  lstart = scanToFlag(fmem, "SCNB_SCALE_FACTOR", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_dihe_vdw_screenings;
  if (lstart >= 0) {
    tmp_dihe_vdw_screenings = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               dihe_parameter_count);
  }

  // Read CHARMM improper dihedral parameters
  std::vector<int> tmp_charmm_impr_atoms;
  std::vector<double> tmp_charmm_impr_stiffnesses;
  std::vector<double> tmp_charmm_impr_phase_angles;
  lstart = scanToFlag(fmem, "CHARMM_NUM_IMPROPERS", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    std::vector<int> tmp_cimpr_count = iAmberPrmtopData(fmem, lstart, 1, 8, 1);
    charmm_impr_term_count = tmp_cimpr_count[0];
    if (charmm_impr_term_count > 0) {
      lstart = scanToFlag(fmem, "CHARMM_IMPROPERS", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_charmm_impr_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               charmm_impr_term_count * 5);
      lstart = scanToFlag(fmem, "CHARMM_NUM_IMPR_TYPES", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      std::vector<int> tmp_cimpr_types = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 1);
      charmm_impr_parameter_count = tmp_cimpr_types[0];
      lstart = scanToFlag(fmem, "CHARMM_IMPROPER_FORCE_CONSTANT", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_charmm_impr_stiffnesses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                     charmm_impr_parameter_count);
      lstart = scanToFlag(fmem, "CHARMM_IMPROPER_PHASE", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_charmm_impr_phase_angles = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                      charmm_impr_parameter_count);
    }
  }

  // Read hydrogen bonding parameters
  lstart = scanToFlag(fmem, "SOLTY", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_solty_info;
  if (lstart >= 0) {
    tmp_solty_info = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, unused_natyp);
  }

  // Read Lennard-Jones coefficients
  const int nljparm = atom_type_count * (atom_type_count + 1) / 2;
  lstart = scanToFlag(fmem, "LENNARD_JONES_ACOEF", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_lj_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         nljparm);
  lstart = scanToFlag(fmem, "LENNARD_JONES_BCOEF", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_lj_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         nljparm);
  lstart = scanToFlag(fmem, "LENNARD_JONES_CCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<double> tmp_lj_c_values;
  if (lstart >= 0) {
    tmp_lj_c_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  std::vector<double> tmp_lj_14_a_values;
  std::vector<double> tmp_lj_14_b_values;
  std::vector<double> tmp_lj_14_c_values;
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_ACOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_a_values = tmp_lj_a_values;
  }
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_BCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_b_values = tmp_lj_b_values;
  }
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_CCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_c_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_c_values = tmp_lj_c_values;
  }

  // Read angle, and dihedral indexing data.  Bond indexing data was read nearer the beginning in
  // order to calculate or confirm the number of separate molecules in the system.
  lstart = scanToFlag(fmem, "BONDS_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_bond_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_term_with_hydrogen * 3);
  lstart = scanToFlag(fmem, "BONDS_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_bond_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_term_without_hydrogen * 3);
  lstart = scanToFlag(fmem, "ANGLES_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_angl_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_term_with_hydrogen * 4);
  lstart = scanToFlag(fmem, "ANGLES_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_angl_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_term_without_hydrogen * 4);
  lstart = scanToFlag(fmem, "DIHEDRALS_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_dihe_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_term_with_hydrogen * 5);
  lstart = scanToFlag(fmem, "DIHEDRALS_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_dihe_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_term_without_hydrogen * 5);

  // Read the excluded atoms list (despite counts being read a very long way up in the file).
  lstart = scanToFlag(fmem, "EXCLUDED_ATOMS_LIST", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_exclusion_list = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         total_exclusions);

  // Read more hydrogen bonding parameters
  std::vector<double> tmp_hbond_a_values;
  std::vector<double> tmp_hbond_b_values;
  std::vector<double> tmp_hbond_cutoffs;
  if (hbond_10_12_parameter_count > 0) {
    lstart = scanToFlag(fmem, "HBOND_ACOEF", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                          hbond_10_12_parameter_count);
    lstart = scanToFlag(fmem, "HBOND_BCOEF", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                          hbond_10_12_parameter_count);
    lstart = scanToFlag(fmem, "HBCUT", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_cutoffs = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                         hbond_10_12_parameter_count);
  }
  expandLennardJonesTables(&tmp_lj_a_values, &tmp_lj_b_values, &tmp_lj_c_values,
                           &tmp_lj_14_a_values, &tmp_lj_14_b_values, &tmp_lj_14_c_values,
                           &tmp_hbond_a_values, &tmp_hbond_b_values, atom_type_count,
                           tmp_nonbonded_parameter_index);

  // Read the atom type names
  lstart = scanToFlag(fmem, "AMBER_ATOM_TYPE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_atom_types = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);

  // Read the tree chain information
  lstart = scanToFlag(fmem, "TREE_CHAIN_CLASSIFICATION", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<char4> tmp_tree_symbols;
  if (lstart >= 0) {
    tmp_tree_symbols = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);
  }
  lstart = scanToFlag(fmem, "JOIN_ARRAY", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_tree_joining_info;
  if (lstart >= 0) {
    tmp_tree_joining_info = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }
  lstart = scanToFlag(fmem, "IROTAT", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_last_rotator_info;
  if (lstart >= 0) {
    tmp_last_rotator_info = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read the GB radius set and other parameters.  These are always allocated, to prepare for
  // applying a radius set as part of the implicit solvent model later.
  std::vector<double> tmp_atomic_pb_radii(atom_count, 0.0);
  std::vector<double> tmp_gb_screening_factors(atom_count, 0.0);
  std::vector<double> tmp_gb_coef(atom_count, 0.0);
  std::vector<int> tmp_neck_gb_indices(atom_count, 0);
  lstart = scanToFlag(fmem, "RADIUS_SET", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    for (int i = tfr.line_limits[lstart]; i < tfr.line_limits[lstart + 1]; i++) {
      pb_radii_set += tfr.text[i];
    }
    lstart = scanToFlag(fmem, "RADII", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_atomic_pb_radii = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
    lstart = scanToFlag(fmem, "SCREEN", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_gb_screening_factors = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read an extra descriptor for polarizability
  lstart = scanToFlag(fmem, "IPOL", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_ipol;
  std::vector<double> tmp_atpol;
  if (lstart >= 0) {
    tmp_ipol = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 1);

    // Proceed to read more information about polarizabilities
    if (tmp_ipol[0] == 1) {
      use_polarization = PolarizationSetting::ON;
      lstart = scanToFlag(fmem, "POLARIZABILITY", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
      tmp_atpol = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
    }
  }

  // Read virtual site frame information
  std::vector<int> vsite_custom_frames;
  std::vector<double> vsite_custom_details;
  lstart = scanToFlag(fmem, "VIRTUAL_SITE_FRAMES", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    vsite_custom_frames = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 0,
                                           virtual_site_count * 6);
  }
  lstart = scanToFlag(fmem, "VIRTUAL_SITE_FRAME_DETAILS", &dfmt,
                      TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    vsite_custom_details = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 0,
                                            virtual_site_count * 3);
  }

  // Read CHARMM CMAP terms and parameters
  std::vector<int> tmp_cmap_atoms;
  std::vector<int> tmp_cmap_surf_dims;
  std::vector<int> tmp_cmap_surf_bounds;
  std::vector<double> tmp_cmap_surfaces;
  std::vector<std::string> cmap_aliases = {"CHARMM_CMAP_", "CMAP_"};
  const int n_alias = cmap_aliases.size();
  for (int i = 0; i < n_alias; i++) {
    const std::string count_flag = cmap_aliases[i] + "COUNT";
    lstart = scanToFlag(fmem, count_flag.c_str(), &dfmt, TopologyRequirement::OPTIONAL, lstart);
    if (lstart >= 0) {
      const std::string res_flag = cmap_aliases[i] + "RESOLUTION";
      const std::string parm_flag = cmap_aliases[i] + "PARAMETER_";
      const std::string index_flag = cmap_aliases[i] + "INDEX";
      const std::vector<int> tmp_cmap_count = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                               2);
      cmap_term_count = tmp_cmap_count[0];
      cmap_surface_count = tmp_cmap_count[1];
      lstart = scanToFlag(fmem, res_flag.c_str(), &dfmt, TopologyRequirement::ESSENTIAL, lstart);
      tmp_cmap_surf_dims = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                            cmap_surface_count);
      tmp_cmap_surf_bounds.resize(cmap_surface_count + 1);
      int cmap_bound_acc = 0;
      for (int i = 0; i < cmap_surface_count; i++) {
        tmp_cmap_surf_bounds[i] = cmap_bound_acc;
        cmap_bound_acc += tmp_cmap_surf_dims[i] * tmp_cmap_surf_dims[i];
        std::string cmap_name = parm_flag;
        cmap_name += (i + 1 < 10) ? "0" + std::to_string(i + 1) : std::to_string(i + 1);
        lstart = scanToFlag(fmem, cmap_name.c_str(), &dfmt, TopologyRequirement::ESSENTIAL,
                            lstart);
        const int npts = tmp_cmap_surf_dims[i] * tmp_cmap_surf_dims[i];

        // The CMAPs are read in row-major format and must be converted to column-major format.
        std::vector<double> tmp_surf = dAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, npts);
        for (int j = 0; j < tmp_cmap_surf_dims[i]; j++) {
          for (int k = 0; k < j; k++) {
            std::swap(tmp_surf[(j * tmp_cmap_surf_dims[i]) + k],
                      tmp_surf[(k * tmp_cmap_surf_dims[i]) + j]);
          }
        }
        tmp_cmap_surfaces.insert(tmp_cmap_surfaces.end(), tmp_surf.begin(), tmp_surf.end());
      }
      tmp_cmap_surf_bounds[cmap_surface_count] = cmap_bound_acc;
      lstart = scanToFlag(fmem, index_flag.c_str(), &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_cmap_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 6 * cmap_term_count);
    }
  }
  CmapAccessories cmap_table = ComputeCmapDerivatives(cmap_surface_count, tmp_cmap_surf_dims,
                                                      tmp_cmap_surf_bounds, tmp_cmap_surfaces);

  // Condense the exclusion list to avoid counting atoms with no actual exclusions as having one
  // "blank atom" exclusion.  In the Amber topology, an atom with no exclusions is customarily
  // listed as having one, and the excluded atom is index 0, which in Fortran is not a valid array
  // index.
  CondensedExclusions cond_excl = processExclusions(tmp_atom_exclusion_counts, tmp_exclusion_list,
                                                    source);

  // Create tables of all valence interactions in preparation for organizing other connections
  // and aspects of the topology
  BasicValenceTable basic_vtable = basicValenceIndexing(atom_count, tmp_bond_atoms_h,
                                                        tmp_bond_atoms_noh, tmp_angl_atoms_h,
                                                        tmp_angl_atoms_noh, tmp_dihe_atoms_h,
                                                        tmp_dihe_atoms_noh);
  CharmmValenceTable charmm_vtable = charmmValenceIndexing(atom_count, tmp_ub_atoms,
                                                           tmp_charmm_impr_atoms, tmp_cmap_atoms);

  // Condense the non-bonded 1:4 scaling factors, and flag dihedrals that are responsible for no
  // such interaction.
  AttenuationParameterSet attn_parm =
    condenseScreeningFactors(basic_vtable, tmp_dihe_elec_screenings, tmp_dihe_vdw_screenings,
                             default_elec14_screening, default_vdw14_screening);
  attenuated_14_type_count = attn_parm.total_14_sets;

  // Create vectors for virtual site frame indexing and dimensions
  VirtualSiteTable vsite_table = listVirtualSites(virtual_site_count, source,
                                                  tmp_masses, basic_vtable, tmp_bond_equilibria,
                                                  tmp_atom_types, tmp_atom_names,
                                                  vsite_custom_frames, vsite_custom_details);

  // Check the largest residue--it may have changed given the virtual site placement
  largest_residue_size = reviewLargestResidue(tmp_residue_limits, largest_residue_size, policy);

  // Elaborate on the bond connections information.  Make arrays of all 1:1 (virtual site to parent
  // atom), 1:2 (bonded atoms), 1:3 (atoms connected by a shortest path of two bonds), and 1:4
  // (atoms connected by a shortest path of three bonds) exclusions, double-counting everything to
  // list every excluded pair interaction, atom by atom.
  Map1234 all_nb_excl = mapExclusions(atom_count, basic_vtable, charmm_vtable, vsite_table);
  checkExclusions(cond_excl, all_nb_excl, source);

  // If the atomic numbers were not read from the topology file itself, infer them from the masses.
  // This can get tricky if there has been hydrogen mass repartitioning, but try and unroll that.
  if (tmp_atomic_numbers.size() == 0) {
    tmp_atomic_numbers = atomicNumbersFromMasses(tmp_masses, tmp_atom_names,
                                                 all_nb_excl.nb12_excl_list,
                                                 all_nb_excl.nb12_excl_bounds, source, policy);
  }

  // Examine dihedral coverge: are all 1:4 interactions covered by some dihedral, with a scaling
  // factor assoicated with those parameters?
  const std::vector<int3> outstanding_14_pairs =
    checkDihedral14Coverage(atom_count, tmp_atomic_numbers, basic_vtable, all_nb_excl, vsite_table,
                            attn_parm, policy);
  inferred_14_attenuations = outstanding_14_pairs.size();
  std::vector<int> tmp_inferred_14_i_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_j_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_param_idx(inferred_14_attenuations);
  for (int i = 0; i < inferred_14_attenuations; i++) {
    tmp_inferred_14_i_atoms[i]   = outstanding_14_pairs[i].x;
    tmp_inferred_14_j_atoms[i]   = outstanding_14_pairs[i].y;
    tmp_inferred_14_param_idx[i] = outstanding_14_pairs[i].z;
  }

  // With the bond connections information at the ready, calculate the number of molecules.  If
  // peridoic box information is present, compare the information generated by the search to the
  // information presented in the topology itself.
  std::vector<int> tmp_molecule_membership, tmp_molecule_limits, tmp_molecule_contents;
  mapMolecules(atom_count, &molecule_count, all_nb_excl, &tmp_molecule_membership,
               &tmp_molecule_limits, &tmp_molecule_contents);
  std::vector<int> tmp_solvent_pointers;
  std::vector<int> tmp_atoms_per_molecule;
  switch(periodic_box_class) {
  case UnitCellType::NONE:
    last_solute_residue = residue_count - 1;
    last_solute_atom = atom_count - 1;
    first_solvent_molecule = molecule_count;
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:

    // Read the demarcations between solute and solvent
    lstart = scanToFlag(fmem, "SOLVENT_POINTERS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_solvent_pointers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 3);
    last_solute_residue = tmp_solvent_pointers[0] - 1;
    last_solute_atom = (last_solute_residue == -1) ?
      -1 : tmp_residue_limits[last_solute_residue + 1] - 1;
    first_solvent_molecule = tmp_solvent_pointers[2] - 1;

    // Read the numbers of atoms per molecule
    lstart = scanToFlag(fmem, "ATOMS_PER_MOLECULE", &dfmt, TopologyRequirement::ESSENTIAL,
                        lstart);
    tmp_atoms_per_molecule = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, molecule_count);
    break;
  }
  for (int i = 0; i < molecule_count; i++) {
    largest_molecule_size = std::max(largest_molecule_size,
                                     tmp_molecule_limits[i + 1] - tmp_molecule_limits[i]);
  }
  
  // Mobile atoms are not stated in the topology, but they are an essential part of the information
  // on how the system will move.  Allocate a bitmask for all atoms to be mobile by default.
  const int mobile_atom_mask_size = roundUp<int>(atom_count / (sizeof(int) * 8), warp_size_int);
  const std::vector<int> tmp_mobile_atoms(mobile_atom_mask_size, -1);

  // Allocate the Hybrid ARRAY-kind objects based on the compiled data
  size_t int_items = roundUp<int>(tmp_desc.size(), warp_size_int) +
                     roundUp(residue_count + 1, warp_size_int) +
                     roundUp(molecule_count + 1, warp_size_int) + mobile_atom_mask_size +
                     23 * roundUp(atom_count, warp_size_int) +
                      7 * roundUp(urey_bradley_term_count, warp_size_int) +
                     10 * roundUp(charmm_impr_term_count, warp_size_int) +
                     12 * roundUp(cmap_term_count, warp_size_int) +
                     roundUp(cmap_surface_count, warp_size_int) +
                      2 * roundUp(cmap_surface_count + 1, warp_size_int) +
                      6 * roundUp(bond_term_count, warp_size_int) +
                      8 * roundUp(angl_term_count, warp_size_int) +
                     11 * roundUp(dihe_term_count, warp_size_int) +
                      6 * roundUp(virtual_site_count, warp_size_int) +
                     roundUp(total_exclusions, warp_size_int) +
                      3 * roundUp(inferred_14_attenuations, warp_size_int) +
                     roundUp<int>(all_nb_excl.nb11_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb12_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb13_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb14_excl_list.size(), warp_size_int);
  size_t double_items =  9 * roundUp(atom_count, warp_size_int) +
                         2 * roundUp(urey_bradley_parameter_count, warp_size_int) +
                         2 * roundUp(charmm_impr_parameter_count, warp_size_int) +
                        20 * roundUp<int>(tmp_cmap_surfaces.size(), warp_size_int) +
                         2 * roundUp(bond_parameter_count, warp_size_int) +
                         2 * roundUp(angl_parameter_count, warp_size_int) +
                         3 * roundUp(dihe_parameter_count, warp_size_int) +
                         3 * roundUp(virtual_site_count, warp_size_int) +
                         8 * roundUp(atom_type_count * atom_type_count, warp_size_int) +
                         2 * roundUp(atom_type_count, warp_size_int) +
                         roundUp(charge_type_count,  warp_size_int) +
                         2 * roundUp(attenuated_14_type_count, warp_size_int) ;
  size_t float_items = double_items;
  size_t char4_items = 4 * roundUp(atom_count, warp_size_int) +
                       roundUp(residue_count, warp_size_int) +
                       2 * roundUp(bond_term_count, warp_size_int) +
                       2 * roundUp(angl_term_count, warp_size_int) +
                       2 * roundUp(dihe_term_count, warp_size_int);
  int_data.resize(int_items);
  double_data.resize(double_items);
  float_data.resize(float_items);
  char4_data.resize(char4_items);

  // Lay out Hybrid POINTER-kind int objects based on the compiled data.  The putHost member
  // function of the Hybrid class has an overloaded version that sets a POINTER-kind object to its
  // target and then fills the object, thus populating the appropriate segment of the target.
  size_t ic = descriptors.putHost(&int_data, tmp_desc, 0, warp_size_zu);
  ic = residue_limits.putHost(&int_data, tmp_residue_limits, ic, warp_size_zu);
  ic = atom_struc_numbers.putHost(&int_data, tmp_atom_struc_numbers, ic, warp_size_zu);
  ic = residue_numbers.putHost(&int_data, tmp_residue_numbers, ic, warp_size_zu);
  ic = molecule_limits.putHost(&int_data, tmp_molecule_limits, ic, warp_size_zu);
  ic = atomic_numbers.putHost(&int_data, tmp_atomic_numbers, ic, warp_size_zu);
  ic = molecule_membership.putHost(&int_data, tmp_molecule_membership, ic, warp_size_zu);
  ic = mobile_atoms.putHost(&int_data, tmp_mobile_atoms, ic, warp_size_zu);
  ic = molecule_contents.putHost(&int_data, tmp_molecule_contents, ic, warp_size_zu);
  ic = urey_bradley_i_atoms.putHost(&int_data, charmm_vtable.ub_i_atoms, ic, warp_size_zu);
  ic = urey_bradley_k_atoms.putHost(&int_data, charmm_vtable.ub_k_atoms, ic, warp_size_zu);
  ic = urey_bradley_parameter_indices.putHost(&int_data, charmm_vtable.ub_param_idx, ic,
                                              warp_size_zu);
  ic = urey_bradley_assigned_atoms.putHost(&int_data, charmm_vtable.ub_assigned_atoms, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_index.putHost(&int_data, charmm_vtable.ub_assigned_index, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_terms.putHost(&int_data, charmm_vtable.ub_assigned_terms, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_bounds.putHost(&int_data, charmm_vtable.ub_assigned_bounds, ic,
                                            warp_size_zu);
  ic = charmm_impr_i_atoms.putHost(&int_data, charmm_vtable.impr_i_atoms, ic, warp_size_zu);
  ic = charmm_impr_j_atoms.putHost(&int_data, charmm_vtable.impr_j_atoms, ic, warp_size_zu);
  ic = charmm_impr_k_atoms.putHost(&int_data, charmm_vtable.impr_k_atoms, ic, warp_size_zu);
  ic = charmm_impr_l_atoms.putHost(&int_data, charmm_vtable.impr_l_atoms, ic, warp_size_zu);
  ic = charmm_impr_parameter_indices.putHost(&int_data, charmm_vtable.impr_param_idx, ic,
                                             warp_size_zu);
  ic = charmm_impr_assigned_atoms.putHost(&int_data, charmm_vtable.impr_assigned_atoms, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_index.putHost(&int_data, charmm_vtable.impr_assigned_index, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_terms.putHost(&int_data, charmm_vtable.impr_assigned_terms, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_bounds.putHost(&int_data, charmm_vtable.impr_assigned_bounds, ic,
                                           warp_size_zu);
  ic = cmap_i_atoms.putHost(&int_data, charmm_vtable.cmap_i_atoms, ic, warp_size_zu);
  ic = cmap_j_atoms.putHost(&int_data, charmm_vtable.cmap_j_atoms, ic, warp_size_zu);
  ic = cmap_k_atoms.putHost(&int_data, charmm_vtable.cmap_k_atoms, ic, warp_size_zu);
  ic = cmap_l_atoms.putHost(&int_data, charmm_vtable.cmap_l_atoms, ic, warp_size_zu);
  ic = cmap_m_atoms.putHost(&int_data, charmm_vtable.cmap_m_atoms, ic, warp_size_zu);
  ic = cmap_surface_dimensions.putHost(&int_data, tmp_cmap_surf_dims, ic, warp_size_zu);
  ic = cmap_surface_bounds.putHost(&int_data, tmp_cmap_surf_bounds, ic, warp_size_zu);
  ic = cmap_patch_bounds.putHost(&int_data, cmap_table.patch_matrix_bounds, ic, warp_size_zu);
  ic = cmap_surface_indices.putHost(&int_data, charmm_vtable.cmap_param_idx, ic, warp_size_zu);
  ic = cmap_assigned_atoms.putHost(&int_data, charmm_vtable.cmap_assigned_atoms, ic, warp_size_zu);
  ic = cmap_assigned_index.putHost(&int_data, charmm_vtable.cmap_assigned_index, ic, warp_size_zu);
  ic = cmap_assigned_terms.putHost(&int_data, charmm_vtable.cmap_assigned_terms, ic, warp_size_zu);
  ic = cmap_assigned_bounds.putHost(&int_data, charmm_vtable.cmap_assigned_bounds, ic,
                                    warp_size_zu);
  ic = bond_i_atoms.putHost(&int_data, basic_vtable.bond_i_atoms, ic, warp_size_zu);
  ic = bond_j_atoms.putHost(&int_data, basic_vtable.bond_j_atoms, ic, warp_size_zu);
  ic = bond_parameter_indices.putHost(&int_data, basic_vtable.bond_param_idx, ic, warp_size_zu);
  ic = bond_assigned_atoms.putHost(&int_data, basic_vtable.bond_assigned_atoms, ic, warp_size_zu);
  ic = bond_assigned_index.putHost(&int_data, basic_vtable.bond_assigned_index, ic, warp_size_zu);
  ic = bond_assigned_terms.putHost(&int_data, basic_vtable.bond_assigned_terms, ic, warp_size_zu);
  ic = bond_assigned_bounds.putHost(&int_data, basic_vtable.bond_assigned_bounds, ic,
                                    warp_size_zu);
  ic = angl_i_atoms.putHost(&int_data, basic_vtable.angl_i_atoms, ic, warp_size_zu);
  ic = angl_j_atoms.putHost(&int_data, basic_vtable.angl_j_atoms, ic, warp_size_zu);
  ic = angl_k_atoms.putHost(&int_data, basic_vtable.angl_k_atoms, ic, warp_size_zu);
  ic = angl_parameter_indices.putHost(&int_data, basic_vtable.angl_param_idx, ic, warp_size_zu);
  ic = angl_assigned_atoms.putHost(&int_data, basic_vtable.angl_assigned_atoms, ic, warp_size_zu);
  ic = angl_assigned_index.putHost(&int_data, basic_vtable.angl_assigned_index, ic, warp_size_zu);
  ic = angl_assigned_terms.putHost(&int_data, basic_vtable.angl_assigned_terms, ic, warp_size_zu);
  ic = angl_assigned_bounds.putHost(&int_data, basic_vtable.angl_assigned_bounds, ic,
                                    warp_size_zu);
  ic = dihe_i_atoms.putHost(&int_data, basic_vtable.dihe_i_atoms, ic, warp_size_zu);
  ic = dihe_j_atoms.putHost(&int_data, basic_vtable.dihe_j_atoms, ic, warp_size_zu);
  ic = dihe_k_atoms.putHost(&int_data, basic_vtable.dihe_k_atoms, ic, warp_size_zu);
  ic = dihe_l_atoms.putHost(&int_data, basic_vtable.dihe_l_atoms, ic, warp_size_zu);
  ic = dihe_parameter_indices.putHost(&int_data, basic_vtable.dihe_param_idx, ic, warp_size_zu);
  ic = dihe14_parameter_indices.putHost(&int_data, attn_parm.dihe14_parameter_indices, ic,
                                        warp_size_zu);
  ic = dihe_assigned_atoms.putHost(&int_data, basic_vtable.dihe_assigned_atoms, ic, warp_size_zu);
  ic = dihe_assigned_index.putHost(&int_data, basic_vtable.dihe_assigned_index, ic, warp_size_zu);
  ic = dihe_assigned_terms.putHost(&int_data, basic_vtable.dihe_assigned_terms, ic, warp_size_zu);
  ic = dihe_assigned_bounds.putHost(&int_data, basic_vtable.dihe_assigned_bounds, ic,
                                    warp_size_zu);
  ic = virtual_site_atoms.putHost(&int_data, vsite_table.vs_atoms, ic, warp_size_zu);
  ic = virtual_site_frame_types.putHost(&int_data, vsite_table.frame_types, ic, warp_size_zu);
  ic = virtual_site_frame1_atoms.putHost(&int_data, vsite_table.frame1_atoms, ic, warp_size_zu);
  ic = virtual_site_frame2_atoms.putHost(&int_data, vsite_table.frame2_atoms, ic, warp_size_zu);
  ic = virtual_site_frame3_atoms.putHost(&int_data, vsite_table.frame3_atoms, ic, warp_size_zu);
  ic = virtual_site_frame4_atoms.putHost(&int_data, vsite_table.frame4_atoms, ic, warp_size_zu);
  ic = charge_indices.putHost(&int_data, tmp_charge_indices, ic, warp_size_zu);
  ic = lennard_jones_indices.putHost(&int_data, tmp_lennard_jones_indices, ic, warp_size_zu);
  ic = atom_exclusion_bounds.putHost(&int_data, cond_excl.atom_excl_bounds, ic, warp_size_zu);
  ic = atom_exclusion_list.putHost(&int_data, cond_excl.atom_excl_list, ic, warp_size_zu);
  ic = nb11_exclusion_bounds.putHost(&int_data, all_nb_excl.nb11_excl_bounds, ic, warp_size_zu);
  ic = nb11_exclusion_list.putHost(&int_data, all_nb_excl.nb11_excl_list, ic, warp_size_zu);
  ic = nb12_exclusion_bounds.putHost(&int_data, all_nb_excl.nb12_excl_bounds, ic, warp_size_zu);
  ic = nb12_exclusion_list.putHost(&int_data, all_nb_excl.nb12_excl_list, ic, warp_size_zu);
  ic = nb13_exclusion_bounds.putHost(&int_data, all_nb_excl.nb13_excl_bounds, ic, warp_size_zu);
  ic = nb13_exclusion_list.putHost(&int_data, all_nb_excl.nb13_excl_list, ic, warp_size_zu);
  ic = nb14_exclusion_bounds.putHost(&int_data, all_nb_excl.nb14_excl_bounds, ic, warp_size_zu);
  ic = nb14_exclusion_list.putHost(&int_data, all_nb_excl.nb14_excl_list, ic, warp_size_zu);
  ic = infr14_i_atoms.putHost(&int_data, tmp_inferred_14_i_atoms, ic, warp_size_zu);
  ic = infr14_j_atoms.putHost(&int_data, tmp_inferred_14_j_atoms, ic, warp_size_zu);
  ic = infr14_parameter_indices.putHost(&int_data, tmp_inferred_14_param_idx, ic, warp_size_zu);
  ic = neck_gb_indices.putHost(&int_data, tmp_neck_gb_indices, ic, warp_size_zu);
  ic = tree_joining_info.putHost(&int_data, tmp_tree_joining_info, ic, warp_size_zu);
  ic = last_rotator_info.putHost(&int_data, tmp_last_rotator_info, ic, warp_size_zu);

  // Do the same for double Hybrid POINTER-kind objects
  size_t dc = atomic_charges.putHost(&double_data, tmp_charges, 0, warp_size_zu);
  dc = atomic_masses.putHost(&double_data, tmp_masses, dc, warp_size_zu);
  std::vector<double> inv_mass(atom_count, 0.0);
  for (int i = 0; i < atom_count; i++) {
    inv_mass[i] = (tmp_masses[i] > constants::tiny) ? 1.0 / tmp_masses[i] : 0.0;
  }
  dc = inverse_atomic_masses.putHost(&double_data, inv_mass, dc, warp_size_zu);
  dc = urey_bradley_stiffnesses.putHost(&double_data, tmp_ub_stiffnesses, dc, warp_size_zu);
  dc = urey_bradley_equilibria.putHost(&double_data, tmp_ub_equilibria, dc, warp_size_zu);
  dc = charmm_impr_stiffnesses.putHost(&double_data, tmp_charmm_impr_stiffnesses, dc,
                                       warp_size_zu);
  dc = charmm_impr_phase_angles.putHost(&double_data, tmp_charmm_impr_phase_angles, dc,
                                        warp_size_zu);
  dc = cmap_surfaces.putHost(&double_data, tmp_cmap_surfaces, dc, warp_size_zu);
  dc = cmap_phi_derivatives.putHost(&double_data, cmap_table.phi_derivatives, dc, warp_size_zu);
  dc = cmap_psi_derivatives.putHost(&double_data, cmap_table.psi_derivatives, dc, warp_size_zu);
  dc = cmap_phi_psi_derivatives.putHost(&double_data, cmap_table.phi_psi_derivatives, dc,
                                        warp_size_zu);
  dc = cmap_patches.putHost(&double_data, cmap_table.patch_matrix_form, dc, warp_size_zu);
  dc = bond_stiffnesses.putHost(&double_data, tmp_bond_stiffnesses, dc, warp_size_zu);
  dc = bond_equilibria.putHost(&double_data, tmp_bond_equilibria, dc, warp_size_zu);
  dc = angl_stiffnesses.putHost(&double_data, tmp_angl_stiffnesses, dc, warp_size_zu);
  dc = angl_equilibria.putHost(&double_data, tmp_angl_equilibria, dc, warp_size_zu);
  dc = dihe_amplitudes.putHost(&double_data, tmp_dihe_amplitudes, dc, warp_size_zu);
  dc = dihe_periodicities.putHost(&double_data, tmp_dihe_periodicities, dc, warp_size_zu);
  dc = dihe_phase_angles.putHost(&double_data, tmp_dihe_phase_angles, dc, warp_size_zu);
  dc = virtual_site_frame_dim1.putHost(&double_data, vsite_table.frame_dim1, dc, warp_size_zu);
  dc = virtual_site_frame_dim2.putHost(&double_data, vsite_table.frame_dim2, dc, warp_size_zu);
  dc = virtual_site_frame_dim3.putHost(&double_data, vsite_table.frame_dim3, dc, warp_size_zu);
  dc = charge_parameters.putHost(&double_data, tmp_charge_parameters, dc, warp_size_zu);
  dc = lj_a_values.putHost(&double_data, tmp_lj_a_values, dc, warp_size_zu);
  dc = lj_b_values.putHost(&double_data, tmp_lj_b_values, dc, warp_size_zu);
  dc = lj_c_values.putHost(&double_data, tmp_lj_c_values, dc, warp_size_zu);
  dc = lj_14_a_values.putHost(&double_data, tmp_lj_14_a_values, dc, warp_size_zu);
  dc = lj_14_b_values.putHost(&double_data, tmp_lj_14_b_values, dc, warp_size_zu);
  dc = lj_14_c_values.putHost(&double_data, tmp_lj_14_c_values, dc, warp_size_zu);
  dc = lj_type_corrections.putHost(&double_data, std::vector<double>(atom_type_count, 0.0), dc,
                                   warp_size_zu);
  dc = attn14_elec_factors.putHost(&double_data, attn_parm.elec_screening_factors, dc,
                                   warp_size_zu);
  dc = attn14_vdw_factors.putHost(&double_data, attn_parm.vdw_screening_factors, dc, warp_size_zu);
  dc = atomic_pb_radii.putHost(&double_data, tmp_atomic_pb_radii, dc, warp_size_zu);
  dc = gb_screening_factors.putHost(&double_data, tmp_gb_screening_factors, dc, warp_size_zu);
  dc = gb_alpha_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);
  dc = gb_beta_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);
  dc = gb_gamma_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);
  dc = solty_info.putHost(&double_data, tmp_solty_info, dc, warp_size_zu);
  dc = hbond_a_values.putHost(&double_data, tmp_hbond_a_values, dc, warp_size_zu);
  dc = hbond_b_values.putHost(&double_data, tmp_hbond_b_values, dc, warp_size_zu);
  dc = hbond_cutoffs.putHost(&double_data, tmp_hbond_cutoffs, dc, warp_size_zu);

  // Do the same for float Hybrid POINTER-kind objects, using the range constructor to make
  // single-precision vectors out of their dobule-precision counterparts before loading the
  // Hybrid objects.
  const std::vector<float>sp_tmp_charges(tmp_charges.begin(), tmp_charges.end());
  size_t fc = sp_atomic_charges.putHost(&float_data, sp_tmp_charges, 0, warp_size_zu);
  const std::vector<float>sp_tmp_masses(tmp_masses.begin(), tmp_masses.end());
  fc = sp_atomic_masses.putHost(&float_data, sp_tmp_masses, fc, warp_size_zu);
  std::vector<float> sp_inv_mass(atom_count, 0.0);
  for (int i = 0; i < atom_count; i++) {
    sp_inv_mass[i] = (tmp_masses[i] > constants::tiny) ? 1.0 / tmp_masses[i] : 0.0;
  }
  fc = sp_inverse_atomic_masses.putHost(&float_data, sp_inv_mass, fc, warp_size_zu);
  const std::vector<float> sp_tmp_ub_stiffnesses(tmp_ub_stiffnesses.begin(),
                                                 tmp_ub_stiffnesses.end());
  fc = sp_urey_bradley_stiffnesses.putHost(&float_data, sp_tmp_ub_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_ub_equilibria(tmp_ub_equilibria.begin(),
                                                tmp_ub_equilibria.end());
  fc = sp_urey_bradley_equilibria.putHost(&float_data, sp_tmp_ub_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_charmm_impr_stiffnesses(tmp_charmm_impr_stiffnesses.begin(),
                                                          tmp_charmm_impr_stiffnesses.end());
  fc = sp_charmm_impr_stiffnesses.putHost(&float_data, sp_tmp_charmm_impr_stiffnesses, fc,
                                          warp_size_zu);
  const std::vector<float> sp_tmp_charmm_impr_phase_angles(tmp_charmm_impr_phase_angles.begin(),
                                                           tmp_charmm_impr_phase_angles.end());
  fc = sp_charmm_impr_phase_angles.putHost(&float_data, sp_tmp_charmm_impr_phase_angles, fc,
                                           warp_size_zu);
  const std::vector<float> sp_tmp_cmap_surfaces(tmp_cmap_surfaces.begin(),
                                                tmp_cmap_surfaces.end());
  fc = sp_cmap_surfaces.putHost(&float_data, sp_tmp_cmap_surfaces, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dphi(cmap_table.phi_derivatives.begin(),
                                            cmap_table.phi_derivatives.end());
  fc = sp_cmap_phi_derivatives.putHost(&float_data, sp_tmp_cmap_dphi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dpsi(cmap_table.psi_derivatives.begin(),
                                            cmap_table.psi_derivatives.end());
  fc = sp_cmap_psi_derivatives.putHost(&float_data, sp_tmp_cmap_dpsi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dphipsi(cmap_table.phi_psi_derivatives.begin(),
                                               cmap_table.phi_psi_derivatives.end());
  fc = sp_cmap_phi_psi_derivatives.putHost(&float_data, sp_tmp_cmap_dphipsi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_patches(cmap_table.patch_matrix_form.begin(),
                                               cmap_table.patch_matrix_form.end());
  fc = sp_cmap_patches.putHost(&float_data, sp_tmp_cmap_patches, fc, warp_size_zu);
  const std::vector<float> sp_tmp_bond_stiffnesses(tmp_bond_stiffnesses.begin(),
                                                   tmp_bond_stiffnesses.end());
  fc = sp_bond_stiffnesses.putHost(&float_data, sp_tmp_bond_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_bond_equilibria(tmp_bond_equilibria.begin(),
                                                  tmp_bond_equilibria.end());
  fc = sp_bond_equilibria.putHost(&float_data, sp_tmp_bond_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_angl_stiffnesses(tmp_angl_stiffnesses.begin(),
                                                   tmp_angl_stiffnesses.end());
  fc = sp_angl_stiffnesses.putHost(&float_data, sp_tmp_angl_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_angl_equilibria(tmp_angl_equilibria.begin(),
                                                  tmp_angl_equilibria.end());
  fc = sp_angl_equilibria.putHost(&float_data, sp_tmp_angl_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_amplitudes(tmp_dihe_amplitudes.begin(),
                                                  tmp_dihe_amplitudes.end());
  fc = sp_dihe_amplitudes.putHost(&float_data, sp_tmp_dihe_amplitudes, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_periodicities(tmp_dihe_periodicities.begin(),
                                                     tmp_dihe_periodicities.end());
  fc = sp_dihe_periodicities.putHost(&float_data, sp_tmp_dihe_periodicities, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_phase_angles(tmp_dihe_phase_angles.begin(),
                                                    tmp_dihe_phase_angles.end());
  fc = sp_dihe_phase_angles.putHost(&float_data, sp_tmp_dihe_phase_angles, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim1(vsite_table.frame_dim1.begin(),
                                         vsite_table.frame_dim1.end());
  fc = sp_virtual_site_frame_dim1.putHost(&float_data, sp_frame_dim1, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim2(vsite_table.frame_dim2.begin(),
                                         vsite_table.frame_dim2.end());
  fc = sp_virtual_site_frame_dim2.putHost(&float_data, sp_frame_dim2, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim3(vsite_table.frame_dim3.begin(),
                                         vsite_table.frame_dim3.end());
  fc = sp_virtual_site_frame_dim3.putHost(&float_data, sp_frame_dim3, fc, warp_size_zu);
  const std::vector<float> sp_tmp_charge_parameters(tmp_charge_parameters.begin(),
                                                    tmp_charge_parameters.end());
  fc = sp_charge_parameters.putHost(&float_data, sp_tmp_charge_parameters, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_a_values(tmp_lj_a_values.begin(), tmp_lj_a_values.end());
  fc = sp_lj_a_values.putHost(&float_data, sp_tmp_lj_a_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_b_values(tmp_lj_b_values.begin(), tmp_lj_b_values.end());
  fc = sp_lj_b_values.putHost(&float_data, sp_tmp_lj_b_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_c_values(tmp_lj_c_values.begin(), tmp_lj_c_values.end());
  fc = sp_lj_c_values.putHost(&float_data, sp_tmp_lj_c_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_a_values(tmp_lj_14_a_values.begin(),
                                                 tmp_lj_14_a_values.end());
  fc = sp_lj_14_a_values.putHost(&float_data, sp_tmp_lj_14_a_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_b_values(tmp_lj_14_b_values.begin(),
                                                 tmp_lj_14_b_values.end());
  fc = sp_lj_14_b_values.putHost(&float_data, sp_tmp_lj_14_b_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_c_values(tmp_lj_14_c_values.begin(),
                                                 tmp_lj_14_c_values.end());
  fc = sp_lj_14_c_values.putHost(&float_data, sp_tmp_lj_14_c_values, fc, warp_size_zu);
  fc = sp_lj_type_corrections.putHost(&float_data, std::vector<float>(atom_type_count, 0.0), fc,
                                      warp_size_zu);
  const std::vector<float> sp_tmp_elec_screening_factors(attn_parm.elec_screening_factors.begin(),
                                                         attn_parm.elec_screening_factors.end());
  fc = sp_attn14_elec_factors.putHost(&float_data, sp_tmp_elec_screening_factors, fc,
                                      warp_size_zu);
  const std::vector<float> sp_tmp_vdw_screening_factors(attn_parm.vdw_screening_factors.begin(),
                                                        attn_parm.vdw_screening_factors.end());
  fc = sp_attn14_vdw_factors.putHost(&float_data, sp_tmp_vdw_screening_factors, fc, warp_size_zu);
  const std::vector<float> sp_tmp_atomic_pb_radii(tmp_atomic_pb_radii.begin(),
                                                  tmp_atomic_pb_radii.end());
  fc = sp_atomic_pb_radii.putHost(&float_data, sp_tmp_atomic_pb_radii, fc, warp_size_zu);
  const std::vector<float> sp_tmp_gb_screening_factors(tmp_gb_screening_factors.begin(),
                                                       tmp_gb_screening_factors.end());
  fc = sp_gb_screening_factors.putHost(&float_data, sp_tmp_gb_screening_factors, fc, warp_size_zu);
  const std::vector<float> sp_tmp_gb_coef(tmp_gb_coef.begin(), tmp_gb_coef.end());
  fc = sp_gb_alpha_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);
  fc = sp_gb_beta_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);
  fc = sp_gb_gamma_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);

  // Do the same for char4 Hybrid POINTER-kind objects
  size_t c4c = atom_names.putHost(&char4_data, tmp_atom_names, 0, warp_size_zu);
  c4c = atom_types.putHost(&char4_data, tmp_atom_types, c4c, warp_size_zu);
  c4c = residue_names.putHost(&char4_data, tmp_residue_names, c4c, warp_size_zu);
  c4c = bond_modifiers.putHost(&char4_data, basic_vtable.bond_mods, c4c, warp_size_zu);
  c4c = angl_modifiers.putHost(&char4_data, basic_vtable.angl_mods, c4c, warp_size_zu);
  c4c = dihe_modifiers.putHost(&char4_data, basic_vtable.dihe_mods, c4c, warp_size_zu);
  c4c = bond_assigned_mods.putHost(&char4_data, basic_vtable.bond_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = angl_assigned_mods.putHost(&char4_data, basic_vtable.angl_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = dihe_assigned_mods.putHost(&char4_data, basic_vtable.dihe_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = tree_symbols.putHost(&char4_data, tmp_tree_symbols, c4c, warp_size_zu);
  
  // The Amber topology read here does not contain overflow names of any sort
  const std::vector<char4> blank_overflow;
  c4c = atom_overflow_names.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
  c4c = atom_overflow_types.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
  c4c = residue_overflow_names.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::rebasePointers() {
  
  // Repair the integer data POINTER-kind Hybrids first.
  descriptors.swapTarget(&int_data);
  residue_limits.swapTarget(&int_data);
  atom_struc_numbers.swapTarget(&int_data);
  residue_numbers.swapTarget(&int_data);
  molecule_limits.swapTarget(&int_data);
  atomic_numbers.swapTarget(&int_data);
  mobile_atoms.swapTarget(&int_data);
  molecule_membership.swapTarget(&int_data);
  molecule_contents.swapTarget(&int_data);
  urey_bradley_i_atoms.swapTarget(&int_data);
  urey_bradley_k_atoms.swapTarget(&int_data);
  urey_bradley_parameter_indices.swapTarget(&int_data);
  urey_bradley_assigned_atoms.swapTarget(&int_data);
  urey_bradley_assigned_index.swapTarget(&int_data);
  urey_bradley_assigned_terms.swapTarget(&int_data);
  urey_bradley_assigned_bounds.swapTarget(&int_data);
  charmm_impr_i_atoms.swapTarget(&int_data);
  charmm_impr_j_atoms.swapTarget(&int_data);
  charmm_impr_k_atoms.swapTarget(&int_data);
  charmm_impr_l_atoms.swapTarget(&int_data);
  charmm_impr_parameter_indices.swapTarget(&int_data);
  charmm_impr_assigned_atoms.swapTarget(&int_data);
  charmm_impr_assigned_index.swapTarget(&int_data);
  charmm_impr_assigned_terms.swapTarget(&int_data);
  charmm_impr_assigned_bounds.swapTarget(&int_data);
  cmap_i_atoms.swapTarget(&int_data);
  cmap_j_atoms.swapTarget(&int_data);
  cmap_k_atoms.swapTarget(&int_data);
  cmap_l_atoms.swapTarget(&int_data);
  cmap_m_atoms.swapTarget(&int_data);
  cmap_surface_dimensions.swapTarget(&int_data);
  cmap_surface_bounds.swapTarget(&int_data);
  cmap_patch_bounds.swapTarget(&int_data);
  cmap_surface_indices.swapTarget(&int_data);
  cmap_assigned_atoms.swapTarget(&int_data);
  cmap_assigned_index.swapTarget(&int_data);
  cmap_assigned_terms.swapTarget(&int_data);
  cmap_assigned_bounds.swapTarget(&int_data);
  bond_i_atoms.swapTarget(&int_data);
  bond_j_atoms.swapTarget(&int_data);
  bond_parameter_indices.swapTarget(&int_data);
  bond_assigned_atoms.swapTarget(&int_data);
  bond_assigned_index.swapTarget(&int_data);
  bond_assigned_terms.swapTarget(&int_data);
  bond_assigned_bounds.swapTarget(&int_data);
  angl_i_atoms.swapTarget(&int_data);
  angl_j_atoms.swapTarget(&int_data);
  angl_k_atoms.swapTarget(&int_data);
  angl_parameter_indices.swapTarget(&int_data);
  angl_assigned_atoms.swapTarget(&int_data);
  angl_assigned_index.swapTarget(&int_data);
  angl_assigned_terms.swapTarget(&int_data);
  angl_assigned_bounds.swapTarget(&int_data);
  dihe_i_atoms.swapTarget(&int_data);
  dihe_j_atoms.swapTarget(&int_data);
  dihe_k_atoms.swapTarget(&int_data);
  dihe_l_atoms.swapTarget(&int_data);
  dihe_parameter_indices.swapTarget(&int_data);
  dihe14_parameter_indices.swapTarget(&int_data);
  dihe_assigned_atoms.swapTarget(&int_data);
  dihe_assigned_index.swapTarget(&int_data);
  dihe_assigned_terms.swapTarget(&int_data);
  dihe_assigned_bounds.swapTarget(&int_data);
  virtual_site_atoms.swapTarget(&int_data);
  virtual_site_frame_types.swapTarget(&int_data);
  virtual_site_frame1_atoms.swapTarget(&int_data);
  virtual_site_frame2_atoms.swapTarget(&int_data);
  virtual_site_frame3_atoms.swapTarget(&int_data);
  virtual_site_frame4_atoms.swapTarget(&int_data);
  charge_indices.swapTarget(&int_data);
  lennard_jones_indices.swapTarget(&int_data);
  atom_exclusion_bounds.swapTarget(&int_data);
  atom_exclusion_list.swapTarget(&int_data);
  nb11_exclusion_bounds.swapTarget(&int_data);
  nb11_exclusion_list.swapTarget(&int_data);
  nb12_exclusion_bounds.swapTarget(&int_data);
  nb12_exclusion_list.swapTarget(&int_data);
  nb13_exclusion_bounds.swapTarget(&int_data);
  nb13_exclusion_list.swapTarget(&int_data);
  nb14_exclusion_bounds.swapTarget(&int_data);
  nb14_exclusion_list.swapTarget(&int_data);
  infr14_i_atoms.swapTarget(&int_data);
  infr14_j_atoms.swapTarget(&int_data);
  infr14_parameter_indices.swapTarget(&int_data);
  neck_gb_indices.swapTarget(&int_data);
  tree_joining_info.swapTarget(&int_data);
  last_rotator_info.swapTarget(&int_data);

  // Repair the double-precision real POINTER-kind Hybrid objects
  atomic_charges.swapTarget(&double_data);
  atomic_masses.swapTarget(&double_data);
  inverse_atomic_masses.swapTarget(&double_data);
  urey_bradley_stiffnesses.swapTarget(&double_data);
  urey_bradley_equilibria.swapTarget(&double_data);
  charmm_impr_stiffnesses.swapTarget(&double_data);
  charmm_impr_phase_angles.swapTarget(&double_data);
  cmap_surfaces.swapTarget(&double_data);
  cmap_phi_derivatives.swapTarget(&double_data);
  cmap_psi_derivatives.swapTarget(&double_data);
  cmap_phi_psi_derivatives.swapTarget(&double_data);
  cmap_patches.swapTarget(&double_data);
  bond_stiffnesses.swapTarget(&double_data);
  bond_equilibria.swapTarget(&double_data);
  angl_stiffnesses.swapTarget(&double_data);
  angl_equilibria.swapTarget(&double_data);
  dihe_amplitudes.swapTarget(&double_data);
  dihe_periodicities.swapTarget(&double_data);
  dihe_phase_angles.swapTarget(&double_data);
  virtual_site_frame_dim1.swapTarget(&double_data);
  virtual_site_frame_dim2.swapTarget(&double_data);
  virtual_site_frame_dim3.swapTarget(&double_data);
  charge_parameters.swapTarget(&double_data);
  lj_a_values.swapTarget(&double_data);
  lj_b_values.swapTarget(&double_data);
  lj_c_values.swapTarget(&double_data);
  lj_14_a_values.swapTarget(&double_data);
  lj_14_b_values.swapTarget(&double_data);
  lj_14_c_values.swapTarget(&double_data);
  lj_type_corrections.swapTarget(&double_data);
  attn14_elec_factors.swapTarget(&double_data);
  attn14_vdw_factors.swapTarget(&double_data);
  atomic_pb_radii.swapTarget(&double_data);
  gb_screening_factors.swapTarget(&double_data);
  gb_alpha_parameters.swapTarget(&double_data);
  gb_beta_parameters.swapTarget(&double_data);
  gb_gamma_parameters.swapTarget(&double_data);
  solty_info.swapTarget(&double_data);
  hbond_a_values.swapTarget(&double_data);
  hbond_b_values.swapTarget(&double_data);
  hbond_cutoffs.swapTarget(&double_data);

  // Repair the single-precision real POINTER-kind Hybrid objects
  sp_atomic_charges.swapTarget(&float_data);
  sp_atomic_masses.swapTarget(&float_data);
  sp_inverse_atomic_masses.swapTarget(&float_data);
  sp_urey_bradley_stiffnesses.swapTarget(&float_data);
  sp_urey_bradley_equilibria.swapTarget(&float_data);
  sp_charmm_impr_stiffnesses.swapTarget(&float_data);
  sp_charmm_impr_phase_angles.swapTarget(&float_data);
  sp_cmap_surfaces.swapTarget(&float_data);
  sp_cmap_phi_derivatives.swapTarget(&float_data);
  sp_cmap_psi_derivatives.swapTarget(&float_data);
  sp_cmap_phi_psi_derivatives.swapTarget(&float_data);
  sp_cmap_patches.swapTarget(&float_data);
  sp_bond_stiffnesses.swapTarget(&float_data);
  sp_bond_equilibria.swapTarget(&float_data);
  sp_angl_stiffnesses.swapTarget(&float_data);
  sp_angl_equilibria.swapTarget(&float_data);
  sp_dihe_amplitudes.swapTarget(&float_data);
  sp_dihe_periodicities.swapTarget(&float_data);
  sp_dihe_phase_angles.swapTarget(&float_data);
  sp_virtual_site_frame_dim1.swapTarget(&float_data);
  sp_virtual_site_frame_dim2.swapTarget(&float_data);
  sp_virtual_site_frame_dim3.swapTarget(&float_data);
  sp_charge_parameters.swapTarget(&float_data);
  sp_lj_a_values.swapTarget(&float_data);
  sp_lj_b_values.swapTarget(&float_data);
  sp_lj_c_values.swapTarget(&float_data);
  sp_lj_14_a_values.swapTarget(&float_data);
  sp_lj_14_b_values.swapTarget(&float_data);
  sp_lj_14_c_values.swapTarget(&float_data);
  sp_lj_type_corrections.swapTarget(&float_data);
  sp_attn14_elec_factors.swapTarget(&float_data);
  sp_attn14_vdw_factors.swapTarget(&float_data);
  sp_atomic_pb_radii.swapTarget(&float_data);
  sp_gb_screening_factors.swapTarget(&float_data);
  sp_gb_alpha_parameters.swapTarget(&float_data);
  sp_gb_beta_parameters.swapTarget(&float_data);
  sp_gb_gamma_parameters.swapTarget(&float_data);

  // Repair the char4 POINTER-kind Hybrid objects
  atom_names.swapTarget(&char4_data);
  atom_types.swapTarget(&char4_data);
  residue_names.swapTarget(&char4_data);
  bond_modifiers.swapTarget(&char4_data);
  angl_modifiers.swapTarget(&char4_data);
  dihe_modifiers.swapTarget(&char4_data);
  bond_assigned_mods.swapTarget(&char4_data);
  angl_assigned_mods.swapTarget(&char4_data);
  dihe_assigned_mods.swapTarget(&char4_data);
  atom_overflow_names.swapTarget(&char4_data);
  atom_overflow_types.swapTarget(&char4_data);
  residue_overflow_names.swapTarget(&char4_data);
  tree_symbols.swapTarget(&char4_data);
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &original) :
    version_stamp{},
    date{original.date},
    title{original.title},
    source{original.source},
    force_fields{original.force_fields},
    atom_count{original.atom_count},
    residue_count{original.residue_count},
    molecule_count{original.molecule_count},
    largest_residue_size{original.largest_residue_size},
    last_solute_residue{original.last_solute_residue},
    last_solute_atom{original.last_solute_atom},
    first_solvent_molecule{original.first_solvent_molecule},
    last_atom_before_cap{original.last_atom_before_cap},
    implicit_copy_count{original.implicit_copy_count},
    largest_molecule_size{original.largest_molecule_size},
    descriptors{original.descriptors},
    residue_limits{original.residue_limits},
    atom_struc_numbers{original.atom_struc_numbers},
    residue_numbers{original.residue_numbers},
    molecule_limits{original.molecule_limits},
    atomic_numbers{original.atomic_numbers},
    mobile_atoms{original.mobile_atoms},
    molecule_membership{original.molecule_membership},
    molecule_contents{original.molecule_contents},
    atomic_charges{original.atomic_charges},
    atomic_masses{original.atomic_masses},
    inverse_atomic_masses{original.inverse_atomic_masses},
    sp_atomic_charges{original.sp_atomic_charges},
    sp_atomic_masses{original.sp_atomic_masses},
    sp_inverse_atomic_masses{original.sp_inverse_atomic_masses},
    atom_names{original.atom_names},
    atom_types{original.atom_types},
    residue_names{original.residue_names},
    urey_bradley_term_count{original.urey_bradley_term_count},
    charmm_impr_term_count{original.charmm_impr_term_count},
    cmap_term_count{original.cmap_term_count},
    urey_bradley_parameter_count{original.urey_bradley_parameter_count},
    charmm_impr_parameter_count{original.charmm_impr_parameter_count},
    cmap_surface_count{original.cmap_surface_count},
    urey_bradley_pert_term_count{original.urey_bradley_pert_term_count},
    charmm_impr_pert_term_count{original.charmm_impr_pert_term_count},
    cmap_pert_term_count{original.cmap_pert_term_count},
    urey_bradleys_in_perturbed_group{original.urey_bradleys_in_perturbed_group},
    charmm_imprs_in_perturbed_group{original.charmm_imprs_in_perturbed_group},
    cmaps_in_perturbed_group{original.cmaps_in_perturbed_group},
    urey_bradley_i_atoms{original.urey_bradley_i_atoms},
    urey_bradley_k_atoms{original.urey_bradley_k_atoms},
    urey_bradley_parameter_indices{original.urey_bradley_parameter_indices},
    urey_bradley_assigned_atoms{original.urey_bradley_assigned_atoms},
    urey_bradley_assigned_index{original.urey_bradley_assigned_index},
    urey_bradley_assigned_terms{original.urey_bradley_assigned_terms},
    urey_bradley_assigned_bounds{original.urey_bradley_assigned_bounds},
    charmm_impr_i_atoms{original.charmm_impr_i_atoms},
    charmm_impr_j_atoms{original.charmm_impr_j_atoms},
    charmm_impr_k_atoms{original.charmm_impr_k_atoms},
    charmm_impr_l_atoms{original.charmm_impr_l_atoms},
    charmm_impr_parameter_indices{original.charmm_impr_parameter_indices},
    charmm_impr_assigned_atoms{original.charmm_impr_assigned_atoms},
    charmm_impr_assigned_index{original.charmm_impr_assigned_index},
    charmm_impr_assigned_terms{original.charmm_impr_assigned_terms},
    charmm_impr_assigned_bounds{original.charmm_impr_assigned_bounds},
    cmap_i_atoms{original.cmap_i_atoms},
    cmap_j_atoms{original.cmap_j_atoms},
    cmap_k_atoms{original.cmap_k_atoms},
    cmap_l_atoms{original.cmap_l_atoms},
    cmap_m_atoms{original.cmap_m_atoms},
    cmap_surface_dimensions{original.cmap_surface_dimensions},
    cmap_surface_bounds{original.cmap_surface_bounds},
    cmap_patch_bounds{original.cmap_patch_bounds},
    cmap_surface_indices{original.cmap_surface_indices},
    cmap_assigned_atoms{original.cmap_assigned_atoms},
    cmap_assigned_index{original.cmap_assigned_index},
    cmap_assigned_terms{original.cmap_assigned_terms},
    cmap_assigned_bounds{original.cmap_assigned_bounds},
    urey_bradley_stiffnesses{original.urey_bradley_stiffnesses},
    urey_bradley_equilibria{original.urey_bradley_equilibria},
    charmm_impr_stiffnesses{original.charmm_impr_stiffnesses},
    charmm_impr_phase_angles{original.charmm_impr_phase_angles},
    cmap_surfaces{original.cmap_surfaces},
    cmap_phi_derivatives{original.cmap_phi_derivatives},
    cmap_psi_derivatives{original.cmap_psi_derivatives},
    cmap_phi_psi_derivatives{original.cmap_phi_psi_derivatives},
    cmap_patches{original.cmap_patches},
    sp_urey_bradley_stiffnesses{original.sp_urey_bradley_stiffnesses},
    sp_urey_bradley_equilibria{original.sp_urey_bradley_equilibria},
    sp_charmm_impr_stiffnesses{original.sp_charmm_impr_stiffnesses},
    sp_charmm_impr_phase_angles{original.sp_charmm_impr_phase_angles},
    sp_cmap_surfaces{original.sp_cmap_surfaces},
    sp_cmap_phi_derivatives{original.sp_cmap_phi_derivatives},
    sp_cmap_psi_derivatives{original.sp_cmap_psi_derivatives},
    sp_cmap_phi_psi_derivatives{original.sp_cmap_phi_psi_derivatives},
    sp_cmap_patches{original.sp_cmap_patches},
    bond_term_with_hydrogen{original.bond_term_with_hydrogen},
    angl_term_with_hydrogen{original.angl_term_with_hydrogen},
    dihe_term_with_hydrogen{original.dihe_term_with_hydrogen},
    bond_term_without_hydrogen{original.bond_term_without_hydrogen},
    angl_term_without_hydrogen{original.angl_term_without_hydrogen},
    dihe_term_without_hydrogen{original.dihe_term_without_hydrogen},
    bond_term_count{original.bond_term_count},
    angl_term_count{original.angl_term_count},
    dihe_term_count{original.dihe_term_count},
    bond_parameter_count{original.bond_parameter_count},
    angl_parameter_count{original.angl_parameter_count},
    dihe_parameter_count{original.dihe_parameter_count},
    bond_perturbation_term_count{original.bond_perturbation_term_count},
    angl_perturbation_term_count{original.angl_perturbation_term_count},
    dihe_perturbation_term_count{original.dihe_perturbation_term_count},
    bonds_in_perturbed_group{original.bonds_in_perturbed_group},
    angls_in_perturbed_group{original.angls_in_perturbed_group},
    dihes_in_perturbed_group{original.dihes_in_perturbed_group},
    bonded_group_count{original.bonded_group_count},
    bond_stiffnesses{original.bond_stiffnesses},
    bond_equilibria{original.bond_equilibria},
    angl_stiffnesses{original.angl_stiffnesses},
    angl_equilibria{original.angl_equilibria},
    dihe_amplitudes{original.dihe_amplitudes},
    dihe_periodicities{original.dihe_periodicities},
    dihe_phase_angles{original.dihe_phase_angles},
    sp_bond_stiffnesses{original.sp_bond_stiffnesses},
    sp_bond_equilibria{original.sp_bond_equilibria},
    sp_angl_stiffnesses{original.sp_angl_stiffnesses},
    sp_angl_equilibria{original.sp_angl_equilibria},
    sp_dihe_amplitudes{original.sp_dihe_amplitudes},
    sp_dihe_periodicities{original.sp_dihe_periodicities},
    sp_dihe_phase_angles{original.sp_dihe_phase_angles},
    bond_i_atoms{original.bond_i_atoms},
    bond_j_atoms{original.bond_j_atoms},
    bond_parameter_indices{original.bond_parameter_indices},
    bond_assigned_atoms{original.bond_assigned_atoms},
    bond_assigned_index{original.bond_assigned_index},
    bond_assigned_terms{original.bond_assigned_terms},
    bond_assigned_bounds{original.bond_assigned_bounds},
    angl_i_atoms{original.angl_i_atoms},
    angl_j_atoms{original.angl_j_atoms},
    angl_k_atoms{original.angl_k_atoms},
    angl_parameter_indices{original.angl_parameter_indices},
    angl_assigned_atoms{original.angl_assigned_atoms},
    angl_assigned_index{original.angl_assigned_index},
    angl_assigned_terms{original.angl_assigned_terms},
    angl_assigned_bounds{original.angl_assigned_bounds},
    dihe_i_atoms{original.dihe_i_atoms},
    dihe_j_atoms{original.dihe_j_atoms},
    dihe_k_atoms{original.dihe_k_atoms},
    dihe_l_atoms{original.dihe_l_atoms},
    dihe_parameter_indices{original.dihe_parameter_indices},
    dihe14_parameter_indices{original.dihe14_parameter_indices},
    dihe_assigned_atoms{original.dihe_assigned_atoms},
    dihe_assigned_index{original.dihe_assigned_index},
    dihe_assigned_terms{original.dihe_assigned_terms},
    dihe_assigned_bounds{original.dihe_assigned_bounds},
    bond_modifiers{original.bond_modifiers},
    angl_modifiers{original.angl_modifiers},
    dihe_modifiers{original.dihe_modifiers},
    bond_assigned_mods{original.bond_assigned_mods},
    angl_assigned_mods{original.angl_assigned_mods},
    dihe_assigned_mods{original.dihe_assigned_mods},
    virtual_site_count{original.virtual_site_count},
    virtual_site_atoms{original.virtual_site_atoms},
    virtual_site_frame_types{original.virtual_site_frame_types},
    virtual_site_frame1_atoms{original.virtual_site_frame1_atoms},
    virtual_site_frame2_atoms{original.virtual_site_frame2_atoms},
    virtual_site_frame3_atoms{original.virtual_site_frame3_atoms},
    virtual_site_frame4_atoms{original.virtual_site_frame4_atoms},
    virtual_site_frame_dim1{original.virtual_site_frame_dim1},
    virtual_site_frame_dim2{original.virtual_site_frame_dim2},
    virtual_site_frame_dim3{original.virtual_site_frame_dim3},
    sp_virtual_site_frame_dim1{original.sp_virtual_site_frame_dim1},
    sp_virtual_site_frame_dim2{original.sp_virtual_site_frame_dim2},
    sp_virtual_site_frame_dim3{original.sp_virtual_site_frame_dim3},
    charge_type_count{original.charge_type_count},
    atom_type_count{original.atom_type_count},
    total_exclusions{original.total_exclusions},
    attenuated_14_type_count{original.attenuated_14_type_count},
    inferred_14_attenuations{original.inferred_14_attenuations},
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    salt_concentration{original.salt_concentration},
    coulomb_constant{original.coulomb_constant},
    pb_radii_set{original.pb_radii_set},
    charge_indices{original.charge_indices},
    lennard_jones_indices{original.lennard_jones_indices},
    atom_exclusion_bounds{original.atom_exclusion_bounds},
    atom_exclusion_list{original.atom_exclusion_list},
    nb11_exclusion_bounds{original.nb11_exclusion_bounds},
    nb11_exclusion_list{original.nb11_exclusion_list},
    nb12_exclusion_bounds{original.nb12_exclusion_bounds},
    nb12_exclusion_list{original.nb12_exclusion_list},
    nb13_exclusion_bounds{original.nb13_exclusion_bounds},
    nb13_exclusion_list{original.nb13_exclusion_list},
    nb14_exclusion_bounds{original.nb14_exclusion_bounds},
    nb14_exclusion_list{original.nb14_exclusion_list},
    infr14_i_atoms{original.infr14_i_atoms},
    infr14_j_atoms{original.infr14_j_atoms},
    infr14_parameter_indices{original.infr14_parameter_indices},
    neck_gb_indices{original.neck_gb_indices},
    charge_parameters{original.charge_parameters},
    lj_a_values{original.lj_a_values},
    lj_b_values{original.lj_b_values},
    lj_c_values{original.lj_c_values},
    lj_14_a_values{original.lj_14_a_values},
    lj_14_b_values{original.lj_14_b_values},
    lj_14_c_values{original.lj_14_c_values},
    lj_type_corrections{original.lj_type_corrections},
    attn14_elec_factors{original.attn14_elec_factors},
    attn14_vdw_factors{original.attn14_vdw_factors},
    atomic_pb_radii{original.atomic_pb_radii},
    gb_screening_factors{original.gb_screening_factors},
    gb_alpha_parameters{original.gb_alpha_parameters},
    gb_beta_parameters{original.gb_beta_parameters},
    gb_gamma_parameters{original.gb_gamma_parameters},
    sp_charge_parameters{original.sp_charge_parameters},
    sp_lj_a_values{original.sp_lj_a_values},
    sp_lj_b_values{original.sp_lj_b_values},
    sp_lj_c_values{original.sp_lj_c_values},
    sp_lj_14_a_values{original.sp_lj_14_a_values},
    sp_lj_14_b_values{original.sp_lj_14_b_values},
    sp_lj_14_c_values{original.sp_lj_14_c_values},
    sp_lj_type_corrections{original.sp_lj_type_corrections},
    sp_attn14_elec_factors{original.sp_attn14_elec_factors},
    sp_attn14_vdw_factors{original.sp_attn14_vdw_factors},
    sp_atomic_pb_radii{original.sp_atomic_pb_radii},
    sp_gb_screening_factors{original.sp_gb_screening_factors},
    sp_gb_alpha_parameters{original.sp_gb_alpha_parameters},
    sp_gb_beta_parameters{original.sp_gb_beta_parameters},
    sp_gb_gamma_parameters{original.sp_gb_gamma_parameters},
    use_bond_constraints{original.use_bond_constraints},
    use_settle{original.use_settle},
    use_perturbation_info{original.use_perturbation_info},
    use_solvent_cap_option{original.use_solvent_cap_option},
    use_polarization{original.use_polarization},
    water_residue_name{original.water_residue_name},
    bond_constraint_mask{original.bond_constraint_mask},
    bond_constraint_omit_mask{original.bond_constraint_omit_mask},
    rigid_water_count{original.rigid_water_count},
    bond_constraint_count{original.bond_constraint_count},
    degrees_of_freedom{original.degrees_of_freedom},
    nonrigid_particle_count{original.nonrigid_particle_count},
    atom_overflow_names{original.atom_overflow_names},
    atom_overflow_types{original.atom_overflow_types},
    residue_overflow_names{original.residue_overflow_names},
    unused_nhparm{original.unused_nhparm},
    unused_nparm{original.unused_nparm},
    unused_natyp{original.unused_natyp},
    hbond_10_12_parameter_count{original.hbond_10_12_parameter_count},
    heavy_bonds_plus_constraints{original.heavy_bonds_plus_constraints},
    heavy_angls_plus_constraints{original.heavy_angls_plus_constraints},
    heavy_dihes_plus_constraints{original.heavy_dihes_plus_constraints},
    tree_joining_info{original.tree_joining_info},
    last_rotator_info{original.last_rotator_info},
    solty_info{original.solty_info},
    hbond_a_values{original.hbond_a_values},
    hbond_b_values{original.hbond_b_values},
    hbond_cutoffs{original.hbond_cutoffs},
    tree_symbols{original.tree_symbols},
    int_data{original.int_data},
    double_data{original.double_data},
    float_data{original.float_data},
    char4_data{original.char4_data}
{
  const int nvst_char = strlen(original.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = original.version_stamp[i];
  }

  // Repair pointers.  It is tedious to write the initializer list in the way that was done, and
  // it adds insult to injury to now have to list out all of the pointers yet again.  However, this
  // obviates the need for additional depth within the topology object to differentiate its trivial
  // and more complex member variables, so best to accept the tedium here for the sake of
  // simplicity elsewhere.
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
AtomGraph& AtomGraph::operator=(const AtomGraph &other) {
  
  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy the elements of the original 
  const int nvst_char = strlen(other.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = other.version_stamp[i];
  }
  date = other.date;
  title = other.title;
  source = other.source;
  force_fields = other.force_fields;
  atom_count = other.atom_count;
  residue_count = other.residue_count;
  molecule_count = other.molecule_count;
  largest_residue_size = other.largest_residue_size;
  last_solute_residue = other.last_solute_residue;
  last_solute_atom = other.last_solute_atom;
  first_solvent_molecule = other.first_solvent_molecule;
  last_atom_before_cap = other.last_atom_before_cap;
  implicit_copy_count = other.implicit_copy_count;
  largest_molecule_size = other.largest_molecule_size;
  descriptors = other.descriptors;
  residue_limits = other.residue_limits;
  atom_struc_numbers = other.atom_struc_numbers;
  residue_numbers = other.residue_numbers;
  molecule_limits = other.molecule_limits;
  atomic_numbers = other.atomic_numbers;
  mobile_atoms = other.mobile_atoms;
  molecule_membership = other.molecule_membership;
  molecule_contents = other.molecule_contents;
  atomic_charges = other.atomic_charges;
  atomic_masses = other.atomic_masses;
  inverse_atomic_masses = other.inverse_atomic_masses;
  sp_atomic_charges = other.sp_atomic_charges;
  sp_atomic_masses = other.sp_atomic_masses;
  sp_inverse_atomic_masses = other.sp_inverse_atomic_masses;
  atom_names = other.atom_names;
  atom_types = other.atom_types;
  residue_names = other.residue_names;
  urey_bradley_term_count = other.urey_bradley_term_count;
  charmm_impr_term_count = other.charmm_impr_term_count;
  cmap_term_count = other.cmap_term_count;
  urey_bradley_parameter_count = other.urey_bradley_parameter_count;
  charmm_impr_parameter_count = other.charmm_impr_parameter_count;
  cmap_surface_count = other.cmap_surface_count;
  urey_bradley_pert_term_count = other.urey_bradley_pert_term_count;
  charmm_impr_pert_term_count = other.charmm_impr_pert_term_count;
  cmap_pert_term_count = other.cmap_pert_term_count;
  urey_bradleys_in_perturbed_group = other.urey_bradleys_in_perturbed_group;
  charmm_imprs_in_perturbed_group = other.charmm_imprs_in_perturbed_group;
  cmaps_in_perturbed_group = other.cmaps_in_perturbed_group;
  urey_bradley_i_atoms = other.urey_bradley_i_atoms;
  urey_bradley_k_atoms = other.urey_bradley_k_atoms;
  urey_bradley_parameter_indices = other.urey_bradley_parameter_indices;
  urey_bradley_assigned_atoms = other.urey_bradley_assigned_atoms;
  urey_bradley_assigned_index = other.urey_bradley_assigned_index;
  urey_bradley_assigned_terms = other.urey_bradley_assigned_terms;
  urey_bradley_assigned_bounds = other.urey_bradley_assigned_bounds;
  charmm_impr_i_atoms = other.charmm_impr_i_atoms;
  charmm_impr_j_atoms = other.charmm_impr_j_atoms;
  charmm_impr_k_atoms = other.charmm_impr_k_atoms;
  charmm_impr_l_atoms = other.charmm_impr_l_atoms;
  charmm_impr_parameter_indices = other.charmm_impr_parameter_indices;
  charmm_impr_assigned_atoms = other.charmm_impr_assigned_atoms;
  charmm_impr_assigned_index = other.charmm_impr_assigned_index;
  charmm_impr_assigned_terms = other.charmm_impr_assigned_terms;
  charmm_impr_assigned_bounds = other.charmm_impr_assigned_bounds;
  cmap_i_atoms = other.cmap_i_atoms;
  cmap_j_atoms = other.cmap_j_atoms;
  cmap_k_atoms = other.cmap_k_atoms;
  cmap_l_atoms = other.cmap_l_atoms;
  cmap_m_atoms = other.cmap_m_atoms;
  cmap_surface_dimensions = other.cmap_surface_dimensions;
  cmap_surface_bounds = other.cmap_surface_bounds;
  cmap_patch_bounds = other.cmap_patch_bounds;
  cmap_surface_indices = other.cmap_surface_indices;
  cmap_assigned_atoms = other.cmap_assigned_atoms;
  cmap_assigned_index = other.cmap_assigned_index;
  cmap_assigned_terms = other.cmap_assigned_terms;
  cmap_assigned_bounds = other.cmap_assigned_bounds;
  urey_bradley_stiffnesses = other.urey_bradley_stiffnesses;
  urey_bradley_equilibria = other.urey_bradley_equilibria;
  charmm_impr_stiffnesses = other.charmm_impr_stiffnesses;
  charmm_impr_phase_angles = other.charmm_impr_phase_angles;
  cmap_surfaces = other.cmap_surfaces;
  cmap_phi_derivatives = other.cmap_phi_derivatives;
  cmap_psi_derivatives = other.cmap_psi_derivatives;
  cmap_phi_psi_derivatives = other.cmap_phi_psi_derivatives;
  cmap_patches = other.cmap_patches;
  sp_urey_bradley_stiffnesses = other.sp_urey_bradley_stiffnesses;
  sp_urey_bradley_equilibria = other.sp_urey_bradley_equilibria;
  sp_charmm_impr_stiffnesses = other.sp_charmm_impr_stiffnesses;
  sp_charmm_impr_phase_angles = other.sp_charmm_impr_phase_angles;
  sp_cmap_surfaces = other.sp_cmap_surfaces;
  sp_cmap_phi_derivatives = other.sp_cmap_phi_derivatives;
  sp_cmap_psi_derivatives = other.sp_cmap_psi_derivatives;
  sp_cmap_phi_psi_derivatives = other.sp_cmap_phi_psi_derivatives;
  sp_cmap_patches = other.sp_cmap_patches;
  bond_term_with_hydrogen = other.bond_term_with_hydrogen;
  angl_term_with_hydrogen = other.angl_term_with_hydrogen;
  dihe_term_with_hydrogen = other.dihe_term_with_hydrogen;
  bond_term_without_hydrogen = other.bond_term_without_hydrogen;
  angl_term_without_hydrogen = other.angl_term_without_hydrogen;
  dihe_term_without_hydrogen = other.dihe_term_without_hydrogen;
  bond_term_count = other.bond_term_count;
  angl_term_count = other.angl_term_count;
  dihe_term_count = other.dihe_term_count;
  bond_parameter_count = other.bond_parameter_count;
  angl_parameter_count = other.angl_parameter_count;
  dihe_parameter_count = other.dihe_parameter_count;
  bond_perturbation_term_count = other.bond_perturbation_term_count;
  angl_perturbation_term_count = other.angl_perturbation_term_count;
  dihe_perturbation_term_count = other.dihe_perturbation_term_count;
  bonds_in_perturbed_group = other.bonds_in_perturbed_group;
  angls_in_perturbed_group = other.angls_in_perturbed_group;
  dihes_in_perturbed_group = other.dihes_in_perturbed_group;
  bonded_group_count = other.bonded_group_count;
  bond_stiffnesses = other.bond_stiffnesses;
  bond_equilibria = other.bond_equilibria;
  angl_stiffnesses = other.angl_stiffnesses;
  angl_equilibria = other.angl_equilibria;
  dihe_amplitudes = other.dihe_amplitudes;
  dihe_periodicities = other.dihe_periodicities;
  dihe_phase_angles = other.dihe_phase_angles;
  sp_bond_stiffnesses = other.sp_bond_stiffnesses;
  sp_bond_equilibria = other.sp_bond_equilibria;
  sp_angl_stiffnesses = other.sp_angl_stiffnesses;
  sp_angl_equilibria = other.sp_angl_equilibria;
  sp_dihe_amplitudes = other.sp_dihe_amplitudes;
  sp_dihe_periodicities = other.sp_dihe_periodicities;
  sp_dihe_phase_angles = other.sp_dihe_phase_angles;
  bond_i_atoms = other.bond_i_atoms;
  bond_j_atoms = other.bond_j_atoms;
  bond_parameter_indices = other.bond_parameter_indices;
  bond_assigned_atoms = other.bond_assigned_atoms;
  bond_assigned_index = other.bond_assigned_index;
  bond_assigned_terms = other.bond_assigned_terms;
  bond_assigned_bounds = other.bond_assigned_bounds;
  angl_i_atoms = other.angl_i_atoms;
  angl_j_atoms = other.angl_j_atoms;
  angl_k_atoms = other.angl_k_atoms;
  angl_parameter_indices = other.angl_parameter_indices;
  angl_assigned_atoms = other.angl_assigned_atoms;
  angl_assigned_index = other.angl_assigned_index;
  angl_assigned_terms = other.angl_assigned_terms;
  angl_assigned_bounds = other.angl_assigned_bounds;
  dihe_i_atoms = other.dihe_i_atoms;
  dihe_j_atoms = other.dihe_j_atoms;
  dihe_k_atoms = other.dihe_k_atoms;
  dihe_l_atoms = other.dihe_l_atoms;
  dihe_parameter_indices = other.dihe_parameter_indices;
  dihe14_parameter_indices = other.dihe14_parameter_indices;
  dihe_assigned_atoms = other.dihe_assigned_atoms;
  dihe_assigned_index = other.dihe_assigned_index;
  dihe_assigned_terms = other.dihe_assigned_terms;
  dihe_assigned_bounds = other.dihe_assigned_bounds;
  bond_modifiers = other.bond_modifiers;
  angl_modifiers = other.angl_modifiers;
  dihe_modifiers = other.dihe_modifiers;
  bond_assigned_mods = other.bond_assigned_mods;
  angl_assigned_mods = other.angl_assigned_mods;
  dihe_assigned_mods = other.dihe_assigned_mods;
  virtual_site_count = other.virtual_site_count;
  virtual_site_atoms = other.virtual_site_atoms;
  virtual_site_frame_types = other.virtual_site_frame_types;
  virtual_site_frame1_atoms = other.virtual_site_frame1_atoms;
  virtual_site_frame2_atoms = other.virtual_site_frame2_atoms;
  virtual_site_frame3_atoms = other.virtual_site_frame3_atoms;
  virtual_site_frame4_atoms = other.virtual_site_frame4_atoms;
  virtual_site_frame_dim1 = other.virtual_site_frame_dim1;
  virtual_site_frame_dim2 = other.virtual_site_frame_dim2;
  virtual_site_frame_dim3 = other.virtual_site_frame_dim3;
  sp_virtual_site_frame_dim1 = other.sp_virtual_site_frame_dim1;
  sp_virtual_site_frame_dim2 = other.sp_virtual_site_frame_dim2;
  sp_virtual_site_frame_dim3 = other.sp_virtual_site_frame_dim3;
  charge_type_count = other.charge_type_count;
  atom_type_count = other.atom_type_count;
  total_exclusions = other.total_exclusions;
  attenuated_14_type_count = other.attenuated_14_type_count;
  inferred_14_attenuations = other.inferred_14_attenuations;
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  salt_concentration = other.salt_concentration;
  coulomb_constant = other.coulomb_constant;
  pb_radii_set = other.pb_radii_set;
  charge_indices = other.charge_indices;
  lennard_jones_indices = other.lennard_jones_indices;
  atom_exclusion_bounds = other.atom_exclusion_bounds;
  atom_exclusion_list = other.atom_exclusion_list;
  nb11_exclusion_bounds = other.nb11_exclusion_bounds;
  nb11_exclusion_list = other.nb11_exclusion_list;
  nb12_exclusion_bounds = other.nb12_exclusion_bounds;
  nb12_exclusion_list = other.nb12_exclusion_list;
  nb13_exclusion_bounds = other.nb13_exclusion_bounds;
  nb13_exclusion_list = other.nb13_exclusion_list;
  nb14_exclusion_bounds = other.nb14_exclusion_bounds;
  nb14_exclusion_list = other.nb14_exclusion_list;
  infr14_i_atoms = other.infr14_i_atoms;
  infr14_j_atoms = other.infr14_j_atoms;
  infr14_parameter_indices = other.infr14_parameter_indices;
  neck_gb_indices = other.neck_gb_indices;
  charge_parameters = other.charge_parameters;
  lj_a_values = other.lj_a_values;
  lj_b_values = other.lj_b_values;
  lj_c_values = other.lj_c_values;
  lj_14_a_values = other.lj_14_a_values;
  lj_14_b_values = other.lj_14_b_values;
  lj_14_c_values = other.lj_14_c_values;
  lj_type_corrections = other.lj_type_corrections;
  attn14_elec_factors = other.attn14_elec_factors;
  attn14_vdw_factors = other.attn14_vdw_factors;
  atomic_pb_radii = other.atomic_pb_radii;
  gb_screening_factors = other.gb_screening_factors;
  gb_alpha_parameters = other.gb_alpha_parameters;
  gb_beta_parameters = other.gb_beta_parameters;
  gb_gamma_parameters = other.gb_gamma_parameters;
  sp_charge_parameters = other.sp_charge_parameters;
  sp_lj_a_values = other.sp_lj_a_values;
  sp_lj_b_values = other.sp_lj_b_values;
  sp_lj_c_values = other.sp_lj_c_values;
  sp_lj_14_a_values = other.sp_lj_14_a_values;
  sp_lj_14_b_values = other.sp_lj_14_b_values;
  sp_lj_14_c_values = other.sp_lj_14_c_values;
  sp_lj_type_corrections = other.sp_lj_type_corrections;
  sp_attn14_elec_factors = other.sp_attn14_elec_factors;
  sp_attn14_vdw_factors = other.sp_attn14_vdw_factors;
  sp_atomic_pb_radii = other.sp_atomic_pb_radii;
  sp_gb_screening_factors = other.sp_gb_screening_factors;
  sp_gb_alpha_parameters = other.sp_gb_alpha_parameters;
  sp_gb_beta_parameters = other.sp_gb_beta_parameters;
  sp_gb_gamma_parameters = other.sp_gb_gamma_parameters;
  use_bond_constraints = other.use_bond_constraints;
  use_settle = other.use_settle;
  use_perturbation_info = other.use_perturbation_info;
  use_solvent_cap_option = other.use_solvent_cap_option;
  use_polarization = other.use_polarization;
  water_residue_name = other.water_residue_name;
  bond_constraint_mask = other.bond_constraint_mask;
  bond_constraint_omit_mask = other.bond_constraint_omit_mask;
  rigid_water_count = other.rigid_water_count;
  bond_constraint_count = other.bond_constraint_count;
  degrees_of_freedom = other.degrees_of_freedom;
  nonrigid_particle_count = other.nonrigid_particle_count;
  atom_overflow_names = other.atom_overflow_names;
  atom_overflow_types = other.atom_overflow_types;
  residue_overflow_names = other.residue_overflow_names;
  unused_nhparm = other.unused_nhparm;
  unused_nparm = other.unused_nparm;
  unused_natyp = other.unused_natyp;
  hbond_10_12_parameter_count = other.hbond_10_12_parameter_count;
  heavy_bonds_plus_constraints = other.heavy_bonds_plus_constraints;
  heavy_angls_plus_constraints = other.heavy_angls_plus_constraints;
  heavy_dihes_plus_constraints = other.heavy_dihes_plus_constraints;
  tree_joining_info = other.tree_joining_info;
  last_rotator_info = other.last_rotator_info;
  solty_info = other.solty_info;
  hbond_a_values = other.hbond_a_values;
  hbond_b_values = other.hbond_b_values;
  hbond_cutoffs = other.hbond_cutoffs;
  tree_symbols = other.tree_symbols;
  int_data = other.int_data;
  double_data = other.double_data;
  float_data = other.float_data;
  char4_data = other.char4_data;
  
  // Repair pointers and return the result, analogous to what happens in the PhaseSpace object
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(AtomGraph &&original) :
    version_stamp{},
    date{original.date},
    title{std::move(original.title)},
    source{std::move(original.source)},
    force_fields{std::move(original.force_fields)},
    atom_count{original.atom_count},
    residue_count{original.residue_count},
    molecule_count{original.molecule_count},
    largest_residue_size{original.largest_residue_size},
    last_solute_residue{original.last_solute_residue},
    last_solute_atom{original.last_solute_atom},
    first_solvent_molecule{original.first_solvent_molecule},
    last_atom_before_cap{original.last_atom_before_cap},
    implicit_copy_count{original.implicit_copy_count},
    largest_molecule_size{original.largest_molecule_size},
    descriptors{std::move(original.descriptors)},
    residue_limits{std::move(original.residue_limits)},
    atom_struc_numbers{std::move(original.atom_struc_numbers)},
    residue_numbers{std::move(original.residue_numbers)},
    molecule_limits{std::move(original.molecule_limits)},
    atomic_numbers{std::move(original.atomic_numbers)},
    mobile_atoms{std::move(original.mobile_atoms)},
    molecule_membership{std::move(original.molecule_membership)},
    molecule_contents{std::move(original.molecule_contents)},
    atomic_charges{std::move(original.atomic_charges)},
    atomic_masses{std::move(original.atomic_masses)},
    inverse_atomic_masses{std::move(original.inverse_atomic_masses)},
    sp_atomic_charges{std::move(original.sp_atomic_charges)},
    sp_atomic_masses{std::move(original.sp_atomic_masses)},
    sp_inverse_atomic_masses{std::move(original.sp_inverse_atomic_masses)},
    atom_names{std::move(original.atom_names)},
    atom_types{std::move(original.atom_types)},
    residue_names{std::move(original.residue_names)},
    urey_bradley_term_count{original.urey_bradley_term_count},
    charmm_impr_term_count{original.charmm_impr_term_count},
    cmap_term_count{original.cmap_term_count},
    urey_bradley_parameter_count{original.urey_bradley_parameter_count},
    charmm_impr_parameter_count{original.charmm_impr_parameter_count},
    cmap_surface_count{original.cmap_surface_count},
    urey_bradley_pert_term_count{original.urey_bradley_pert_term_count},
    charmm_impr_pert_term_count{original.charmm_impr_pert_term_count},
    cmap_pert_term_count{original.cmap_pert_term_count},
    urey_bradleys_in_perturbed_group{original.urey_bradleys_in_perturbed_group},
    charmm_imprs_in_perturbed_group{original.charmm_imprs_in_perturbed_group},
    cmaps_in_perturbed_group{original.cmaps_in_perturbed_group},
    urey_bradley_i_atoms{std::move(original.urey_bradley_i_atoms)},
    urey_bradley_k_atoms{std::move(original.urey_bradley_k_atoms)},
    urey_bradley_parameter_indices{std::move(original.urey_bradley_parameter_indices)},
    urey_bradley_assigned_atoms{std::move(original.urey_bradley_assigned_atoms)},
    urey_bradley_assigned_index{std::move(original.urey_bradley_assigned_index)},
    urey_bradley_assigned_terms{std::move(original.urey_bradley_assigned_terms)},
    urey_bradley_assigned_bounds{std::move(original.urey_bradley_assigned_bounds)},
    charmm_impr_i_atoms{std::move(original.charmm_impr_i_atoms)},
    charmm_impr_j_atoms{std::move(original.charmm_impr_j_atoms)},
    charmm_impr_k_atoms{std::move(original.charmm_impr_k_atoms)},
    charmm_impr_l_atoms{std::move(original.charmm_impr_l_atoms)},
    charmm_impr_parameter_indices{std::move(original.charmm_impr_parameter_indices)},
    charmm_impr_assigned_atoms{std::move(original.charmm_impr_assigned_atoms)},
    charmm_impr_assigned_index{std::move(original.charmm_impr_assigned_index)},
    charmm_impr_assigned_terms{std::move(original.charmm_impr_assigned_terms)},
    charmm_impr_assigned_bounds{std::move(original.charmm_impr_assigned_bounds)},
    cmap_i_atoms{std::move(original.cmap_i_atoms)},
    cmap_j_atoms{std::move(original.cmap_j_atoms)},
    cmap_k_atoms{std::move(original.cmap_k_atoms)},
    cmap_l_atoms{std::move(original.cmap_l_atoms)},
    cmap_m_atoms{std::move(original.cmap_m_atoms)},
    cmap_surface_dimensions{std::move(original.cmap_surface_dimensions)},
    cmap_surface_bounds{std::move(original.cmap_surface_bounds)},
    cmap_patch_bounds{std::move(original.cmap_patch_bounds)},
    cmap_surface_indices{std::move(original.cmap_surface_indices)},
    cmap_assigned_atoms{std::move(original.cmap_assigned_atoms)},
    cmap_assigned_index{std::move(original.cmap_assigned_index)},
    cmap_assigned_terms{std::move(original.cmap_assigned_terms)},
    cmap_assigned_bounds{std::move(original.cmap_assigned_bounds)},
    urey_bradley_stiffnesses{std::move(original.urey_bradley_stiffnesses)},
    urey_bradley_equilibria{std::move(original.urey_bradley_equilibria)},
    charmm_impr_stiffnesses{std::move(original.charmm_impr_stiffnesses)},
    charmm_impr_phase_angles{std::move(original.charmm_impr_phase_angles)},
    cmap_surfaces{std::move(original.cmap_surfaces)},
    cmap_phi_derivatives{std::move(original.cmap_phi_derivatives)},
    cmap_psi_derivatives{std::move(original.cmap_psi_derivatives)},
    cmap_phi_psi_derivatives{std::move(original.cmap_phi_psi_derivatives)},
    cmap_patches{std::move(original.cmap_patches)},
    sp_urey_bradley_stiffnesses{std::move(original.sp_urey_bradley_stiffnesses)},
    sp_urey_bradley_equilibria{std::move(original.sp_urey_bradley_equilibria)},
    sp_charmm_impr_stiffnesses{std::move(original.sp_charmm_impr_stiffnesses)},
    sp_charmm_impr_phase_angles{std::move(original.sp_charmm_impr_phase_angles)},
    sp_cmap_surfaces{std::move(original.sp_cmap_surfaces)},
    sp_cmap_phi_derivatives{std::move(original.sp_cmap_phi_derivatives)},
    sp_cmap_psi_derivatives{std::move(original.sp_cmap_psi_derivatives)},
    sp_cmap_phi_psi_derivatives{std::move(original.sp_cmap_phi_psi_derivatives)},
    sp_cmap_patches{std::move(original.sp_cmap_patches)},
    bond_term_with_hydrogen{original.bond_term_with_hydrogen},
    angl_term_with_hydrogen{original.angl_term_with_hydrogen},
    dihe_term_with_hydrogen{original.dihe_term_with_hydrogen},
    bond_term_without_hydrogen{original.bond_term_without_hydrogen},
    angl_term_without_hydrogen{original.angl_term_without_hydrogen},
    dihe_term_without_hydrogen{original.dihe_term_without_hydrogen},
    bond_term_count{original.bond_term_count},
    angl_term_count{original.angl_term_count},
    dihe_term_count{original.dihe_term_count},
    bond_parameter_count{original.bond_parameter_count},
    angl_parameter_count{original.angl_parameter_count},
    dihe_parameter_count{original.dihe_parameter_count},
    bond_perturbation_term_count{original.bond_perturbation_term_count},
    angl_perturbation_term_count{original.angl_perturbation_term_count},
    dihe_perturbation_term_count{original.dihe_perturbation_term_count},
    bonds_in_perturbed_group{original.bonds_in_perturbed_group},
    angls_in_perturbed_group{original.angls_in_perturbed_group},
    dihes_in_perturbed_group{original.dihes_in_perturbed_group},
    bonded_group_count{original.bonded_group_count},
    bond_stiffnesses{std::move(original.bond_stiffnesses)},
    bond_equilibria{std::move(original.bond_equilibria)},
    angl_stiffnesses{std::move(original.angl_stiffnesses)},
    angl_equilibria{std::move(original.angl_equilibria)},
    dihe_amplitudes{std::move(original.dihe_amplitudes)},
    dihe_periodicities{std::move(original.dihe_periodicities)},
    dihe_phase_angles{std::move(original.dihe_phase_angles)},
    sp_bond_stiffnesses{std::move(original.sp_bond_stiffnesses)},
    sp_bond_equilibria{std::move(original.sp_bond_equilibria)},
    sp_angl_stiffnesses{std::move(original.sp_angl_stiffnesses)},
    sp_angl_equilibria{std::move(original.sp_angl_equilibria)},
    sp_dihe_amplitudes{std::move(original.sp_dihe_amplitudes)},
    sp_dihe_periodicities{std::move(original.sp_dihe_periodicities)},
    sp_dihe_phase_angles{std::move(original.sp_dihe_phase_angles)},
    bond_i_atoms{std::move(original.bond_i_atoms)},
    bond_j_atoms{std::move(original.bond_j_atoms)},
    bond_parameter_indices{std::move(original.bond_parameter_indices)},
    bond_assigned_atoms{std::move(original.bond_assigned_atoms)},
    bond_assigned_index{std::move(original.bond_assigned_index)},
    bond_assigned_terms{std::move(original.bond_assigned_terms)},
    bond_assigned_bounds{std::move(original.bond_assigned_bounds)},
    angl_i_atoms{std::move(original.angl_i_atoms)},
    angl_j_atoms{std::move(original.angl_j_atoms)},
    angl_k_atoms{std::move(original.angl_k_atoms)},
    angl_parameter_indices{std::move(original.angl_parameter_indices)},
    angl_assigned_atoms{std::move(original.angl_assigned_atoms)},
    angl_assigned_index{std::move(original.angl_assigned_index)},
    angl_assigned_terms{std::move(original.angl_assigned_terms)},
    angl_assigned_bounds{std::move(original.angl_assigned_bounds)},
    dihe_i_atoms{std::move(original.dihe_i_atoms)},
    dihe_j_atoms{std::move(original.dihe_j_atoms)},
    dihe_k_atoms{std::move(original.dihe_k_atoms)},
    dihe_l_atoms{std::move(original.dihe_l_atoms)},
    dihe_parameter_indices{std::move(original.dihe_parameter_indices)},
    dihe14_parameter_indices{std::move(original.dihe14_parameter_indices)},
    dihe_assigned_atoms{std::move(original.dihe_assigned_atoms)},
    dihe_assigned_index{std::move(original.dihe_assigned_index)},
    dihe_assigned_terms{std::move(original.dihe_assigned_terms)},
    dihe_assigned_bounds{std::move(original.dihe_assigned_bounds)},
    bond_modifiers{std::move(original.bond_modifiers)},
    angl_modifiers{std::move(original.angl_modifiers)},
    dihe_modifiers{std::move(original.dihe_modifiers)},
    bond_assigned_mods{std::move(original.bond_assigned_mods)},
    angl_assigned_mods{std::move(original.angl_assigned_mods)},
    dihe_assigned_mods{std::move(original.dihe_assigned_mods)},
    virtual_site_count{original.virtual_site_count},
    virtual_site_atoms{std::move(original.virtual_site_atoms)},
    virtual_site_frame_types{std::move(original.virtual_site_frame_types)},
    virtual_site_frame1_atoms{std::move(original.virtual_site_frame1_atoms)},
    virtual_site_frame2_atoms{std::move(original.virtual_site_frame2_atoms)},
    virtual_site_frame3_atoms{std::move(original.virtual_site_frame3_atoms)},
    virtual_site_frame4_atoms{std::move(original.virtual_site_frame4_atoms)},
    virtual_site_frame_dim1{std::move(original.virtual_site_frame_dim1)},
    virtual_site_frame_dim2{std::move(original.virtual_site_frame_dim2)},
    virtual_site_frame_dim3{std::move(original.virtual_site_frame_dim3)},
    sp_virtual_site_frame_dim1{std::move(original.sp_virtual_site_frame_dim1)},
    sp_virtual_site_frame_dim2{std::move(original.sp_virtual_site_frame_dim2)},
    sp_virtual_site_frame_dim3{std::move(original.sp_virtual_site_frame_dim3)},
    charge_type_count{original.charge_type_count},
    atom_type_count{original.atom_type_count},
    total_exclusions{original.total_exclusions},
    attenuated_14_type_count{original.attenuated_14_type_count},
    inferred_14_attenuations{original.inferred_14_attenuations},
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    salt_concentration{original.salt_concentration},
    coulomb_constant{original.coulomb_constant},
    pb_radii_set{std::move(original.pb_radii_set)},
    charge_indices{std::move(original.charge_indices)},
    lennard_jones_indices{std::move(original.lennard_jones_indices)},
    atom_exclusion_bounds{std::move(original.atom_exclusion_bounds)},
    atom_exclusion_list{std::move(original.atom_exclusion_list)},
    nb11_exclusion_bounds{std::move(original.nb11_exclusion_bounds)},
    nb11_exclusion_list{std::move(original.nb11_exclusion_list)},
    nb12_exclusion_bounds{std::move(original.nb12_exclusion_bounds)},
    nb12_exclusion_list{std::move(original.nb12_exclusion_list)},
    nb13_exclusion_bounds{std::move(original.nb13_exclusion_bounds)},
    nb13_exclusion_list{std::move(original.nb13_exclusion_list)},
    nb14_exclusion_bounds{std::move(original.nb14_exclusion_bounds)},
    nb14_exclusion_list{std::move(original.nb14_exclusion_list)},
    infr14_i_atoms{std::move(original.infr14_i_atoms)},
    infr14_j_atoms{std::move(original.infr14_j_atoms)},
    infr14_parameter_indices{std::move(original.infr14_parameter_indices)},
    neck_gb_indices{std::move(original.neck_gb_indices)},
    charge_parameters{std::move(original.charge_parameters)},
    lj_a_values{std::move(original.lj_a_values)},
    lj_b_values{std::move(original.lj_b_values)},
    lj_c_values{std::move(original.lj_c_values)},
    lj_14_a_values{std::move(original.lj_14_a_values)},
    lj_14_b_values{std::move(original.lj_14_b_values)},
    lj_14_c_values{std::move(original.lj_14_c_values)},
    lj_type_corrections{std::move(original.lj_type_corrections)},
    attn14_elec_factors{std::move(original.attn14_elec_factors)},
    attn14_vdw_factors{std::move(original.attn14_vdw_factors)},
    atomic_pb_radii{std::move(original.atomic_pb_radii)},
    gb_screening_factors{std::move(original.gb_screening_factors)},
    gb_alpha_parameters{std::move(original.gb_alpha_parameters)},
    gb_beta_parameters{std::move(original.gb_beta_parameters)},
    gb_gamma_parameters{std::move(original.gb_gamma_parameters)},
    sp_charge_parameters{std::move(original.sp_charge_parameters)},
    sp_lj_a_values{std::move(original.sp_lj_a_values)},
    sp_lj_b_values{std::move(original.sp_lj_b_values)},
    sp_lj_c_values{std::move(original.sp_lj_c_values)},
    sp_lj_14_a_values{std::move(original.sp_lj_14_a_values)},
    sp_lj_14_b_values{std::move(original.sp_lj_14_b_values)},
    sp_lj_14_c_values{std::move(original.sp_lj_14_c_values)},
    sp_lj_type_corrections{std::move(original.sp_lj_type_corrections)},
    sp_attn14_elec_factors{std::move(original.sp_attn14_elec_factors)},
    sp_attn14_vdw_factors{std::move(original.sp_attn14_vdw_factors)},
    sp_atomic_pb_radii{std::move(original.sp_atomic_pb_radii)},
    sp_gb_screening_factors{std::move(original.sp_gb_screening_factors)},
    sp_gb_alpha_parameters{std::move(original.sp_gb_alpha_parameters)},
    sp_gb_beta_parameters{std::move(original.sp_gb_beta_parameters)},
    sp_gb_gamma_parameters{std::move(original.sp_gb_gamma_parameters)},
    use_bond_constraints{original.use_bond_constraints},
    use_settle{original.use_settle},
    use_perturbation_info{original.use_perturbation_info},
    use_solvent_cap_option{original.use_solvent_cap_option},
    use_polarization{original.use_polarization},
    water_residue_name{original.water_residue_name},
    bond_constraint_mask{std::move(original.bond_constraint_mask)},
    bond_constraint_omit_mask{std::move(original.bond_constraint_omit_mask)},
    rigid_water_count{original.rigid_water_count},
    bond_constraint_count{original.bond_constraint_count},
    degrees_of_freedom{original.degrees_of_freedom},
    nonrigid_particle_count{original.nonrigid_particle_count},
    atom_overflow_names{std::move(original.atom_overflow_names)},
    atom_overflow_types{std::move(original.atom_overflow_types)},
    residue_overflow_names{std::move(original.residue_overflow_names)},
    unused_nhparm{original.unused_nhparm},
    unused_nparm{original.unused_nparm},
    unused_natyp{original.unused_natyp},
    hbond_10_12_parameter_count{std::move(original.hbond_10_12_parameter_count)},
    heavy_bonds_plus_constraints{std::move(original.heavy_bonds_plus_constraints)},
    heavy_angls_plus_constraints{std::move(original.heavy_angls_plus_constraints)},
    heavy_dihes_plus_constraints{std::move(original.heavy_dihes_plus_constraints)},
    tree_joining_info{std::move(original.tree_joining_info)},
    last_rotator_info{std::move(original.last_rotator_info)},
    solty_info{std::move(original.solty_info)},
    hbond_a_values{std::move(original.hbond_a_values)},
    hbond_b_values{std::move(original.hbond_b_values)},
    hbond_cutoffs{std::move(original.hbond_cutoffs)},
    tree_symbols{std::move(original.tree_symbols)},
    int_data{std::move(original.int_data)},
    double_data{std::move(original.double_data)},
    float_data{std::move(original.float_data)},
    char4_data{std::move(original.char4_data)}
{
  // Copy the last item, the explicit character version stamp
  const int nvst_char = strlen(original.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = original.version_stamp[i];
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraph& AtomGraph::operator=(AtomGraph &&other) {
  
  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy or move the elements of the original, as appropriate
  const int nvst_char = strlen(other.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = other.version_stamp[i];
  }
  date = other.date;
  title = std::move(other.title);
  source = std::move(other.source);
  force_fields = std::move(other.force_fields);
  atom_count = other.atom_count;
  residue_count = other.residue_count;
  molecule_count = other.molecule_count;
  largest_residue_size = other.largest_residue_size;
  last_solute_residue = other.last_solute_residue;
  last_solute_atom = other.last_solute_atom;
  first_solvent_molecule = other.first_solvent_molecule;
  last_atom_before_cap = other.last_atom_before_cap;
  implicit_copy_count = other.implicit_copy_count;
  largest_molecule_size = other.largest_molecule_size;
  descriptors = std::move(other.descriptors);
  residue_limits = std::move(other.residue_limits);
  atom_struc_numbers = std::move(other.atom_struc_numbers);
  residue_numbers = std::move(other.residue_numbers);
  molecule_limits = std::move(other.molecule_limits);
  atomic_numbers = std::move(other.atomic_numbers);
  mobile_atoms = std::move(other.mobile_atoms);
  molecule_membership = std::move(other.molecule_membership);
  molecule_contents = std::move(other.molecule_contents);
  atomic_charges = std::move(other.atomic_charges);
  atomic_masses = std::move(other.atomic_masses);
  inverse_atomic_masses = std::move(other.inverse_atomic_masses);
  sp_atomic_charges = std::move(other.sp_atomic_charges);
  sp_atomic_masses = std::move(other.sp_atomic_masses);
  sp_inverse_atomic_masses = std::move(other.sp_inverse_atomic_masses);
  atom_names = std::move(other.atom_names);
  atom_types = std::move(other.atom_types);
  residue_names = std::move(other.residue_names);
  urey_bradley_term_count = other.urey_bradley_term_count;
  charmm_impr_term_count = other.charmm_impr_term_count;
  cmap_term_count = other.cmap_term_count;
  urey_bradley_parameter_count = other.urey_bradley_parameter_count;
  charmm_impr_parameter_count = other.charmm_impr_parameter_count;
  cmap_surface_count = other.cmap_surface_count;
  urey_bradley_pert_term_count = other.urey_bradley_pert_term_count;
  charmm_impr_pert_term_count = other.charmm_impr_pert_term_count;
  cmap_pert_term_count = other.cmap_pert_term_count;
  urey_bradleys_in_perturbed_group = other.urey_bradleys_in_perturbed_group;
  charmm_imprs_in_perturbed_group = other.charmm_imprs_in_perturbed_group;
  cmaps_in_perturbed_group = other.cmaps_in_perturbed_group;
  urey_bradley_i_atoms = std::move(other.urey_bradley_i_atoms);
  urey_bradley_k_atoms = std::move(other.urey_bradley_k_atoms);
  urey_bradley_parameter_indices = std::move(other.urey_bradley_parameter_indices);
  urey_bradley_assigned_atoms = std::move(other.urey_bradley_assigned_atoms);
  urey_bradley_assigned_index = std::move(other.urey_bradley_assigned_index);
  urey_bradley_assigned_terms = std::move(other.urey_bradley_assigned_terms);
  urey_bradley_assigned_bounds = std::move(other.urey_bradley_assigned_bounds);
  charmm_impr_i_atoms = std::move(other.charmm_impr_i_atoms);
  charmm_impr_j_atoms = std::move(other.charmm_impr_j_atoms);
  charmm_impr_k_atoms = std::move(other.charmm_impr_k_atoms);
  charmm_impr_l_atoms = std::move(other.charmm_impr_l_atoms);
  charmm_impr_parameter_indices = std::move(other.charmm_impr_parameter_indices);
  charmm_impr_assigned_atoms = std::move(other.charmm_impr_assigned_atoms);
  charmm_impr_assigned_index = std::move(other.charmm_impr_assigned_index);
  charmm_impr_assigned_terms = std::move(other.charmm_impr_assigned_terms);
  charmm_impr_assigned_bounds = std::move(other.charmm_impr_assigned_bounds);
  cmap_i_atoms = std::move(other.cmap_i_atoms);
  cmap_j_atoms = std::move(other.cmap_j_atoms);
  cmap_k_atoms = std::move(other.cmap_k_atoms);
  cmap_l_atoms = std::move(other.cmap_l_atoms);
  cmap_m_atoms = std::move(other.cmap_m_atoms);
  cmap_surface_dimensions = std::move(other.cmap_surface_dimensions);
  cmap_surface_bounds = std::move(other.cmap_surface_bounds);
  cmap_patch_bounds = std::move(other.cmap_patch_bounds);
  cmap_surface_indices = std::move(other.cmap_surface_indices);
  cmap_assigned_atoms = std::move(other.cmap_assigned_atoms);
  cmap_assigned_index = std::move(other.cmap_assigned_index);
  cmap_assigned_terms = std::move(other.cmap_assigned_terms);
  cmap_assigned_bounds = std::move(other.cmap_assigned_bounds);
  urey_bradley_stiffnesses = std::move(other.urey_bradley_stiffnesses);
  urey_bradley_equilibria = std::move(other.urey_bradley_equilibria);
  charmm_impr_stiffnesses = std::move(other.charmm_impr_stiffnesses);
  charmm_impr_phase_angles = std::move(other.charmm_impr_phase_angles);
  cmap_surfaces = std::move(other.cmap_surfaces);
  cmap_phi_derivatives = std::move(other.cmap_phi_derivatives);
  cmap_psi_derivatives = std::move(other.cmap_psi_derivatives);
  cmap_phi_psi_derivatives = std::move(other.cmap_phi_psi_derivatives);
  cmap_patches = std::move(other.cmap_patches);
  sp_urey_bradley_stiffnesses = std::move(other.sp_urey_bradley_stiffnesses);
  sp_urey_bradley_equilibria = std::move(other.sp_urey_bradley_equilibria);
  sp_charmm_impr_stiffnesses = std::move(other.sp_charmm_impr_stiffnesses);
  sp_charmm_impr_phase_angles = std::move(other.sp_charmm_impr_phase_angles);
  sp_cmap_surfaces = std::move(other.sp_cmap_surfaces);
  sp_cmap_phi_derivatives = std::move(other.sp_cmap_phi_derivatives);
  sp_cmap_psi_derivatives = std::move(other.sp_cmap_psi_derivatives);
  sp_cmap_phi_psi_derivatives = std::move(other.sp_cmap_phi_psi_derivatives);
  sp_cmap_patches = std::move(other.sp_cmap_patches);
  bond_term_with_hydrogen = other.bond_term_with_hydrogen;
  angl_term_with_hydrogen = other.angl_term_with_hydrogen;
  dihe_term_with_hydrogen = other.dihe_term_with_hydrogen;
  bond_term_without_hydrogen = other.bond_term_without_hydrogen;
  angl_term_without_hydrogen = other.angl_term_without_hydrogen;
  dihe_term_without_hydrogen = other.dihe_term_without_hydrogen;
  bond_term_count = other.bond_term_count;
  angl_term_count = other.angl_term_count;
  dihe_term_count = other.dihe_term_count;
  bond_parameter_count = other.bond_parameter_count;
  angl_parameter_count = other.angl_parameter_count;
  dihe_parameter_count = other.dihe_parameter_count;
  bond_perturbation_term_count = other.bond_perturbation_term_count;
  angl_perturbation_term_count = other.angl_perturbation_term_count;
  dihe_perturbation_term_count = other.dihe_perturbation_term_count;
  bonds_in_perturbed_group = other.bonds_in_perturbed_group;
  angls_in_perturbed_group = other.angls_in_perturbed_group;
  dihes_in_perturbed_group = other.dihes_in_perturbed_group;
  bonded_group_count = other.bonded_group_count;
  bond_stiffnesses = std::move(other.bond_stiffnesses);
  bond_equilibria = std::move(other.bond_equilibria);
  angl_stiffnesses = std::move(other.angl_stiffnesses);
  angl_equilibria = std::move(other.angl_equilibria);
  dihe_amplitudes = std::move(other.dihe_amplitudes);
  dihe_periodicities = std::move(other.dihe_periodicities);
  dihe_phase_angles = std::move(other.dihe_phase_angles);
  sp_bond_stiffnesses = std::move(other.sp_bond_stiffnesses);
  sp_bond_equilibria = std::move(other.sp_bond_equilibria);
  sp_angl_stiffnesses = std::move(other.sp_angl_stiffnesses);
  sp_angl_equilibria = std::move(other.sp_angl_equilibria);
  sp_dihe_amplitudes = std::move(other.sp_dihe_amplitudes);
  sp_dihe_periodicities = std::move(other.sp_dihe_periodicities);
  sp_dihe_phase_angles = std::move(other.sp_dihe_phase_angles);
  bond_i_atoms = std::move(other.bond_i_atoms);
  bond_j_atoms = std::move(other.bond_j_atoms);
  bond_parameter_indices = std::move(other.bond_parameter_indices);
  bond_assigned_atoms = std::move(other.bond_assigned_atoms);
  bond_assigned_index = std::move(other.bond_assigned_index);
  bond_assigned_terms = std::move(other.bond_assigned_terms);
  bond_assigned_bounds = std::move(other.bond_assigned_bounds);
  angl_i_atoms = std::move(other.angl_i_atoms);
  angl_j_atoms = std::move(other.angl_j_atoms);
  angl_k_atoms = std::move(other.angl_k_atoms);
  angl_parameter_indices = std::move(other.angl_parameter_indices);
  angl_assigned_atoms = std::move(other.angl_assigned_atoms);
  angl_assigned_index = std::move(other.angl_assigned_index);
  angl_assigned_terms = std::move(other.angl_assigned_terms);
  angl_assigned_bounds = std::move(other.angl_assigned_bounds);
  dihe_i_atoms = std::move(other.dihe_i_atoms);
  dihe_j_atoms = std::move(other.dihe_j_atoms);
  dihe_k_atoms = std::move(other.dihe_k_atoms);
  dihe_l_atoms = std::move(other.dihe_l_atoms);
  dihe_parameter_indices = std::move(other.dihe_parameter_indices);
  dihe14_parameter_indices = std::move(other.dihe14_parameter_indices);
  dihe_assigned_atoms = std::move(other.dihe_assigned_atoms);
  dihe_assigned_index = std::move(other.dihe_assigned_index);
  dihe_assigned_terms = std::move(other.dihe_assigned_terms);
  dihe_assigned_bounds = std::move(other.dihe_assigned_bounds);
  bond_modifiers = std::move(other.bond_modifiers);
  angl_modifiers = std::move(other.angl_modifiers);
  dihe_modifiers = std::move(other.dihe_modifiers);
  bond_assigned_mods = std::move(other.bond_assigned_mods);
  angl_assigned_mods = std::move(other.angl_assigned_mods);
  dihe_assigned_mods = std::move(other.dihe_assigned_mods);
  virtual_site_count = std::move(other.virtual_site_count);
  virtual_site_atoms = std::move(other.virtual_site_atoms);
  virtual_site_frame_types = std::move(other.virtual_site_frame_types);
  virtual_site_frame1_atoms = std::move(other.virtual_site_frame1_atoms);
  virtual_site_frame2_atoms = std::move(other.virtual_site_frame2_atoms);
  virtual_site_frame3_atoms = std::move(other.virtual_site_frame3_atoms);
  virtual_site_frame4_atoms = std::move(other.virtual_site_frame4_atoms);
  virtual_site_frame_dim1 = std::move(other.virtual_site_frame_dim1);
  virtual_site_frame_dim2 = std::move(other.virtual_site_frame_dim2);
  virtual_site_frame_dim3 = std::move(other.virtual_site_frame_dim3);
  sp_virtual_site_frame_dim1 = std::move(other.sp_virtual_site_frame_dim1);
  sp_virtual_site_frame_dim2 = std::move(other.sp_virtual_site_frame_dim2);
  sp_virtual_site_frame_dim3 = std::move(other.sp_virtual_site_frame_dim3);
  charge_type_count = other.charge_type_count;
  atom_type_count = other.atom_type_count;
  total_exclusions = other.total_exclusions;
  attenuated_14_type_count = other.attenuated_14_type_count;
  inferred_14_attenuations = other.inferred_14_attenuations;
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  salt_concentration = other.salt_concentration;
  coulomb_constant = other.coulomb_constant;
  pb_radii_set = std::move(other.pb_radii_set);
  charge_indices = std::move(other.charge_indices);
  lennard_jones_indices = std::move(other.lennard_jones_indices);
  atom_exclusion_bounds = std::move(other.atom_exclusion_bounds);
  atom_exclusion_list = std::move(other.atom_exclusion_list);
  nb11_exclusion_bounds = std::move(other.nb11_exclusion_bounds);
  nb11_exclusion_list = std::move(other.nb11_exclusion_list);
  nb12_exclusion_bounds = std::move(other.nb12_exclusion_bounds);
  nb12_exclusion_list = std::move(other.nb12_exclusion_list);
  nb13_exclusion_bounds = std::move(other.nb13_exclusion_bounds);
  nb13_exclusion_list = std::move(other.nb13_exclusion_list);
  nb14_exclusion_bounds = std::move(other.nb14_exclusion_bounds);
  nb14_exclusion_list = std::move(other.nb14_exclusion_list);
  infr14_i_atoms = std::move(other.infr14_i_atoms);
  infr14_j_atoms = std::move(other.infr14_j_atoms);
  infr14_parameter_indices = std::move(other.infr14_parameter_indices);
  neck_gb_indices = std::move(other.neck_gb_indices);
  charge_parameters = std::move(other.charge_parameters);
  lj_a_values = std::move(other.lj_a_values);
  lj_b_values = std::move(other.lj_b_values);
  lj_c_values = std::move(other.lj_c_values);
  lj_14_a_values = std::move(other.lj_14_a_values);
  lj_14_b_values = std::move(other.lj_14_b_values);
  lj_14_c_values = std::move(other.lj_14_c_values);
  lj_type_corrections = std::move(other.lj_type_corrections);
  attn14_elec_factors = std::move(other.attn14_elec_factors);
  attn14_vdw_factors = std::move(other.attn14_vdw_factors);
  atomic_pb_radii = std::move(other.atomic_pb_radii);
  gb_screening_factors = std::move(other.gb_screening_factors);
  gb_alpha_parameters = std::move(other.gb_alpha_parameters);
  gb_beta_parameters = std::move(other.gb_beta_parameters);
  gb_gamma_parameters = std::move(other.gb_gamma_parameters);
  sp_charge_parameters = std::move(other.sp_charge_parameters);
  sp_lj_a_values = std::move(other.sp_lj_a_values);
  sp_lj_b_values = std::move(other.sp_lj_b_values);
  sp_lj_c_values = std::move(other.sp_lj_c_values);
  sp_lj_14_a_values = std::move(other.sp_lj_14_a_values);
  sp_lj_14_b_values = std::move(other.sp_lj_14_b_values);
  sp_lj_14_c_values = std::move(other.sp_lj_14_c_values);
  sp_lj_type_corrections = std::move(other.sp_lj_type_corrections);
  sp_attn14_elec_factors = std::move(other.sp_attn14_elec_factors);
  sp_attn14_vdw_factors = std::move(other.sp_attn14_vdw_factors);
  sp_atomic_pb_radii = std::move(other.sp_atomic_pb_radii);
  sp_gb_screening_factors = std::move(other.sp_gb_screening_factors);
  sp_gb_alpha_parameters = std::move(other.sp_gb_alpha_parameters);
  sp_gb_beta_parameters = std::move(other.sp_gb_beta_parameters);
  sp_gb_gamma_parameters = std::move(other.sp_gb_gamma_parameters);
  use_bond_constraints = std::move(other.use_bond_constraints);
  use_settle = other.use_settle;
  use_perturbation_info = other.use_perturbation_info;
  use_solvent_cap_option = other.use_solvent_cap_option;
  use_polarization = other.use_polarization;
  water_residue_name = other.water_residue_name;
  bond_constraint_mask = std::move(other.bond_constraint_mask);
  bond_constraint_omit_mask = std::move(other.bond_constraint_omit_mask);
  rigid_water_count = other.rigid_water_count;
  bond_constraint_count = other.bond_constraint_count;
  degrees_of_freedom = other.degrees_of_freedom;
  nonrigid_particle_count = other.nonrigid_particle_count;
  atom_overflow_names = std::move(other.atom_overflow_names);
  atom_overflow_types = std::move(other.atom_overflow_types);
  residue_overflow_names = std::move(other.residue_overflow_names);
  unused_nhparm = other.unused_nhparm;
  unused_nparm = other.unused_nparm;
  unused_natyp = other.unused_natyp;
  hbond_10_12_parameter_count = other.hbond_10_12_parameter_count;
  heavy_bonds_plus_constraints = std::move(other.heavy_bonds_plus_constraints);
  heavy_angls_plus_constraints = std::move(other.heavy_angls_plus_constraints);
  heavy_dihes_plus_constraints = std::move(other.heavy_dihes_plus_constraints);
  tree_joining_info = std::move(other.tree_joining_info);
  last_rotator_info = std::move(other.last_rotator_info);
  solty_info = std::move(other.solty_info);
  hbond_a_values = std::move(other.hbond_a_values);
  hbond_b_values = std::move(other.hbond_b_values);
  hbond_cutoffs = std::move(other.hbond_cutoffs);
  tree_symbols = std::move(other.tree_symbols);
  int_data = std::move(other.int_data);
  double_data = std::move(other.double_data);
  float_data = std::move(other.float_data);
  char4_data = std::move(other.char4_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getFileName() const {
  return source;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueCount() const {
  return residue_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getMoleculeCount() const {
  return molecule_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLargestResidueSize() const {
  return largest_residue_size;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLastSoluteResidue() const {
  return last_solute_residue;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLastSoluteAtom() const {
  return last_solute_atom;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getFirstSolventMolecule() const {
  return first_solvent_molecule;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getFirstSolventAtom() const {

  // The first solvent molecule may be outside of the system bounds, indicating that there is, in
  // fact, no solvent.  Trap that case and return -1.
  if (first_solvent_molecule < 0 || first_solvent_molecule >= molecule_count) {
    return -1;
  }
  return molecule_contents.readHost(molecule_limits.readHost(first_solvent_molecule));
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLargestMoleculeSize() const {
  return largest_molecule_size;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDescriptor(const TopologyDescriptor choice) const {
  return descriptors.readHost(static_cast<ulint>(choice));
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDescriptor(const SanderDescriptor choice) const {
  return descriptors.readHost(static_cast<ulint>(choice));
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueLimits() const {
  return residue_limits.readHost();
}

//-------------------------------------------------------------------------------------------------
int2 AtomGraph::getResidueLimits(const int index) const {
  int2 tmp = {residue_limits.readHost(index), residue_limits.readHost(index + 1)};
  return tmp;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueIndex() const {

  // Fill in an entire vector based on the residue limits array
  std::vector<int> result(atom_count);
  for (int i = 0; i < residue_count; i++) {
    const int llim = residue_limits.readHost(i);
    const int hlim = residue_limits.readHost(i + 1);
    for (int j = llim; j < hlim; j++) {
      result[j] = i;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueIndex(const int atom_index) const {

  // This will still happen on the fly, rather than storing a long list of numbers.  It's simply
  // not common, or critical in performant code, to access the residue number.
  int lguess = 0;
  int hguess = residue_count;
  const int* limits_ptr = residue_limits.data();
  while (lguess < hguess - 1) {

    // Choose a residue intermediate between the lower and upper bounds
    const int mguess = lguess + ((hguess - lguess) / 2);
    if (atom_index >= limits_ptr[mguess + 1]) {
      lguess = mguess;
    }
    else if (atom_index < limits_ptr[mguess]) {
      hguess = mguess;
    }
    else {
      return mguess;
    }
  }
  return lguess;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueNumber() const {
  return residue_numbers.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getResidueNumber(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot report residue numbers "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "getResidueNumber");
  }
  return residue_numbers.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getResidueNumber(const int index) const {
  return residue_numbers.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeLimits() const {
  return molecule_limits.readHost();
}

//-------------------------------------------------------------------------------------------------
int2 AtomGraph::getMoleculeLimits(const int index) const {
  int2 tmp = {molecule_limits.readHost(index), molecule_limits.readHost(index + 1)};
  return tmp;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getParticlesPerMolecule() const {
  std::vector<int> atoms_per_molecule = molecule_limits.readHost();
  for (int i = 0; i < molecule_count; i++) {
    atoms_per_molecule[i] = atoms_per_molecule[i + 1] - atoms_per_molecule[i];
  }
  atoms_per_molecule.pop_back();
  return atoms_per_molecule;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getParticlesPerMolecule(const int index) const {
  return molecule_limits.readHost(index + 1) - molecule_limits.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomicNumber() const {
  return atomic_numbers.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomicNumber(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot report atomic numbers "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "getAtomicNumber");
  }
  return atomic_numbers.readHost(low_index, high_index - low_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAtomicNumber(const int index) const {
  return atomic_numbers.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> AtomGraph::getAtomMobility() const {
  return getAtomMobility(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> AtomGraph::getAtomMobility(const int low_index, const int high_index) const {

  // Range check as this will use the pointer
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot report atom mobility "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "getAtomMobility");
  }
  std::vector<bool> mobiles(high_index - low_index, true);
  const int* m_ptr = mobile_atoms.data();
  const int int_bits = sizeof(int) * 8;
  for (int i = low_index; i < high_index; i++) {
    const int access_index = i / int_bits;
    mobiles[i - low_index] = ((static_cast<uint>(m_ptr[access_index]) <<
                               (i - (access_index * int_bits))) & 0x1);
  }
  return mobiles;
}

//-------------------------------------------------------------------------------------------------
bool AtomGraph::getAtomMobility(const int index) const {
  const int int_bits = sizeof(int) * 8;
  const int access_index = index / int_bits;
  const uint m_val = static_cast<uint>(mobile_atoms.readHost(access_index));
  return ((m_val << (index - (access_index * int_bits))) & 0x1);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> AtomGraph::getAtomMobilityMask() const {
  return getAtomMobilityMask(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> AtomGraph::getAtomMobilityMask(const int low_index, const int high_index) const {
  const int uint_bits = sizeof(uint) * 8;
  std::vector<uint> result((high_index - low_index + uint_bits - 1) / uint_bits, 0U);
  int result_pos = 0;
  int result_bit = 0;
  const int* mobility_ptr = mobile_atoms.data();
  int mask_pos = low_index / uint_bits;
  int mask_bit = low_index - mask_pos * uint_bits;
  for (int i = low_index; i < high_index; i++) {
    result[result_pos] |= (((static_cast<uint>(mobility_ptr[mask_pos]) >> mask_bit) & 0x1) <<
                           result_bit);
    result_bit++;
    if (result_bit == uint_bits) {
      result_bit = 0;
      result_pos++;
    }
    mask_bit++;
    if (mask_bit == uint_bits) {
      mask_bit = 0;
      mask_pos++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const MobilitySetting movement) {

  // Just feed into the next, most general case.  The cost of the unneccessary bounds check is
  // trivial compared to the simplicity this affords in the code.
  modifyAtomMobility(0, atom_count, movement);
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int low_index, const int high_index,
                                   const MobilitySetting movement) {

  // Range check as this will use the pointer
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot change atom mobility "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "modifyAtomMobility");
  }
  const int int_bits = sizeof(int) * 8;
  int* m_ptr = mobile_atoms.data();
  for (int i = low_index; i < high_index; i++) {
    const int access_index = i / int_bits;
    const int bshift = i - (access_index * int_bits);
    const uint orig_mask = static_cast<uint>(m_ptr[access_index]);
    switch (movement) {
    case MobilitySetting::OFF:
      m_ptr[access_index] = (orig_mask & (~(0x1 << bshift)));
      break;
    case MobilitySetting::ON:
      m_ptr[access_index] = (orig_mask | (0x1 << bshift));
      break;
    case MobilitySetting::TOGGLE:
      m_ptr[access_index] = (orig_mask ^ (0x1 << bshift));
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int index, const MobilitySetting movement) {
  const int int_bits = sizeof(int) * 8;
  const int access_index = index / int_bits;
  const int m_val = mobile_atoms.readHost(access_index);
  const int bshift = index - (access_index * int_bits);
  switch (movement) {
  case MobilitySetting::OFF:
    mobile_atoms.putHost(access_index, m_val & (~(0x1 << bshift)));
    break;
  case MobilitySetting::ON:
    mobile_atoms.putHost(access_index, m_val | (0x1 << bshift));
    break;
  case MobilitySetting::TOGGLE:
    mobile_atoms.putHost(access_index, m_val ^ (0x1 << bshift));
    break;
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeMembership() const {
  return molecule_membership.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeMembership(const int low_index,
                                                  const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid molecule range " + std::to_string(low_index) + " to " +
          std::to_string(high_index) + ".", "AtomGraph", "getMoleculeMemberhip");
  }
  return molecule_membership.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getMoleculeMembership(const int index) const {
  return molecule_membership.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeContents() const {
  return molecule_contents.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getMoleculeContents(const int index) const {
  return molecule_contents.readHost(molecule_limits.readHost(index),
                                    molecule_limits.readHost(index + 1));
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomName() const {
  return atom_names.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomName(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          ".", "AtomGraph", "getAtomName");
  }
  return atom_names.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getAtomName(const int index) const {
  return atom_names.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomType() const {
  return atom_types.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getAtomType(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          ".", "AtomGraph", "getAtomType");
  }
  return atom_types.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getAtomType(const int index) const {
  return atom_types.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::getResidueName() const {
  return residue_names.readHost();
}

//-------------------------------------------------------------------------------------------------
char4 AtomGraph::getResidueName(const int index) const {
  return residue_names.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getUreyBradleyTermCount() const {
  return urey_bradley_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCharmmImprTermCount() const {
  return charmm_impr_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapTermCount() const {
  return cmap_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getUreyBradleyParameterCount() const {
  return urey_bradley_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCharmmImprParameterCount() const {
  return charmm_impr_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapSurfaceCount() const {
  return cmap_surface_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getCmapDimension(const int index) const {
  return cmap_surface_dimensions.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondTermCount() const {
  return bond_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAngleTermCount() const {
  return angl_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDihedralTermCount() const {
  return dihe_term_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondParameterCount() const {
  return bond_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAngleParameterCount() const {
  return angl_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDihedralParameterCount() const {
  return dihe_parameter_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteCount() const {
  return virtual_site_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::findVirtualSites() const {
  return virtual_site_atoms.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> AtomGraph::findVirtualSites(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid virtual site range " + std::to_string(low_index) + " to " +
          std::to_string(high_index) + ".", "AtomGraph", "findVirtualSites");
  }
  std::vector<int2> result;
  const int* vstmp = virtual_site_atoms.data();
  const int nvs = virtual_site_count;
  for (int i = 0; i < nvs; i++) {
    if (vstmp[i] >= low_index && vstmp[i] < high_index) {
      result.push_back({vstmp[i], i});
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::findVirtualSites(const int index) const {
  return virtual_site_atoms.readHost(index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteFrameType(const int index) const {
  return virtual_site_frame_types.readHost(index);  
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getVirtualSiteFrameAtom(const int index, const int nfrm) const {
  switch (nfrm) {
  case 1:
    return virtual_site_frame1_atoms.readHost(index);
  case 2:
    return virtual_site_frame2_atoms.readHost(index);
  case 3:
    return virtual_site_frame3_atoms.readHost(index);
  case 4:
    return virtual_site_frame4_atoms.readHost(index);
  default:
    rtErr("Virtual sites cannot have a frame atom number " + std::to_string(nfrm) + ".",
          "AtomGraph", "getVirtualSiteFrameAtom");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getChargeTypeCount() const {
  return charge_type_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getAtomTypeCount() const {
  return atom_type_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getTotalExclusions() const {
  return total_exclusions;
}

//-------------------------------------------------------------------------------------------------
UnitCellType AtomGraph::getUnitCellType() const {
  return periodic_box_class;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel AtomGraph::getImplicitSolventModel() const {
  return gb_style;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getDielectricConstant() const {
  return dielectric_constant;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getSaltConcentration() const {
  return salt_concentration;
}

//-------------------------------------------------------------------------------------------------
double AtomGraph::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
std::string AtomGraph::getPBRadiiSet() const {
  return pb_radii_set;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getRigidWaterCount() const {
  return rigid_water_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getBondConstraintCount() const {
  return bond_constraint_count;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getDegreesOfFreedom() const {
  return degrees_of_freedom;
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getNonrigidParticleCount() const {
  return nonrigid_particle_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getChargeIndex() const {
  return charge_indices.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getChargeIndex(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          ".", "AtomGraph", "getChargeIndex");
  }
  return charge_indices.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getChargeIndex(const int index) const {
  return charge_indices.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getLennardJonesIndex() const {
  return lennard_jones_indices.readHost();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getLennardJonesIndex(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          ".", "AtomGraph", "getLennardJonesIndex");
  }
  return lennard_jones_indices.readHost(low_index, high_index);
}

//-------------------------------------------------------------------------------------------------
int AtomGraph::getLennardJonesIndex(const int index) const {
  return lennard_jones_indices.readHost(index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getAtomExclusions(const int index) const {
  
  // Assemble 1:1, then 1:2, 1:3, and finally 1:4 exclusions.  They will all get returned.
  std::vector<int> result = getNonbonded11Exclusions(index);
  std::vector<int> result2 = getNonbonded12Exclusions(index);
  std::vector<int> result3 = getNonbonded13Exclusions(index);
  std::vector<int> result4 = getNonbonded14Exclusions(index);
  result.insert(result.end(), result2.begin(), result2.end());
  result.insert(result.end(), result3.begin(), result3.end());
  result.insert(result.end(), result4.begin(), result4.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded11Exclusions(const int index) const {
  std::vector<int> result;
  const int nb11_low = nb11_exclusion_bounds.readHost(index);
  const int nb11_high = nb11_exclusion_bounds.readHost(index + 1);
  const int* nb11_list = nb11_exclusion_list.data();
  for (int i = nb11_low; i < nb11_high; i++) {
    result.push_back(nb11_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded12Exclusions(const int index) const {
  std::vector<int> result;
  const int nb12_low = nb12_exclusion_bounds.readHost(index);
  const int nb12_high = nb12_exclusion_bounds.readHost(index + 1);
  const int* nb12_list = nb12_exclusion_list.data();
  for (int i = nb12_low; i < nb12_high; i++) {
    result.push_back(nb12_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded13Exclusions(const int index) const {
  std::vector<int> result;
  const int nb13_low = nb13_exclusion_bounds.readHost(index);
  const int nb13_high = nb13_exclusion_bounds.readHost(index + 1);
  const int* nb13_list = nb13_exclusion_list.data();
  for (int i = nb13_low; i < nb13_high; i++) {
    result.push_back(nb13_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomGraph::getNonbonded14Exclusions(const int index) const {
  std::vector<int> result;
  const int nb14_low = nb14_exclusion_bounds.readHost(index);
  const int nb14_high = nb14_exclusion_bounds.readHost(index + 1);
  const int* nb14_list = nb14_exclusion_list.data();
  for (int i = nb14_low; i < nb14_high; i++) {
    result.push_back(nb14_list[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowAtomName(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(atom_overflow_names.data(), atom_overflow_names.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(atom_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowAtomType(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(atom_overflow_types.data(),
                      atom_overflow_types.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(atom_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowResidueName(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(residue_overflow_names.data(), residue_overflow_names.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(residue_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ValenceKit<double> AtomGraph::getDoublePrecisionValenceKit(HybridTargetLevel tier) const {
  return ValenceKit<double>(bond_term_count, angl_term_count, dihe_term_count,
                            bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                            inferred_14_attenuations, attenuated_14_type_count,
                            urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                            urey_bradley_parameter_count, charmm_impr_parameter_count,
                            cmap_surface_count, bond_stiffnesses.data(tier),
                            bond_equilibria.data(tier), angl_stiffnesses.data(tier),
                            angl_equilibria.data(tier), dihe_amplitudes.data(tier),
                            dihe_periodicities.data(tier), dihe_phase_angles.data(tier),
                            attn14_elec_factors.data(tier), attn14_vdw_factors.data(tier),
                            bond_i_atoms.data(tier), bond_j_atoms.data(tier),
                            bond_parameter_indices.data(tier), bond_modifiers.data(tier),
                            angl_i_atoms.data(tier), angl_j_atoms.data(tier),
                            angl_k_atoms.data(tier), angl_parameter_indices.data(tier),
                            angl_modifiers.data(tier), dihe_i_atoms.data(tier),
                            dihe_j_atoms.data(tier), dihe_k_atoms.data(tier),
                            dihe_l_atoms.data(tier), dihe_parameter_indices.data(tier),
                            dihe14_parameter_indices.data(tier), dihe_modifiers.data(tier),
                            infr14_i_atoms.data(tier), infr14_j_atoms.data(tier),
                            infr14_parameter_indices.data(tier), urey_bradley_i_atoms.data(tier),
                            urey_bradley_k_atoms.data(tier),
                            urey_bradley_parameter_indices.data(tier),
                            charmm_impr_i_atoms.data(tier), charmm_impr_j_atoms.data(tier),
                            charmm_impr_k_atoms.data(tier), charmm_impr_l_atoms.data(tier),
                            charmm_impr_parameter_indices.data(tier), cmap_i_atoms.data(tier),
                            cmap_j_atoms.data(tier), cmap_k_atoms.data(tier),
                            cmap_l_atoms.data(tier), cmap_m_atoms.data(tier),
                            cmap_surface_dimensions.data(tier), cmap_surface_bounds.data(tier),
                            cmap_patch_bounds.data(tier), cmap_surface_indices.data(tier),
                            urey_bradley_stiffnesses.data(tier),
                            urey_bradley_equilibria.data(tier), charmm_impr_stiffnesses.data(tier),
                            charmm_impr_phase_angles.data(tier), cmap_surfaces.data(tier),
                            cmap_phi_derivatives.data(tier), cmap_psi_derivatives.data(tier),
                            cmap_phi_psi_derivatives.data(tier), cmap_patches.data(tier),
                            bond_assigned_atoms.data(tier), bond_assigned_index.data(tier),
                            bond_assigned_terms.data(tier), bond_assigned_bounds.data(tier),
                            angl_assigned_atoms.data(tier), angl_assigned_index.data(tier),
                            angl_assigned_terms.data(tier), angl_assigned_bounds.data(tier),
                            dihe_assigned_atoms.data(tier), dihe_assigned_index.data(tier),
                            dihe_assigned_terms.data(tier), dihe_assigned_bounds.data(tier),
                            urey_bradley_assigned_atoms.data(tier),
                            urey_bradley_assigned_index.data(tier),
                            urey_bradley_assigned_terms.data(tier),
                            urey_bradley_assigned_bounds.data(tier),
                            charmm_impr_assigned_atoms.data(tier),
                            charmm_impr_assigned_index.data(tier),
                            charmm_impr_assigned_terms.data(tier),
                            charmm_impr_assigned_bounds.data(tier), cmap_assigned_atoms.data(tier),
                            cmap_assigned_index.data(tier), cmap_assigned_terms.data(tier),
                            cmap_assigned_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
ValenceKit<float> AtomGraph::getSinglePrecisionValenceKit(HybridTargetLevel tier) const {
  return ValenceKit<float>(bond_term_count, angl_term_count, dihe_term_count,
                           bond_parameter_count, angl_parameter_count, dihe_parameter_count,
                           inferred_14_attenuations, attenuated_14_type_count,
                           urey_bradley_term_count, charmm_impr_term_count, cmap_term_count,
                           urey_bradley_parameter_count, charmm_impr_parameter_count,
                           cmap_surface_count, sp_bond_stiffnesses.data(tier),
                           sp_bond_equilibria.data(tier), sp_angl_stiffnesses.data(tier),
                           sp_angl_equilibria.data(tier), sp_dihe_amplitudes.data(tier),
                           sp_dihe_periodicities.data(tier), sp_dihe_phase_angles.data(tier),
                           sp_attn14_elec_factors.data(tier), sp_attn14_vdw_factors.data(tier),
                           bond_i_atoms.data(tier), bond_j_atoms.data(tier),
                           bond_parameter_indices.data(tier), bond_modifiers.data(tier),
                           angl_i_atoms.data(tier), angl_j_atoms.data(tier),
                           angl_k_atoms.data(tier), angl_parameter_indices.data(tier),
                           angl_modifiers.data(tier), dihe_i_atoms.data(tier),
                           dihe_j_atoms.data(tier), dihe_k_atoms.data(tier),
                           dihe_l_atoms.data(tier), dihe_parameter_indices.data(tier),
                           dihe14_parameter_indices.data(tier), dihe_modifiers.data(tier),
                           infr14_i_atoms.data(tier), infr14_j_atoms.data(tier),
                           infr14_parameter_indices.data(tier), urey_bradley_i_atoms.data(tier),
                           urey_bradley_k_atoms.data(tier),
                           urey_bradley_parameter_indices.data(tier),
                           charmm_impr_i_atoms.data(tier), charmm_impr_j_atoms.data(tier),
                           charmm_impr_k_atoms.data(tier), charmm_impr_l_atoms.data(tier),
                           charmm_impr_parameter_indices.data(tier), cmap_i_atoms.data(tier),
                           cmap_j_atoms.data(tier), cmap_k_atoms.data(tier),
                           cmap_l_atoms.data(tier), cmap_m_atoms.data(tier),
                           cmap_surface_dimensions.data(tier), cmap_surface_bounds.data(tier),
                           cmap_patch_bounds.data(tier), cmap_surface_indices.data(tier),
                           sp_urey_bradley_stiffnesses.data(tier),
                           sp_urey_bradley_equilibria.data(tier),
                           sp_charmm_impr_stiffnesses.data(tier),
                           sp_charmm_impr_phase_angles.data(tier), sp_cmap_surfaces.data(tier),
                           sp_cmap_phi_derivatives.data(tier), sp_cmap_psi_derivatives.data(tier),
                           sp_cmap_phi_psi_derivatives.data(tier), sp_cmap_patches.data(tier),
                           bond_assigned_atoms.data(tier), bond_assigned_index.data(tier),
                           bond_assigned_terms.data(tier), bond_assigned_bounds.data(tier),
                           angl_assigned_atoms.data(tier), angl_assigned_index.data(tier),
                           angl_assigned_terms.data(tier), angl_assigned_bounds.data(tier),
                           dihe_assigned_atoms.data(tier), dihe_assigned_index.data(tier),
                           dihe_assigned_terms.data(tier), dihe_assigned_bounds.data(tier),
                           urey_bradley_assigned_atoms.data(tier),
                           urey_bradley_assigned_index.data(tier),
                           urey_bradley_assigned_terms.data(tier),
                           urey_bradley_assigned_bounds.data(tier),
                           charmm_impr_assigned_atoms.data(tier),
                           charmm_impr_assigned_index.data(tier),
                           charmm_impr_assigned_terms.data(tier),
                           charmm_impr_assigned_bounds.data(tier), cmap_assigned_atoms.data(tier),
                           cmap_assigned_index.data(tier), cmap_assigned_terms.data(tier),
                           cmap_assigned_bounds.data(tier));
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<double>
AtomGraph::getDoublePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return NonbondedKit<double>(atom_count, charge_type_count, atom_type_count, coulomb_constant,
                              atomic_charges.data(tier), charge_indices.data(tier),
                              lennard_jones_indices.data(tier), charge_parameters.data(tier),
                              lj_a_values.data(tier), lj_b_values.data(tier),
                              lj_c_values.data(tier), lj_14_a_values.data(tier),
                              lj_14_b_values.data(tier), lj_14_c_values.data(tier),
                              nb11_exclusion_list.data(tier), nb11_exclusion_bounds.data(tier),
                              nb12_exclusion_list.data(tier), nb12_exclusion_bounds.data(tier),
                              nb13_exclusion_list.data(tier), nb13_exclusion_bounds.data(tier),
                              nb14_exclusion_list.data(tier), nb14_exclusion_bounds.data(tier),
                              lj_type_corrections.data(tier));
}

//-------------------------------------------------------------------------------------------------
NonbondedKit<float> AtomGraph::getSinglePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  return NonbondedKit<float>(atom_count, charge_type_count, atom_type_count, coulomb_constant,
                             sp_atomic_charges.data(tier), charge_indices.data(tier),
                             lennard_jones_indices.data(tier), sp_charge_parameters.data(tier),
                             sp_lj_a_values.data(tier), sp_lj_b_values.data(tier),
                             sp_lj_c_values.data(tier), sp_lj_14_a_values.data(tier),
                             sp_lj_14_b_values.data(tier), sp_lj_14_c_values.data(tier),
                             nb11_exclusion_list.data(tier), nb11_exclusion_bounds.data(tier),
                             nb12_exclusion_list.data(tier), nb12_exclusion_bounds.data(tier),
                             nb13_exclusion_list.data(tier), nb13_exclusion_bounds.data(tier),
                             nb14_exclusion_list.data(tier), nb14_exclusion_bounds.data(tier),
                             sp_lj_type_corrections.data(tier));
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<double>
AtomGraph::getDoublePrecisionImplicitSolventKit(const HybridTargetLevel tier) const {
  return ImplicitSolventKit<double>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                    neck_gb_indices.data(tier), atomic_pb_radii.data(tier),
                                    gb_screening_factors.data(tier),
                                    gb_alpha_parameters.data(tier), gb_beta_parameters.data(tier),
                                    gb_gamma_parameters.data(tier));
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventKit<float>
AtomGraph::getSinglePrecisionImplicitSolventKit(const HybridTargetLevel tier) const {
  return ImplicitSolventKit<float>(atom_count, gb_style, dielectric_constant, salt_concentration,
                                   neck_gb_indices.data(tier), sp_atomic_pb_radii.data(tier),
                                   sp_gb_screening_factors.data(tier),
                                   sp_gb_alpha_parameters.data(tier),
                                   sp_gb_beta_parameters.data(tier),
                                   sp_gb_gamma_parameters.data(tier));
}

//-------------------------------------------------------------------------------------------------
ChemicalDetailsKit AtomGraph::getChemicalDetailsKit(HybridTargetLevel tier) const {
  return ChemicalDetailsKit(atom_count, residue_count, molecule_count, atom_names.data(tier),
                            residue_names.data(tier), atom_types.data(tier),
                            atomic_numbers.data(tier), residue_limits.data(tier),
                            atom_struc_numbers.data(tier), residue_numbers.data(tier),
                            molecule_membership.data(tier), molecule_contents.data(tier),
                            molecule_limits.data(tier));
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<double>
AtomGraph::getDoublePrecisionVirtualSiteKit(const HybridTargetLevel tier) const {
  return VirtualSiteKit<double>(virtual_site_count, virtual_site_atoms.data(tier),
                                virtual_site_frame_types.data(tier),
                                virtual_site_frame1_atoms.data(tier),
                                virtual_site_frame2_atoms.data(tier),
                                virtual_site_frame3_atoms.data(tier),
                                virtual_site_frame4_atoms.data(tier),
                                virtual_site_frame_dim1.data(tier),
                                virtual_site_frame_dim2.data(tier),
                                virtual_site_frame_dim3.data(tier));
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKit<float>
AtomGraph::getSinglePrecisionVirtualSiteKit(const HybridTargetLevel tier) const {
  return VirtualSiteKit<float>(virtual_site_count, virtual_site_atoms.data(tier),
                               virtual_site_frame_types.data(tier),
                               virtual_site_frame1_atoms.data(tier),
                               virtual_site_frame2_atoms.data(tier),
                               virtual_site_frame3_atoms.data(tier),
                               virtual_site_frame4_atoms.data(tier),
                               sp_virtual_site_frame_dim1.data(tier),
                               sp_virtual_site_frame_dim2.data(tier),
                               sp_virtual_site_frame_dim3.data(tier));
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void AtomGraph::upload() {
  int_data.upload();
  double_data.upload();
  float_data.upload();
  char4_data.upload();
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::download() {
  int_data.download();
  double_data.download();
  float_data.download();
  char4_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void AtomGraph::setSource(const std::string &new_source) {
  source = new_source;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                        const double dielectric_in, const double saltcon_in,
                                        const AtomicRadiusSet radii_set,
                                        const ExceptionResponse policy) {
  gb_style = igb_in;

  // Trap GB use cases with Parse radii (these are much older, for Poisson-Boltzmann calculations)
  switch (radii_set) {
  case AtomicRadiusSet::BONDI:
  case AtomicRadiusSet::AMBER6:
  case AtomicRadiusSet::MBONDI:
  case AtomicRadiusSet::MBONDI2:
  case AtomicRadiusSet::MBONDI3:
    break;
  case AtomicRadiusSet::PARSE:

    // The only use of Parse radii is PB calculations, which are currently not supported in OMNI.
    // Nonetheless, the radius set is included for completeness.  Just make sure that no one is
    // using it for Generalized Born calculations.
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Parse radii are not intended for Generalized Born calculations and cannot "
              "therefore be applied to topology " + source + ".", "AtomGraph",
              "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("Parse radii are not intended for Generalized Born calculations and should not "
               "therefore be applied to topology " + source + ".", "AtomGraph",
               "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  case AtomicRadiusSet::NONE:
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      {
        // Radii are not being applied, so check that nonzero radii are at least present
        bool nonzero_radii_found = false;
        const double* radii_ptr = atomic_pb_radii.data();
        for (int i = 0; i < atom_count; i++) {
          nonzero_radii_found = (nonzero_radii_found || radii_ptr[i] > constants::tiny);
        }
        if (nonzero_radii_found == false) {

          // This is probably always an error, but use the policy switch in case someone is
          // trying something clever.
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("No nonzero radii were found in topology file " + source + ", or applied "
                  "when setting a Generalized Born implicit solvent model.", "AtomGraph",
                  "setImplicitSolventModel");
          case ExceptionResponse::WARN:
            rtWarn("No nonzero radii were found in topology file " + source + ", or applied "
                   "when setting a Generalized Born implicit solvent model.  This is likely to "
                   "cause problems later.", "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
      }
    }
  }

  // Check for offending radii, specifically those that might not fit within the bounds of the
  // "neck" GB models.  If such radii are found, recursively call this function and set the
  // appropriate (Bondi or mBondi3) radii and screening parameters.
  if ((gb_style == ImplicitSolventModel::NECK_GB    && radii_set != AtomicRadiusSet::BONDI) ||
      (gb_style == ImplicitSolventModel::NECK_GB_II && radii_set != AtomicRadiusSet::MBONDI3)) {
    const double* radii_ptr = atomic_pb_radii.data();
    bool bad_radii_found = false;
    for (int i = 0; i < atom_count; i++) {
      bad_radii_found = (bad_radii_found || radii_ptr [i] < 1.0 || radii_ptr[i] > 2.0);
    }
    if (bad_radii_found) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
              "Anstroms or larger than 2.0 Angstroms, which were found in the topology.",
              "AtomGraph", "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
               "Anstroms or larger than 2.0 Angstroms, which were found in the topology.",
               "AtomGraph", "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, AtomicRadiusSet::BONDI, policy);
      return;
    }
  }

  double* alpha_ptr = gb_alpha_parameters.data();
  double* beta_ptr  = gb_beta_parameters.data();
  double* gamma_ptr = gb_gamma_parameters.data();
  double* screen_ptr = gb_screening_factors.data();
  const int* znum_ptr = atomic_numbers.data();

  // Set the dielectric constant and salt concentration, if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    dielectric_constant = dielectric_in;
    salt_concentration = saltcon_in;
    break;
  }

  // Set the radius set and screening parameters (a radius set of NONE will leave these
  // parameters unchanged in the topology)
  const int* nb12_bounds_ptr = nb12_exclusion_bounds.data();
  const int* nb12_excl_ptr = nb12_exclusion_list.data();
  const int* nb13_bounds_ptr = nb13_exclusion_bounds.data();
  const int* nb13_excl_ptr = nb13_exclusion_list.data();
  const char4* atom_type_ptr = atom_types.data();
  const char4* atom_name_ptr = atom_names.data();
  const char4* res_name_ptr = residue_names.data();
  const double* mass_ptr = atomic_masses.data();
  const std::vector<int> atom_residue_idx = getResidueIndex();
  for (int i = 0; i < atom_count; i++) {

    // Impart the radius set 
    double atom_rad = atomic_pb_radii.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.2;
        if (nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] > 0) {
          const int bonded_atom_idx = nb12_excl_ptr[nb12_bounds_ptr[i]];
          if (radii_set == AtomicRadiusSet::AMBER6 || radii_set == AtomicRadiusSet::MBONDI) {
            switch (znum_ptr[bonded_atom_idx]) {
            case 1:
              if (strncmpCased(char4ToString(atom_type_ptr[bonded_atom_idx]).c_str(), "HW", 2,
                               CaseSensitivity::NO)) {
                atom_rad = 0.8;
              }
              break;
            case 6:
              atom_rad = 1.3;
              break;
            case 7:
              if (radii_set == AtomicRadiusSet::MBONDI) {
                atom_rad = 1.3;
              }
              break;
            case 8:
            case 16:
              atom_rad = 0.8;
              break;
            default:
              break;
            }
          }
          else if (radii_set == AtomicRadiusSet::MBONDI2 ||
                   radii_set == AtomicRadiusSet::MBONDI3) {
            if (znum_ptr[bonded_atom_idx] == 7) {
              atom_rad = 1.3;
              const char4 atmc4 = atom_name_ptr[i];
              if (radii_set == AtomicRadiusSet::MBONDI3 &&
                  char4ToString(res_name_ptr[atom_residue_idx[bonded_atom_idx]]) == "ARG " &&
                  (atmc4.x == 'H' && (atmc4.y == 'H' || atmc4.y == 'E'))) {
                atom_rad = 1.17;
              }
            }
          }
        }
        else {
          switch (policy) {
          case ExceptionResponse::DIE:
            break;
          case ExceptionResponse::WARN:
            rtErr("Unbonded hydrogen atom " + char4ToString(atom_name_ptr[i]) + " detected in "
                  "residue " + char4ToString(res_name_ptr[atom_residue_idx[i]]) + ".  This "
                  "hydrogen's radius will keep the default Bondi value of 1.2 Angstroms.\n",
                  "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
        break;
      case 6:
        {
          // Detect repartitioning in hydrogen masses when calculating the atom's true mass
          double atomi_mass = mass_ptr[i];
          for (int j = nb12_bounds_ptr[i]; j < nb12_bounds_ptr[i + 1]; j++) {
            const int neighbor_atom_idx = nb12_excl_ptr[j];
            if (znum_ptr[neighbor_atom_idx] == 1 && mass_ptr[neighbor_atom_idx] > 1.008) {
              atomi_mass += mass_ptr[neighbor_atom_idx] - 1.008;
            }
          }

          // Assuming all masses are based on the natural abundances, identify unified carbon
          // atoms as having (mass of C + mass of one or more H).
          const char4 attc4 = atom_type_ptr[i];
          if ((attc4.x == 'C' && attc4.y == '1' && atomi_mass > 13.018) ||
              (attc4.x == 'C' && attc4.y == '2' && atomi_mass > 14.026) ||
              (attc4.x == 'C' && attc4.y == '3' && atomi_mass > 15.034)) {

            // United atom carbon radius
            atom_rad = 2.2;
          }
          else {

            // Standard carbon radius
            atom_rad = 1.7;
          }
        }
        break;
      case 7:
        atom_rad = 1.55;
        break;
      case 8:
        atom_rad = 1.5;
        if (radii_set == AtomicRadiusSet::MBONDI3) {

          // Adjust carboxylic oxygens on proteins (side-chains as well as C-termini)
          const std::string resc4 = char4ToString(res_name_ptr[atom_residue_idx[i]]);
          const std::string atmc4 = char4ToString(atom_name_ptr[i]);
          if (((resc4 == "ASP " || resc4 == "AS4 ") && (atmc4 == "OD1 " || atmc4 == "OD2 ")) ||
              ((resc4 == "GLU " || resc4 == "GL4 ") && (atmc4 == "OE1 " || atmc4 == "OE2 ")) ||
              atmc4 == "OXT ") {
            atom_rad = 1.4;
          }
          if (znum_ptr[i] == 8 && nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] == 1 &&
              znum_ptr[nb12_excl_ptr[nb12_bounds_ptr[i]]] == 6) {
            for (int j = nb13_bounds_ptr[i]; j < nb13_bounds_ptr[i + 1]; j++) {
              if (znum_ptr[nb13_excl_ptr[j]] == 8) {
                atom_rad = 1.4;
              }
            }
          }
        }
        break;
      case 9:
        atom_rad = 1.5;
        break;
      case 14:
        atom_rad = 2.1;
        break;
      case 15:
        atom_rad = 1.85;
        break;
      case 16:
        atom_rad = 1.8;
        break;
      case 17:
        atom_rad = 1.7;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::PARSE:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.00;
        break;
      case 6:
        atom_rad = 1.7;
        break;
      case 7:
        atom_rad = 1.5;
        break;
      case 8:
        atom_rad = 1.4;
        break;
      case 16:
        atom_rad = 1.85;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::NONE:
      break;
    }
    atomic_pb_radii.putHost(atom_rad, i);

    // Impart the GB screening parameters
    double atom_screen = gb_screening_factors.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_screen = 0.85;
        break;
      case 6:
        atom_screen = 0.72;
        break;
      case 7:
        atom_screen = 0.79;
        break;
      case 8:
        atom_screen = 0.85;
        break;
      case 9:
        atom_screen = 0.88;
        break;
      case 15:
        atom_screen = 0.86;
        break;
      case 16:
        atom_screen = 0.96;
        break;
      default:
        atom_screen = 0.8;
        break;
      }
    case AtomicRadiusSet::PARSE:
    case AtomicRadiusSet::NONE:
      break;
    }
    gb_screening_factors.putHost(atom_screen, i);
  }

  // Set alpha, beta, gamma, and screening factors if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = 0.0;
      beta_ptr[i]  = 0.0;
      gamma_ptr[i] = 0.0;
    }
    break;
  case ImplicitSolventModel::OBC_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_i_alpha;
      beta_ptr[i]  = gb_obc_i_beta;
      gamma_ptr[i] = gb_obc_i_gamma;
    }
    break;
  case ImplicitSolventModel::OBC_GB_II:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_ii_alpha;
      beta_ptr[i]  = gb_obc_ii_beta;
      gamma_ptr[i] = gb_obc_ii_gamma;
    }
    break;
  case ImplicitSolventModel::NECK_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_neck_i_alpha;
      beta_ptr[i]  = gb_neck_i_beta;
      gamma_ptr[i] = gb_neck_i_gamma;
      switch(znum_ptr[i]) {
      case 1:
        screen_ptr[i] = gb_neck_i_screen_h;
        break;
      case 6:
        screen_ptr[i] = gb_neck_i_screen_c;
        break;
      case 7:
        screen_ptr[i] = gb_neck_i_screen_n;
        break;
      case 8:
        screen_ptr[i] = gb_neck_i_screen_o;
        break;
      case 16:
        screen_ptr[i] = gb_neck_i_screen_s;
        break;
      default:
        screen_ptr[i] = gb_neck_i_screen_default;
        break;
      }
    }
    break;
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      switch(znum_ptr[i]) {
      case 1:
        alpha_ptr[i]  = gb_neck_ii_alpha_h;
        beta_ptr[i]   = gb_neck_ii_beta_h;
        gamma_ptr[i]  = gb_neck_ii_gamma_h;
        screen_ptr[i] = gb_neck_ii_screen_h;
        break;
      case 6:
        alpha_ptr[i]  = gb_neck_ii_alpha_c;
        beta_ptr[i]   = gb_neck_ii_beta_c;
        gamma_ptr[i]  = gb_neck_ii_gamma_c;
        screen_ptr[i] = gb_neck_ii_screen_c;
        break;
      case 7:
        alpha_ptr[i]  = gb_neck_ii_alpha_n;
        beta_ptr[i]   = gb_neck_ii_beta_n;
        gamma_ptr[i]  = gb_neck_ii_gamma_n;
        screen_ptr[i] = gb_neck_ii_screen_n;
        break;
      case 8:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_o;
        break;
      case 15:
        alpha_ptr[i]  = gb_neck_ii_alpha_p;
        beta_ptr[i]   = gb_neck_ii_beta_p;
        gamma_ptr[i]  = gb_neck_ii_gamma_p;
        screen_ptr[i] = gb_neck_ii_screen_p;
        break;
      case 16:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_s;
        break;
      default:
        alpha_ptr[i]  = 1.0;
        beta_ptr[i]   = 0.8;
        gamma_ptr[i]  = 4.85;
        screen_ptr[i] = 0.5;
        break;
      }
    }
    break;
  }

  // Note the radius set in the topology, if it has indeed changed
  if (pb_radii_set.size() == 0 || radii_set != AtomicRadiusSet::NONE) {
    pb_radii_set = getAtomicRadiusSetName(radii_set);
  }

  // Compute the neck GB indices based on the baseline atomic PB radii.  These values must later
  // be checked against the available table size.
  int* neck_idx_ptr = neck_gb_indices.data();
  double* radii_ptr = atomic_pb_radii.data();
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      neck_idx_ptr[i] = static_cast<int>(((radii_ptr[i] - 1.0) * 20.0) + 0.5);
    }
    break;
  }
}

} // namespace topology
} // namespace omni
