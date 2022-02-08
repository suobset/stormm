#include "Math/rounding.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_enumerators.h"
#include "atomgraph_synthesis.h"

namespace omni {
namespace topology {

using cuda::HybridKind;
using math::roundUp;
using math::maxAbsoluteDifference;
using math::maxValue;
using math::minValue;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::sum;
using parse::realToString;
using testing::Approx;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using topology::ValenceKit;
  
//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<int> topology_indices_in) :
    topology_count{static_cast<int>(topologies_in.size())},
    system_count{static_cast<int>(topology_indices_in.size())},
    total_atoms{0}, total_virtual_sites{0}, total_bond_terms{0}, total_angl_terms{0},
    total_dihe_terms{0}, total_ubrd_terms{0}, total_cimp_terms{0}, total_cmap_terms{0},
    total_atom_types{0}, total_charge_types{0}, total_bond_params{0}, total_angl_params{0},
    total_dihe_params{0}, total_ubrd_params{0}, total_cimp_params{0}, total_cmap_surfaces{0},
    periodic_box_class{UnitCellType::NONE}, gb_style{ImplicitSolventModel::NONE},
    dielectric_constant{1.0}, salt_concentration{0.0}, coulomb_constant{accepted_coulomb_constant},
    use_bond_constraints{ShakeSetting::OFF}, use_settle{SettleSetting::OFF},
    water_residue_name{' ', ' ', ' ', ' '}, pb_radii_sets{},
    topologies{topologies_in},
    topology_indices{HybridKind::POINTER, "tpsyn_indices"},
    atom_counts{HybridKind::POINTER, "typsyn_atom_counts"},
    residue_counts{HybridKind::POINTER, "typsyn_res_counts"},
    molecule_counts{HybridKind::POINTER, "typsyn_mol_counts"},
    largest_residue_sizes{HybridKind::POINTER, "typsyn_max_res"},
    last_solute_residues{HybridKind::POINTER, "typsyn_last_sol_res"},
    last_solute_atoms{HybridKind::POINTER, "typsyn_last_sol_atm"},
    first_solvent_molecules{HybridKind::POINTER, "typsyn_1st_solv_mol"},
    ubrd_term_counts{HybridKind::POINTER, "typsyn_ubrd_counts"},
    cimp_term_counts{HybridKind::POINTER, "typsyn_cimp_counts"},
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
    atom_bit_offsets{HybridKind::POINTER, "typsyn_abit_offsets"},
    residue_offsets{HybridKind::POINTER, "typsyn_res_offsets"},
    molecule_offsets{HybridKind::POINTER, "typsyn_res_offsets"},
    ubrd_term_offsets{HybridKind::POINTER, "typsyn_ubrd_offset"},
    cimp_term_offsets{HybridKind::POINTER, "typsyn_cimp_offset"},
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
    residue_names{HybridKind::ARRAY, "tpsyn_res_names"},
    ubrd_stiffnesses{HybridKind::ARRAY, "tpsyn_ub_stiff"},
    ubrd_equilibria{HybridKind::ARRAY, "tpsyn_ub_equil"},
    cimp_stiffnesses{HybridKind::ARRAY, "tpsyn_cimp_stiff"},
    cimp_phase_angles{HybridKind::ARRAY, "tpsyn_cimp_equil"},
    cmap_surface_dimensions{HybridKind::ARRAY, "tpsyn_cmap_dims"},
    cmap_surface_bounds{HybridKind::ARRAY, "tpsyn_cmap_bounds"},
    cmap_surfaces{HybridKind::ARRAY, "tpsyn_cmap_surf"},
    sp_ubrd_stiffnesses{HybridKind::ARRAY, "tpsyn_ub_stiff_sp"},
    sp_ubrd_equilibria{HybridKind::ARRAY, "tpsyn_ub_equil_sp"},
    sp_cimp_stiffnesses{HybridKind::ARRAY, "tpsyn_cimp_stiff_sp"},
    sp_cimp_phase_angles{HybridKind::ARRAY, "tpsyn_cimp_equil_sp"},
    sp_cmap_surfaces{HybridKind::ARRAY, "tpsyn_cmap_surf_sp"},
    bond_stiffnesses{HybridKind::ARRAY, "tpsyn_bondk"},
    bond_equilibria{HybridKind::ARRAY, "tpsyn_bondl0"},
    angl_stiffnesses{HybridKind::ARRAY, "tpsyn_anglk"},
    angl_equilibria{HybridKind::ARRAY, "tpsyn_anglt0"},
    dihe_amplitudes{HybridKind::ARRAY, "tpsyn_dihek"},
    dihe_periodicities{HybridKind::ARRAY, "tpsyn_dihen"},
    dihe_phase_angles{HybridKind::ARRAY, "tpsyn_dihepsi"},
    sp_bond_stiffnesses{HybridKind::ARRAY, "tpsyn_bondk_sp"},
    sp_bond_equilibria{HybridKind::ARRAY, "tpsyn_bondl0_sp"},
    sp_angl_stiffnesses{HybridKind::ARRAY, "tpsyn_anglk_sp"},
    sp_angl_equilibria{HybridKind::ARRAY, "tpsyn_anglt0_sp"},
    sp_dihe_amplitudes{HybridKind::ARRAY, "tpsyn_dihek_sp"},
    sp_dihe_periodicities{HybridKind::ARRAY, "tpsyn_dihen_sp"},
    sp_dihe_phase_angles{HybridKind::ARRAY, "tpsyn_dihepsi_sp"},
    ubrd_i_atoms{HybridKind::ARRAY, "tpsyn_ubrd_i"},
    ubrd_k_atoms{HybridKind::ARRAY, "tpsyn_ubrd_k"},
    ubrd_param_idx{HybridKind::ARRAY, "tpsyn_ubrd_idx"},
    cimp_i_atoms{HybridKind::ARRAY, "tpsyn_cimp_i"},
    cimp_j_atoms{HybridKind::ARRAY, "tpsyn_cimp_j"},
    cimp_k_atoms{HybridKind::ARRAY, "tpsyn_cimp_k"},
    cimp_l_atoms{HybridKind::ARRAY, "tpsyn_cimp_l"},
    cimp_param_idx{HybridKind::ARRAY, "tpsyn_cimp_idx"},
    cmap_i_atoms{HybridKind::ARRAY, "tpsyn_cmap_i"},
    cmap_j_atoms{HybridKind::ARRAY, "tpsyn_cmap_j"},
    cmap_k_atoms{HybridKind::ARRAY, "tpsyn_cmap_k"},
    cmap_l_atoms{HybridKind::ARRAY, "tpsyn_cmap_l"},
    cmap_m_atoms{HybridKind::ARRAY, "tpsyn_cmap_m"},
    cmap_param_idx{HybridKind::ARRAY, "tpsyn_cmap_idx"},
    bond_i_atoms{HybridKind::ARRAY, "tpsyn_bond_i"},
    bond_j_atoms{HybridKind::ARRAY, "tpsyn_bond_j"},
    bond_param_idx{HybridKind::ARRAY, "tpsyn_bond_idx"},
    angl_i_atoms{HybridKind::ARRAY, "tpsyn_angl_i"},
    angl_j_atoms{HybridKind::ARRAY, "tpsyn_angl_j"},
    angl_k_atoms{HybridKind::ARRAY, "tpsyn_angl_k"},
    angl_param_idx{HybridKind::ARRAY, "tpsyn_angl_idx"},
    dihe_i_atoms{HybridKind::ARRAY, "tpsyn_dihe_i"},
    dihe_j_atoms{HybridKind::ARRAY, "tpsyn_dihe_j"},
    dihe_k_atoms{HybridKind::ARRAY, "tpsyn_dihe_k"},
    dihe_l_atoms{HybridKind::ARRAY, "tpsyn_dihe_l"},
    dihe_param_idx{HybridKind::ARRAY, "tpsyn_dihe_idx"},
    charge_indices{HybridKind::ARRAY, "tpsyn_q_idx"},
    lennard_jones_indices{HybridKind::ARRAY, "tpsyn_lj_idx"},
    nmr_initial_steps{HybridKind::ARRAY, "tpsyn_nmr_init_step"},
    nmr_final_steps{HybridKind::ARRAY, "tpsyn_nmr_final_step"},
    nmr_increments{HybridKind::ARRAY, "tpsyn_nmr_inc"},
    nmr_k_initial_values{HybridKind::ARRAY, "tpsyn_k_init"},
    nmr_r_initial_values{HybridKind::ARRAY, "tpsyn_r_init"},
    nmr_k_final_values{HybridKind::ARRAY, "tpsyn_k_final"},
    nmr_r_final_values{HybridKind::ARRAY, "tpsyn_r_final"},
    sp_nmr_k_initial_values{HybridKind::ARRAY, "tpsyn_k_init_sp"},
    sp_nmr_r_initial_values{HybridKind::ARRAY, "tpsyn_r_init_sp"},
    sp_nmr_k_final_values{HybridKind::ARRAY, "tpsyn_k_final_sp"},
    sp_nmr_r_final_values{HybridKind::ARRAY, "tpsyn_r_final_sp"},
    atom_imports{HybridKind::ARRAY, "tpsyn_atom_imports"},
    vwu_instruction_sets{HybridKind::ARRAY, "tpsyn_vwu_insr_sets"},
    bond_instructions{HybridKind::ARRAY, "tpsyn_bond_insr"},
    angl_instructions{HybridKind::ARRAY, "tpsyn_angl_insr"},
    dihe_instructions{HybridKind::ARRAY, "tpsyn_dihe_insr"},
    cmap_instructions{HybridKind::ARRAY, "tpsyn_cmap_insr"},
    nmr2_instructions{HybridKind::ARRAY, "tpsyn_nmr2_insr"},
    nmr3_instructions{HybridKind::ARRAY, "tpsyn_nmr3_insr"},
    nmr4_instructions{HybridKind::ARRAY, "tpsyn_nmr4_insr"}
{
  // Setup and memory layout
  const std::vector<int> topology_index_rebase = checkTopologyList(topology_indices_in);
  checkCommonSettings();
  buildAtomAndTermArrays(topology_indices_in, topology_index_rebase);

  // Condense valence and charge parameters into compact tables of unique values
  condenseParameterTables();
  
  // The unique Lennard-Jones parameters are hard to map out.  To the degree that there are
  // unique parameters in each system, the Lennard-Jones tables might as well be block matrices.
  // Keep a list of the A, B, and possibly C coefficients of each Lennard-Jones type interacting
  // with all other types.
  extendLJMatrices();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
AtomGraphSynthesis::checkTopologyList(const std::vector<int> &topology_indices_in) {

  // Check that no system indexes a topology outside of the given list
  if (maxValue(topology_indices_in) >= topology_count || minValue(topology_indices_in) < 0) {
    rtErr("One or more systems references a topology index outside of the list supplied.",
          "AtomGraphSynthesis");
  }

  // Check that each topology in the supplied list is being used.  Roll that result into an
  // analysis of the uniqueness of each topology.
  std::vector<bool> topology_unique(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    topology_unique[topology_indices_in[i]] = true;
  }
  int n_unused_topology = 0;
  for (int i = 0; i < topology_count; i++) {
    n_unused_topology += (topology_unique[i] == false);
  }
  if (n_unused_topology > 0) {
    rtWarn("Out of " + std::to_string(topology_count) + " topologies, " +
           std::to_string(n_unused_topology) + " are not referenced by any systems in this "
           "synthesis.", "AtomGraphSynthesis");
  }

  // Check that all topologies are, in fact, unique.  Compact the list if necessary and update
  // the system topology indexing.
  for (int i = 0; i < topology_count; i++) {
    if (topology_unique[i]) {
      const int topi_natom = topologies[i]->getAtomCount();
      for (int j = i + 1; j < topology_count; j++) {
        if (topologies[j]->getAtomCount() == topi_natom &&
            topologies[i]->getFileName() == topologies[j]->getFileName()) {
          topology_unique[j] = false;
        }
      }
    }
  }
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
          "AtomGraphSynthesis");
  }
  if (system_count == 0) {
    rtErr("No systems making use of any of the " + std::to_string(topology_count) +
          " topologies were detected.", "AtomGraphSynthesis");
  }
  return topology_index_rebase;
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
  for (int i = 1; i < topology_count; i++) {
    if (topologies[i]->getImplicitSolventModel() != gb_style) {
      rtErr("All topologies must have a consistent implicit solvent model setting.  The first, "
            "set to " + getImplicitSolventModelName(gb_style) + ", does not agree with topology " +
            std::to_string(i + 1) + ", using " +
            getImplicitSolventModelName(topologies[i]->getImplicitSolventModel()) + ".",
            "AtomGraphSynthesis", "checkCommonSettings");
    }
    if (appr_dielcon.test(topologies[i]->getDielectricConstant()) == false) {
      rtErr("All topologies must have a consistent dielectric constant, but values of " +
            realToString(dielectric_constant) + " and " +
            realToString(topologies[i]->getDielectricConstant()) + " were found."
            "AtomGraphSynthesis", "checkCommonSettings");
    }
    if (appr_saltcon.test(topologies[i]->getSaltConcentration()) == false) {
      rtErr("All systems must use the same salt concentration.  Values of " +
            realToString(salt_concentration) + " and " +
            realToString(topologies[i]->getSaltConcentration()) + " were found.",
            "AtomGraphSynthesis", "checkCommonSettings");
    }
    if (appr_coulomb.test(topologies[i]->getCoulombConstant()) == false) {
      rtErr("All topologies must make use of the same definitions of physical constants.  "
            "Coulomb's constant is defined as " + realToString(coulomb_constant, 6) +
            " in the first topology and " + realToString(topologies[i]->getCoulombConstant()) +
            " in topology " + std::to_string(i + 1) + ".", "AtomGraphSynthesis",
            "checkCommonSettings");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::buildAtomAndTermArrays(const std::vector<int> &topology_indices_in,
                                                const std::vector<int> &topology_index_rebase) {

  // Allocate memory and set POINTER-kind arrays for the small packets of data
  const int padded_system_count = roundUp(system_count, warp_size_int);
  int_system_data.resize(33 * padded_system_count);
  int pivot = 0;
  topology_indices.setPointer(&int_system_data, pivot, system_count);
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
  nb_exclusion_offsets.setPointer(&int_system_data, pivot, system_count);
  
  // Load the topology indexing first
  for (int i = 0; i < system_count; i++) {
    topology_indices.putHost(topology_index_rebase[topology_indices_in[i]], i);
  }

  // Loop over all systems, fill in the above details, and compute the sizes of various arrays
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
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(i)];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    atom_counts.putHost(cdk.natom, i);
    total_atoms += cdk.natom;
    virtual_site_counts.putHost(ag_ptr->getVirtualSiteCount(), i);
    total_virtual_sites += vsk.nsite;
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
    vste_offset += roundUp(ag_ptr->getVirtualSiteCount(), warp_size_int);
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
  }

  // Allocate detailed arrays for each descriptor, then collate all topologies.  This
  // fills the "atom and residue details" arrays listed in atomgraph_synthesis.h.
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
    const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
    const int synth_bit_base = atom_bit_offsets.readHost(i);
    const int synth_atom_base = atom_offsets.readHost(i);
    const int synth_residue_base = residue_offsets.readHost(i);
    const int synth_molecule_base = molecule_offsets.readHost(i);
    for (int j = 0; j < cdk.nres + 1; j++) {
      residue_limits_ptr[synth_residue_base + j] = cdk.res_limits[j] + synth_atom_base;
    }
    for (int j = 0; j < cdk.natom; j++) {
      atom_struc_numbers_ptr[synth_atom_base + j] = cdk.atom_numbers[j];
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
      residue_numbers_ptr[synth_residue_base + j] = cdk.res_numbers[j];
      residue_names_ptr[synth_residue_base + j].x = cdk.res_names[j].x;
      residue_names_ptr[synth_residue_base + j].y = cdk.res_names[j].y;
      residue_names_ptr[synth_residue_base + j].z = cdk.res_names[j].z;
      residue_names_ptr[synth_residue_base + j].w = cdk.res_names[j].w;
    }
    for (int j = 0; j < cdk.nmol + 1; j++) {
      molecule_limits_ptr[synth_molecule_base + j] = cdk.mol_limits[j] + synth_atom_base;
    }
    const int uint_bits = sizeof(uint) * 8;
    mobile_atoms.putHost(ag_ptr->getAtomMobilityMask(), synth_bit_base,
                         (cdk.natom + uint_bits - 1) / uint_bits);
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
  ubrd_i_atoms.resize(ubrd_offset);
  ubrd_k_atoms.resize(ubrd_offset);
  cimp_i_atoms.resize(cimp_offset);
  cimp_j_atoms.resize(cimp_offset);
  cimp_k_atoms.resize(cimp_offset);
  cimp_l_atoms.resize(cimp_offset);
  cmap_i_atoms.resize(cmap_offset);
  cmap_j_atoms.resize(cmap_offset);
  cmap_k_atoms.resize(cmap_offset);
  cmap_l_atoms.resize(cmap_offset);
  cmap_m_atoms.resize(cmap_offset);
  bond_i_atoms.resize(bond_offset);
  bond_j_atoms.resize(bond_offset);
  angl_i_atoms.resize(angl_offset);
  angl_j_atoms.resize(angl_offset);
  angl_k_atoms.resize(angl_offset);
  dihe_i_atoms.resize(dihe_offset);
  dihe_j_atoms.resize(dihe_offset);
  dihe_k_atoms.resize(dihe_offset);
  dihe_l_atoms.resize(dihe_offset);
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
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::condenseParameterTables() {

  // Compute the numbers of unique parameters.  Take the opportunity to compute offsets (starting
  // bounds) for various sets of terms.
  int bond_offset = 0;
  int angl_offset = 0;
  int dihe_offset = 0;
  int ubrd_offset = 0;
  int cimp_offset = 0;
  int cmap_offset = 0;
  int chrg_offset = 0;
  std::vector<int> topology_bond_table_offsets(topology_count);
  std::vector<int> topology_angl_table_offsets(topology_count);
  std::vector<int> topology_dihe_table_offsets(topology_count);
  std::vector<int> topology_ubrd_table_offsets(topology_count);
  std::vector<int> topology_cimp_table_offsets(topology_count);
  std::vector<int> topology_cmap_table_offsets(topology_count);
  std::vector<int> topology_chrg_table_offsets(topology_count);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* ag_ptr = topologies[i];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
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
    bond_offset += vk.nbond;
    angl_offset += vk.nangl;
    dihe_offset += vk.ndihe;
    ubrd_offset += vk.nubrd;
    cimp_offset += vk.ncimp;
    cmap_offset += vk.ncmap;
    chrg_offset += nbk.natom;
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
  std::vector<double> filtered_chrg;
  std::vector<float> sp_filtered_chrg;
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
  int n_unique_bond = 0;
  int n_unique_angl = 0;
  int n_unique_dihe = 0;
  int n_unique_ubrd = 0;
  int n_unique_cimp = 0;
  int n_unique_cmap = 0;
  int n_unique_chrg = 0;
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    const ChemicalDetailsKit i_cdk     = iag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> i_nbk   = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk      = iag_ptr->getDoublePrecisionValenceKit();
    const NonbondedKit<float> i_nbk_sp = iag_ptr->getSinglePrecisionNonbondedKit();
    const ValenceKit<float> i_vk_sp    = iag_ptr->getSinglePrecisionValenceKit();

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

    // Seek out unique angle parameters
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
        const NonbondedKit<float> k_nbk_sp = kag_ptr->getSinglePrecisionNonbondedKit();
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
  }

  // Record synthesis parameters found thus far.  Pad the relevant ARRAY-kind Hybrid objects in
  // case some developer tries to have all lanes of a warp access an element.
  total_bond_params   = n_unique_bond;
  total_angl_params   = n_unique_angl;
  total_dihe_params   = n_unique_dihe;
  total_ubrd_params   = n_unique_ubrd;
  total_cimp_params   = n_unique_cimp;
  total_cmap_surfaces = n_unique_cmap;
  total_charge_types  = n_unique_chrg;
  bond_stiffnesses.resize(roundUp(total_bond_params, warp_size_int));
  bond_equilibria.resize(roundUp(total_bond_params, warp_size_int));
  angl_stiffnesses.resize(roundUp(total_angl_params, warp_size_int));
  angl_equilibria.resize(roundUp(total_angl_params, warp_size_int));
  dihe_amplitudes.resize(roundUp(total_dihe_params, warp_size_int));
  dihe_periodicities.resize(roundUp(total_dihe_params, warp_size_int));
  dihe_phase_angles.resize(roundUp(total_dihe_params, warp_size_int));
  sp_bond_stiffnesses.resize(roundUp(total_bond_params, warp_size_int));
  sp_bond_equilibria.resize(roundUp(total_bond_params, warp_size_int));
  sp_angl_stiffnesses.resize(roundUp(total_angl_params, warp_size_int));
  sp_angl_equilibria.resize(roundUp(total_angl_params, warp_size_int));
  sp_dihe_amplitudes.resize(roundUp(total_dihe_params, warp_size_int));
  sp_dihe_periodicities.resize(roundUp(total_dihe_params, warp_size_int));
  sp_dihe_phase_angles.resize(roundUp(total_dihe_params, warp_size_int));
  ubrd_stiffnesses.resize(roundUp(total_ubrd_params, warp_size_int));
  ubrd_equilibria.resize(roundUp(total_ubrd_params, warp_size_int));
  cimp_stiffnesses.resize(roundUp(total_cimp_params, warp_size_int));
  cimp_stiffnesses.resize(roundUp(total_cimp_params, warp_size_int));
  cmap_surface_dimensions.resize(roundUp(total_cmap_surfaces, warp_size_int));
  cmap_surface_bounds.resize(roundUp(total_cmap_surfaces + 1, warp_size_int));
  cmap_surfaces.resize(roundUp(filtered_cmap_surf.size(), warp_size_zu));
  sp_ubrd_stiffnesses.resize(roundUp(total_ubrd_params, warp_size_int));
  sp_ubrd_equilibria.resize(roundUp(total_ubrd_params, warp_size_int));
  sp_cimp_stiffnesses.resize(roundUp(total_cimp_params, warp_size_int));
  sp_cimp_stiffnesses.resize(roundUp(total_cimp_params, warp_size_int));
  sp_cmap_surfaces.resize(roundUp(sp_filtered_cmap_surf.size(), warp_size_zu));
  bond_stiffnesses.putHost(filtered_bond_keq);
  bond_equilibria.putHost(filtered_bond_leq);
  angl_stiffnesses.putHost(filtered_angl_keq);
  angl_equilibria.putHost(filtered_angl_theta);
  dihe_amplitudes.putHost(filtered_dihe_amp);
  dihe_periodicities.putHost(filtered_dihe_freq);
  dihe_phase_angles.putHost(filtered_dihe_phi);
  sp_bond_stiffnesses.putHost(sp_filtered_bond_keq);
  sp_bond_equilibria.putHost(sp_filtered_bond_leq);
  sp_angl_stiffnesses.putHost(sp_filtered_angl_keq);
  sp_angl_equilibria.putHost(sp_filtered_angl_theta);
  sp_dihe_amplitudes.putHost(sp_filtered_dihe_amp);
  sp_dihe_periodicities.putHost(sp_filtered_dihe_freq);
  sp_dihe_phase_angles.putHost(sp_filtered_dihe_phi);
  ubrd_stiffnesses.putHost(filtered_ubrd_keq);
  ubrd_equilibria.putHost(filtered_ubrd_leq);
  cimp_stiffnesses.putHost(filtered_cimp_keq);
  cimp_phase_angles.putHost(filtered_cimp_phi);
  cmap_surface_dimensions.putHost(filtered_cmap_dim);
  cmap_surface_bounds.putHost(filtered_cmap_surf_bounds);
  cmap_surfaces.putHost(filtered_cmap_surf);
  sp_ubrd_stiffnesses.putHost(sp_filtered_ubrd_keq);
  sp_ubrd_equilibria.putHost(sp_filtered_ubrd_leq);
  sp_cimp_stiffnesses.putHost(sp_filtered_cimp_keq);
  sp_cimp_phase_angles.putHost(sp_filtered_cimp_phi);
  sp_cmap_surfaces.putHost(sp_filtered_cmap_surf);

  // Lay out the tables for valence and charge parameter indexing within the synthesis
  const int scm1 = system_count - 1;
  ubrd_param_idx.resize(ubrd_term_offsets.readHost(scm1) + ubrd_term_counts.readHost(scm1));
  cimp_param_idx.resize(cimp_term_offsets.readHost(scm1) + cimp_term_counts.readHost(scm1));
  cmap_param_idx.resize(cmap_term_offsets.readHost(scm1) + cmap_term_counts.readHost(scm1));
  bond_param_idx.resize(bond_term_offsets.readHost(scm1) + bond_term_counts.readHost(scm1));
  angl_param_idx.resize(angl_term_offsets.readHost(scm1) + angl_term_counts.readHost(scm1));
  dihe_param_idx.resize(dihe_term_offsets.readHost(scm1) + dihe_term_counts.readHost(scm1));
  charge_indices.resize(atom_offsets.readHost(scm1) + atom_counts.readHost(scm1));
                        
  // Loop back over all systems and copy the known mapping of individual topologies to the
  // synthesis as a whole
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[topology_indices.readHost(i)];
    const NonbondedKit<double> i_nbk   = iag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> i_vk      = iag_ptr->getDoublePrecisionValenceKit();
    const int tp_index = topology_indices.readHost(i);
    const int ag_bond_table_offset = topology_bond_table_offsets[tp_index];
    const int ag_angl_table_offset = topology_angl_table_offsets[tp_index];
    const int ag_dihe_table_offset = topology_dihe_table_offsets[tp_index];
    const int ag_ubrd_table_offset = topology_ubrd_table_offsets[tp_index];
    const int ag_cimp_table_offset = topology_cimp_table_offsets[tp_index];
    const int ag_cmap_table_offset = topology_cmap_table_offsets[tp_index];
    const int ag_chrg_table_offset = topology_chrg_table_offsets[tp_index];
    const int synth_bond_table_offset = bond_term_offsets.readHost(i);
    const int synth_angl_table_offset = angl_term_offsets.readHost(i);
    const int synth_dihe_table_offset = dihe_term_offsets.readHost(i);
    const int synth_ubrd_table_offset = ubrd_term_offsets.readHost(i);
    const int synth_cimp_table_offset = cimp_term_offsets.readHost(i);
    const int synth_cmap_table_offset = cmap_term_offsets.readHost(i);
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
  }
}
  
//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::extendLJMatrices() {

  // Search the diagonal entries first.  Use approximate comparisons with tolerances of 1.0e-5 in
  // an effort to accommodate very large Lennard-Jones A coefficients (~10^7).  Make a vector of
  // all unique diagonal entries.
  std::vector<double3> diag_entries;
  std::vector<int> lj_type_bounds(topology_count + 1, 0);
  for (int i = 0; i < topology_count; i++) {
    lj_type_bounds[i + 1] = lj_type_bounds[i] + topologies[i]->getAtomTypeCount();
  }
  const int total_types = lj_type_bounds[topology_count];
  std::vector<bool> lj_coverage(total_types, false);
  std::vector<int> lj_diagtypes(total_types);
  int n_unique_diags = 0;
  for (int i = 0; i < topology_count; i++) {
    const NonbondedKit<double> i_nbk = topologies[i]->getDoublePrecisionNonbondedKit();
    for (int j = lj_type_bounds[i]; j < lj_type_bounds[i + 1]; j++) {
      if (lj_coverage[j]) {
        continue;
      }
      const size_t matpos = (j - lj_type_bounds[i]) * (i_nbk.n_lj_types + 1);
      const Approx ent_a(i_nbk.lja_coeff[matpos], 1.0e-5);
      const Approx ent_b(i_nbk.ljb_coeff[matpos], 1.0e-5);
      const Approx ent_c(i_nbk.ljc_coeff[matpos], 1.0e-5);
      for (int k = i; k < topology_count; k++) {
        const int mstart = (k == i) ? j : lj_type_bounds[k];
        const NonbondedKit<double> k_nbk = topologies[k]->getDoublePrecisionNonbondedKit();
        for (int m = mstart; m < lj_type_bounds[k + 1]; m++) {
          const size_t testpos = (m - lj_type_bounds[k]) * (k_nbk.n_lj_types + 1);
          if (lj_coverage[m] == false && ent_a.test(k_nbk.lja_coeff[testpos]) &&
              ent_b.test(k_nbk.ljb_coeff[testpos]) && ent_c.test(k_nbk.ljc_coeff[testpos])) {
            lj_coverage[m] = true;
            lj_diagtypes[m] = n_unique_diags;
          }
        }
      }
      diag_entries.push_back({i_nbk.lja_coeff[matpos], i_nbk.ljb_coeff[matpos],
                              i_nbk.ljc_coeff[matpos]});
      n_unique_diags++;
    }
  }

  // CHECK
#if 0
  printf("Global LJ IDs (%d topologies) = [\n", topology_count);
  for (int i = 0; i < topology_count; i++) {
    printf("%4d -> %4d    ", topologies[i]->getAtomTypeCount(),
           lj_type_bounds[i + 1] - lj_type_bounds[i]);
  }
  printf("\n];\n");
  for (int i = 0; i < n_unique_diags; i++) {
    for (int j = 0; j < topology_count; j++) {
      if (lj_type_bounds[j + 1] - lj_type_bounds[j] > i) {
        printf("%4d ", lj_diagtypes[lj_type_bounds[j] + i]);
      }
      else {
        printf("     ");
      }
    }
    printf("\n");
  }
  printf("];\n");
#endif
  // END CHECK
  
  // Loop back over each topology and find all cross-terms involving the various diagonals.
  // Fill out a map of the numbers of unique cross-terms.
  std::vector<double3> lj_entries;
  std::vector<int> lj_entry_bounds(1, 0);
  std::vector<int2> lj_entry_map;
  int pos = 0;
  int n_lj_entry = 0;
  for (int i = 0; i < n_unique_diags; i++) {
    for (int j = i; j < n_unique_diags; j++) {

      // Each element p of the map stores the row and column for the Lennard-Jones table entries
      // catalogged between entries p and (p + 1) of the bounds array.
      lj_entry_map.push_back({i, j});
      for (int k = 0; k < topology_count; k++) {
        const NonbondedKit<double> k_nbk = topologies[k]->getDoublePrecisionNonbondedKit();
        for (int m = lj_type_bounds[k]; m < lj_type_bounds[k + 1]; m++) {
          if (lj_diagtypes[m] == i) {
            for (int n = lj_type_bounds[k]; n < lj_type_bounds[k + 1]; n++) {
              if (lj_diagtypes[n] == j) {
                const int matpos = (m - lj_type_bounds[k]) +
                                   ((n - lj_type_bounds[k]) * k_nbk.n_lj_types);
                const Approx ent_a(k_nbk.lja_coeff[matpos], 1.0e-5);
                const Approx ent_b(k_nbk.ljb_coeff[matpos], 1.0e-5);
                const Approx ent_c(k_nbk.ljc_coeff[matpos], 1.0e-5);
                bool found = false;
                for (int p = lj_entry_bounds[pos]; p < n_lj_entry; p++) {
                  found = (found || (ent_a.test(lj_entries[p].x) && ent_b.test(lj_entries[p].y) &&
                                     ent_c.test(lj_entries[p].z)));
                }
                if (found == false) {
                  lj_entries.push_back({k_nbk.lja_coeff[matpos], k_nbk.ljb_coeff[matpos],
                                        k_nbk.ljc_coeff[matpos]});
                  n_lj_entry++;
                }
              }
            }
          }
        }
      }
      lj_entry_bounds.push_back(n_lj_entry);
      pos++;
    }
  }

  // CHECK
#if 0
  printf("matrix_ideas = [\n");
  int jj = 0;
  for (int i = 0; i < n_unique_diags; i++) {
    for (int j = i; j < n_unique_diags; j++) {
      printf(" %d", lj_entry_bounds[jj + 1] - lj_entry_bounds[jj]);
      jj++;
    }
    printf("\n");
  }
  printf("];\n");
#endif
  // END CHECK
  
  // Loop back over the compiled matrix and find cases where a single type pair interaction
  // can have more than a single set of Lennard-Jones A, B, and C coefficients.  Split those
  // N unique sets of values across ceil(sqrt(N)) different sub-types of each atom.
  pos = 0;
  std::vector<int> n_type_copy(n_unique_diags + 1, 0);
  for (int i = 0; i < n_unique_diags; i++) {
    for (int j = i; j < n_unique_diags; j++) {
      const int n_takes = lj_entry_bounds[pos + 1] - lj_entry_bounds[pos];
      switch (n_takes) {
      case 1:
        n_type_copy[i] = std::max(n_type_copy[i], 1);
        break;
      case 4:
        n_type_copy[i] = std::max(n_type_copy[i], 2);
        break;
      case 9:
        n_type_copy[i] = std::max(n_type_copy[i], 3);
        break;
      default:
        n_type_copy[i] = std::max(n_type_copy[i], static_cast<int>(ceil(sqrt(n_takes))));
        break;
      }
      pos++;
    }
  }
  prefixSumInPlace<int>(&n_type_copy, PrefixSumType::EXCLUSIVE, "ExtendLJMatrices");

  // CHECK
#if 0
  printf("LJ types = [\n");
  for (int i = 0; i < n_unique_diags; i++) {
    printf("  %12.4lf %12.4lf %12.4lf  %2d\n", diag_entries[i].x, diag_entries[i].y,
           diag_entries[i].z, n_type_copy[i]);
  }
  printf("];\n");
#endif
  // END CHECK
  
  // Expand the matrices of Lennard-Jones atom types
  for (int i = 0; i < n_unique_diags; i++) {
    for (int j = n_type_copy[i]; j < n_type_copy[i + 1]; j++) {

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
  
} // namespace topology
} // namespace omni
