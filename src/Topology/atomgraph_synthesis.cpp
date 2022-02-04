#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "atomgraph_synthesis.h"
#include "atomgraph_abstracts.h"

namespace omni {
namespace topology {

using cuda::HybridKind;
using math::roundUp;
using math::maxAbsoluteDifference;
using math::maxValue;
using math::minValue;
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
    urey_bradley_parameter_indices{HybridKind::ARRAY, "tpsyn_urey_parm"},
    charmm_impr_parameter_indices{HybridKind::ARRAY, "tpsyn_cimp_parm"},
    cmap_surface_indices{HybridKind::ARRAY, "tpsyn_cmap_parm"},
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
    sp_nmr_k_final_values{HybridKind::ARRAY, "tpsyn_"},
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
  // Check that no system indexes a topology outside of the given list
  if (maxValue(topology_indices) >= topology_count || minValue(topology_indices) < 0) {
    rtErr("One or more systems references a topology index outside of the list supplied.",
          "AtomGraphSynthesis");
  }

  // Check that each topology in the supplied list is being used.  Roll that result into an
  // analysis of the uniqueness of each topology.
  std::vector<bool> topology_unique(topology_count, false);
  for (int i = 0; i < system_count; i++) {
    topology_unique[topology_indices.readHost(i)] = true;
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
    if (topology_unique[i]) {
      topology_index_rebase[i] = n_unique_top;
      topologies[n_unique_top] = topologies[i];
      n_unique_top++;
    }
  }
  topology_count = n_unique_top;
  
  // Allocate memory and set POINTER-kind arrays for the small packets of data
  const int padded_system_count = roundUp(system_count, warp_size_int);
  int_system_data.resize(32 * padded_system_count);
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
  last_solute_atoms.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  first_solvent_molecules.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  ubrd_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  cimp_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  cmap_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  bond_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  angl_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  dihe_term_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  virtual_site_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  atom_type_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  total_exclusion_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  rigid_water_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  bond_constraint_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  degrees_of_freedom.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  nonrigid_particle_counts.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  atom_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  residue_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  molecule_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  ubrd_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  cimp_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  cmap_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  bond_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  angl_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  dihe_term_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  virtual_site_offsets.setPointer(&int_system_data, 0, system_count);
  pivot += padded_system_count;
  nb_exclusion_offsets.setPointer(&int_system_data, 0, system_count);

  // Load the topology indexing first
  for (int i = 0; i < system_count; i++) {
    topology_indices.putHost(topology_index_rebase[topology_indices_in[i]], i);
  }

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
    const AtomGraph* ag_ptr = topologies[topology_indices.readHost(i)];
    const ChemicalDetailsKit cdk     = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk   = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk      = ag_ptr->getDoublePrecisionValenceKit();
    total_atoms += cdk.natom;
    atom_counts.putHost(cdk.natom, i);
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
    virtual_site_offsets.putHost(vste_offset, i);
    vste_offset += roundUp(ag_ptr->getVirtualSiteCount(), warp_size_int);
    nb_exclusion_offsets.putHost(excl_offset, i);
    excl_offset += roundUp(ag_ptr->getTotalExclusions(), warp_size_int);
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
  urey_bradley_parameter_indices.resize(ubrd_offset);
  charmm_impr_parameter_indices.resize(cimp_offset);
  cmap_surface_indices.resize(cmap_offset);
  bond_parameter_indices.resize(bond_offset);
  angl_parameter_indices.resize(angl_offset);
  dihe_parameter_indices.resize(dihe_offset);

  // Compute the numbers of unique parameters.  Take the opportunity to compute offsets (starting
  // bounds) for various sets of terms.
  int max_unique_atom = 0;
  int max_unique_bond = 0;
  int max_unique_angl = 0;
  int max_unique_dihe = 0;
  int max_unique_ubrd = 0;
  int max_unique_cimp = 0;
  int max_unique_cmap = 0;
  int max_unique_chrg = 0;
  int max_unique_atyp = 0;
  std::vector<int> topology_atom_offsets(topology_count);
  std::vector<int> topology_bond_table_offsets(topology_count);
  std::vector<int> topology_angl_table_offsets(topology_count);
  std::vector<int> topology_dihe_table_offsets(topology_count);
  std::vector<int> topology_ubrd_table_offsets(topology_count);
  std::vector<int> topology_cimp_table_offsets(topology_count);
  std::vector<int> topology_cmap_table_offsets(topology_count);
  std::vector<int> topology_chrg_table_offsets(topology_count);
  std::vector<int> topology_atyp_table_offsets(topology_count);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* ag_ptr = topologies[i];
    const ChemicalDetailsKit cdk   = ag_ptr->getChemicalDetailsKit();
    const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
    const ValenceKit<double> vk    = ag_ptr->getDoublePrecisionValenceKit();
    topology_atom_offsets[i] = max_unique_atom;
    topology_bond_table_offsets[i] = max_unique_bond;
    topology_angl_table_offsets[i] = max_unique_angl;
    topology_dihe_table_offsets[i] = max_unique_dihe;
    topology_ubrd_table_offsets[i] = max_unique_ubrd;
    topology_cimp_table_offsets[i] = max_unique_cimp;
    topology_cmap_table_offsets[i] = max_unique_cmap;    
    topology_chrg_table_offsets[i] = max_unique_chrg;
    topology_atyp_table_offsets[i] = max_unique_atyp;
    max_unique_atom += cdk.natom;
    max_unique_bond += vk.nbond_param;
    max_unique_angl += vk.nangl_param;
    max_unique_dihe += vk.ndihe_param;
    max_unique_ubrd += vk.nubrd_param;
    max_unique_cimp += vk.ncimp_param;
    max_unique_cmap += vk.ncmap_surf;
    max_unique_chrg += nbk.n_q_types;
    max_unique_atyp += nbk.n_lj_types;
  }

  // Pre-compute some quantities relating to CMAPs that will help distinguish these gargantuan
  // "parameters" in the inner loops that follow.
  std::vector<double> cmap_parameter_sums(max_unique_cmap);
  for (int i = 0; i < topology_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    const ValenceKit<double> i_vk = iag_ptr->getDoublePrecisionValenceKit();
    for (int j = 0; j < i_vk.ncmap_surf; j++) {
      double cm_sum = 0.0;
      const int cmj_llim = i_vk.cmap_surf_bounds[j];
      const int cmj_hlim = cmj_llim + (i_vk.cmap_dim[j] * i_vk.cmap_dim[j]);
      for (int k = cmj_llim; k < cmj_hlim; k++) {
        cm_sum += i_vk.cmap_surf[k];
      }
      cmap_parameter_sums[topology_cmap_table_offsets[i] + j] = cm_sum;
    }
  }
  
  // Create lists of unique parameters for the valence and non-bonded calculations.
  std::vector<int> bond_synthesis_index(max_unique_bond, false);
  std::vector<int> angl_synthesis_index(max_unique_angl, false);
  std::vector<int> dihe_synthesis_index(max_unique_dihe, false);
  std::vector<int> ubrd_synthesis_index(max_unique_ubrd, false);
  std::vector<int> cimp_synthesis_index(max_unique_cimp, false);
  std::vector<int> cmap_synthesis_index(max_unique_cmap, false);
  std::vector<int> chrg_synthesis_index(max_unique_chrg, false);
  std::vector<int> atyp_synthesis_index(max_unique_atyp, false);
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
  int n_unique_atyp = 0;
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
              ij_bond_keq.test(k_vk.bond_keq[m]) == false &&
              ij_bond_leq.test(k_vk.bond_leq[m]) == false) {
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
              ij_angl_keq.test(k_vk.angl_keq[m]) == false &&
              ij_angl_theta.test(k_vk.angl_theta[m]) == false) {
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
              ij_dihe_amp.test(k_vk.dihe_amp[m]) == false &&
              ij_dihe_freq.test(k_vk.dihe_freq[m]) == false &&
              ij_dihe_phi.test(k_vk.dihe_phi[m]) == false) {
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
              ij_ubrd_keq.test(k_vk.ubrd_keq[m]) == false &&
              ij_ubrd_leq.test(k_vk.ubrd_leq[m]) == false) {
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
              ij_cimp_keq.test(k_vk.cimp_keq[m]) == false &&
              ij_cimp_phi.test(k_vk.cimp_phi[m]) == false) {
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
              ij_cmap_sum.test(cmap_parameter_sums[topology_cmap_table_offsets[k] + m]) == false &&
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

  // With the unique parameters enumerated and maps leading from parameters in any individual
  // system into the unified arrays within the synthesis, make collated arrays for each atom and
  // energy term.
  
}

} // namespace topology
} // namespace omni
