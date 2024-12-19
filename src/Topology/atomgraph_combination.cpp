#include <ctime>
#include <unistd.h>
#include "copyright.h"
#include "atomgraph.h"
#include "lennard_jones_analysis.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const int n_a, const AtomGraph &ag_b, const int n_b,
                     const ExceptionResponse policy) :
    AtomGraph()
{
  snprintf(version_stamp, 16, "STORMM 0.1");
  std::time_t raw_time = std::time(nullptr);
  std::tm* current_time = std::localtime(&raw_time);
  date = *current_time;
  title = ag_a.title;
  source = "Combination of ";
  if (n_a > 1) {
    source += std::to_string(n_a) + "x " + ag_a.source;
  }
  else {
    source += ag_a.source;
  }
  source += " and ";
  if (n_b > 1) {
    source += std::to_string(n_b) + "x " + ag_b.source;
  }
  else {
    source += ag_b.source;
  }
  force_fields = ag_a.force_fields;
  force_fields.insert(force_fields.end(), ag_b.force_fields.begin(), ag_b.force_fields.end());
  atom_count = (n_a * ag_a.atom_count) + (n_b * ag_b.atom_count);
  residue_count = (n_a * ag_a.residue_count) + (n_b * ag_b.residue_count);
  molecule_count  = (n_a * ag_a.molecule_count) + (n_b * ag_b.molecule_count);
  largest_residue_size = std::max(ag_a.largest_residue_size, ag_b.largest_residue_size);

  // After this point, it must be decided how the topologies are to be combined.  The procedure
  // will be to place solute atoms before solvent atoms.
  last_solute_residue = (n_a * ag_a.last_solute_residue) + (n_b * ag_b.last_solute_residue);
  last_solute_atom = (n_a * ag_a.last_solute_atom) + (n_b * ag_b.last_solute_atom);
  first_solvent_molecule = (n_a * ag_a.first_solvent_molecule) +
                           (n_b * ag_b.first_solvent_molecule);
  last_atom_before_cap = (n_a * ag_a.last_atom_before_cap) + (n_b * ag_b.last_atom_before_cap);
  implicit_copy_count = 0;
  largest_molecule_size = std::max(ag_a.largest_molecule_size, ag_b.largest_molecule_size);
  water_residue_size = std::max(ag_a.water_residue_size, ag_b.water_residue_size);
  water_residue_count = (n_a * ag_a.water_residue_count) + (n_b * ag_b.water_residue_count);
  unconstrained_dof = (n_a * ag_a.unconstrained_dof) + (n_b * ag_b.unconstrained_dof);
  constrained_dof = (n_a * ag_a.constrained_dof) + (n_b * ag_b.constrained_dof);

  // Form consensus tables of bond, angle, dihedral, Urey-Bradley, CHARMM improper, and CMAP terms.
  // Try to preserve the order of parameter indices in each topology as much as possible.
  std::vector<int2> sysb_bond_parm_idx_map;
  const ValenceKit<double> vk_a = ag_a.getDoublePrecisionValenceKit();
  const ValenceKit<double> vk_b = ag_b.getDoublePrecisionValenceKit();
  ParameterUnion<double> uni_bonds(vk_a.bond_keq, vk_a.bond_leq, vk_a.nbond_param,
                                   vk_b.bond_keq, vk_b.bond_leq, vk_b.nbond_param);
  ParameterUnion<double> uni_angls(vk_a.angl_keq, vk_a.angl_theta, vk_a.nangl_param,
                                   vk_b.angl_keq, vk_b.angl_theta, vk_b.nangl_param);
  ParameterUnion<double> uni_dihes(vk_a.dihe_amp, vk_a.dihe_freq, vk_a.dihe_phi, vk_a.nangl_param,
                                   vk_b.dihe_amp, vk_b.dihe_freq, vk_b.dihe_phi, vk_b.nangl_param);
  ParameterUnion<double> uni_hbonds(ag_a.hbond_a_values.data(), ag_a.hbond_b_values.data(),
                                    ag_a.hbond_10_12_parameter_count, ag_b.hbond_a_values.data(),
                                    ag_b.hbond_b_values.data(), ag_b.hbond_10_12_parameter_count);
  const NonbondedKit<double> nbk_a = ag_a.getDoublePrecisionNonbondedKit();
  const NonbondedKit<double> nbk_b = ag_b.getDoublePrecisionNonbondedKit();
  const std::vector<std::vector<char4>> atyp_list_a = ag_a.getAtomTypeNameTable();
  const std::vector<std::vector<char4>> atyp_list_b = ag_b.getAtomTypeNameTable();
  LennardJonesAnalysis uni_ljtab(nbk_a, atyp_list_a);
  uni_ljtab.addSet(nbk_b, atyp_list_b);

  // CHECK
  printf("A: %s\n", ag_a.getFileName().c_str());
  printf("B: %s\n", ag_b.getFileName().c_str());
  printf("There are %2d + %2d = %2d Lennard-Jones types\n", nbk_a.n_lj_types, nbk_b.n_lj_types,
         uni_ljtab.getLJTypeCount());
  // END CHECK
  
  // The arrays that enter into the topology must be constructed from scratch in order to feed them
  // to the loadHybridArrays() function used by other constructors.
  const int n_desc = static_cast<int>(TopologyDescriptor::N_VALUES);
  std::vector<int> tmp_desc(n_desc);
  for (int i = 0; i < n_desc; i++) {
    switch (static_cast<TopologyDescriptor>(i)) {
    case TopologyDescriptor::ATOM_COUNT:
      tmp_desc[i] = atom_count;
      break;
    case TopologyDescriptor::ATOM_TYPE_COUNT:

      // Draw upon the consensus Lennard-Jones tables for the two topologies.
      
      break;
    case TopologyDescriptor::BONDS_WITH_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.bond_term_with_hydrogen) + (n_b * ag_b.bond_term_with_hydrogen);
      break;
    case TopologyDescriptor::BONDS_WITHOUT_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.bond_term_without_hydrogen) +
                    (n_b * ag_b.bond_term_without_hydrogen);
      break;
    case TopologyDescriptor::ANGLES_WITH_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.angl_term_with_hydrogen) + (n_b * ag_b.angl_term_with_hydrogen);
      break;
    case TopologyDescriptor::ANGLES_WITHOUT_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.angl_term_without_hydrogen) +
                    (n_b * ag_b.angl_term_without_hydrogen);
      break;
    case TopologyDescriptor::DIHEDRALS_WITH_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.dihe_term_with_hydrogen) + (n_b * ag_b.dihe_term_with_hydrogen);
      break;
    case TopologyDescriptor::DIHEDRALS_WITHOUT_HYDROGEN:
      tmp_desc[i] = (n_a * ag_a.dihe_term_without_hydrogen) +
                    (n_b * ag_b.dihe_term_without_hydrogen);
      break;
    case TopologyDescriptor::NHPARM_UNUSED:
      tmp_desc[i] = (n_a * ag_a.unused_nhparm) + (n_b * ag_b.unused_nhparm);
      break;
    case TopologyDescriptor::ADDLES_CREATED:
      tmp_desc[i] = (n_a * ag_a.unused_nparm) + (n_b * ag_b.unused_nparm);
      break;
    case TopologyDescriptor::TOTAL_EXCLUDED_ATOMS:
      tmp_desc[i] = (n_a * ag_a.total_exclusions) + (n_b * ag_b.total_exclusions);
      break;
    case TopologyDescriptor::RESIDUE_COUNT:
      tmp_desc[i] = (n_a * ag_a.residue_count) + (n_b * ag_b.residue_count);
      break;
    case TopologyDescriptor::NBONA_UNUSED:
      tmp_desc[i] = (n_a * ag_a.heavy_bonds_plus_constraints) +
                    (n_b * ag_b.heavy_bonds_plus_constraints);
      break;
    case TopologyDescriptor::NTHETA_UNUSED:
      tmp_desc[i] = (n_a * ag_a.heavy_angls_plus_constraints) +
                    (n_b * ag_b.heavy_angls_plus_constraints);
      break;
    case TopologyDescriptor::NPHIA_UNUSED:
      tmp_desc[i] = (n_a * ag_a.heavy_dihes_plus_constraints) +
                    (n_b * ag_b.heavy_dihes_plus_constraints);
      break;
    case TopologyDescriptor::BOND_TYPE_COUNT:

      // Draw upon the consensus bond parameter tables.
      tmp_desc[i] = uni_bonds.getUniqueParameterCount();
      break;
    case TopologyDescriptor::ANGLE_TYPE_COUNT:

      // Draw upon the consensus angle parameter tables.
      tmp_desc[i] = uni_angls.getUniqueParameterCount();
      break;
    case TopologyDescriptor::DIHEDRAL_TYPE_COUNT:

      // Draw upon the consensus dihedral parameter tables.
      tmp_desc[i] = uni_dihes.getUniqueParameterCount();
      break;
    case TopologyDescriptor::NATYP_UNUSED:
      tmp_desc[i] = (n_a * ag_a.unused_natyp) + (n_b * ag_b.unused_natyp);
      break;
    case TopologyDescriptor::NPHB_UNUSED:

      // Draw upon the consensus hydrogen bonding parameter tables
      tmp_desc[i] = uni_hbonds.getUniqueParameterCount();
      break;
    case TopologyDescriptor::PERTURBATION:
      tmp_desc[i] = (ag_a.use_perturbation_info == PerturbationSetting::ON ||
                     ag_b.use_perturbation_info == PerturbationSetting::ON);
      break;
    case TopologyDescriptor::BOND_PERTURBATIONS:
      tmp_desc[i] = (n_a * ag_a.bond_perturbation_term_count) +
                    (n_b * ag_b.bond_perturbation_term_count);
      break;
    case TopologyDescriptor::ANGLE_PERTURBATIONS:
      tmp_desc[i] = (n_a * ag_a.angl_perturbation_term_count) +
                    (n_b * ag_b.angl_perturbation_term_count);
      break;
    case TopologyDescriptor::DIHEDRAL_PERTURBATIONS:
      tmp_desc[i] = (n_a * ag_a.dihe_perturbation_term_count) +
                    (n_b * ag_b.dihe_perturbation_term_count);
      break;
    case TopologyDescriptor::BONDS_IN_PERTURBED_GROUP:
      tmp_desc[i] = (n_a * ag_a.bonds_in_perturbed_group) + (n_b * ag_b.bonds_in_perturbed_group);
      break;
    case TopologyDescriptor::ANGLES_IN_PERTURBED_GROUP:
      tmp_desc[i] = (n_a * ag_a.angls_in_perturbed_group) + (n_b * ag_b.angls_in_perturbed_group);
      break;
    case TopologyDescriptor::DIHEDRALS_IN_PERTURBED_GROUP:
      tmp_desc[i] = (n_a * ag_a.dihes_in_perturbed_group) + (n_b * ag_b.dihes_in_perturbed_group);
      break;
    case TopologyDescriptor::BOX_TYPE_INDEX:

      // Type promotion: while it could be infeasible to combine certain systems' coordinates,
      // assume that the system with the "higher order" unit cell dictates the unit cell of the
      // combined system.
      int a_cls_idx, b_cls_idx;
      switch (ag_a.periodic_box_class) {
      case UnitCellType::NONE:
        a_cls_idx = 0;
        break;
      case UnitCellType::ORTHORHOMBIC:
        a_cls_idx = 1;
        break;
      case UnitCellType::TRICLINIC:
        a_cls_idx = 2;
        break;
      }
      switch (ag_b.periodic_box_class) {
      case UnitCellType::NONE:
        b_cls_idx = 0;
        break;
      case UnitCellType::ORTHORHOMBIC:
        b_cls_idx = 1;
        break;
      case UnitCellType::TRICLINIC:
        b_cls_idx = 2;
        break;
      }
      tmp_desc[i] = std::max(a_cls_idx, b_cls_idx);
      break;
    case TopologyDescriptor::ATOM_COUNT_LARGEST_RESIDUE:
      tmp_desc[i] = std::max(ag_a.largest_residue_size, ag_b.largest_residue_size);
      break;
    case TopologyDescriptor::CAP:
      tmp_desc[i] = (ag_a.use_solvent_cap_option == SolventCapSetting::ON ||
                     ag_b.use_solvent_cap_option == SolventCapSetting::ON);
      break;
    case TopologyDescriptor::EXTRA_POINT_COUNT:
      tmp_desc[i] = (n_a * ag_a.virtual_site_count) + (n_b * ag_b.virtual_site_count);
      break;
    case TopologyDescriptor::PIMD_SLICE_COUNT:
      tmp_desc[i] = std::max(ag_a.implicit_copy_count, ag_b.implicit_copy_count);
      break;
    case TopologyDescriptor::N_VALUES:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const AtomGraph &ag_b, const int n_b,
                     const ExceptionResponse policy) :
    AtomGraph(ag_a, 1, ag_b, n_b, policy)
{}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const AtomGraph &ag_b,
                     const ExceptionResponse policy) :
    AtomGraph(ag_a, 1, ag_b, 1, policy)
{
}
  
} // namespace topology
} // namespace stormm
