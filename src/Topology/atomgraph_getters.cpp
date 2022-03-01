#include <cmath>
#include <cstdio>
#include <climits>
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"

namespace omni {
namespace topology {

using card::HybridTargetLevel;
using math::findBin;

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
  return findBin(residue_limits.data(), atom_index, residue_count);
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
ValenceKit<double> AtomGraph::getDoublePrecisionValenceKit(HybridTargetLevel tier) const {
  return ValenceKit<double>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
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
  return ValenceKit<float>(atom_count, bond_term_count, angl_term_count, dihe_term_count,
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

} // namespace topology
} // namespace omni
