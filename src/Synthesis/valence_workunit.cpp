#include <algorithm>
#include "Math/statistics.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Topology/topology_util.h"
#include "valence_workunit.h"

// CHECK
#include "Parsing/parse.h"
// END CHECK

namespace omni {
namespace synthesis {

using math::DataOrder;
using math::locateValue;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::reduceUniqueValues;
using topology::VirtualSiteKind;
using topology::markAffectorAtoms;
using topology::extractBoundedListEntries;
using topology::writeAtomList;
  
//-------------------------------------------------------------------------------------------------
ValenceDelegator::ValenceDelegator(const AtomGraph *ag_in, const RestraintApparatus *ra_in) :
    atom_count{ag_in->getAtomCount()},
    first_unassigned_atom{0}, max_presence_allocation{2},
    bond_affector_list{}, bond_affector_bounds{}, angl_affector_list{}, angl_affector_bounds{},
    dihe_affector_list{}, dihe_affector_bounds{}, ubrd_affector_list{}, ubrd_affector_bounds{},
    cimp_affector_list{}, cimp_affector_bounds{}, cmap_affector_list{}, cmap_affector_bounds{},
    infr_affector_list{}, infr_affector_bounds{}, vste_affector_list{}, vste_affector_bounds{},
    cnst_affector_list{}, cnst_affector_bounds{}, sett_affector_list{}, sett_affector_bounds{},
    work_unit_assignment_count{}, work_unit_presence{}, assigned_update_work_units{},
    ag_pointer{ag_in}, ra_pointer{ra_in}
{
  // Get relevant abstracts
  const ValenceKit<double> vk = ag_in->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_in->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_in->getDoublePrecisionConstraintKit();
  const RestraintApparatusDpReader rar = ra_in->dpData();
  
  // Allocate and fill the arrays
  allocate();
  fillAffectorArrays(vk, vsk, cnk, rar);
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getAtomAssignmentCount(const int atom_index) const {
  return work_unit_assignment_count[atom_index];
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getFirstUnassignedAtom() const {
  return first_unassigned_atom;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getBondAffectors(const int atom_index) const {
  return extractBoundedListEntries(bond_affector_list, bond_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getAngleAffectors(const int atom_index) const {
  return extractBoundedListEntries(angl_affector_list, angl_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getDihedralAffectors(const int atom_index) const {
  return extractBoundedListEntries(dihe_affector_list, dihe_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getUreyBradleyAffectors(const int atom_index) const {
  return extractBoundedListEntries(ubrd_affector_list, ubrd_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getCharmmImproperAffectors(const int atom_index) const {
  return extractBoundedListEntries(cimp_affector_list, cimp_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getCmapAffectors(const int atom_index) const {
  return extractBoundedListEntries(cmap_affector_list, cmap_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getInferred14Affectors(const int atom_index) const {
  return extractBoundedListEntries(infr_affector_list, infr_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getPositionalRestraintAffectors(const int atom_index) const {
  return extractBoundedListEntries(rposn_affector_list, rposn_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getDistanceRestraintAffectors(const int atom_index) const {
  return extractBoundedListEntries(rbond_affector_list, rbond_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getAngleRestraintAffectors(const int atom_index) const {
  return extractBoundedListEntries(rangl_affector_list, rangl_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getDihedralRestraintAffectors(const int atom_index) const {
  return extractBoundedListEntries(rdihe_affector_list, rdihe_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getSettleGroupAffectors(const int atom_index) const {
  return extractBoundedListEntries(sett_affector_list, sett_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getConstraintGroupAffectors(const int atom_index) const {
  return extractBoundedListEntries(cnst_affector_list, cnst_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getVirtualSiteAffectors(const int atom_index) const {
  return extractBoundedListEntries(vste_affector_list, vste_affector_bounds, atom_index);
}

//-------------------------------------------------------------------------------------------------
int ValenceDelegator::getUpdateWorkUnit(const int atom_index) const {
  return assigned_update_work_units[atom_index];
}
  
//-------------------------------------------------------------------------------------------------
const AtomGraph* ValenceDelegator::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* ValenceDelegator::getRestraintApparatusPointer() const {
  return ra_pointer;
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::checkPresence(const int atom_index, const int vwu_index) const {
  bool result = false;
  const int offset = atom_index * max_presence_allocation;
  for (int i = 0; i < work_unit_assignment_count[atom_index]; i++) {
    result = (result || work_unit_presence[offset + i] == vwu_index);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markAtomAddition(const int vwu_index, const int atom_index) {
  const int atom_wua = work_unit_assignment_count[atom_index];

  // If necessary, resize and repaginate the record of which work units each atom takes part in
  if (atom_wua == max_presence_allocation) {
    work_unit_presence.resize(atom_count * (max_presence_allocation + 1));
    for (int i = atom_count - 1; i >= 0; i++) {
      const int new_offset = (max_presence_allocation + 1) * i;
      const int old_offset = max_presence_allocation * i;
      for (int j = 0; j < max_presence_allocation; j++) {
        work_unit_presence[new_offset + j] = work_unit_presence[old_offset + j];
      }
    }
    max_presence_allocation += 1;
  }

  // Mark the new work unit to which the indexed atom has been assigned
  work_unit_presence[(max_presence_allocation * atom_index) + atom_wua] = vwu_index;
  work_unit_assignment_count[atom_index] = atom_wua + 1;

  // Update the first unassigned atom, if appropriate
  while (first_unassigned_atom < atom_count &&
         work_unit_assignment_count[first_unassigned_atom] > 0) {
    first_unassigned_atom += 1;
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markBondAddition(const int vwu_index, const int bond_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markAngleAddition(const int vwu_index, const int angl_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markDihedralAddition(const int vwu_index, const int dihe_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markUreyBradleyAddition(const int vwu_index, const int ubrd_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markCharmmImproperAddition(const int vwu_index, const int cimp_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markCmapAddition(const int vwu_index, const int cmap_term_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markPositionalRestraintAddition(const int vwu_index,
                                                       const int posn_rstr_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markDistanceRestraintAddition(const int vwu_index,
                                                     const int dist_rstr_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markAngleRestraintAddition(const int vwu_index, const int angl_rstr_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markDihedralRestraintAddition(const int vwu_index,
                                                     const int dihe_rstr_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markConstraintGroupAddition(const int vwu_index,
                                                   const int cnst_group_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markSettleGroupAddition(const int vwu_index, const int sett_group_index) {

}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::markVirtualSiteAddition(const int vwu_index, const int vsite_index) {

}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setUpdateWorkUnit(const int atom_index, const int vwu_index) {
  assigned_update_work_units[atom_index] = vwu_index;
}

//-------------------------------------------------------------------------------------------------
void ValenceDelegator::allocate() {
 
  // Allocate memory as needed
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  const RestraintApparatusDpReader rar = ra_pointer->dpData();
  bond_affector_list.resize(2 * vk.nbond);
  bond_affector_bounds.resize(atom_count + 1, 0);
  angl_affector_list.resize(3 * vk.nangl);
  angl_affector_bounds.resize(atom_count + 1, 0);
  dihe_affector_list.resize(4 * vk.ndihe);
  dihe_affector_bounds.resize(atom_count + 1, 0);
  ubrd_affector_list.resize(2 * vk.nubrd);
  ubrd_affector_bounds.resize(atom_count + 1, 0);
  cimp_affector_list.resize(4 * vk.ncimp);
  cimp_affector_bounds.resize(atom_count + 1, 0);
  cmap_affector_list.resize(5 * vk.ncmap);
  cmap_affector_bounds.resize(atom_count + 1, 0);
  infr_affector_list.resize(2 * vk.ninfr14);
  infr_affector_bounds.resize(atom_count + 1, 0);
  rposn_affector_list.resize(rar.nposn);
  rposn_affector_bounds.resize(atom_count + 1, 0);
  rbond_affector_list.resize(2 * rar.nbond);
  rbond_affector_bounds.resize(atom_count + 1, 0);
  rangl_affector_list.resize(3 * rar.nangl);
  rangl_affector_bounds.resize(atom_count + 1, 0);
  rdihe_affector_list.resize(4 * rar.ndihe);
  rdihe_affector_bounds.resize(atom_count + 1, 0);
  vste_affector_list.resize(5 * vsk.nsite);
  vste_affector_bounds.resize(atom_count + 1, 0);
  sett_affector_list.resize(3 * cnk.nsettle);
  sett_affector_bounds.resize(atom_count + 1, 0);
  cnst_affector_list.resize(cnk.group_bounds[cnk.ngroup], -1);
  cnst_affector_bounds.resize(atom_count + 1, 0);  
  work_unit_assignment_count.resize(atom_count, 0);
  work_unit_presence.resize(max_presence_allocation * atom_count);
  assigned_update_work_units.resize(atom_count, -1);
}
  
//-------------------------------------------------------------------------------------------------
void ValenceDelegator::fillAffectorArrays(const ValenceKit<double> &vk,
                                          const VirtualSiteKit<double> &vsk,
                                          const ConstraintKit<double> &cnk,
                                          const RestraintApparatusDpReader &rar) {

  // Pass through the topology, filling out the valence term affector arrays and the virtual site
  // frame atom arrays.
  markAffectorAtoms(&bond_affector_bounds, &bond_affector_list, vk.nbond, vk.bond_i_atoms,
                    vk.bond_j_atoms);
  markAffectorAtoms(&angl_affector_bounds, &angl_affector_list, vk.nangl, vk.angl_i_atoms,
                    vk.angl_j_atoms, vk.angl_k_atoms);
  markAffectorAtoms(&dihe_affector_bounds, &dihe_affector_list, vk.ndihe, vk.dihe_i_atoms,
                    vk.dihe_j_atoms, vk.dihe_k_atoms, vk.dihe_l_atoms);
  markAffectorAtoms(&ubrd_affector_bounds, &ubrd_affector_list, vk.nubrd, vk.ubrd_i_atoms,
                    vk.ubrd_k_atoms);
  markAffectorAtoms(&cimp_affector_bounds, &cimp_affector_list, vk.ncimp, vk.cimp_i_atoms,
                    vk.cimp_j_atoms, vk.cimp_k_atoms, vk.cimp_l_atoms);
  markAffectorAtoms(&cmap_affector_bounds, &cmap_affector_list, vk.ncmap, vk.cmap_i_atoms,
                    vk.cmap_j_atoms, vk.cmap_k_atoms, vk.cmap_l_atoms, vk.cmap_m_atoms);
  markAffectorAtoms(&infr_affector_bounds, &infr_affector_list, vk.ninfr14, vk.infr14_i_atoms,
                    vk.infr14_l_atoms);
  markAffectorAtoms(&sett_affector_bounds, &sett_affector_list, cnk.nsettle, cnk.settle_ox_atoms,
                    cnk.settle_h1_atoms, cnk.settle_h2_atoms);
  markAffectorAtoms(&rposn_affector_bounds, &rposn_affector_list, rar.nposn, rar.rposn_atoms);
  markAffectorAtoms(&rbond_affector_bounds, &rbond_affector_list, rar.nbond, rar.rbond_i_atoms,
                    rar.rbond_j_atoms);
  markAffectorAtoms(&rangl_affector_bounds, &rangl_affector_list, rar.nangl, rar.rangl_i_atoms,
                    rar.rangl_j_atoms, rar.rangl_k_atoms);
  markAffectorAtoms(&rdihe_affector_bounds, &rdihe_affector_list, rar.ndihe, rar.rdihe_i_atoms,
                    rar.rdihe_j_atoms, rar.rdihe_k_atoms, rar.rdihe_l_atoms);

  // Virtual sites can contain -1 in the atom arrays and therefore cannot go through the
  // general procedure.  Best to handle them with a filter over the frame type, anyway.
  for (int pos = 0; pos < vsk.nsite; pos++) {

    // There will be a minimum of two frame atoms, but the virtual site itself is part of the
    // "term" that puts forces on particles.  The non-bonded forces thus far determined to act
    // on the virtual site trasmit onto as many as four frame atoms.
    vste_affector_bounds[vsk.vs_atoms[pos]] += 1;
    vste_affector_bounds[vsk.frame1_idx[pos]] += 1;
    vste_affector_bounds[vsk.frame2_idx[pos]] += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      break;
    case VirtualSiteKind::FIXED_4:
      vste_affector_bounds[vsk.frame3_idx[pos]] += 1;
      vste_affector_bounds[vsk.frame4_idx[pos]] += 1;
      break;
    }
  }

  // Constraint groups are not of any pre-defined size and therefore cannot go through the
  // general procedure.
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      cnst_affector_bounds[cnk.group_list[j]] += 1;
    }
  }
  prefixSumInPlace<int>(&vste_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");
  prefixSumInPlace<int>(&cnst_affector_bounds, PrefixSumType::EXCLUSIVE, "ValenceDelegator");

  // Finish up the virtual site mapping
  for (int pos = 0; pos < vsk.nsite; pos++) {
    const int virtual_atom = vsk.vs_atoms[pos];
    const int parent_atom = vsk.frame1_idx[pos];
    const int frame2_atom = vsk.frame2_idx[pos];
    const int list_vs_idx = vste_affector_bounds[virtual_atom];
    const int list_p_idx  = vste_affector_bounds[parent_atom];
    const int list_f2_idx = vste_affector_bounds[frame2_atom];
    vste_affector_list[list_vs_idx] = pos;
    vste_affector_list[list_p_idx ] = pos;
    vste_affector_list[list_f2_idx] = pos;
    vste_affector_bounds[virtual_atom] += 1;
    vste_affector_bounds[parent_atom]  += 1;
    vste_affector_bounds[frame2_atom]  += 1;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vsk.vs_param_idx[pos]])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom = vsk.frame3_idx[pos];
        const int frame4_atom = vsk.frame4_idx[pos];
        const int list_f3_idx = vste_affector_bounds[frame3_atom];
        const int list_f4_idx = vste_affector_bounds[frame4_atom];
        vste_affector_list[list_f3_idx] = pos;
        vste_affector_list[list_f4_idx] = pos;
        vste_affector_bounds[frame3_atom] += 1;
        vste_affector_bounds[frame4_atom] += 1;
      }
      break;
    }
  }

  // Finish up the constraint group mapping
  for (int pos = 0; pos < cnk.ngroup; pos++) {
    for (int j = cnk.group_bounds[pos]; j < cnk.group_bounds[pos + 1]; j++) {
      const int c_atom = cnk.group_list[j];
      const int list_c_idx = cnst_affector_bounds[c_atom];
      cnst_affector_list[list_c_idx] = pos;
      cnst_affector_bounds[c_atom] += 1;
    }
  }

  // Rewind the final prefix sums after using them to populate the virtual site and constraint
  // affector lists.
  const int natom = atom_count;
  for (int i = natom; i > 0; i--) {
    vste_affector_bounds[i]  = vste_affector_bounds[i - 1];
    cnst_affector_bounds[i]  = cnst_affector_bounds[i - 1];
  }
  vste_affector_bounds[0] = 0;
  cnst_affector_bounds[0] = 0;
}
  
//-------------------------------------------------------------------------------------------------
ValenceWorkUnit::ValenceWorkUnit(ValenceDelegator *vdel_in, const int list_index_in,
                                 const int seed_atom_in, const int max_atoms_in) :
    atom_count{0}, bond_term_count{0}, angl_term_count{0}, dihe_term_count{0}, ubrd_term_count{0},
    cimp_term_count{0}, cmap_term_count{0}, rposn_term_count{0}, rbond_term_count{0},
    rangl_term_count{0}, rdihe_term_count{0}, cnst_group_count{0}, sett_group_count{0},
    vste_count{0}, list_index{list_index_in}, min_atom_index{-1}, max_atom_index{-1},
    atom_limit{max_atoms_in}, atom_import_list{}, bond_term_list{}, angl_term_list{},
    dihe_term_list{}, ubrd_term_list{}, cimp_term_list{}, cmap_term_list{}, bond_i_atoms{},
    bond_j_atoms{}, angl_i_atoms{}, angl_j_atoms{}, angl_k_atoms{}, dihe_i_atoms{}, dihe_j_atoms{},
    dihe_k_atoms{}, dihe_l_atoms{}, ubrd_i_atoms{}, ubrd_k_atoms{}, cimp_i_atoms{}, cimp_j_atoms{},
    cimp_k_atoms{}, cimp_l_atoms{}, cmap_i_atoms{}, cmap_j_atoms{}, cmap_k_atoms{}, cmap_l_atoms{},
    cmap_m_atoms{}, rposn_term_list{}, rbond_term_list{}, rangl_term_list{}, rdihe_term_list{},
    rposn_atoms{}, rbond_i_atoms{}, rbond_j_atoms{}, rangl_i_atoms{}, rangl_j_atoms{},
    rangl_k_atoms{}, rdihe_i_atoms{}, rdihe_j_atoms{}, rdihe_k_atoms{}, rdihe_l_atoms{},
    cnst_group_list{}, sett_group_list{}, cnst_group_atoms{}, cnst_group_bounds{}, sett_ox_atoms{},
    sett_h1_atoms{}, sett_h2_atoms{}, virtual_site_list{}, vsite_atoms{}, vsite_parent_atoms{},
    vsite_frame2_atoms{}, vsite_frame3_atoms{}, vsite_frame4_atoms{}, vdel_pointer{vdel_in},
    ag_pointer{vdel_in->getTopologyPointer()}, ra_pointer{vdel_in->getRestraintApparatusPointer()}
{
  // Check the atom bounds
  if (atom_limit < minimum_valence_work_unit_atoms ||
      atom_limit > maximum_valence_work_unit_atoms) {
    const std::string err_msg = (atom_limit < minimum_valence_work_unit_atoms) ?
                                "is too small to make effective work units." :
                                "could lead to work units that do not fit in HPC resources.";
    rtErr("The maximum allowed number of atoms should be between " +
          std::to_string(minimum_valence_work_unit_atoms) + " and " +
          std::to_string(maximum_valence_work_unit_atoms) + ".  A value of " +
          std::to_string(atom_limit) + err_msg, "ValenceWorkUnit");
  }

  // Reserve space for atoms
  atom_import_list.reserve(atom_limit);
  
  // Unpack the original topology
  ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  
  // Starting with the seed atom, branch out and add new atoms until approaching the maximum
  // allowed number of atoms.  If the atom is part of some molecule that can all be accommodated
  // by the one work unit, include all atoms.
  std::vector<int> candidate_additions(1, seed_atom_in);
  std::vector<int> growth_points;
  growth_points.reserve(32);
  candidate_additions.reserve(32);
  const int mol_idx = cdk.mol_home[seed_atom_in];
  if (cdk.mol_limits[mol_idx + 1] - cdk.mol_limits[mol_idx] < atom_limit) {
    for (int i = cdk.mol_limits[mol_idx]; i < cdk.mol_limits[mol_idx + 1]; i++) {
      addNewAtom(i);
    }
    candidate_additions.resize(0);
    const int fua_atom = vdel_pointer->getFirstUnassignedAtom();
    if (fua_atom < cdk.natom) {
      candidate_additions.push_back(fua_atom);
    }
  }

  // Proceed to add atoms layer by layer in the current molecule, or if that has been totally
  // covered, jump to the next molecule and continue the process.
  int ncandidate = candidate_additions.size();
  while (ncandidate > 0 && atom_count + ncandidate < atom_limit) {

    // Transfer candidate atoms to become growth points.  In construction, valence work units
    // do not overlap.  Atoms included in some previous work unit do not become candidates.
    growth_points.resize(0);
    for (int i = 0; i < ncandidate; i++) {
      addNewAtom(candidate_additions[i]);
      growth_points.push_back(candidate_additions[i]);
    }
    candidate_additions.resize(0);

    // Loop over the growth points and determine new candidate atoms.  During construction,
    // valence work units do not overlap.  Atoms included in some previous work unit do not
    // become new candidates.
    const int ngrow = growth_points.size();
    for (int i = 0; i < ngrow; i++) {
      const int grow_atom = growth_points[i];
      for (int j = nbk.nb12_bounds[grow_atom]; j < nbk.nb12_bounds[grow_atom + 1]; j++) {
        if (vdel_pointer->getAtomAssignmentCount(nbk.nb12x[j]) == 0) {
          candidate_additions.push_back(nbk.nb12x[j]);
        }
      }
    }
    ncandidate = candidate_additions.size();

    // If no candidate molecules have yet been found, try jumping to the first unassigned
    // atom.  In all likelihood, this will be on another molecule.  That will be the seed for the
    // next round of additions.
    if (ncandidate == 0) {
      const int fua_atom = vdel_pointer->getFirstUnassignedAtom();
      if (fua_atom < cdk.natom) {
        candidate_additions.push_back(fua_atom);
        ncandidate = 1;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getListIndex() const {
  return list_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMinAtomIndex() const {
  return min_atom_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMaxAtomIndex() const {
  return max_atom_index;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMaxAtoms() const {
  return atom_limit;
}

//-------------------------------------------------------------------------------------------------
ValenceDelegator* ValenceWorkUnit::getDelegatorPointer() {
  return vdel_pointer;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ValenceWorkUnit::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* ValenceWorkUnit::getRestraintApparatusPointer() const {
  return ra_pointer;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getImportedAtom(const int index) const {
  return atom_import_list[index];
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::setListIndex(const int list_index_in) {
  list_index = list_index_in;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::setAtomLimit(const int new_limit) {
  if (new_limit < atom_count) {
    rtErr("The atom limit cannot be set below the number of atoms currently in a work unit.  "
          "ValenceWorkUnit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + " contains " + std::to_string(atom_count) + " atoms and "
          "cannot reduce its limit to " + std::to_string(new_limit) + ".", "ValenceWorkUnit",
          "setAtomLimit");
  }
  atom_limit = new_limit;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAtom(const int atom_index) {

  // Add the new atom to the list of atom imports.
  atom_import_list.push_back(atom_index);
  atom_count += 1;

  // Check the appropriate boxes in the delegator.
  vdel_pointer->markAtomAddition(list_index, atom_index);  
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewBondTerm(const int bond_term_index) {
  bond_term_list.push_back(bond_term_index);
  vdel_pointer->markBondAddition(list_index, bond_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAngleTerm(const int angl_term_index) {
  angl_term_list.push_back(angl_term_index);
  vdel_pointer->markAngleAddition(list_index, angl_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewDihedralTerm(const int dihe_term_index) {
  dihe_term_list.push_back(dihe_term_index);
  vdel_pointer->markDihedralAddition(list_index, dihe_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewUreyBradleyTerm(const int ubrd_term_index) {
  ubrd_term_list.push_back(ubrd_term_index);
  vdel_pointer->markUreyBradleyAddition(list_index, ubrd_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewCharmmImproperTerm(const int cimp_term_index) {
  cimp_term_list.push_back(cimp_term_index);
  vdel_pointer->markCharmmImproperAddition(list_index, cimp_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewCmapTerm(const int cmap_term_index) {
  cmap_term_list.push_back(cmap_term_index);
  vdel_pointer->markCmapAddition(list_index, cmap_term_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewPositionalRestraint(const int posn_rstr_index) {
  rposn_term_list.push_back(posn_rstr_index);
  vdel_pointer->markPositionalRestraintAddition(list_index, posn_rstr_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewDistanceRestraint(const int dist_rstr_index) {
  rbond_term_list.push_back(dist_rstr_index);
  vdel_pointer->markDistanceRestraintAddition(list_index, dist_rstr_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAngleRestraint(const int angl_rstr_index) {
  rangl_term_list.push_back(angl_rstr_index);
  vdel_pointer->markAngleRestraintAddition(list_index, angl_rstr_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewDihedralRestraint(const int dihe_rstr_index) {
  rdihe_term_list.push_back(dihe_rstr_index);
  vdel_pointer->markDihedralRestraintAddition(list_index, dihe_rstr_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewConstraintGroup(const int cnst_group_index) {
  cnst_group_list.push_back(cnst_group_index);
  vdel_pointer->markConstraintGroupAddition(list_index, cnst_group_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewSettleGroup(const int sett_group_index) {
  sett_group_list.push_back(sett_group_index);
  vdel_pointer->markSettleGroupAddition(list_index, sett_group_index);
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewVirtualSite(const int vsite_index) {
  virtual_site_list.push_back(vsite_index);
  vdel_pointer->markVirtualSiteAddition(list_index, vsite_index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceWorkUnit::findForcePartners(const int atom_idx,
                                                    const ValenceKit<double> &vk,
                                                    const RestraintApparatusDpReader &rar,
                                                    const VirtualSiteKit<double> &vsk,
                                                    const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(64);

  // Form a new call stack from the input.  Check that the call stack has not grown too large,
  // as this might indicate an infinite loop and some nonsensical feature of the topology.
  std::vector<int> updated_caller_stack(caller_stack);
  updated_caller_stack.push_back(atom_idx);
  if (updated_caller_stack.size() > max_atom_search_stacks) {
    const std::string atom_list = writeAtomList(updated_caller_stack,
                                                ag_pointer->getChemicalDetailsKit());
    rtErr("Call stack for finding force partners is too deep.  This likely indicates that, "
          "somehow, a virtual site in the topology is including itself as one of its frame atoms "
          "or some similar paradox.  Calls include atoms " + atom_list + ".", "findForcePartners");
  }
  
  // Push the atom itself onto the list first.  Even if there are no interactions, having the
  // atom is, by definition, critical to evaluating its movement.
  result.push_back(atom_idx);

  // Puhs atoms relevant to all valence force field terms onto the list.
  const std::vector<int> relevant_bonds = vdel_pointer->getBondAffectors(atom_idx);
  const int nbond = relevant_bonds.size();
  for (int i = 0; i < nbond; i++) {
    const int bond_term_index = relevant_bonds[i];
    result.push_back(vk.bond_i_atoms[bond_term_index]);
    result.push_back(vk.bond_j_atoms[bond_term_index]);
  }
  const std::vector<int> relevant_angls = vdel_pointer->getAngleAffectors(atom_idx);
  const int nangl = relevant_angls.size();
  for (int i = 0; i < nangl; i++) {
    const int angl_term_index = relevant_angls[i];
    result.push_back(vk.angl_i_atoms[angl_term_index]);
    result.push_back(vk.angl_j_atoms[angl_term_index]);
    result.push_back(vk.angl_k_atoms[angl_term_index]);
  }
  const std::vector<int> relevant_dihes = vdel_pointer->getDihedralAffectors(atom_idx);
  const int ndihe = relevant_dihes.size();
  for (int i = 0; i < ndihe; i++) {
    const int dihe_term_index = relevant_dihes[i];
    result.push_back(vk.dihe_i_atoms[dihe_term_index]);
    result.push_back(vk.dihe_j_atoms[dihe_term_index]);
    result.push_back(vk.dihe_k_atoms[dihe_term_index]);
    result.push_back(vk.dihe_l_atoms[dihe_term_index]);
  }
  const std::vector<int> relevant_ubrds = vdel_pointer->getUreyBradleyAffectors(atom_idx);
  const int nubrd = relevant_ubrds.size();
  for (int i = 0; i < nubrd; i++) {
    const int ubrd_term_index = relevant_ubrds[i];
    result.push_back(vk.ubrd_i_atoms[ubrd_term_index]);
    result.push_back(vk.ubrd_k_atoms[ubrd_term_index]);
  }
  const std::vector<int> relevant_cimps = vdel_pointer->getCharmmImproperAffectors(atom_idx);
  const int ncimp = relevant_cimps.size();
  for (int i = 0; i < ncimp; i++) {
    const int cimp_term_index = relevant_cimps[i];
    result.push_back(vk.cimp_i_atoms[cimp_term_index]);
    result.push_back(vk.cimp_j_atoms[cimp_term_index]);
    result.push_back(vk.cimp_k_atoms[cimp_term_index]);
    result.push_back(vk.cimp_l_atoms[cimp_term_index]);
  }
  const std::vector<int> relevant_cmaps = vdel_pointer->getCmapAffectors(atom_idx);
  const int ncmap = relevant_cmaps.size();
  for (int i = 0; i < ncmap; i++) {
    const int cmap_term_index = relevant_cmaps[i];
    result.push_back(vk.cmap_i_atoms[cmap_term_index]);
    result.push_back(vk.cmap_j_atoms[cmap_term_index]);
    result.push_back(vk.cmap_k_atoms[cmap_term_index]);
    result.push_back(vk.cmap_l_atoms[cmap_term_index]);
    result.push_back(vk.cmap_m_atoms[cmap_term_index]);
  }

  // Most 1:4 attenuated interactions are implicitly covered by dihedral interactions to which
  // they can be linked.  Most virtual sites will have additional 1:4 attenuated interactions
  // that are not covered by dihedral interactions, as the virtual sites do not typically
  // participate in dihedral terms.  Push these onto the list.
  const std::vector<int> relevant_infrs = vdel_pointer->getInferred14Affectors(atom_idx);
  const int ninfr = relevant_infrs.size();
  for (int i = 0; i < ninfr; i++) {
    const int infr_term_index = relevant_infrs[i];
    result.push_back(vk.infr14_i_atoms[infr_term_index]);
    result.push_back(vk.infr14_l_atoms[infr_term_index]);
  }

  // Push atoms relevant to NMR restraints onto the list.
  const std::vector<int> relevant_rposn = vdel_pointer->getPositionalRestraintAffectors(atom_idx);
  const int nrposn = relevant_rposn.size();
  for (int i = 0; i < nrposn; i++) {
    const int rposn_term_index = relevant_rposn[i];
    result.push_back(rar.rposn_atoms[rposn_term_index]);
  }
  const std::vector<int> relevant_rbond = vdel_pointer->getDistanceRestraintAffectors(atom_idx);
  const int nrbond = relevant_rbond.size();
  for (int i = 0; i < nrbond; i++) {
    const int rbond_term_index = relevant_rbond[i];
    result.push_back(rar.rbond_i_atoms[rbond_term_index]);
    result.push_back(rar.rbond_j_atoms[rbond_term_index]);
  }
  const std::vector<int> relevant_rangl = vdel_pointer->getAngleRestraintAffectors(atom_idx);
  const int nrangl = relevant_rangl.size();
  for (int i = 0; i < nrangl; i++) {
    const int rangl_term_index = relevant_rangl[i];
    result.push_back(rar.rangl_i_atoms[rangl_term_index]);
    result.push_back(rar.rangl_j_atoms[rangl_term_index]);
    result.push_back(rar.rangl_k_atoms[rangl_term_index]);
  }
  const std::vector<int> relevant_rdihe = vdel_pointer->getDihedralRestraintAffectors(atom_idx);
  const int nrdihe = relevant_rdihe.size();
  for (int i = 0; i < nrdihe; i++) {
    const int rdihe_term_index = relevant_rdihe[i];
    result.push_back(rar.rdihe_i_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_j_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_k_atoms[rdihe_term_index]);
    result.push_back(rar.rdihe_l_atoms[rdihe_term_index]);
  }

  // Virtual sites transfer the non-bonded forces they have accumuated (including contributions
  // from attenuated 1:4 non-bonded interactions, possibly standard valence terms if the force
  // field is very abstract) to their frame atoms.  The virtual site, all of its frame atoms, and
  // any atoms making interactions with the virtual site must therefore be included.
  const std::vector<int> relevant_vstes = vdel_pointer->getVirtualSiteAffectors(atom_idx);
  const int nvsite = relevant_vstes.size();
  for (int i = 0; i < nvsite; i++) {
    const int vste_index = relevant_vstes[i];

    // If the virtual site is atom_idx, then the primary and inferred 1:4 interactions have already
    // been counted in the work above.  This part of the step can be skipped, and since the atom
    // itself has already been added to the list of relevant atoms, it is not necessary to add it
    // again.  If the virtual site is not atom_idx, this step will cover adding the virtual site
    // itself to the list of relevant atoms.
    const int vs_idx = vsk.vs_atoms[vste_index];
    if (vs_idx != atom_idx) {
      const std::vector<int> tmpv = findForcePartners(vs_idx, vk, rar, vsk, updated_caller_stack);
      result.insert(result.end(), tmpv.begin(), tmpv.end());
    }

    // Regardless of whether atoms relevant to additional interactions to the virtual site were
    // included, the frame atoms must be added to the list of relevant atoms so that it is known
    // how forces from the virtual site will be split up, some of which will land upon atom_idx
    // if atom_idx was not the virtual site itself.
    result.push_back(vsk.frame1_idx[vste_index]);
    result.push_back(vsk.frame2_idx[vste_index]);
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[vste_index])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::NONE:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      result.push_back(vsk.frame3_idx[vste_index]);
      break;
    case VirtualSiteKind::FIXED_4:
      result.push_back(vsk.frame3_idx[vste_index]);
      result.push_back(vsk.frame4_idx[vste_index]);
      break;
    }
  }
    
  // Sort the list of required atoms and prune duplicate entries
  reduceUniqueValues(&result);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceWorkUnit::findMovementPartners(const int atom_idx, const ConstraintKit<double> &cnk,
                                      const VirtualSiteKit<double> &vsk,
                                      const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(64);
  
  // Form a new call stack from the input.  Check that the call stack has not grown too large,
  // as this might indicate an infinite loop and some nonsensical feature of the topology.
  std::vector<int> updated_caller_stack(caller_stack);
  updated_caller_stack.push_back(atom_idx);
  if (updated_caller_stack.size() > max_atom_search_stacks) {
    const std::string atom_list = writeAtomList(updated_caller_stack,
                                                ag_pointer->getChemicalDetailsKit());
    rtErr("Call stack for finding movement partners is too deep.  This likely indicates that, "
          "somehow, a constraint group and virtual site in the topology are interacting in a "
          "way that creates a paradox.  Calls include atoms " + atom_list + ".",
          "findMovementPartners");
  }

  // The movement partners must necessarily include the atom itself, even if there are no
  // constraint or virtual sites affecting it.
  result.push_back(atom_idx);

  // If the atom is a virtual site, then all of its frame atoms are likewise required to move
  // before its new position can be determined.  Furthermore, any constraint groups which they
  // are a part of must be included.  Use a recursive call to sweep up the frame atoms'
  // constraint groups.
  if (ag_pointer->getAtomicNumber(atom_idx) == 0) {

    // Find the precise virtual site frame index to include only those atoms
    const std::vector<int> relevant_vstes = vdel_pointer->getVirtualSiteAffectors(atom_idx);
    const int nvsite = relevant_vstes.size();
    for (int i = 0; i < nvsite; i++) {
      const int vste_idx = relevant_vstes[i];
      if (vsk.vs_atoms[vste_idx] == atom_idx) {
        const std::vector<int> tmp1 = findMovementPartners(vsk.frame1_idx[vste_idx], cnk, vsk);
        result.insert(result.end(), tmp1.begin(), tmp1.end());
        const std::vector<int> tmp2 = findMovementPartners(vsk.frame2_idx[vste_idx], cnk, vsk);
        result.insert(result.end(), tmp2.begin(), tmp2.end());
        switch (static_cast<VirtualSiteKind>(vsk.vs_types[vste_idx])) {
        case VirtualSiteKind::FLEX_2:
        case VirtualSiteKind::FIXED_2:
        case VirtualSiteKind::NONE:
          break;
        case VirtualSiteKind::FLEX_3:
        case VirtualSiteKind::FIXED_3:
        case VirtualSiteKind::FAD_3:
        case VirtualSiteKind::OUT_3:
          {
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx], cnk, vsk);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
          }
          break;
        case VirtualSiteKind::FIXED_4:
          {
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx], cnk, vsk);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
            const std::vector<int> tmp4 = findMovementPartners(vsk.frame4_idx[vste_idx], cnk, vsk);
            result.insert(result.end(), tmp4.begin(), tmp4.end());
          }
          break;
        }
      }
    }
  }
  
  // Any constraint groups that affect the atom must have all of their atoms moved along
  // with the atom.
  const std::vector<int> relevant_setts = vdel_pointer->getSettleGroupAffectors(atom_idx);
  const int nsett = relevant_setts.size();
  for (int i = 0; i < nsett; i++) {
    const int sett_index = relevant_setts[i];
    result.push_back(cnk.settle_ox_atoms[sett_index]);
    result.push_back(cnk.settle_h1_atoms[sett_index]);
    result.push_back(cnk.settle_h2_atoms[sett_index]);
  }
  const std::vector<int> relevant_cnsts = vdel_pointer->getConstraintGroupAffectors(atom_idx);
  const int ncnst = relevant_cnsts.size();
  for (int i = 0; i < ncnst; i++) {
    const int cnst_index = relevant_cnsts[i];
    for (int j = cnk.group_bounds[cnst_index]; j < cnk.group_bounds[cnst_index + 1]; j++) {
      result.push_back(cnk.group_list[j]);
    }
  }

  // Sort the list of required atoms and prune duplicate entries
  reduceUniqueValues(&result);  
  return result;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::assignUpdateTasks(const ValenceKit<double> &vk,
                                        const RestraintApparatusDpReader &rar,
                                        const ConstraintKit<double> &cnk,
                                        const VirtualSiteKit<double> &vsk) {
  std::vector<int> required_atoms(64);
  std::vector<int> atom_requests;
  atom_requests.reserve(32);
  std::vector<int> atom_request_bounds(atom_count + 1, 0);
  atom_update_list.reserve(atom_count);
  for (int i = 0; i < atom_count; i++) {
    const int atomi = atom_import_list[i];

    // Cycle if this atom is already scheduled for update by some other work unit
    if (vdel_pointer->getUpdateWorkUnit(atomi) >= 0) {
      atom_request_bounds[i + 1] = atom_requests.size();
      continue;
    }

    // Accumulate a list of atoms to check for.  Scan valence terms, restraints, constraints, and
    // virtual sites.
    required_atoms.resize(0);

    // Accumulate a list of all atoms which must be moved in order to know the correct position
    // of this one.
    std::vector<int> move_req = findMovementPartners(atomi, cnk, vsk);
    const size_t nmove = move_req.size();
    for (size_t j = 0; j < nmove; j++) {
      const std::vector<int> tmpj = findForcePartners(move_req[j], vk, rar, vsk);
      required_atoms.insert(required_atoms.end(), tmpj.begin(), tmpj.end());
    }
    reduceUniqueValues(&required_atoms);

    // Check that each of the required atoms is present in this work unit.  Assemble a list of
    // required atoms if it is not.
    const size_t nreq = required_atoms.size();
    for (size_t j = 0; j < nreq; j++) {
      if (vdel_pointer->checkPresence(required_atoms[j], list_index) == false) {
        atom_requests.push_back(required_atoms[j]);
      }
    }

    // Record the request bounds for the ith atom.  The atom requests array will likely contain
    // duplicates, but these bounds permit an analysis of what is needed in order to make update
    // of the ith atom possible.
    atom_request_bounds[i + 1] = atom_requests.size();
  }

  // Get a list of the unique atom requests
  const std::vector<int> unique_requests = reduceUniqueValues(atom_requests);
  
  // If the number of unique requested atoms can fit within the import size limits of this work
  // unit, mark all atoms in the current list for movement and update and add the requests to the
  // list of imports.
  const int n_atoms_needed = unique_requests.size();
  if (atom_count + n_atoms_needed < atom_limit) {
    for (int i = 0; i < atom_count; i++) {
      if (vdel_pointer->setUpdateWorkUnit(atom_import_list[i], list_index)) {
        atom_update_list.push_back(atom_import_list[i]);
      }
    }
    for (int i = 0; i < n_atoms_needed; i++) {
      addNewAtom(unique_requests[i]);
    }
  }
  else {

    // Scan over all atoms and set this unit as the update work unit for those that have no atom
    // requests (this means that the present import list is sufficient to prepare all force
    // computations and other movements in order to move and update these atoms).
    for (int i = 0; i < atom_count; i++) {
      if (atom_request_bounds[i + 1] == atom_request_bounds[i] &&
          vdel_pointer->setUpdateWorkUnit(atom_import_list[i], list_index)) {
        atom_update_list.push_back(atom_import_list[i]);
      }
    }

    // Work with whatever remaining space to add as many atoms as possible to the list of
    // moving atoms.  Use the array of unique requests to determine which are the most important.
    int addition_threshold = 1;
    std::vector<bool> atom_request_fulfilled(atom_request_bounds[atom_count], false);
    const int all_request_count = atom_request_bounds[atom_count];
    std::vector<int2> unique_rc(n_atoms_needed);
    for (int i = 0; i < n_atoms_needed; i++) {
      unique_rc[i].x = 0;
      unique_rc[i].y = unique_requests[i];
    }
    for (int i = 0; i < all_request_count; i++) {
      unique_rc[locateValue(unique_requests, atom_requests[i], DataOrder::ASCENDING)].x += 1;
    }

    // CHECK
    printf("There are %2d atoms yet to import in work unit %2d:\n", n_atoms_needed, list_index);
    printf("   ");
    for (int i = 0; i < n_atoms_needed; i++) {
      printf(" %4d", unique_rc[i].y);
    }
    printf("\n   ");
    for (int i = 0; i < n_atoms_needed; i++) {
      printf(" %4d", unique_rc[i].x);
    }
    printf("\n");
    // END CHECK

  }
}
  
//-------------------------------------------------------------------------------------------------
std::vector<ValenceWorkUnit> buildValenceWorkUnits(const AtomGraph *ag,
                                                   const RestraintApparatus *ra,
                                                   const int max_atoms_per_vwu) {
  ValenceDelegator vdel(ag, ra);
  return buildValenceWorkUnits(&vdel, max_atoms_per_vwu);
}

//-------------------------------------------------------------------------------------------------
std::vector<ValenceWorkUnit> buildValenceWorkUnits(ValenceDelegator *vdel,
                                                   const int max_atoms_per_vwu) {
  std::vector<ValenceWorkUnit> result;

  // Spread atoms over all work units
  const AtomGraph *ag_ptr = vdel->getTopologyPointer();
  const RestraintApparatus *ra_ptr = vdel->getRestraintApparatusPointer();
  while (vdel->getFirstUnassignedAtom() < ag_ptr->getAtomCount()) {
    const int n_units = result.size();
    result.emplace_back(vdel, n_units, vdel->getFirstUnassignedAtom(), max_atoms_per_vwu);
  }
  
  // Shift atoms between work units in an effort to conserve bonding within any particular work
  // unit.

  // Assign one and only one unit to log the moves to each atom, including updating the positions
  // of constrained atoms and virtual sites.  By construction, the ValenceWorkUnit objects will
  // have a modicum of extra space to store additional atoms.  In order to determine the position
  // of an atom that is part of a constraint group, the work unit must possess and move all atoms
  // in the constraint group.  Likewise, in order to place a virtual site, the work unit must have
  // current locations of all frame atoms for the virtual site in their updated positions.  In
  // order to move any atom, all forces on that atom must be known.  The hierarchy is therefore:
  //
  // Atom to move and log
  //  - Atoms that contribute directly to forces on the atom to move
  //  - Other atoms in constraint group
  //     - Atoms that contribute to forces on the constraint group atoms
  //  - Frame atoms, if the atom is a virtual site
  //     - Atoms that contribute to forces on the frame atoms
  //
  // Atom A can contribute to forces on some other atom B by being part of a valence term, NMR
  // restraint, or a virtual site frame that also involves B.  The hierarchy cannot chain
  // indefinitely, as only one cycle of force computation (including transmission from a virtual
  // site) can occur before moving atoms.
  const ValenceKit<double> vk = ag_ptr->getDoublePrecisionValenceKit();
  const RestraintApparatusDpReader rar = ra_ptr->dpData();
  const ConstraintKit<double> cnk = ag_ptr->getDoublePrecisionConstraintKit();
  const VirtualSiteKit<double> vsk = ag_ptr->getDoublePrecisionVirtualSiteKit();
  const int vwu_count = result.size();
  for (int i = 0; i < vwu_count; i++) {
    const int natom = result[i].getAtomCount();
    result[i].assignUpdateTasks(vk, rar, cnk, vsk);
  }  
  
  return result;
}

} // namespace synthesis
} // namespace omni
