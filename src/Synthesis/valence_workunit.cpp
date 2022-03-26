#include "Math/summation.h"
#include "Topology/topology_util.h"
#include "valence_workunit.h"

namespace omni {
namespace synthesis {

using math::prefixSumInPlace;
using math::PrefixSumType;
using topology::VirtualSiteKind;
using topology::markAffectorAtoms;
using topology::extractBoundedListEntries;
  
//-------------------------------------------------------------------------------------------------
ValenceDelegator::ValenceDelegator(const AtomGraph *ag_in, const RestraintApparatus *ra_in) :
    atom_count{ag_in->getAtomCount()},
    first_unassigned_atom{0}, max_presence_allocation{2},
    bond_affector_list{}, bond_affector_bounds{}, angl_affector_list{}, angl_affector_bounds{},
    dihe_affector_list{}, dihe_affector_bounds{}, ubrd_affector_list{}, ubrd_affector_bounds{},
    cimp_affector_list{}, cimp_affector_bounds{}, cmap_affector_list{}, cmap_affector_bounds{},
    vste_affector_list{}, vste_affector_bounds{}, cnst_affector_list{}, cnst_affector_bounds{},
    sett_affector_list{}, sett_affector_bounds{}, work_unit_assignment_count{},
    work_unit_presence{}, ag_pointer{ag_in}, ra_pointer{ra_in}
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
  rposn_affector_list.resize(rar.nposn);
  rposn_affector_bounds.resize(atom_count + 1, 0);
  rbond_affector_list.resize(2 * rar.nbond);
  rbond_affector_bounds.resize(atom_count + 1, 0);
  rangl_affector_list.resize(3 * rar.nangl);
  rangl_affector_bounds.resize(atom_count + 1, 0);
  rdihe_affector_list.resize(4 * rar.ndihe);
  rdihe_affector_bounds.resize(atom_count + 1, 0);
  vste_affector_list.resize(4 * vsk.nsite);
  vste_affector_bounds.resize(atom_count + 1, 0);
  sett_affector_list.resize(3 * cnk.nsettle);
  sett_affector_bounds.resize(atom_count + 1, 0);
  cnst_affector_list.resize(cnk.group_bounds[cnk.ngroup], -1);
  cnst_affector_bounds.resize(atom_count + 1, 0);  
  work_unit_assignment_count.resize(atom_count, 0);
  work_unit_presence.resize(max_presence_allocation * atom_count);
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
    const int parent_atom = vsk.frame1_idx[pos];
    const int frame2_atom = vsk.frame2_idx[pos];
    const int list_p_idx  = vste_affector_bounds[parent_atom];
    const int list_f2_idx = vste_affector_bounds[frame2_atom];
    vste_affector_list[list_p_idx ] = pos;
    vste_affector_list[list_f2_idx] = pos;
    vste_affector_bounds[parent_atom] += 1;
    vste_affector_bounds[frame2_atom] += 1;
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
ValenceWorkUnit::ValenceWorkUnit(const AtomGraph *ag_in, const RestraintApparatus *ra_in,
                                 ValenceDelegator *vdel_in, const int list_index_in,
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
    ag_pointer{ag_in}, ra_pointer{ra_in}
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
void ValenceWorkUnit::assignUpdateTasks(const ValenceKit<double> &vk,
                                        const RestraintApparatusDpReader &rar,
                                        const ConstraintKit<double> &cnk,
                                        const VirtualSiteKit<double> &vsk) {
  std::vector<int> required_atoms(64);
  for (int i = 0; i < atom_count; i++) {

    // Accumulate a list of atoms to check for.  Scan valence terms, restraints, constraints, and
    // virtual sites.
    required_atoms.resize(0);
    const int atomi = atom_import_list[i];
    const std::vector<int> relevant_bonds = vdel_pointer->getBondAffectors(atomi);
    const int nbond = relevant_bonds.size();
    for (int j = 0; j < nbond; j++) {
      const int bond_term_index = relevant_bonds[j];
      required_atoms.push_back(vk.bond_i_atoms[bond_term_index]);
      required_atoms.push_back(vk.bond_j_atoms[bond_term_index]);
    }
    const std::vector<int> relevant_angls = vdel_pointer->getAngleAffectors(atomi);
    const int nangl = relevant_angls.size();
    for (int j = 0; j < nangl; j++) {
      const int angl_term_index = relevant_angls[j];
      required_atoms.push_back(vk.angl_i_atoms[angl_term_index]);
      required_atoms.push_back(vk.angl_j_atoms[angl_term_index]);
      required_atoms.push_back(vk.angl_k_atoms[angl_term_index]);
    }
    const std::vector<int> relevant_dihes = vdel_pointer->getDihedralAffectors(atomi);
    const int ndihe = relevant_dihes.size();
    for (int j = 0; j < ndihe; j++) {
      const int dihe_term_index = relevant_dihes[j];
      required_atoms.push_back(vk.dihe_i_atoms[dihe_term_index]);
      required_atoms.push_back(vk.dihe_j_atoms[dihe_term_index]);
      required_atoms.push_back(vk.dihe_k_atoms[dihe_term_index]);
      required_atoms.push_back(vk.dihe_l_atoms[dihe_term_index]);
    }
    const std::vector<int> relevant_ubrds = vdel_pointer->getUreyBradleyAffectors(atomi);
    const int nubrd = relevant_ubrds.size();
    for (int j = 0; j < nubrd; j++) {
      const int ubrd_term_index = relevant_ubrds[j];
      required_atoms.push_back(vk.ubrd_i_atoms[ubrd_term_index]);
      required_atoms.push_back(vk.ubrd_k_atoms[ubrd_term_index]);
    }
    const std::vector<int> relevant_cimps = vdel_pointer->getCharmmImproperAffectors(atomi);
    const int ncimp = relevant_cimps.size();
    for (int j = 0; j < ncimp; j++) {
      const int cimp_term_index = relevant_cimps[j];
      required_atoms.push_back(vk.cimp_i_atoms[cimp_term_index]);
      required_atoms.push_back(vk.cimp_j_atoms[cimp_term_index]);
      required_atoms.push_back(vk.cimp_k_atoms[cimp_term_index]);
      required_atoms.push_back(vk.cimp_l_atoms[cimp_term_index]);
    }
    const std::vector<int> relevant_cmaps = vdel_pointer->getCmapAffectors(atomi);
    const int ncmap = relevant_cmaps.size();
    for (int j = 0; j < ncmap; j++) {
      const int cmap_term_index = relevant_cmaps[j];
      required_atoms.push_back(vk.cmap_i_atoms[cmap_term_index]);
      required_atoms.push_back(vk.cmap_j_atoms[cmap_term_index]);
      required_atoms.push_back(vk.cmap_k_atoms[cmap_term_index]);
      required_atoms.push_back(vk.cmap_l_atoms[cmap_term_index]);
      required_atoms.push_back(vk.cmap_m_atoms[cmap_term_index]);
    }
    const std::vector<int> relevant_rposn = vdel_pointer->getPositionalRestraintAffectors(atomi);
    const int nrposn = relevant_rposn.size();
    for (int j = 0; j < nrposn; j++) {
      const int rposn_term_index = relevant_rposn[j];
      required_atoms.push_back(rar.rposn_atoms[rposn_term_index]);
    }
    const std::vector<int> relevant_rbond = vdel_pointer->getDistanceRestraintAffectors(atomi);
    const int nrbond = relevant_rbond.size();
    for (int j = 0; j < nrbond; j++) {
      const int rbond_term_index = relevant_rbond[j];
      required_atoms.push_back(rar.rbond_i_atoms[rbond_term_index]);
      required_atoms.push_back(rar.rbond_j_atoms[rbond_term_index]);
    }
    const std::vector<int> relevant_rangl = vdel_pointer->getAngleRestraintAffectors(atomi);
    const int nrangl = relevant_rangl.size();
    for (int j = 0; j < nrangl; j++) {
      const int rangl_term_index = relevant_rangl[j];
      required_atoms.push_back(rar.rangl_i_atoms[rangl_term_index]);
      required_atoms.push_back(rar.rangl_j_atoms[rangl_term_index]);
      required_atoms.push_back(rar.rangl_k_atoms[rangl_term_index]);
    }
    const std::vector<int> relevant_rdihe = vdel_pointer->getDihedralRestraintAffectors(atomi);
    const int nrdihe = relevant_rdihe.size();
    for (int j = 0; j < nrdihe; j++) {
      const int rdihe_term_index = relevant_rdihe[j];
      required_atoms.push_back(rar.rdihe_i_atoms[rdihe_term_index]);
      required_atoms.push_back(rar.rdihe_j_atoms[rdihe_term_index]);
      required_atoms.push_back(rar.rdihe_k_atoms[rdihe_term_index]);
      required_atoms.push_back(rar.rdihe_l_atoms[rdihe_term_index]);
    }
    const std::vector<int> relevant_setts = vdel_pointer->getSettleGroupAffectors(atomi);
    const int nsett = relevant_setts.size();
    for (int j = 0; j < nsett; j++) {
      const int sett_index = relevant_setts[j];
      required_atoms.push_back(cnk.settle_ox_atoms[sett_index]);
      required_atoms.push_back(cnk.settle_h1_atoms[sett_index]);
      required_atoms.push_back(cnk.settle_h2_atoms[sett_index]);
    }
    const std::vector<int> relevant_cnsts = vdel_pointer->getConstraintGroupAffectors(atomi);
    const int ncnst = relevant_cnsts.size();
    for (int j = 0; j < ncnst; j++) {
      const int cnst_index = relevant_cnsts[j];
      for (int k = cnk.group_bounds[cnst_index]; k < cnk.group_bounds[cnst_index + 1]; k++) {
        required_atoms.push_back(cnk.group_list[k]);
      }
    }
    const std::vector<int> relevant_vstes = vdel_pointer->getVirtualSiteAffectors(atomi);
    const int nvsite = relevant_vstes.size();
    for (int j = 0; j < nvsite; j++) {
      const int vste_index = relevant_vstes[j];
      required_atoms.push_back(vsk.frame1_idx[vste_index]);
      required_atoms.push_back(vsk.frame2_idx[vste_index]);
      switch (static_cast<VirtualSiteKind>(vsk.vs_types[vste_index])) {
      case VirtualSiteKind::FLEX_2:
      case VirtualSiteKind::FIXED_2:
      case VirtualSiteKind::NONE:
        break;
      case VirtualSiteKind::FLEX_3:
      case VirtualSiteKind::FIXED_3:
      case VirtualSiteKind::FAD_3:
      case VirtualSiteKind::OUT_3:
        required_atoms.push_back(vsk.frame3_idx[vste_index]);
        break;
      case VirtualSiteKind::FIXED_4:
        required_atoms.push_back(vsk.frame3_idx[vste_index]);
        required_atoms.push_back(vsk.frame4_idx[vste_index]);
        break;
      }
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
std::vector<ValenceWorkUnit> buildValenceWorkUnits(const AtomGraph *ag,
                                                   const RestraintApparatus *ra,
                                                   const int max_atoms_per_vwu) {
  ValenceDelegator vdel(ag, ra);
  std::vector<ValenceWorkUnit> result;
  
  // Spread atoms over all work units
  while (vdel.getFirstUnassignedAtom() < ag->getAtomCount()) {
    const int n_units = result.size();
    result.emplace_back(ag, ra, &vdel, n_units, vdel.getFirstUnassignedAtom(), max_atoms_per_vwu);
  }

  // CHECK
  for (size_t i = 0; i < result.size(); i++) {
    printf("Work unit %4zu has %3d atoms:\n", i, result[i].getAtomCount());
    int k = 0;
    for (int j = 0; j < result[i].getAtomCount(); j++) {
      printf(" %5d", result[i].getImportedAtom(j));
      k++;
      if (k == 16) {
        printf("\n");
        k = 0;
      }
    }
    if (k > 0) {
      printf("\n");
    }
  }
  // END CHECK
  
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
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const RestraintApparatusDpReader rar = ra->dpData();
  const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  const int vwu_count = result.size();
  for (int i = 0; i < vwu_count; i++) {
    const int natom = result[i].getAtomCount();
    result[i].assignUpdateTasks(vk, rar, cnk, vsk);
  }
  
  
  return result;
}

} // namespace synthesis
} // namespace omni
