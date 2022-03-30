#include <algorithm>
#include "Math/statistics.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/topology_util.h"
#include "valence_workunit.h"

namespace omni {
namespace synthesis {

using math::DataOrder;
using math::locateValue;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::reduceUniqueValues;
using math::roundUp;
using parse::char4ToString;
using topology::extractBoundedListEntries;
using topology::markAffectorAtoms;
using topology::TorsionKind;
using topology::VirtualSiteKind;
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
std::vector<int> ValenceDelegator::getBondAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(bond_affector_list, bond_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getAngleAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(angl_affector_list, angl_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDihedralAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(dihe_affector_list, dihe_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getUreyBradleyAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(ubrd_affector_list, ubrd_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getCharmmImproperAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cimp_affector_list, cimp_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceDelegator::getCmapAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cmap_affector_list, cmap_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getInferred14Affectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(infr_affector_list, infr_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getPositionalRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rposn_affector_list, rposn_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDistanceRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rbond_affector_list, rbond_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getAngleRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rangl_affector_list, rangl_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getDihedralRestraintAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(rdihe_affector_list, rdihe_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getVirtualSiteAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(vste_affector_list, vste_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getSettleGroupAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(sett_affector_list, sett_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::getConstraintGroupAffectors(const std::vector<int> &atom_indices) const {
  return extractBoundedListEntries(cnst_affector_list, cnst_affector_bounds, atom_indices);
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ValenceDelegator::findMovementPartners(const int atom_idx,
                                       const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(64);
  const ConstraintKit<double> cnk  = ag_pointer->getDoublePrecisionConstraintKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  
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
    const std::vector<int> relevant_vstes = extractBoundedListEntries(vste_affector_list,
                                                                      vste_affector_bounds,
                                                                      atom_idx);
    const int nvsite = relevant_vstes.size();
    for (int i = 0; i < nvsite; i++) {
      const int vste_idx = relevant_vstes[i];
      if (vsk.vs_atoms[vste_idx] == atom_idx) {
        const std::vector<int> tmp1 = findMovementPartners(vsk.frame1_idx[vste_idx],
                                                           updated_caller_stack);
        result.insert(result.end(), tmp1.begin(), tmp1.end());
        const std::vector<int> tmp2 = findMovementPartners(vsk.frame2_idx[vste_idx],
                                                           updated_caller_stack);
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
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
          }
          break;
        case VirtualSiteKind::FIXED_4:
          {
            const std::vector<int> tmp3 = findMovementPartners(vsk.frame3_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp3.begin(), tmp3.end());
            const std::vector<int> tmp4 = findMovementPartners(vsk.frame4_idx[vste_idx],
                                                               updated_caller_stack);
            result.insert(result.end(), tmp4.begin(), tmp4.end());
          }
          break;
        }
      }
    }
  }
  
  // Any constraint groups that affect the atom must have all of their atoms moved along
  // with the atom.
  const std::vector<int> relevant_setts = extractBoundedListEntries(sett_affector_list,
                                                                    sett_affector_bounds,
                                                                    atom_idx);
  const int nsett = relevant_setts.size();
  for (int i = 0; i < nsett; i++) {
    const int sett_index = relevant_setts[i];
    result.push_back(cnk.settle_ox_atoms[sett_index]);
    result.push_back(cnk.settle_h1_atoms[sett_index]);
    result.push_back(cnk.settle_h2_atoms[sett_index]);
  }
  const std::vector<int> relevant_cnsts = extractBoundedListEntries(cnst_affector_list,
                                                                    cnst_affector_bounds,
                                                                    atom_idx);
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
std::vector<int> ValenceDelegator::findForcePartners(const int atom_idx,
                                                     const std::vector<int> &caller_stack) const {
  std::vector<int> result;
  result.reserve(64);
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const RestraintApparatusDpReader rar = ra_pointer->dpData();
  
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

  // Push atoms relevant to all valence force field terms onto the list.
  const std::vector<int> relevant_bonds = extractBoundedListEntries(bond_affector_list,
                                                                    bond_affector_bounds,
                                                                    atom_idx);
  const int nbond = relevant_bonds.size();
  for (int i = 0; i < nbond; i++) {
    const int bond_term_index = relevant_bonds[i];
    result.push_back(vk.bond_i_atoms[bond_term_index]);
    result.push_back(vk.bond_j_atoms[bond_term_index]);
  }
  const std::vector<int> relevant_angls = extractBoundedListEntries(angl_affector_list,
                                                                    angl_affector_bounds,
                                                                    atom_idx);
  const int nangl = relevant_angls.size();
  for (int i = 0; i < nangl; i++) {
    const int angl_term_index = relevant_angls[i];
    result.push_back(vk.angl_i_atoms[angl_term_index]);
    result.push_back(vk.angl_j_atoms[angl_term_index]);
    result.push_back(vk.angl_k_atoms[angl_term_index]);
  }
  const std::vector<int> relevant_dihes = extractBoundedListEntries(dihe_affector_list,
                                                                    dihe_affector_bounds,
                                                                    atom_idx);
  const int ndihe = relevant_dihes.size();
  for (int i = 0; i < ndihe; i++) {
    const int dihe_term_index = relevant_dihes[i];
    result.push_back(vk.dihe_i_atoms[dihe_term_index]);
    result.push_back(vk.dihe_j_atoms[dihe_term_index]);
    result.push_back(vk.dihe_k_atoms[dihe_term_index]);
    result.push_back(vk.dihe_l_atoms[dihe_term_index]);
  }
  const std::vector<int> relevant_ubrds = extractBoundedListEntries(ubrd_affector_list,
                                                                    ubrd_affector_bounds,
                                                                    atom_idx);
  const int nubrd = relevant_ubrds.size();
  for (int i = 0; i < nubrd; i++) {
    const int ubrd_term_index = relevant_ubrds[i];
    result.push_back(vk.ubrd_i_atoms[ubrd_term_index]);
    result.push_back(vk.ubrd_k_atoms[ubrd_term_index]);
  }
  const std::vector<int> relevant_cimps = extractBoundedListEntries(cimp_affector_list,
                                                                    cimp_affector_bounds,
                                                                    atom_idx);
  const int ncimp = relevant_cimps.size();
  for (int i = 0; i < ncimp; i++) {
    const int cimp_term_index = relevant_cimps[i];
    result.push_back(vk.cimp_i_atoms[cimp_term_index]);
    result.push_back(vk.cimp_j_atoms[cimp_term_index]);
    result.push_back(vk.cimp_k_atoms[cimp_term_index]);
    result.push_back(vk.cimp_l_atoms[cimp_term_index]);
  }
  const std::vector<int> relevant_cmaps = extractBoundedListEntries(cmap_affector_list,
                                                                    cmap_affector_bounds,
                                                                    atom_idx);
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
  const std::vector<int> relevant_infrs = extractBoundedListEntries(infr_affector_list,
                                                                    infr_affector_bounds,
                                                                    atom_idx);
  const int ninfr = relevant_infrs.size();
  for (int i = 0; i < ninfr; i++) {
    const int infr_term_index = relevant_infrs[i];
    result.push_back(vk.infr14_i_atoms[infr_term_index]);
    result.push_back(vk.infr14_l_atoms[infr_term_index]);
  }

  // Push atoms relevant to NMR restraints onto the list.
  const std::vector<int> relevant_rposn = extractBoundedListEntries(rposn_affector_list,
                                                                    rposn_affector_bounds,
                                                                    atom_idx);
  const int nrposn = relevant_rposn.size();
  for (int i = 0; i < nrposn; i++) {
    const int rposn_term_index = relevant_rposn[i];
    result.push_back(rar.rposn_atoms[rposn_term_index]);
  }
  const std::vector<int> relevant_rbond = extractBoundedListEntries(rbond_affector_list,
                                                                    rbond_affector_bounds,
                                                                    atom_idx);
  const int nrbond = relevant_rbond.size();
  for (int i = 0; i < nrbond; i++) {
    const int rbond_term_index = relevant_rbond[i];
    result.push_back(rar.rbond_i_atoms[rbond_term_index]);
    result.push_back(rar.rbond_j_atoms[rbond_term_index]);
  }
  const std::vector<int> relevant_rangl = extractBoundedListEntries(rangl_affector_list,
                                                                    rangl_affector_bounds,
                                                                    atom_idx);
  const int nrangl = relevant_rangl.size();
  for (int i = 0; i < nrangl; i++) {
    const int rangl_term_index = relevant_rangl[i];
    result.push_back(rar.rangl_i_atoms[rangl_term_index]);
    result.push_back(rar.rangl_j_atoms[rangl_term_index]);
    result.push_back(rar.rangl_k_atoms[rangl_term_index]);
  }
  const std::vector<int> relevant_rdihe = extractBoundedListEntries(rdihe_affector_list,
                                                                    rdihe_affector_bounds,
                                                                    atom_idx);
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
  const std::vector<int> relevant_vstes = extractBoundedListEntries(vste_affector_list,
                                                                    vste_affector_bounds,
                                                                    atom_idx);
  const int nvsite = relevant_vstes.size();
  for (int i = 0; i < nvsite; i++) {
    const int vste_index = relevant_vstes[i];

    // If the virtual site is atom_idx, then the primary and inferred 1:4 interactions have already
    // been counted in the work above.  This part of the step can be skipped, and since the atom
    // itself has already been added to the list of relevant atoms, it is not necessary to add it
    // again.  If the virtual site is not atom_idx, this step will cover adding the virtual site
    // itself, and all of the particles that participate in valence or restraint forces with it, to
    // the list of relevant atoms.
    const int vs_idx = vsk.vs_atoms[vste_index];
    if (vs_idx != atom_idx) {
      const std::vector<int> tmpv = findForcePartners(vs_idx, updated_caller_stack);
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
std::vector<int> ValenceDelegator::getUpdateDependencies(const int atom_index) const {
  std::vector<int> result;
  std::vector<int> move_req = findMovementPartners(atom_index);
  const size_t nmove = move_req.size();
  for (size_t j = 0; j < nmove; j++) {

    // Each of the atoms in move_req will be included in the tmpj vector, so no need to
    // explicitly add the contents of move_req to the list.
    const std::vector<int> tmpj = findForcePartners(move_req[j]);
    result.insert(result.end(), tmpj.begin(), tmpj.end());
  }
  reduceUniqueValues(&result);
  return result;
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
    for (int i = atom_count - 1; i >= 0; i--) {
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
}

//-------------------------------------------------------------------------------------------------
bool ValenceDelegator::setUpdateWorkUnit(const int atom_index, const int vwu_index) {
  if (assigned_update_work_units[atom_index] >= 0) {
    return false;
  }
  assigned_update_work_units[atom_index] = vwu_index;

  // Update the first unassigned atom, if appropriate
  while (first_unassigned_atom < atom_count &&
         work_unit_assignment_count[first_unassigned_atom] > 0) {
    first_unassigned_atom += 1;
  }
  return true;
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
    imported_atom_count{0}, moved_atom_count{0}, updated_atom_count{0}, bond_term_count{0},
    angl_term_count{0}, dihe_term_count{0}, ubrd_term_count{0}, cimp_term_count{0},
    cdhe_term_count{0}, cmap_term_count{0}, infr14_term_count{0}, rposn_term_count{0},
    rbond_term_count{0}, rangl_term_count{0}, rdihe_term_count{0}, cnst_group_count{0},
    sett_group_count{0}, vste_count{0}, list_index{list_index_in}, min_atom_index{-1},
    max_atom_index{-1}, atom_limit{max_atoms_in}, atom_import_list{}, bond_term_list{},
    angl_term_list{}, dihe_term_list{}, ubrd_term_list{}, cimp_term_list{}, cmap_term_list{},
    infr14_term_list{}, bond_i_atoms{}, bond_j_atoms{}, angl_i_atoms{}, angl_j_atoms{},
    angl_k_atoms{}, dihe_i_atoms{}, dihe_j_atoms{}, dihe_k_atoms{}, dihe_l_atoms{}, ubrd_i_atoms{},
    ubrd_k_atoms{}, cimp_i_atoms{}, cimp_j_atoms{}, cimp_k_atoms{}, cimp_l_atoms{}, cmap_i_atoms{},
    cmap_j_atoms{}, cmap_k_atoms{}, cmap_l_atoms{}, cmap_m_atoms{}, infr14_i_atoms{},
    infr14_l_atoms{}, cdhe_term_list{}, cdhe_is_cimp{}, cdhe_i_atoms{}, cdhe_j_atoms{},
    cdhe_k_atoms{}, cdhe_l_atoms{}, rposn_term_list{}, rbond_term_list{}, rangl_term_list{},
    rdihe_term_list{}, rposn_atoms{}, rbond_i_atoms{}, rbond_j_atoms{}, rangl_i_atoms{},
    rangl_j_atoms{}, rangl_k_atoms{}, rdihe_i_atoms{}, rdihe_j_atoms{}, rdihe_k_atoms{},
    rdihe_l_atoms{}, cnst_group_list{}, sett_group_list{}, cnst_group_atoms{}, cnst_group_bounds{},
    sett_ox_atoms{}, sett_h1_atoms{}, sett_h2_atoms{}, virtual_site_list{}, vsite_atoms{},
    vsite_frame1_atoms{}, vsite_frame2_atoms{}, vsite_frame3_atoms{}, vsite_frame4_atoms{},
    vdel_pointer{vdel_in}, ag_pointer{vdel_in->getTopologyPointer()},
    ra_pointer{vdel_in->getRestraintApparatusPointer()}
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
  int fua_atom = vdel_pointer->getFirstUnassignedAtom();
  if (fua_atom != seed_atom_in && fua_atom < cdk.natom) {
    candidate_additions.push_back(fua_atom);
  }
  
  // Proceed to add atoms layer by layer in the current molecule, or if that has been totally
  // covered, jump to the next molecule and continue the process.
  int ncandidate = candidate_additions.size();
  bool capacity_reached = false;
  while (ncandidate > 0 && capacity_reached == false) {

    // A single successful addition is needed to prove that the capacity has not yet been reached.
    capacity_reached = true;
    
    // Transfer candidate atoms to become growth points.  In construction, valence work units
    // do not overlap.  Atoms included in some previous work unit do not become candidates.
    growth_points.resize(0);
    for (int i = 0; i < ncandidate; i++) {

      // This is an atom that the work unit will be responsible for updating.  Check the number
      // of its dependencies and add them to the import list, if there is room.  Mark the atom
      // as one for the work unit to update.
      const std::vector<int> up_deps = vdel_pointer->getUpdateDependencies(candidate_additions[i]);
      const int ndeps = up_deps.size();
      std::vector<bool> incl_up_deps(ndeps, false);
      int n_new_atoms = 0;
      for (int j = 0; j < ndeps; j++) {
        if (vdel_pointer->checkPresence(up_deps[j], list_index) == false) {
          n_new_atoms++;
          incl_up_deps[j] = true;
        }
      }
      if (imported_atom_count + n_new_atoms < atom_limit) {

        // The candidate atom will be part of its own dependencies list.
        for (int j = 0; j < ndeps; j++) {
          if (incl_up_deps[j]) {
            addNewAtomImport(up_deps[j]);
          }
        }
        growth_points.push_back(candidate_additions[i]);
        addNewAtomUpdate(candidate_additions[i]);
        capacity_reached = false;
      }
    }
    candidate_additions.resize(0);
    
    // Loop over the growth points and determine new candidate atoms.  During construction,
    // valence work units do not overlap.  Atoms included in some previous work unit do not
    // become new candidates.
    const int ngrow = growth_points.size();
    for (int i = 0; i < ngrow; i++) {
      const int grow_atom = growth_points[i];
      for (int j = nbk.nb12_bounds[grow_atom]; j < nbk.nb12_bounds[grow_atom + 1]; j++) {
        if (vdel_pointer->getUpdateWorkUnit(nbk.nb12x[j]) == -1) {
          candidate_additions.push_back(nbk.nb12x[j]);
        }
      }
    }
    reduceUniqueValues(&candidate_additions);
    ncandidate = candidate_additions.size();

    // If no candidate molecules have yet been found, try jumping to the first unassigned
    // atom.  In all likelihood, this will be on another molecule.  That will be the seed for the
    // next round of additions.
    if (ncandidate == 0) {
      fua_atom = vdel_pointer->getFirstUnassignedAtom();
      if (fua_atom < cdk.natom) {
        candidate_additions.push_back(fua_atom);
        ncandidate = 1;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getImportedAtomCount() const {
  return imported_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getMovedAtomCount() const {
  return moved_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ValenceWorkUnit::getUpdatedAtomCount() const {
  return updated_atom_count;
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
  if (new_limit < imported_atom_count) {
    rtErr("The atom limit cannot be set below the number of atoms currently in a work unit.  "
          "ValenceWorkUnit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + " contains " + std::to_string(imported_atom_count) +
          " atoms and cannot reduce its limit to " + std::to_string(new_limit) + ".",
          "ValenceWorkUnit", "setAtomLimit");
  }
  atom_limit = new_limit;
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAtomImport(const int atom_index) {

  // Check that the atom is not already in the work unit
  if (vdel_pointer->checkPresence(atom_index, list_index)) {
    return;
  }
  
  // Add the new atom to the list of atom imports.
  atom_import_list.push_back(atom_index);
  imported_atom_count += 1;

  // Check the appropriate boxes in the delegator.
  vdel_pointer->markAtomAddition(list_index, atom_index);  
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::addNewAtomUpdate(const int atom_index) {

  // Check that the atom is on the import list
  if (vdel_pointer->checkPresence(atom_index, list_index) == false) {
    rtErr("Atom with topological index " + std::to_string(atom_index) + " is not imported by "
          "ValenceWorkUnit " + std::to_string(list_index) + ".", "ValenceWorkUnit",
          "addNewAtomUpdate");
  }
  
  // Check the appropriate boxes in the delegator.
  if (vdel_pointer->setUpdateWorkUnit(atom_index, list_index)) {
    atom_update_list.push_back(atom_index);
    updated_atom_count += 1;
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::makeAtomMoveList() {
  std::vector<int> required_atoms;
  for (int i = 0; i < updated_atom_count; i++) {

    // The list of movement partners comprises the particle itself.
    const std::vector<int> tmpv = vdel_pointer->findMovementPartners(atom_update_list[i]);
    required_atoms.insert(required_atoms.end(), tmpv.begin(), tmpv.end());
  }
  atom_move_list = reduceUniqueValues(required_atoms);
  moved_atom_count = atom_move_list.size();
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::sortAtomSets() {
  std::sort(atom_import_list.begin(), atom_import_list.end(), [](int a, int b) { return a < b; });
  std::sort(atom_move_list.begin(), atom_move_list.end(), [](int a, int b) { return a < b; });
  std::sort(atom_update_list.begin(), atom_update_list.end(), [](int a, int b) { return a < b; });

  // Check the sizing of each list
  if (static_cast<int>(atom_import_list.size()) != imported_atom_count) {
    rtErr("Atom import count (" + std::to_string(imported_atom_count) + ") does not reflect the "
          "true length of the imported atom list (" + std::to_string(atom_import_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");
  }
  if (static_cast<int>(atom_move_list.size()) != moved_atom_count) {
    rtErr("Moving atom count (" + std::to_string(moved_atom_count) + ") does not reflect the "
          "true length of the atom move list (" + std::to_string(atom_move_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");    
  }
  if (static_cast<int>(atom_update_list.size()) != updated_atom_count) {
    rtErr("Atom update count (" + std::to_string(updated_atom_count) + ") does not reflect the "
          "true length of the atom move list (" + std::to_string(atom_update_list.size()) +
          ") in work unit " + std::to_string(list_index) + " serving topology " +
          ag_pointer->getFileName() + ".", "ValenceWorkUnit", "sortAtomSets");    
  }
  
  // Check the atom import list to ensure that all entries are unique.
  for (int i = 1; i < imported_atom_count; i++) {
    if (atom_import_list[i] == atom_import_list[i - 1]) {
      const int ires = ag_pointer->getResidueIndex(atom_import_list[i]);
      rtErr("A duplicate entry is present in the atom import list of work unit " +
            std::to_string(list_index) + " serving topology " + ag_pointer->getFileName() +
            ".  The topological atom index is " + std::to_string(atom_import_list[i]) + ", name " +
            char4ToString(ag_pointer->getAtomName(atom_import_list[i])) + ", residue " +
            char4ToString(ag_pointer->getResidueName(ires)) + " " +
            std::to_string(ag_pointer->getResidueNumber(atom_import_list[i])) + ".",
            "ValenceWorkUnit", "sortAtomSets");
    }
  }

  // Check the movement list to ensure that each atom is present in the import list.
  for (int i = 0; i < moved_atom_count; i++) {
    if (locateValue(atom_import_list, atom_move_list[i],
                    DataOrder::ASCENDING) == imported_atom_count) {
      rtErr("Atom index " + std::to_string(atom_move_list[i]) + "(" +
            char4ToString(ag_pointer->getAtomName(atom_move_list[i])) + ") in topology " +
            ag_pointer->getFileName() + " is scheduled to be moved by work unit " +
            std::to_string(list_index) + " but not present in the import list of " +
            std::to_string(imported_atom_count) + " atoms.", "ValenceWorkUnit", "sortAtomSets");
    }
  }
  for (int i = 0; i < updated_atom_count; i++) {
    if (locateValue(atom_move_list, atom_update_list[i],
                    DataOrder::ASCENDING) == moved_atom_count) {
      rtErr("Atom index " + std::to_string(atom_update_list[i]) + "(" +
            char4ToString(ag_pointer->getAtomName(atom_update_list[i])) + ") in topology " +
            ag_pointer->getFileName() + " is scheduled to be updated and logged by work unit " +
            std::to_string(list_index) + " but not present in the movement list of " +
            std::to_string(moved_atom_count) + " atoms.", "ValenceWorkUnit", "sortAtomSets");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ValenceWorkUnit::logActivities() {

  // Extract information from the topology
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag_pointer->getDoublePrecisionVirtualSiteKit();
  const ConstraintKit<double> cnk = ag_pointer->getDoublePrecisionConstraintKit();
  const RestraintApparatusDpReader rar = ra_pointer->dpData();

  // Make a "straight table" based on the minimum and maximum imported atom indices.  Most work
  // units will involve a confined sequence of atoms from within the topology.  While it might
  // be inefficient to have each work unit allocate a list mapping the entire topology to one of
  // its internal atoms (a bounded memory requirement in that ValenceWorkUnits are constructed
  // serially, but memory allocation effort still growing as the number of work units times the
  // number of atoms), it is feasible to take an offset and then allocate an array to map the
  // relevant part of the topology to the list of atoms that each work unit does hold.
  const int mapping_offset = atom_import_list[0];
  const int map_length = atom_import_list[imported_atom_count - 1] + 1 - mapping_offset;
  std::vector<int> import_map(map_length, -1);
  for (int i = 0; i < imported_atom_count; i++) {
    import_map[atom_import_list[i] - mapping_offset] = i;
  }
  
  // Detail valence interactions for this work unit
  const DataOrder order = DataOrder::ASCENDING;
  const std::vector<int> relevant_bonds = vdel_pointer->getBondAffectors(atom_move_list);
  const size_t nbond = relevant_bonds.size();
  for (size_t pos = 0; pos < nbond; pos++) {
    const int bond_idx = relevant_bonds[pos];
    bond_term_list.push_back(bond_idx);
    bond_i_atoms.push_back(import_map[vk.bond_i_atoms[bond_idx] - mapping_offset]);
    bond_j_atoms.push_back(import_map[vk.bond_j_atoms[bond_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_angls = vdel_pointer->getAngleAffectors(atom_move_list);
  const size_t nangl = relevant_angls.size();
  for (size_t pos = 0; pos < nangl; pos++) {
    const int angl_idx = relevant_angls[pos];
    angl_term_list.push_back(angl_idx);
    angl_i_atoms.push_back(import_map[vk.angl_i_atoms[angl_idx] - mapping_offset]);
    angl_j_atoms.push_back(import_map[vk.angl_j_atoms[angl_idx] - mapping_offset]);
    angl_k_atoms.push_back(import_map[vk.angl_k_atoms[angl_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_dihes = vdel_pointer->getDihedralAffectors(atom_move_list);
  const size_t ndihe = relevant_dihes.size();
  for (size_t pos = 0; pos < ndihe; pos++) {
    const int dihe_idx = relevant_dihes[pos];
    dihe_term_list.push_back(dihe_idx);
    dihe_i_atoms.push_back(import_map[vk.dihe_i_atoms[dihe_idx] - mapping_offset]);
    dihe_j_atoms.push_back(import_map[vk.dihe_j_atoms[dihe_idx] - mapping_offset]);
    dihe_k_atoms.push_back(import_map[vk.dihe_k_atoms[dihe_idx] - mapping_offset]);
    dihe_l_atoms.push_back(import_map[vk.dihe_l_atoms[dihe_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_ubrds = vdel_pointer->getUreyBradleyAffectors(atom_move_list);
  const size_t nubrd = relevant_ubrds.size();
  for (size_t pos = 0; pos < nubrd; pos++) {
    const int ubrd_idx = relevant_ubrds[pos];
    ubrd_term_list.push_back(ubrd_idx);
    ubrd_i_atoms.push_back(import_map[vk.ubrd_i_atoms[ubrd_idx] - mapping_offset]);
    ubrd_k_atoms.push_back(import_map[vk.ubrd_k_atoms[ubrd_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_cimps = vdel_pointer->getCharmmImproperAffectors(atom_move_list);
  const size_t ncimp = relevant_cimps.size();
  for (size_t pos = 0; pos < ncimp; pos++) {
    const int cimp_idx = relevant_cimps[pos];
    cimp_term_list.push_back(cimp_idx);
    cimp_i_atoms.push_back(import_map[vk.cimp_i_atoms[cimp_idx] - mapping_offset]);
    cimp_j_atoms.push_back(import_map[vk.cimp_j_atoms[cimp_idx] - mapping_offset]);
    cimp_k_atoms.push_back(import_map[vk.cimp_k_atoms[cimp_idx] - mapping_offset]);
    cimp_l_atoms.push_back(import_map[vk.cimp_l_atoms[cimp_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_cmaps = vdel_pointer->getCmapAffectors(atom_move_list);
  const size_t ncmap = relevant_cmaps.size();
  for (size_t pos = 0; pos < ncmap; pos++) {
    const int cmap_idx = relevant_cmaps[pos];
    cmap_term_list.push_back(cmap_idx);
    cmap_i_atoms.push_back(import_map[vk.cmap_i_atoms[cmap_idx] - mapping_offset]);
    cmap_j_atoms.push_back(import_map[vk.cmap_j_atoms[cmap_idx] - mapping_offset]);
    cmap_k_atoms.push_back(import_map[vk.cmap_k_atoms[cmap_idx] - mapping_offset]);
    cmap_l_atoms.push_back(import_map[vk.cmap_l_atoms[cmap_idx] - mapping_offset]);
    cmap_m_atoms.push_back(import_map[vk.cmap_m_atoms[cmap_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_infrs = vdel_pointer->getInferred14Affectors(atom_move_list);
  const size_t ninfr = relevant_infrs.size();
  for (size_t pos = 0; pos < ninfr; pos++) {
    const int infr_idx = relevant_infrs[pos];
    infr14_term_list.push_back(infr_idx);
    infr14_i_atoms.push_back(import_map[vk.infr14_i_atoms[infr_idx] - mapping_offset]);
    infr14_l_atoms.push_back(import_map[vk.infr14_l_atoms[infr_idx] - mapping_offset]);
  }

  // Detail the restraint tersm for this work unit to manage
  const std::vector<int> relevant_rposns =
    vdel_pointer->getPositionalRestraintAffectors(atom_move_list);
  const size_t nrposn = relevant_rposns.size();
  for (size_t pos = 0; pos < nrposn; pos++) {
    const int rposn_idx = relevant_rposns[pos];
    rposn_term_list.push_back(rposn_idx);
    rposn_atoms.push_back(import_map[rar.rposn_atoms[rposn_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_rbonds =
    vdel_pointer->getDistanceRestraintAffectors(atom_move_list);
  const size_t nrbond = relevant_rbonds.size();
  for (size_t pos = 0; pos < nrbond; pos++) {
    const int rbond_idx = relevant_rbonds[pos];
    rbond_term_list.push_back(rbond_idx);
    rbond_i_atoms.push_back(import_map[rar.rbond_i_atoms[rbond_idx] - mapping_offset]);
    rbond_j_atoms.push_back(import_map[rar.rbond_j_atoms[rbond_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_rangls =
    vdel_pointer->getAngleRestraintAffectors(atom_move_list);
  const size_t nrangl = relevant_rangls.size();
  for (size_t pos = 0; pos < nrangl; pos++) {
    const int rangl_idx = relevant_rangls[pos];
    rangl_term_list.push_back(rangl_idx);
    rangl_i_atoms.push_back(import_map[rar.rangl_i_atoms[rangl_idx] - mapping_offset]);
    rangl_j_atoms.push_back(import_map[rar.rangl_j_atoms[rangl_idx] - mapping_offset]);
    rangl_k_atoms.push_back(import_map[rar.rangl_k_atoms[rangl_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_rdihes =
    vdel_pointer->getDihedralRestraintAffectors(atom_move_list);
  const size_t nrdihe = relevant_rdihes.size();
  for (size_t pos = 0; pos < nrdihe; pos++) {
    const int rdihe_idx = relevant_rdihes[pos];
    rdihe_term_list.push_back(rdihe_idx);
    rdihe_i_atoms.push_back(import_map[rar.rdihe_i_atoms[rdihe_idx] - mapping_offset]);
    rdihe_j_atoms.push_back(import_map[rar.rdihe_j_atoms[rdihe_idx] - mapping_offset]);
    rdihe_k_atoms.push_back(import_map[rar.rdihe_k_atoms[rdihe_idx] - mapping_offset]);
    rdihe_l_atoms.push_back(import_map[rar.rdihe_l_atoms[rdihe_idx] - mapping_offset]);
  }

  // Detail virtual sites for thie work unit to manage
  const std::vector<int> relevant_vstes = vdel_pointer->getVirtualSiteAffectors(atom_move_list);
  const size_t nvste = relevant_vstes.size();
  for (size_t pos = 0; pos < nvste; pos++) {
    const int vste_idx = relevant_vstes[pos];
    virtual_site_list.push_back(vste_idx);
    vsite_atoms.push_back(import_map[vsk.vs_atoms[vste_idx] - mapping_offset]);
    vsite_frame1_atoms.push_back(import_map[vsk.frame1_idx[vste_idx] - mapping_offset]);
    vsite_frame2_atoms.push_back(import_map[vsk.frame2_idx[vste_idx] - mapping_offset]);
    vsite_frame3_atoms.push_back(import_map[vsk.frame3_idx[vste_idx] - mapping_offset]);
    vsite_frame4_atoms.push_back(import_map[vsk.frame4_idx[vste_idx] - mapping_offset]);
  }

  // Detail the constraint groups for this work unit
  const std::vector<int> relevant_setts = vdel_pointer->getSettleGroupAffectors(atom_move_list);
  const size_t nsett = relevant_setts.size();
  for (size_t pos = 0; pos < nsett; pos++) {
    const int sett_idx = relevant_setts[pos];
    sett_group_list.push_back(sett_idx);
    sett_ox_atoms.push_back(import_map[cnk.settle_ox_atoms[sett_idx] - mapping_offset]);
    sett_h1_atoms.push_back(import_map[cnk.settle_h1_atoms[sett_idx] - mapping_offset]);
    sett_h2_atoms.push_back(import_map[cnk.settle_h2_atoms[sett_idx] - mapping_offset]);
  }
  const std::vector<int> relevant_cnsts =
    vdel_pointer->getConstraintGroupAffectors(atom_move_list);
  const size_t ncnst = relevant_cnsts.size();
  cnst_group_bounds.push_back(0);
  for (size_t pos = 0; pos < ncnst; pos++) {
    const int cnst_idx = relevant_cnsts[pos];
    cnst_group_list.push_back(cnst_idx);
    for (int i = cnk.group_bounds[cnst_idx]; i < cnk.group_bounds[cnst_idx + 1]; i++) {
      cnst_group_atoms.push_back(import_map[cnk.group_list[i] - mapping_offset]);
    }
    cnst_group_bounds.push_back(cnst_group_atoms.size());
  }

  // Log the counts of each task
  bond_term_count   = bond_term_list.size();
  angl_term_count   = angl_term_list.size();
  dihe_term_count   = dihe_term_list.size();
  ubrd_term_count   = ubrd_term_list.size();
  cimp_term_count   = cimp_term_list.size();
  cmap_term_count   = cmap_term_list.size();
  infr14_term_count = infr14_term_list.size();
  rposn_term_count  = rposn_term_list.size();
  rbond_term_count  = rbond_term_list.size();
  rangl_term_count  = rangl_term_list.size();
  rdihe_term_count  = rdihe_term_list.size();
  cnst_group_count  = cnst_group_list.size();
  sett_group_count  = sett_group_list.size();
  vste_count        = virtual_site_list.size();

  // Compute the number of composite dihedrals: the dihedrals already subsume a great deal of the
  // 1:4 attenuated interactions (probably the entirety of these interactions if the system has no
  // virtual sites).  However, many dihedrals have additional cosine terms overlaid on the same
  // four atoms.  Bundle these terms in pairs, if possible, and eat the blank calculation if there
  // is no secondary dihedral to pair in any one case.  For a typical simulation with an Amber
  // force field, this will cut down on the total number of dihedral computations by about 40%,
  // and bundling the 1:4 interactions will likewise reduce the memory traffic.  The overall
  // reduction in memory bandwidth for reading instructions (relative to an approach that reads
  // each dihedral and 1:4 term separately, with a 64-bit instruction for each dihedral and a
  // 32-bit instruction for each 1:4 term) is roughly 30%, on top of other optimizations that
  // streamline the information access.  The composite dihedrals will replace the lists of standard
  // dihedral instructions and CHARMM impropers, fusing them into one.
  std::vector<int2> dihe_presence_bounds(imported_atom_count);
  std::vector<bool> presence_found(imported_atom_count, false);
  std::vector<bool> is_improper(dihe_term_count, false);
  for (int pos = 0; pos < dihe_term_count; pos++) {
    const int term_idx = dihe_term_list[pos];
    switch (static_cast<TorsionKind>(vk.dihe_modifiers[term_idx].w)) {
    case TorsionKind::PROPER_NO_14:
    case TorsionKind::PROPER:
      break;
    case TorsionKind::IMPROPER_NO_14:
    case TorsionKind::IMPROPER:
      is_improper[pos] = true;
    }
    const int i_atom = dihe_i_atoms[pos];
    const int j_atom = dihe_j_atoms[pos];
    const int k_atom = dihe_k_atoms[pos];
    const int l_atom = dihe_l_atoms[pos];
    if (presence_found[i_atom]) {
      dihe_presence_bounds[i_atom].y = pos;
    }
    else {
      dihe_presence_bounds[i_atom].x = pos;
      presence_found[i_atom] = true;
    }
    if (presence_found[j_atom]) {
      dihe_presence_bounds[j_atom].y = pos;
    }
    else {
      dihe_presence_bounds[j_atom].x = pos;
      presence_found[j_atom] = true;
    }
    if (presence_found[k_atom]) {
      dihe_presence_bounds[k_atom].y = pos;
    }
    else {
      dihe_presence_bounds[k_atom].x = pos;
      presence_found[k_atom] = true;
    }
    if (presence_found[l_atom]) {
      dihe_presence_bounds[l_atom].y = pos;
    }
    else {
      dihe_presence_bounds[l_atom].x = pos;
      presence_found[l_atom] = true;
    }
  }

  // Find all dihedrals that can be paired when acting on the same atoms
  std::vector<bool> coverage(dihe_term_count, false);
  for (int pos = 0; pos < dihe_term_count; pos++) {
    if (coverage[pos] || is_improper[pos]) {
      continue;
    }
    const int i_atom = dihe_i_atoms[pos];
    const int j_atom = dihe_j_atoms[pos];
    const int k_atom = dihe_k_atoms[pos];
    const int l_atom = dihe_l_atoms[pos];
    int min_np = std::max(dihe_presence_bounds[i_atom].x, dihe_presence_bounds[j_atom].x);
    min_np = std::max(min_np, dihe_presence_bounds[k_atom].x);
    min_np = std::max(min_np, dihe_presence_bounds[l_atom].x);
    min_np = std::max(min_np, pos + 1);
    int max_np = std::min(dihe_presence_bounds[i_atom].y, dihe_presence_bounds[j_atom].y);
    max_np = std::min(max_np, dihe_presence_bounds[k_atom].y);
    max_np = std::min(max_np, dihe_presence_bounds[l_atom].y);
    for (int npos = min_np; npos < min_np; npos++) {
      if (coverage[npos] || is_improper[npos]) {
        continue;
      }
      if ((dihe_i_atoms[npos] == dihe_i_atoms[pos] && dihe_j_atoms[npos] == dihe_j_atoms[pos] &&
           dihe_k_atoms[npos] == dihe_k_atoms[pos] && dihe_l_atoms[npos] == dihe_l_atoms[pos]) ||
          (dihe_i_atoms[npos] == dihe_l_atoms[pos] && dihe_j_atoms[npos] == dihe_k_atoms[pos] &&
           dihe_k_atoms[npos] == dihe_j_atoms[pos] && dihe_l_atoms[npos] == dihe_i_atoms[pos])) {

        // This dihedral applies to the same atoms.  Make a composite dihedral, then bail out.
        cdhe_term_list.push_back({dihe_term_list[pos], dihe_term_list[npos]});
        cdhe_is_cimp.push_back(false);
        cdhe_i_atoms.push_back(dihe_i_atoms[pos]);
        cdhe_j_atoms.push_back(dihe_j_atoms[pos]);
        cdhe_k_atoms.push_back(dihe_k_atoms[pos]);
        cdhe_l_atoms.push_back(dihe_l_atoms[pos]);
        coverage[pos] = true;
        coverage[npos] = true;
        break;
      }
    }
  }

  // Add the remaining dihedrals to the composite dihedral list as singlets
  for (int pos = 0; pos < dihe_term_count; pos++) {
    if (coverage[pos]) {
      continue;
    }
    cdhe_term_list.push_back({dihe_term_list[pos], -1});
    cdhe_is_cimp.push_back(false);
    cdhe_i_atoms.push_back(dihe_i_atoms[pos]);
    cdhe_j_atoms.push_back(dihe_j_atoms[pos]);
    cdhe_k_atoms.push_back(dihe_k_atoms[pos]);
    cdhe_l_atoms.push_back(dihe_l_atoms[pos]);
    coverage[pos] = true;
  }

  // Add CHARMM improper dihedrals to the composite list as more singlet terms
  for (int pos = 0; pos < cimp_term_count; pos++) {
    cdhe_term_list.push_back({cimp_term_list[pos], -1});
    cdhe_is_cimp.push_back(true);
    cdhe_i_atoms.push_back(cimp_i_atoms[pos]);
    cdhe_j_atoms.push_back(cimp_j_atoms[pos]);
    cdhe_k_atoms.push_back(cimp_k_atoms[pos]);
    cdhe_l_atoms.push_back(cimp_l_atoms[pos]);
    coverage[pos] = true;
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceWorkUnit::getAtomImportList() const {
  return atom_import_list;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> ValenceWorkUnit::getAtomManipulationMasks() const {
  const int nbits = static_cast<int>(sizeof(uint)) * 8;
  const int n_segment = roundUp((imported_atom_count + nbits - 1) / nbits, warp_size_int);
  std::vector<uint2> result(n_segment, {0U, 0U});
  for (int i = 0; i < moved_atom_count; i++) {
    const int seg_idx = (i / nbits);
    const int bit_idx = i - (seg_idx * nbits);
    result[seg_idx].x |= (0x1 << bit_idx);
  }
  for (int i = 0; i < updated_atom_count; i++) {
    const int seg_idx = (i / nbits);
    const int bit_idx = i - (seg_idx * nbits);
    result[seg_idx].y |= (0x1 << bit_idx);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ValenceWorkUnit::getTaskCounts() const {
  std::vector<int> result(static_cast<int>(VwuTask::TOTAL_TASKS));
  result[static_cast<int>(VwuTask::BOND)]   = bond_term_count;
  result[static_cast<int>(VwuTask::ANGL)]   = angl_term_count;
  result[static_cast<int>(VwuTask::DIHE)]   = dihe_term_count;
  result[static_cast<int>(VwuTask::UBRD)]   = ubrd_term_count;
  result[static_cast<int>(VwuTask::CIMP)]   = cimp_term_count;
  result[static_cast<int>(VwuTask::CDHE)]   = cdhe_term_count;
  result[static_cast<int>(VwuTask::CMAP)]   = cmap_term_count;
  result[static_cast<int>(VwuTask::INFR14)] = infr14_term_count;
  result[static_cast<int>(VwuTask::RPOSN)]  = rposn_term_count;
  result[static_cast<int>(VwuTask::RBOND)]  = rbond_term_count;
  result[static_cast<int>(VwuTask::RANGL)]  = rangl_term_count;
  result[static_cast<int>(VwuTask::RDIHE)]  = rdihe_term_count;
  result[static_cast<int>(VwuTask::CGROUP)] = cnst_group_count;
  result[static_cast<int>(VwuTask::SETTLE)] = sett_group_count;
  result[static_cast<int>(VwuTask::VSITE)]  = vste_count;
  return result;
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
  const AtomGraph *ag_ptr = vdel->getTopologyPointer();
  const RestraintApparatus *ra_ptr = vdel->getRestraintApparatusPointer();
  while (vdel->getFirstUnassignedAtom() < ag_ptr->getAtomCount()) {
    const int n_units = result.size();
    result.emplace_back(vdel, n_units, vdel->getFirstUnassignedAtom(), max_atoms_per_vwu);
  }

  // Loop once more over the update list and construct the atom movement list.  The movement list
  // is a superset of all the atoms that the work unit will update, and the import list is a
  // superset of the movement list.  This completes the atom sets that define any work unit.
  // Sort each of the sets.
  const int nvwu = result.size();
  for (int i = 0; i < nvwu; i++) {
    result[i].makeAtomMoveList();
    result[i].sortAtomSets();
  }
  
  // With the atom update assignments of each work unit known and the import list furnishing any
  // additional required dependencies (halo atoms), the task is now to loop over all updates and
  // trace back the force terms that must be computed.
  for (int i = 0; i < nvwu; i++) {
    result[i].logActivities();
  }
  

  return result;
}

} // namespace synthesis
} // namespace omni
