#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "bounded_restraint.h"
#include "restraint_apparatus.h"
#include "restraint_builder.h"

namespace omni {
namespace restraints {

using chemistry::MaskTraversalMode;
using topology::TorsionKind;
using topology::ValenceKit;

//-------------------------------------------------------------------------------------------------
bool restraintTopologyChecks(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                             const AtomMask &mask) {
  const natom_expected = ag.getAtomCount();
  if (cframe.natom != natom_expected) {
    rtErr("The coordinates must match the input topology.  Atom counts of " +
          std::to_string(natom_expected) + " and " + std::to_string(cframe.natom) + " differ.",
          "applyPositionalRestraints");
  }
  if (natom_expected != mask.ag_pointer->getAtomCount) {
    rtErr("The atom mask must match the input topology.  Atom counts of " +
          std::to_string(natom_expected) + " and " +
          std::to_string(mask.ag_pointer->getAtomCount()) + " differ.",
          "applyPositionalRestraints");
  }
}
  
//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                          const CoordinateFrameReader &reference_cframe, const AtomMask &mask,
                          double displacement_penalty, double displacement_onset,
                          double displacement_plateau, double proximity_penalty,
                          double proximity_onset, double proximity_plateau) {
  restraintTopologyChecks(ag, cframe, mask);
  
  // Loop over all atoms in the mask and create restraints
  const int nmasked = mask.getMaskedAtomCount();
  std::vector<BoundedRestraint> result(nmasked);
  const std::vector<int> masked_atoms = mask.getMaskedAtomList();
  for (size_t i = 0; i < nmasked; i++) {
    result[i] = BoundedRestraint(masked_atoms[i], ag, cframe, proximity_penalty,
                                 displacement_penalty, proximity_plateau, proximity_onset,
                                 displacement_onset, displacement_plateau);
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyChiralityInversionRestraints(const AtomGraph &ag, const CoordinateFrameReader &cframe,
                                  const std::vector<int> chiral_atoms) {
  restraintTopologyChecks(ag, cframe, mask);

  // Loop over the assigned dihedrals of all atoms in the topology.  If there are multiple
  // dihedrals impacting the same set of atoms, they will have been assigned to the same atom, or
  // if in reverse to one alternative atom.  Search all of these dihedrals and find all unique
  // dihedrals such that all four atoms are heavy atoms (non-hydrogen) and neither of the central
  // two atoms are one of the chiral centers slated for inversion.
  const int nchiral = chiral_atoms.size();
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  int max_dihe_assignments = 0;
  for (int i = 0; i < vk.natom; i++) {
    max_dihe_assignments = std::max(max_dihe_assignments,
                                    vk.dihe_asgn_bounds[i + 1] - vk.dihe_asgn_bounds[i]);
  }
  std::vector<bool> coverage(max_dihe_assignments);
  std::vector<BoundedRestraint> result;
  for (int i = 0; i < vk.natom; i++) {
    for (int j = vk.dihe_asgn_bounds[i]; j < vk.dihe_asgn_bounds[i + 1]; j++) {
      coverage[j - vk.dihe_asgn_bounds[i]] = false;
    }
    for (int j = vk.dihe_asgn_bounds[i]; j < vk.dihe_asgn_bounds[i + 1]; j++) {
      if (coverage[j]) {
        continue;
      }

      // Check that the dihedral is not an improper.
      const int term_index =  vk.dihe_asgn_terms[j];
      const TorsionKind tk = static_cast<TorsionKind>(vk.dihe_modifiers[term_index].w);
      if (tk == TorsionKind::IMPROPER || tk == TorsionKind::IMPROPER_NO_14) {
        continue;
      }

      // Mark this combination of four atoms as representative of any similar combinations.  When
      // checking the second of two central atoms, give priority to dihedrals with the lower atom
      // index in the controlling position.
      const int atom_i = dihe_asgn_atoms[3 * j];
      const int atom_j = dihe_asgn_atoms[(3 * j) + 1];
      const int atom_k = i;
      const int atom_l = dihe_asgn_atoms[(3 * j) + 2];
      for (int k = j; k < vk.dihe_asgn_bounds[i + 1]; k++) {
        if (dihe_asgn_atoms[3 * k] == atom_i && dihe_asgn_atoms[(3 * k) + 1] == atom_j &&
            dihe_asgn_atoms[(3 * k) + 2] == atom_l) {
          coverage[k - vk.dihe_asgn_bounds[i]] = true;
        }
      }
      bool unique = true;
      if (atom_k > atom_j) {
        for (int k = vk.dihe_asgn_bounds[atom_j]; k < vk.dihe_asgn_bounds[atom_j + 1]; k++) {
          if (dihe_asgn_atoms[3 * k] == atom_l && dihe_asgn_atoms[(3 * k) + 1] == atom_k &&
              dihe_asgn_atoms[(3 * k) + 2] == atom_i) {
            unique = false;
          }
        }
      }
      if (unique) {
        result.push_back(atom_i, atom_j, atom_k, atom_l,         
      }
    }
  }
}

} // namespace restraints
} // namespace omni
