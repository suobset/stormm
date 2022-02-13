#include "Constants/symbol_values.h"
#include "Restraints/bounded_restraint.h"
#include "Structure/local_arrangement.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "bounded_restraint.h"
#include "restraint_apparatus.h"
#include "restraint_builder.h"

namespace omni {
namespace restraints {

using chemistry::MaskTraversalMode;
using restraints::BoundedRestraint;
using structure::dihedral_angle;
using structure::imageValue;
using structure::ImagingMethod;
using symbols::twopi;
using topology::TorsionKind;
using topology::ValenceKit;

//-------------------------------------------------------------------------------------------------
void restraintTopologyChecks(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                             const AtomMask &mask) {
  const int natom_expected = ag->getAtomCount();
  if (cframe.natom != natom_expected) {
    rtErr("The coordinates must match the input topology.  Atom counts of " +
          std::to_string(natom_expected) + " and " + std::to_string(cframe.natom) + " differ.",
          "applyPositionalRestraints");
  }
  if (natom_expected != mask.getAtomGraphPointer()->getAtomCount()) {
    rtErr("The atom mask must match the input topology.  Atom counts of " +
          std::to_string(natom_expected) + " and " +
          std::to_string(mask.getAtomGraphPointer()->getAtomCount()) + " differ.",
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
  std::vector<BoundedRestraint> result(nmasked, BoundedRestraint(ag));
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
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                       const AtomMask &mask, const double penalty,
                       const double flat_bottom_half_width, const double harmonic_maximum) {
  restraintTopologyChecks(ag, cframe, mask);

  // Loop over the assigned dihedrals of all atoms in the topology.  If there are multiple
  // dihedrals impacting the same set of atoms, they will have been assigned to the same atom, or
  // if in reverse to one alternative atom.  Search all of these dihedrals and find all unique
  // dihedrals such that all four atoms are heavy atoms (non-hydrogen) and neither of the central
  // two atoms are one of the chiral centers slated for inversion.
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const std::vector<bool> holding_mask = mask.getMask();
  const double harmonic_width = sqrt(harmonic_maximum / penalty);
  std::vector<bool> coverage(vk.ndihe, false);
  std::vector<BoundedRestraint> result;
  for (int i = 0; i < vk.natom; i++) {
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
      const int atom_i = vk.dihe_asgn_atoms[3 * j];
      const int atom_j = vk.dihe_asgn_atoms[(3 * j) + 1];
      const int atom_k = i;
      const int atom_l = vk.dihe_asgn_atoms[(3 * j) + 2];
      for (int k = j; k < vk.dihe_asgn_bounds[i + 1]; k++) {
        if (vk.dihe_asgn_atoms[3 * k] == atom_i && vk.dihe_asgn_atoms[(3 * k) + 1] == atom_j &&
            vk.dihe_asgn_atoms[(3 * k) + 2] == atom_l) {
          coverage[k] = true;
        }
      }
      if (atom_k < atom_j) {
        for (int k = vk.dihe_asgn_bounds[atom_j]; k < vk.dihe_asgn_bounds[atom_j + 1]; k++) {
          if (vk.dihe_asgn_atoms[3 * k] == atom_l && vk.dihe_asgn_atoms[(3 * k) + 1] == atom_k &&
              vk.dihe_asgn_atoms[(3 * k) + 2] == atom_i) {
            coverage[k] = true;
          }
        }
      }
      const double current_value = dihedral_angle(atom_i, atom_j, atom_k, atom_l, cframe);
      const double r1 = imageValue(current_value - flat_bottom_half_width - harmonic_width, twopi,
                                   ImagingMethod::MINIMUM_IMAGE);
      const double r2 = imageValue(current_value - flat_bottom_half_width, twopi,
                                   ImagingMethod::MINIMUM_IMAGE);
      const double r3 = imageValue(current_value + flat_bottom_half_width, twopi,
                                   ImagingMethod::MINIMUM_IMAGE);
      const double r4 = imageValue(current_value + flat_bottom_half_width + harmonic_width, twopi,
                                   ImagingMethod::MINIMUM_IMAGE);
      result.push_back(BoundedRestraint(atom_i, atom_j, atom_k, atom_l, ag, penalty, penalty,
                                        r1, r2, r3, r4));
    }
  }
  return result;
}

} // namespace restraints
} // namespace omni
