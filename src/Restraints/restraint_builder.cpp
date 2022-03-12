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
FlatBottomPlan::FlatBottomPlan() :
    activation_step{0}, k2{0.0}, k3{0.3}, r1{0.0}, r2{0.0}, r3{0.0}, r4{0.0}
{}

//-------------------------------------------------------------------------------------------------
FlatBottomPlan::FlatBottomPlan(const double k_in, const double r_in) :
    activation_step{0}, k2{k_in}, k3{k_in}, r1{r_in - 10.0}, r2{r_in}, r3{r_in}, r4{r_in + 10.0}
{}

//-------------------------------------------------------------------------------------------------
FlatBottomPlan::FlatBottomPlan(const double k_in, const double r2_in, const double r3_in) :
    activation_step{0}, k2{k_in}, k3{k_in}, r1{r2_in - 10.0}, r2{r2_in}, r3{r3_in},
    r4{r3_in + 10.0}
{}

//-------------------------------------------------------------------------------------------------
FlatBottomPlan::FlatBottomPlan(const double k2_in, const double k3_in, const double r2_in,
                               const double r3_in) :
    activation_step{0}, k2{k2_in}, k3{k3_in}, r1{r2_in - 10.0}, r2{r2_in}, r3{r3_in},
    r4{r3_in + 10.0}
{}

//-------------------------------------------------------------------------------------------------
FlatBottomPlan::FlatBottomPlan(const double k2_in, const double k3_in, const double r1_in,
                               const double r2_in, const double r3_in, const double r4_in,
                               const int step_in) :
    activation_step{step_in}, k2{k2_in}, k3{k3_in}, r1{r1_in}, r2{r2_in}, r3{r3_in}, r4{r4_in}
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint applyDistanceRestraint(const AtomGraph *ag, const int atom_i, const int atom_j,
                                        const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final) {
  return BoundedRestraint(atom_i, atom_j, ag, fb_init.activation_step,
                          fb_final.activation_step, fb_init.k2, fb_init.k3, fb_init.r1, fb_init.r2,
                          fb_init.r3, fb_init.r4, fb_final.k2, fb_final.k3, fb_final.r1,
                          fb_final.r2, fb_final.r3, fb_final.r4);
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint applyAngleRestraint(const AtomGraph *ag, const int atom_i, const int atom_j,
                                     const int atom_k, const FlatBottomPlan fb_init,
                                     const FlatBottomPlan fb_final) {
  return BoundedRestraint(atom_i, atom_j, atom_k, ag, fb_init.activation_step,
                          fb_final.activation_step, fb_init.k2, fb_init.k3, fb_init.r1, fb_init.r2,
                          fb_init.r3, fb_init.r4, fb_final.k2, fb_final.k3, fb_final.r1,
                          fb_final.r2, fb_final.r3, fb_final.r4);
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint applyDihedralRestraint(const AtomGraph *ag, const int atom_i, const int atom_j,
                                        const int atom_k, const int atom_l,
                                        const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final) {
  return BoundedRestraint(atom_i, atom_j, atom_k, atom_l, ag, fb_init.activation_step,
                          fb_final.activation_step, fb_init.k2, fb_init.k3, fb_init.r1, fb_init.r2,
                          fb_init.r3, fb_init.r4, fb_final.k2, fb_final.k3, fb_final.r1,
                          fb_final.r2, fb_final.r3, fb_final.r4);
}

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
                          const double displacement_penalty, const double displacement_onset,
                          const double displacement_plateau, const double proximity_penalty,
                          const double proximity_onset, const double proximity_plateau) {
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

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyHydrogenBondingPreventors(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                               const AtomMask &mask, const double penalty,
                               const double flat_bottom_half_width,
                               const double harmonic_maximum) {
  restraintTopologyChecks(ag, cframe, mask);
  
  // Seek out all proton donors and acceptors: donors are defined as electronegative atoms N, O,
  // S, or P which have a hydrogen attached to them, while acceptors are electronegative atoms N,
  // O, S, or P that have a lone pair available in the Lewis Structure.
  
}

} // namespace restraints
} // namespace omni
