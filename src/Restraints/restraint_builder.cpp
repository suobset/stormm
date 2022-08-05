#include "Constants/symbol_values.h"
#include "DataTypes/stormm_vector_types.h"
#include "Restraints/bounded_restraint.h"
#include "Structure/local_arrangement.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "bounded_restraint.h"
#include "restraint_apparatus.h"
#include "restraint_builder.h"

namespace stormm {
namespace restraints {

using chemistry::ChemicalFeatures;
using chemistry::MaskTraversalMode;
using restraints::BoundedRestraint;
using structure::dihedral_angle;
using structure::imageValue;
using structure::ImagingMethod;
using symbols::twopi;
using topology::TorsionKind;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
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
                                        const FlatBottomPlan fb_init) {
  return applyDistanceRestraint(ag, atom_i, atom_j, fb_init, fb_init);
}

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
                                     const int atom_k, const FlatBottomPlan fb_init) {
  return applyAngleRestraint(ag, atom_i, atom_j, atom_k, fb_init, fb_init);
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
                                        const FlatBottomPlan fb_init) {
  return applyDihedralRestraint(ag, atom_i, atom_j, atom_k, atom_l, fb_init, fb_init);
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
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const AtomMask &mask, const double displacement_penalty,
                          const double displacement_onset, const double displacement_plateau,
                          const double proximity_penalty, const double proximity_onset,
                          const double proximity_plateau) {
  restraintTopologyChecks(ag, ref_cfr, mask);
  
  // Loop over all atoms in the mask and create restraints
  return applyPositionalRestraints(ag, ref_cfr,  mask.getMaskedAtomList(), displacement_penalty,
                                   displacement_onset, displacement_plateau, proximity_penalty,
                                   proximity_onset, proximity_plateau);
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const std::string &mask, const double displacement_penalty,
                          const double displacement_onset, const double displacement_plateau,
                          const double proximity_penalty, const double proximity_onset,
                          const double proximity_plateau) {
  ChemicalFeatures chemfe(ag, ref_cfr);
  return applyPositionalRestraints(ag, ref_cfr, AtomMask(mask, ag, &chemfe, ref_cfr),
                                   displacement_penalty, displacement_onset, displacement_plateau,
                                   proximity_penalty, proximity_onset, proximity_plateau);
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const std::vector<int> &masked_atoms, const double displacement_penalty,
                          const double displacement_onset, const double displacement_plateau,
                          const double proximity_penalty, const double proximity_onset,
                          const double proximity_plateau) {
  const int nmasked = masked_atoms.size();
  std::vector<BoundedRestraint> result;
  result.reserve(nmasked);
  for (size_t i = 0; i < nmasked; i++) {
    const int atom_i = masked_atoms[i];
    const double3 target = { ref_cfr.xcrd[atom_i], ref_cfr.ycrd[atom_i], ref_cfr.zcrd[atom_i] };
    result.emplace_back(atom_i, ag, proximity_penalty, displacement_penalty, proximity_plateau,
                        proximity_onset, displacement_onset, displacement_plateau, target);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                       const AtomMask &mask, const double penalty,
                       const double flat_bottom_half_width, const double harmonic_maximum) {
  restraintTopologyChecks(ag, cfr, mask);

  // Loop over the assigned dihedrals of all atoms in the topology.  If there are multiple
  // dihedrals impacting the same set of atoms, they will have been assigned to the same atom, or
  // if in reverse to one alternative atom.  Search all of these dihedrals and find all unique
  // dihedrals such that all four atoms are heavy atoms (non-hydrogen) and neither of the central
  // two atoms are one of the chiral centers slated for inversion.
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
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
      // index in the controlling position.  Skip dihedrals containing hydrogens.
      const int atom_i = vk.dihe_asgn_atoms[3 * j];
      const int atom_j = vk.dihe_asgn_atoms[(3 * j) + 1];
      const int atom_k = i;
      const int atom_l = vk.dihe_asgn_atoms[(3 * j) + 2];
      if (cdk.z_numbers[atom_i] <= 1 || cdk.z_numbers[atom_j] <= 1 || cdk.z_numbers[atom_k] <= 1 ||
          cdk.z_numbers[atom_l] <= 1) {
        continue;
      }
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
      const double current_value = dihedral_angle(atom_i, atom_j, atom_k, atom_l, cfr);
      const double r1 = current_value - flat_bottom_half_width - harmonic_width;
      const double r2 = current_value - flat_bottom_half_width;
      const double r3 = current_value + flat_bottom_half_width;
      const double r4 = current_value + flat_bottom_half_width + harmonic_width;
      result.push_back(BoundedRestraint(atom_i, atom_j, atom_k, atom_l, ag, penalty, penalty,
                                        r1, r2, r3, r4));
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                       const std::string &mask, const double penalty,
                       const double flat_bottom_half_width, const double harmonic_maximum) {
  ChemicalFeatures chemfe(ag, cfr);
  return applyHoldingRestraints(ag, cfr, AtomMask(mask, ag, &chemfe, cfr), penalty,
                                flat_bottom_half_width, harmonic_maximum);
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyHydrogenBondPreventors(const AtomGraph *ag, const ChemicalFeatures &chemfe,
                            const double penalty, const double approach_point) {
  if (chemfe.getTopologyPointer() != ag) {
    rtErr("AtomGraph pointers do not match.", "applyHydrogenBondingPreventors");
  }

  // Prepare the flat-bottom well template based on the input and a generous assumption
  // about the maximum size of a molecule
  FlatBottomPlan fbhw_template(penalty, 0.0, 0.0, approach_point, 1000.0, 1100.0, 0);
  
  // Get lists of all proton donors and acceptors
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  std::vector<int> donors = chemfe.getHydrogenBondDonorList();
  std::vector<int> acceptors = chemfe.getHydrogenBondAcceptorList();

  // Apply restraints such that both donor and acceptor are on the same molecule and not
  // so close that they might be separated by a 1:4 interaction.
  std::vector<BoundedRestraint> result;
  int ndonor = donors.size();
  int nacceptor = acceptors.size();

  // The fact that this algorithm will ultimately test ndonor x nacceptor atoms is potentially
  // harmful to performance.  Remove atoms which cannot be of interest.
  if (cdk.nmol > 1) {
    std::vector<int> heavy_atom_content(cdk.nmol, 0);
    for (int i = 0; i < cdk.natom; i++) {
      if (cdk.z_numbers[i] > 4) {
        heavy_atom_content[cdk.mol_home[i]] += 1;
      }
    }
    int j = 0;
    for (int i = 0; i < ndonor; i++) {
      donors[j] = donors[i];
      j += (heavy_atom_content[cdk.mol_home[donors[i]]] > 4);
    }
    donors.resize(j);
    j = 0;
    for (int i = 0; i < nacceptor; i++) {
      acceptors[j] = acceptors[i];
      j += (heavy_atom_content[cdk.mol_home[acceptors[i]]] > 4);
    }
    acceptors.resize(j);
  }

  // Run through the (trimmed) lists to add restraints that prevent internal hydrogen bonding.
  // For systems with many distinct molecules each having donors and / or acceptors, this search
  // is still far less efficient than it could be, but the problem has been reduced to edge cases.
  for (int i = 0; i < ndonor; i++) {
    const int donor_atom = donors[i];
    for (int j = 0; j < nacceptor; j++) {
      const int acceptor_atom = acceptors[j];

      // Do not apply restraints to self-interactions.  Only apply restraints to atoms on the
      // same molecule.
      if (donor_atom == acceptor_atom ||
          cdk.mol_home[donor_atom] != cdk.mol_home[acceptor_atom]) {
        continue;
      }
      bool excluded = false;
      for (int k = nbk.nb12_bounds[donor_atom]; k < nbk.nb12_bounds[donor_atom]; k++) {
        excluded = (excluded || nbk.nb12x[k] == acceptor_atom);
      }
      if (excluded == false) {
        for (int k = nbk.nb13_bounds[donor_atom]; k < nbk.nb13_bounds[donor_atom]; k++) {
          excluded = (excluded || nbk.nb13x[k] == acceptor_atom);
        }
      }
      if (excluded == false) {
        for (int k = nbk.nb14_bounds[donor_atom]; k < nbk.nb14_bounds[donor_atom]; k++) {
          excluded = (excluded || nbk.nb14x[k] == acceptor_atom);
        }
      }
      if (excluded == false) {

        // Use the push_back method without pre-allocation, as it is not clear from the start
        // how large the final array will need to be.
        result.push_back(applyDistanceRestraint(ag, donor_atom, acceptor_atom, fbhw_template));
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyHydrogenBondPreventors(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                            const double penalty, const double approach_point) {
  ChemicalFeatures chemfe(ag, cfr);
  return applyHydrogenBondPreventors(ag, chemfe, penalty, approach_point);
}

} // namespace restraints
} // namespace stormm
