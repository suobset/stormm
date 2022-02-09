#include "Chemistry/atommask.h"
#include "bounded_restraint.h"

namespace omni {
namespace restraints {

using chemistry::AtomMask;
using chemistry::MaskInputMode;

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const std::string &mask_l_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in, const double3 init_ref_crd_in,
                                   const double3 final_ref_crd_in) :
    atom_i{-1}, atom_j{-1}, atom_k{-1}, atom_l{-1}, order{0},
    initial_step{init_step_in},
    final_step{final_step_in},
    initial_keq{init_k2_in, init_k3_in},
    initial_r{init_r1_in, init_r2_in, init_r3_in, init_r4_in},
    final_keq{final_k2_in, final_k3_in},
    final_r{final_r1_in, final_r2_in, final_r3_in, final_r4_in},
    init_center{init_ref_crd_in},
    final_center{final_ref_crd_in},
    ag_pointer{ag_in}
{
  // At least one atom is required
  if (mask_i_in.size() == 0LLU) {
    rtErr("At least one atom must be be specified in order to create a bounded restraint.",
          "BoundedRestraint");
  }

  // Obtain the atoms to which this restraint applies based on the input atom masks.
  std::vector<AtomMask> all_masks = {
    AtomMask(mask_i_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom I mask"),
    AtomMask(mask_j_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom J mask"),
    AtomMask(mask_k_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom K mask"),
    AtomMask(mask_l_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom : mask") };
  std::vector<int> n_mask_atom(4);
  bool running = true;
  for (int i = 0; i < 4; i++) {
    n_mask_atom[i] = all_masks[i].getMaskedAtomCount();
    if (n_mask_atom[i] > 1) {
      rtErr("Atom mask \"" + all_masks[i].getInputText() + "\" specifies " +
            std::to_string(n_mask_atom[i]) + " atoms in topology " + ag_pointer->getFileName() +
            ".  The mask must specify one and only one atom.", "BoundedRestraint");
    }
    running = (running && (n_mask_atom[i] == 1));
    order += running;
  }

  // Report if the first atom mask is invalid
  if (n_mask_atom[0] == 1) {
    atom_i = all_masks[0].getMaskedAtomList()[0];
  }
  else {
    rtErr("Atom mask \"" + mask_i_in + "\" specifies " + std::to_string(n_mask_atom[0]) +
          " atoms in topology " + ag_pointer->getFileName() + ".  The mask must specify one and "
          "only one atom.", "BoundedRestraint");
  }

  // Set other atoms
  if (order > 1) {
    atom_j = n_mask_atom[1];
  }
  if (order > 2) {
    atom_k = n_mask_atom[2];
  }
  if (order > 3) {
    atom_l = n_mask_atom[3];
  }
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, std::string(""), ag_in, chemfe, cfr,
                     init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in, init_r2_in,
                     init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in, final_r2_in,
                     final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, std::string(""), std::string(""), ag_in, chemfe, cfr,
                     init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in, init_r2_in,
                     init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in, final_r2_in,
                     final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in, const double3 init_ref_crd_in,
                                   const double3 final_ref_crd_in) :
    BoundedRestraint(mask_i_in, std::string(""), std::string(""), std::string(""), ag_in, chemfe,
                     cfr, init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in,
                     init_r2_in, init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in,
                     final_r2_in, final_r3_in, final_r4_in, init_ref_crd_in, final_ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const std::string &mask_l_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, mask_l_in, ag_in, chemfe, cfr, 0, 0,
                     k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in,
                     r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, std::string(""), ag_in, chemfe, cfr, 0, 0,
                     k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in,
                     r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, std::string(""), std::string(""), ag_in, chemfe, cfr, 0,
                     0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in,
                     r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures *chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in,
                                   const double3 ref_crd_in) :
    BoundedRestraint(mask_i_in, std::string(""), std::string(""), std::string(""), ag_in, chemfe,
                     cfr, 0, 0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, ref_crd_in, ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const int atom_l_in, const AtomGraph *ag_in,
                                   const int init_step_in, const int final_step_in,
                                   const double init_k2_in, const double init_k3_in,
                                   const double init_r1_in, const double init_r2_in,
                                   const double init_r3_in, const double init_r4_in,
                                   const double final_k2_in, const double final_k3_in,
                                   const double final_r1_in, const double final_r2_in,
                                   const double final_r3_in, const double final_r4_in,
                                   const double3 init_ref_crd_in, const double3 final_ref_crd_in) :
    atom_i{atom_i_in}, atom_j{atom_j_in}, atom_k{atom_k_in}, atom_l{atom_l_in}, order{0},
    initial_step{init_step_in},
    final_step{final_step_in},
    initial_keq{init_k2_in, init_k3_in},
    initial_r{init_r1_in, init_r2_in, init_r3_in, init_r4_in},
    final_keq{final_k2_in, final_k3_in},
    final_r{final_r1_in, final_r2_in, final_r3_in, final_r4_in},
    init_center{init_ref_crd_in},
    final_center{final_ref_crd_in},
    ag_pointer{ag_in}
{
  // At least one atom is required
  if (atom_i < 0) {
    rtErr("At least one atom must be be specified in order to create a bounded restraint.",
          "BoundedRestraint");
  }

  // Get the order of the restraint
  order = (atom_i >= 0) + (atom_j >= 0) + (atom_k >= 0) + (atom_l >= 0);

  // Confirm that all atom selections are valid
  const int natom = ag_pointer->getAtomCount();
  if (atom_i >= natom || (order > 1 && atom_j >= natom) || (order > 2 && atom_k >= natom) ||
      (order > 3 && atom_l >= natom)) {
    std::string atom_number_list = std::to_string(atom_i);
    if (order > 1) {
      atom_number_list += ", " + std::to_string(atom_j);
    }
    if (order > 2) {
      atom_number_list += ", " + std::to_string(atom_k);
    }
    if (order > 3) {
      atom_number_list += ", " + std::to_string(atom_l);
    }
    rtErr("The atom selection " + atom_number_list + " does not fall within the atom indexing of "
          "topology " + ag_pointer->getFileName() + ".", "BoundedRestraint");
  }
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const AtomGraph *ag_in, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, -1, ag_in, init_step_in, final_step_in,
                     init_k2_in, init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in,
                     final_k2_in, final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in,
                                   const AtomGraph *ag_in, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, -1, -1, ag_in, init_step_in, final_step_in,
                     init_k2_in, init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in,
                     final_k2_in, final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const AtomGraph *ag_in,
                                   const int init_step_in, const int final_step_in,
                                   const double init_k2_in, const double init_k3_in,
                                   const double init_r1_in, const double init_r2_in,
                                   const double init_r3_in, const double init_r4_in,
                                   const double final_k2_in, const double final_k3_in,
                                   const double final_r1_in, const double final_r2_in,
                                   const double final_r3_in, const double final_r4_in,
                                   const double3 init_ref_crd_in, const double3 final_ref_crd_in) :
    BoundedRestraint(atom_i_in, -1, -1, -1, ag_in, init_step_in, final_step_in, init_k2_in,
                     init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in, final_k2_in,
                     final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in,
                     init_ref_crd_in, final_ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const int atom_l_in, const AtomGraph *ag_in, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, atom_l_in, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const AtomGraph *ag_in, const double k2_in, const double k3_in,
                                   const double r1_in, const double r2_in, const double r3_in,
                                   const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in,
                                   const AtomGraph *ag_in, const double k2_in, const double k3_in,
                                   const double r1_in, const double r2_in, const double r3_in,
                                   const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const AtomGraph *ag_in, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in, double3 ref_crd_in) :
    BoundedRestraint(atom_i_in, -1, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_index, const AtomGraph *ag_in,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in, const int refr_index) :
  BoundedRestraint(atom_index, -1, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in,
                   k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 })
{
  // Make a bounds check on the input before setting the positional restraint center
  if (refr_index < 0 || refr_index >= cfr.natom) {
    rtErr("The requested atom index " + std::to_string(refr_index) + " was invalid for a "
          "coordinate frame with " + std::to_string(cfr.natom) + " atoms.", "BoundedRestraint");
  }
  init_center = { cfr.xcrd[refr_index], cfr.ycrd[refr_index], cfr.zcrd[refr_index] };
  final_center = init_center;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getAtomIndex(const int restrained_atom_number) const {
  switch (restrained_atom_number) {
  case 1:
    return atom_i;
  case 2:
    return atom_j;
  case 3:
    return atom_k;
  case 4:
    return atom_l;
  default:
    rtErr("Valid atoms include 1, 2, 3, and 4, not " + std::to_string(restrained_atom_number) +
          ".", "BoundedRestraint", "getAtomIndex");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getOrder() const {
  return order;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getInitialStep() const {
  return initial_step;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getFinalStep() const {
  return final_step;
}

//-------------------------------------------------------------------------------------------------
double2 BoundedRestraint::getInitialStiffness() const {
  return initial_keq;
}

//-------------------------------------------------------------------------------------------------
double2 BoundedRestraint::getFinalStiffness() const {
  return final_keq;
}

//-------------------------------------------------------------------------------------------------
double4 BoundedRestraint::getInitialDisplacements() const {
  return initial_r;
}

//-------------------------------------------------------------------------------------------------
double4 BoundedRestraint::getFinalDisplacements() const {
  return final_r;
}

//-------------------------------------------------------------------------------------------------
double3 BoundedRestraint::getInitialTargetSite() const {

  // Check the order of the restraint 
  if (order != 1) {
    rtErr("The target site is meaningless outside of positional restraints on individual atoms.",
          "BoundedRestraint", "getInitialTargetSite");
  }
  return init_center;
}

//-------------------------------------------------------------------------------------------------
double3 BoundedRestraint::getFinalTargetSite() const {

  // Check the order of the restraint 
  if (order != 1) {
    rtErr("The target site is meaningless outside of positional restraints on individual atoms.",
          "BoundedRestraint", "getFinalTargetSite");
  }
  return final_center;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* BoundedRestraint::getTopologyPointer() const {
  return ag_pointer;
}

} // namespace restraints
} // namespace omni
