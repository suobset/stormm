#include <algorithm>
#include "Constants/scaling.h"
#include "DataTypes/mixed_types.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Topology/atomgraph_analysis.h"
#include "UnitTesting/approx.h"
#include "chemical_features.h"
#include "indigo.h"

namespace omni {
namespace chemistry {

using card::HybridKind;
using topology::TorsionKind;
using math::crossProduct;
using math::dot;
using math::numberSeriesToBitMask;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::project;
using math::roundUp;
using testing::Approx;
using topology::selectRotatingAtoms;
using trajectory::CoordinateFrameReader;

//-------------------------------------------------------------------------------------------------
BondedNode::BondedNode() :
    previous_atom_index{-1}, atom_index{-1}, layer_index{-1}, root_bond_order{0.0},
    branch_count{0}, branch_atoms{nullptr}, rings_completed{0}
{}

//-------------------------------------------------------------------------------------------------
void BondedNode::setBranchPointer(std::vector<int> *vi, const size_t pos,
                                  const size_t max_branches) {

  // Check that the vector has enough data to accept the pointer
  if (vi->size() < (pos + 1) * max_branches) {
    rtErr("Storage vector does not have sufficient space.", "BondedNode");
  }
  branch_atoms = &(vi->data()[pos * max_branches]);
}

//-------------------------------------------------------------------------------------------------
void BondedNode::addToTree(const int previous_in, const int current_atom, const int current_layer,
                           const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk) {
  previous_atom_index = previous_in;
  atom_index = current_atom;
  layer_index = current_layer;
  const int excl_start = nbk.nb12_bounds[current_atom];
  const int link_candidate_count = nbk.nb12_bounds[current_atom + 1] - excl_start;
  int j = 0;
  for (int i = 0; i < link_candidate_count; i++) {
    const int candidate_atom = nbk.nb12x[i + excl_start];
    if (candidate_atom != previous_in && cdk.z_numbers[candidate_atom] != 0) {
      branch_atoms[j] = candidate_atom;
      j++;
    }
  }
  branch_count = j;
}

//-------------------------------------------------------------------------------------------------
void BondedNode::addBondOrder(const ValenceKit<double> &vk, const Hybrid<double> &bond_orders) {
  for (int i = vk.bond_asgn_bounds[atom_index]; i < vk.bond_asgn_bounds[atom_index + 1]; i++) {
    if (vk.bond_asgn_atoms[i] == previous_atom_index) {
      root_bond_order = bond_orders.readHost(vk.bond_asgn_terms[i]);
      return;
    }
  }
  for (int i = vk.bond_asgn_bounds[previous_atom_index];
       i < vk.bond_asgn_bounds[previous_atom_index + 1]; i++) {
    if (vk.bond_asgn_atoms[i] == atom_index) {
      root_bond_order = bond_orders.readHost(vk.bond_asgn_terms[i]);
      return;
    }
  }
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getPreviousAtom() const {
  return previous_atom_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getAtom() const {
  return atom_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getLayer() const {
  return layer_index;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getBranchCount() const {
  return branch_count;
}

//-------------------------------------------------------------------------------------------------
int BondedNode::getBranchAtom(const int index) const {
  return branch_atoms[index];
}

//-------------------------------------------------------------------------------------------------
int BondedNode::findBranchAtom(const int search_index) const {
  for (int i = 0; i < branch_count; i++) {
    if (branch_atoms[i] == search_index) {
      return i;
    }
  }
  rtErr("The atom with topology index " + std::to_string(search_index + 1) + " is no branch of " +
        "atom " + std::to_string(atom_index + 1) + ".", "BondedNode", "findBranchAtom");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
uint BondedNode::getRingCompletion(const int branch_index) const {
  return ((rings_completed >> branch_index) & 0x1);
}

//-------------------------------------------------------------------------------------------------
double BondedNode::getRootBondOrder() const {
  return root_bond_order;
}
  
//-------------------------------------------------------------------------------------------------
void BondedNode::setRingCompletion(const int branch_index) {
  rings_completed |= (0x1 << branch_index);
}

//-------------------------------------------------------------------------------------------------
void BondedNode::wipeRingCompletion() {
  rings_completed = 0U;
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures() :
    atom_count{0}, planar_atom_count{0}, ring_count{0}, fused_ring_count{0},
    twistable_ring_count{0}, conjugated_group_count{0}, aromatic_group_count{0},
    polar_hydrogen_count{0}, hbond_donor_count{0}, hbond_acceptor_count{0}, chiral_center_count{0},
    rotatable_bond_count{0}, double_bond_count{0}, triple_bond_count{0}, max_ring_size{0},
    temperature{0.0}, rotating_groups_mapped{false},
    planar_centers{HybridKind::POINTER, "chemfe_planarity"},
    ring_inclusion{HybridKind::ARRAY, "chemfe_rings"},
    ring_atom_bounds{HybridKind::POINTER, "chemfe_ring_bounds"},
    ring_atoms{HybridKind::POINTER, "chemfe_ring_atoms"},
    aromatic_group_bounds{HybridKind::POINTER, "chemfe_arom_bounds"},
    aromatic_pi_electrons{HybridKind::POINTER, "chemfe_pi_elec"},
    aromatic_groups{HybridKind::POINTER, "chemfe_arom_groups"},
    polar_hydrogens{HybridKind::POINTER, "chemfe_polar_h"},
    hydrogen_bond_donors{HybridKind::POINTER, "chemfe_hbond_donor"},
    hydrogen_bond_acceptors{HybridKind::POINTER, "chemfe_hbond_accpt"},
    chiral_centers{HybridKind::POINTER, "chemfe_chirals"},
    rotatable_groups{HybridKind::POINTER, "chemfe_rotators"},
    rotatable_group_bounds{HybridKind::POINTER, "chemfe_rotator_bounds"},
    formal_charges{HybridKind::POINTER, "chemfe_formal_charges"},
    bond_orders{HybridKind::POINTER, "chemfe_bond_orders"},
    free_electrons{HybridKind::POINTER, "chemfe_free_e"},
    int_data{HybridKind::ARRAY, "chemfe_int"},
    double_data{HybridKind::ARRAY, "chemfe_double"},
    ag_pointer{nullptr}
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrameReader &cfr,
                                   const MapRotatableGroups map_group_in,
                                   const double temperature_in) :
    atom_count{ag_in->getAtomCount()}, planar_atom_count{0}, ring_count{0}, fused_ring_count{0},
    twistable_ring_count{0}, conjugated_group_count{0}, aromatic_group_count{0},
    polar_hydrogen_count{0}, hbond_donor_count{0}, hbond_acceptor_count{0}, chiral_center_count{0},
    rotatable_bond_count{0}, double_bond_count{0}, triple_bond_count{0},
    max_ring_size{8 * sizeof(ullint)}, temperature{temperature_in}, rotating_groups_mapped{false},
    planar_centers{HybridKind::POINTER, "chemfe_planarity"},
    ring_inclusion{HybridKind::ARRAY, "chemfe_rings"},
    ring_atom_bounds{HybridKind::POINTER, "chemfe_ring_bounds"},
    ring_atoms{HybridKind::POINTER, "chemfe_ring_atoms"},
    aromatic_group_bounds{HybridKind::POINTER, "chemfe_arom_bounds"},
    aromatic_pi_electrons{HybridKind::POINTER, "chemfe_pi_elec"},
    aromatic_groups{HybridKind::POINTER, "chemfe_arom_groups"},
    polar_hydrogens{HybridKind::POINTER, "chemfe_polar_h"},
    hydrogen_bond_donors{HybridKind::POINTER, "chemfe_hbond_donor"},
    hydrogen_bond_acceptors{HybridKind::POINTER, "chemfe_hbond_accpt"},
    chiral_centers{HybridKind::POINTER, "chemfe_chirals"},
    rotatable_groups{HybridKind::POINTER, "chemfe_rotators"},
    rotatable_group_bounds{HybridKind::POINTER, "chemfe_rotator_bounds"},
    formal_charges{HybridKind::POINTER, "chemfe_formal_charges"},
    bond_orders{HybridKind::POINTER, "chemfe_bond_orders"},
    free_electrons{HybridKind::POINTER, "chemfe_free_e"},
    int_data{HybridKind::ARRAY, "chemfe_int"},
    double_data{HybridKind::ARRAY, "chemfe_double"},
    ag_pointer{ag_in}
{
  // Obtain abstracts from the topology
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  
  // Detect improper torsion terms and the atoms that make up their centers.  This is a first pass
  // at what is inhibited from rotating.
  std::vector<int> tmp_planar_centers = findPlanarAtoms(vk);
  
  // Detect rings by the bonding pattern
  std::vector<ullint> tmp_ring_inclusion(atom_count, 0LLU);
  std::vector<int> tmp_ring_atom_bounds(1, 0);
  std::vector<int> tmp_ring_atoms;
  traceTopologicalRings(nbk, cdk, &tmp_ring_inclusion, &tmp_ring_atoms, &tmp_ring_atom_bounds);
  ring_inclusion.resize(atom_count);
  ring_inclusion.putHost(tmp_ring_inclusion);
  
  // Allocate the double-precision real data that will result from Lewis structure determination
  const int padded_atom_count = roundUp(atom_count, warp_size_int);
  const int ndbl = 2 * padded_atom_count + roundUp(vk.nbond, warp_size_int);
  double_data.resize(ndbl);
  formal_charges.setPointer(&double_data, 0, padded_atom_count);
  free_electrons.setPointer(&double_data, padded_atom_count, padded_atom_count);
  bond_orders.setPointer(&double_data, 2 * padded_atom_count, roundUp(vk.nbond, warp_size_int));
  
  // Draw a Lewis structure and record the results.  The formal_charges, bond orders, and
  // free_electrons arrays can now be allocated.  Lewis structures will be drawn for each unique
  // molecule and copied otherwise, but that requires a list of all unique molecules.
  drawLewisStructures(vk, nbk, cdk);

  // CHECK
#if 0
  if (cdk.mol_limits[1] == 1185) {
    printf("Formal charges = [\n");
    for (int i = 0; i < cdk.mol_limits[1]; i++) {
      if (fabs(formal_charges.readHost(i)) < 1.0e-4) {
        continue;
      }
      const char4 atom_name = cdk.atom_names[i];
      const char4 res_name  = cdk.res_names[ag_pointer->getResidueIndex(i)];
      printf("  %c%c%c%c  %c%c%c%c  %9.4lf\n", atom_name.x, atom_name.y, atom_name.z, atom_name.w,
             res_name.x, res_name.y, res_name.z, res_name.w, formal_charges.readHost(i));
    }
    printf("];\n");
  }
#endif
  // END CHECK
  
  // Mark all aromatic atoms after dissecting numerous details about atoms in rings.
  std::vector<int> tmp_aromatic_group_bounds(1, 0);
  std::vector<int> tmp_aromatic_pi_electrons;
  std::vector<int> tmp_aromatic_groups;
  findAromaticGroups(cdk, vk, tmp_ring_atoms, tmp_ring_atom_bounds, &tmp_aromatic_group_bounds,
                     &tmp_aromatic_pi_electrons, &tmp_aromatic_groups);

  // Mark polar hydrogens and hydrogen bond donors and acceptors
  std::vector<int> tmp_polar_hydrogens;
  std::vector<int> tmp_hydrogen_bond_donors;
  std::vector<int> tmp_hydrogen_bond_acceptors;
  findHydrogenBondElements(nbk, cdk, &tmp_polar_hydrogens, &tmp_hydrogen_bond_donors,
                           &tmp_hydrogen_bond_acceptors);
  polar_hydrogen_count = tmp_polar_hydrogens.size();
  hbond_donor_count = tmp_hydrogen_bond_donors.size();
  hbond_acceptor_count = tmp_hydrogen_bond_acceptors.size();
  
  // Find chiral centers
  const std::vector<int> tmp_chiral_centers = findChiralCenters(nbk, vk, cdk, cfr);
  chiral_center_count = tmp_chiral_centers.size();

  // Find rotatable bonds
  std::vector<int> tmp_rotatable_groups;
  std::vector<int> tmp_rotatable_group_bounds;
  switch (map_group_in) {
  case MapRotatableGroups::YES:
    findRotatableBonds(vk, cdk, nbk, tmp_ring_atoms, tmp_ring_atom_bounds, &tmp_rotatable_groups,
                       &tmp_rotatable_group_bounds);
    rotatable_bond_count = static_cast<int>(tmp_rotatable_group_bounds.size()) - 1;
    break;
  case MapRotatableGroups::NO:
    break;
  }
  
  // Store the integer results
  const size_t nint = roundUp(tmp_planar_centers.size(), warp_size_zu) +
                      roundUp(tmp_ring_atom_bounds.size(), warp_size_zu) +
                      roundUp(tmp_ring_atoms.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_group_bounds.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_pi_electrons.size(), warp_size_zu) +
                      roundUp(tmp_aromatic_groups.size(), warp_size_zu) +
                      roundUp(static_cast<size_t>(chiral_center_count), warp_size_zu) +
                      roundUp(tmp_rotatable_groups.size(), warp_size_zu) +
                      roundUp(tmp_rotatable_group_bounds.size(), warp_size_zu) +
                      roundUp(tmp_polar_hydrogens.size(), warp_size_zu) +
                      roundUp(tmp_hydrogen_bond_donors.size(), warp_size_zu) +
                      roundUp(tmp_hydrogen_bond_acceptors.size(), warp_size_zu);
  int_data.resize(nint);
  size_t ic = planar_centers.putHost(&int_data, tmp_planar_centers, 0, warp_size_zu);
  ic = ring_atom_bounds.putHost(&int_data, tmp_ring_atom_bounds, ic, warp_size_zu);
  ic = ring_atoms.putHost(&int_data, tmp_ring_atoms, ic, warp_size_zu);
  ic = aromatic_group_bounds.putHost(&int_data, tmp_aromatic_group_bounds, ic, warp_size_zu);
  ic = aromatic_pi_electrons.putHost(&int_data, tmp_aromatic_pi_electrons, ic, warp_size_zu);
  ic = aromatic_groups.putHost(&int_data, tmp_aromatic_groups, ic, warp_size_zu);
  ic = polar_hydrogens.putHost(&int_data, tmp_polar_hydrogens, ic, warp_size_zu);
  ic = hydrogen_bond_donors.putHost(&int_data, tmp_hydrogen_bond_donors, ic, warp_size_zu);
  ic = hydrogen_bond_acceptors.putHost(&int_data, tmp_hydrogen_bond_acceptors, ic, warp_size_zu);
  ic = chiral_centers.putHost(&int_data, tmp_chiral_centers, ic, warp_size_zu);
  ic = rotatable_groups.putHost(&int_data, tmp_rotatable_groups, ic, warp_size_zu);
  ic = rotatable_group_bounds.putHost(&int_data, tmp_rotatable_group_bounds, ic, warp_size_zu);
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrame &cf,
                                   const MapRotatableGroups map_group_in,
                                   const double temperature_in) :
  ChemicalFeatures(ag_in, cf.data(), map_group_in, temperature_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const AtomGraph *ag_in, const PhaseSpace &ps,
                                   const MapRotatableGroups map_group_in,
                                   const double temperature_in) :
  ChemicalFeatures(ag_in, CoordinateFrameReader(ps), map_group_in, temperature_in)
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(const ChemicalFeatures &original) :
    atom_count{original.atom_count},
    planar_atom_count{original.planar_atom_count},
    ring_count{original.ring_count},
    fused_ring_count{original.fused_ring_count},
    twistable_ring_count{original.twistable_ring_count},
    conjugated_group_count{original.conjugated_group_count},
    aromatic_group_count{original.aromatic_group_count},
    polar_hydrogen_count{original.polar_hydrogen_count},
    hbond_donor_count{original.hbond_donor_count},
    hbond_acceptor_count{original.hbond_acceptor_count},
    chiral_center_count{original.chiral_center_count},
    rotatable_bond_count{original.rotatable_bond_count},
    double_bond_count{original.double_bond_count},
    triple_bond_count{original.triple_bond_count},
    max_ring_size{original.max_ring_size},
    temperature{original.temperature},
    rotating_groups_mapped{original.rotating_groups_mapped},
    planar_centers{original.planar_centers},
    ring_inclusion{original.ring_inclusion},
    ring_atom_bounds{original.ring_atom_bounds},
    ring_atoms{original.ring_atoms},
    aromatic_group_bounds{original.aromatic_group_bounds},
    aromatic_pi_electrons{original.aromatic_pi_electrons},
    aromatic_groups{original.aromatic_groups},
    polar_hydrogens{original.polar_hydrogens},
    hydrogen_bond_donors{original.hydrogen_bond_donors},
    hydrogen_bond_acceptors{original.hydrogen_bond_acceptors},
    chiral_centers{original.chiral_centers},
    rotatable_groups{original.rotatable_groups},
    rotatable_group_bounds{original.rotatable_group_bounds},
    formal_charges{original.formal_charges},
    bond_orders{original.bond_orders},
    free_electrons{original.free_electrons},
    int_data{original.int_data},
    double_data{original.double_data},
    ag_pointer{original.ag_pointer}
{
  repairPointers();
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures::ChemicalFeatures(ChemicalFeatures &&original) :
    atom_count{original.atom_count},
    planar_atom_count{original.planar_atom_count},
    ring_count{original.ring_count},
    fused_ring_count{original.fused_ring_count},
    twistable_ring_count{original.twistable_ring_count},
    conjugated_group_count{original.conjugated_group_count},
    aromatic_group_count{original.aromatic_group_count},
    polar_hydrogen_count{original.polar_hydrogen_count},
    hbond_donor_count{original.hbond_donor_count},
    hbond_acceptor_count{original.hbond_acceptor_count},
    chiral_center_count{original.chiral_center_count},
    rotatable_bond_count{original.rotatable_bond_count},
    double_bond_count{original.double_bond_count},
    triple_bond_count{original.triple_bond_count},
    max_ring_size{original.max_ring_size},
    temperature{original.temperature},
    rotating_groups_mapped{original.rotating_groups_mapped},
    planar_centers{std::move(original.planar_centers)},
    ring_inclusion{std::move(original.ring_inclusion)},
    ring_atom_bounds{std::move(original.ring_atom_bounds)},
    ring_atoms{std::move(original.ring_atoms)},
    aromatic_group_bounds{std::move(original.aromatic_group_bounds)},
    aromatic_pi_electrons{std::move(original.aromatic_pi_electrons)},
    aromatic_groups{std::move(original.aromatic_groups)},
    polar_hydrogens{std::move(original.polar_hydrogens)},
    hydrogen_bond_donors{std::move(original.hydrogen_bond_donors)},
    hydrogen_bond_acceptors{std::move(original.hydrogen_bond_acceptors)},
    chiral_centers{std::move(original.chiral_centers)},
    rotatable_groups{std::move(original.rotatable_groups)},
    rotatable_group_bounds{std::move(original.rotatable_group_bounds)},
    formal_charges{std::move(original.formal_charges)},
    bond_orders{std::move(original.bond_orders)},
    free_electrons{std::move(original.free_electrons)},
    int_data{std::move(original.int_data)},
    double_data{std::move(original.double_data)},
    ag_pointer{std::move(original.ag_pointer)}
{}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures& ChemicalFeatures::operator=(const ChemicalFeatures &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy elements of the other object
  atom_count = other.atom_count;
  planar_atom_count = other.planar_atom_count;
  ring_count = other.ring_count;
  fused_ring_count = other.fused_ring_count;
  twistable_ring_count = other.twistable_ring_count;
  conjugated_group_count = other.conjugated_group_count;
  aromatic_group_count = other.aromatic_group_count;
  polar_hydrogen_count = other.polar_hydrogen_count;
  hbond_donor_count = other.hbond_donor_count;
  hbond_acceptor_count = other.hbond_acceptor_count;
  chiral_center_count = other.chiral_center_count;
  rotatable_bond_count = other.rotatable_bond_count;
  double_bond_count = other.double_bond_count;
  triple_bond_count = other.triple_bond_count;
  max_ring_size = other.max_ring_size;
  temperature = other.temperature;
  rotating_groups_mapped = other.rotating_groups_mapped;
  planar_centers = other.planar_centers;
  ring_inclusion = other.ring_inclusion;
  ring_atom_bounds = other.ring_atom_bounds;
  ring_atoms = other.ring_atoms;
  aromatic_group_bounds = other.aromatic_group_bounds;
  aromatic_pi_electrons = other.aromatic_pi_electrons;
  aromatic_groups = other.aromatic_groups;
  polar_hydrogens = other.polar_hydrogens;
  hydrogen_bond_donors = other.hydrogen_bond_donors;
  hydrogen_bond_acceptors = other.hydrogen_bond_acceptors;
  chiral_centers = other.chiral_centers;
  rotatable_groups = other.rotatable_groups;
  rotatable_group_bounds = other.rotatable_group_bounds;
  formal_charges = other.formal_charges;
  bond_orders = other.bond_orders;
  free_electrons = other.free_electrons;
  int_data = other.int_data;
  double_data = other.double_data;
  ag_pointer = other.ag_pointer;
  repairPointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures& ChemicalFeatures::operator=(ChemicalFeatures &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy elements of the other object
  atom_count = other.atom_count;
  planar_atom_count = other.planar_atom_count;
  ring_count = other.ring_count;
  fused_ring_count = other.fused_ring_count;
  twistable_ring_count = other.twistable_ring_count;
  conjugated_group_count = other.conjugated_group_count;
  polar_hydrogen_count = other.polar_hydrogen_count;
  hbond_donor_count = other.hbond_donor_count;
  hbond_acceptor_count = other.hbond_acceptor_count;
  aromatic_group_count = other.aromatic_group_count;
  chiral_center_count = other.chiral_center_count;
  rotatable_bond_count = other.rotatable_bond_count;
  double_bond_count = other.double_bond_count;
  triple_bond_count = other.triple_bond_count;
  max_ring_size = other.max_ring_size;
  temperature = other.temperature;
  rotating_groups_mapped = other.rotating_groups_mapped;
  planar_centers = std::move(other.planar_centers);
  ring_inclusion = std::move(other.ring_inclusion);
  ring_atom_bounds = std::move(other.ring_atom_bounds);
  ring_atoms = std::move(other.ring_atoms);
  aromatic_group_bounds = std::move(other.aromatic_group_bounds);
  aromatic_pi_electrons = std::move(other.aromatic_pi_electrons);
  aromatic_groups = std::move(other.aromatic_groups);
  polar_hydrogens = std::move(other.polar_hydrogens);
  hydrogen_bond_donors = std::move(other.hydrogen_bond_donors);
  hydrogen_bond_acceptors = std::move(other.hydrogen_bond_acceptors);
  chiral_centers = std::move(other.chiral_centers);
  rotatable_groups = std::move(other.rotatable_groups);
  rotatable_group_bounds = std::move(other.rotatable_group_bounds);
  formal_charges = std::move(other.formal_charges);
  bond_orders = std::move(other.bond_orders);
  free_electrons = std::move(other.free_electrons);
  int_data = std::move(other.int_data);
  double_data = std::move(other.double_data);
  ag_pointer = std::move(other.ag_pointer);
  return *this;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::findPlanarAtoms(const ValenceKit<double> &vk) const {
  std::vector<int> result;
  for (int i = 0; i < vk.ndihe; i++) {
    const TorsionKind ikind = static_cast<TorsionKind>(vk.dihe_modifiers[i].w);
    if (ikind == TorsionKind::IMPROPER_NO_14 || ikind == TorsionKind::IMPROPER) {

      // The third atom of an improper is the central atom, about which a plane is enforced
      result.push_back(vk.dihe_k_atoms[i]);
    }
  }
  for (int i = 0; i < vk.ncimp; i++) {
    result.push_back(vk.cimp_k_atoms[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::traceTopologicalRings(const NonbondedKit<double> &nbk,
                                             const ChemicalDetailsKit &cdk,
                                             std::vector<ullint> *tmp_ring_inclusion,
                                             std::vector<int> *tmp_ring_atoms,
                                             std::vector<int> *tmp_ring_atom_bounds) {
  int max_branches = 0;
  for (int i = 0; i < atom_count; i++) {
    max_branches = std::max(nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i], max_branches);
  }
  const int max_molecule = ag_pointer->getLargestMoleculeSize();
  std::vector<int> all_branch_atoms(max_branches * max_molecule);
  std::vector<BondedNode> links(atom_count);
  std::vector<int> tree_positions(atom_count, -1);
  for (int i = 0; i < max_molecule; i++) {
    links[i].setBranchPointer(&all_branch_atoms, i, max_branches);
  }

  // Loop over all atoms, proceeding to explore all available bonds, until all atoms have been
  // either used as the start of chain / ring exploration or have been included in the exploration
  // initiated by some other atom.
  std::vector<bool> atom_touched(atom_count, false);
  for (int i = 0; i < atom_count; i++) {
    if (cdk.z_numbers[i] == 0) {
      atom_touched[i] = true;
    }
    if (atom_touched[i]) {
      continue;
    }
    
    // Initiate the chain
    atom_touched[i] = true;
    links[0].addToTree(-1, i, 0, nbk, cdk);
    tree_positions[i] = 0;
    int node_count = 1;
    int current_layer = 1;
    int layer_llim = 0;
    int layer_hlim = 1;
    while (layer_hlim > layer_llim) {
      const int next_layer_llim = node_count;
      for (int j = layer_llim; j < layer_hlim; j++) {
        const int j_atom = links[j].getAtom();
        const int j_branch_count = links[j].getBranchCount();
        for (int k = 0; k < j_branch_count; k++) {
          const int k_atom = links[j].getBranchAtom(k);

          // Check the status of the next atom: has it already been incorporated into this chain
          // or an earlier one?  If so, this signifies the completion of some loop.  Determine that
          // loop based on the histories of the current atom and the atom it touches.  Otherwise,
          // add the next atom to the chain.
          if (atom_touched[k_atom] && cdk.z_numbers[k_atom] != 0) {
            markRingAtoms(j_atom, k_atom, tree_positions, node_count, &links, tmp_ring_inclusion,
                          tmp_ring_atoms, tmp_ring_atom_bounds, cdk);
          }
          else {
            atom_touched[k_atom] = true;
            links[node_count].addToTree(j_atom, k_atom, current_layer, nbk, cdk);
            tree_positions[k_atom] = node_count;
            node_count++;
          }
        }
      }
      layer_llim = next_layer_llim;
      layer_hlim = node_count;
      current_layer++;
    }
    
    // Wipe any rings that were found, to allow future trees to evaluate new rings
    for (int j = 0; j < layer_hlim; j++) {
      links[j].wipeRingCompletion();
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::markRingAtoms(const int j_atom, const int k_atom,
                                     const std::vector<int> &tree_positions, const int node_count,
                                     std::vector<BondedNode> *links,
                                     std::vector<ullint> *tmp_ring_inclusion,
                                     std::vector<int> *tmp_ring_atoms,
                                     std::vector<int> *tmp_ring_atom_bounds,
                                     const ChemicalDetailsKit &cdk) {

  // Check that the ring has not already been marked
  BondedNode* linkdata = links->data();
  const int jpos = tree_positions[j_atom];
  const int kpos = tree_positions[k_atom];
  const int kbranch_in_j_atom = linkdata[jpos].findBranchAtom(k_atom);
  const int jbranch_in_k_atom = linkdata[kpos].findBranchAtom(j_atom);
  const bool j_marked = (linkdata[jpos].getRingCompletion(kbranch_in_j_atom));
  const bool k_marked = (linkdata[kpos].getRingCompletion(jbranch_in_k_atom));
  if (j_marked != k_marked) {
    rtErr("Inconsistency detected in ring formation for atoms " + std::to_string(j_atom + 1) +
          " and " + std::to_string(k_atom + 1) + ".", "ChemicalFeatures", "markRingAtoms");
  }
  if (j_marked) {
    return;
  }

  // Search backwards from each atom to find the point at which they reach a common atom
  std::vector<int> j_history, k_history;
  j_history.push_back(j_atom);
  k_history.push_back(k_atom);
  int j_length = 1;
  int k_length = 1;
  bool common_point_found = false;
  bool jgrow = true;
  bool kgrow = true;
  int jtrack = jpos;
  int ktrack = kpos;
  int k_pivot, j_pivot;
  while (common_point_found == false) {
    const int jprev = linkdata[jtrack].getPreviousAtom();
    if (jprev >= 0) {
      jtrack = tree_positions[jprev];
    }
    jgrow = (jgrow && (jprev >= 0));
    if (jgrow) {
      j_history.push_back(jprev);
      j_length++;
      for (int k = 0; k < k_length; k++) {
        if (jprev == k_history[k]) {
          common_point_found = true;
          j_pivot = j_length;
          k_pivot = k + 1;
          break;
        }
      }
    }
    if (common_point_found == false) {
      const int kprev = linkdata[ktrack].getPreviousAtom();
      if (kprev >= 0) {
        ktrack = tree_positions[kprev];
      }
      kgrow = (kgrow && (kprev >= 0));
      if (kgrow) {
        k_history.push_back(kprev);
        k_length++;
        for (int j = 0; j < j_length; j++) {
          if (kprev == j_history[j]) {
            common_point_found = true;
            j_pivot = j + 1;
            k_pivot = k_length;
            break;
          }
        }
      }
    }
    if (jgrow == false && kgrow == false) {
      common_point_found = true;
    }
  }
  
  // Set the appropriate bit on the various ring atoms.  Store the atoms of this ring in a
  // growing list.  This needs to be done with push_back as the ring computation is somewhat
  // involved and would be inefficient to repeat.
  const int ring_size = j_pivot + k_pivot - 1;
  if (ring_size < max_ring_size && ring_size >= 3) {
    
    // Step through the J history until the pivot point, then backwards through the K history
    // starting at the pivot point.  This will list all atoms in the ring.    
    int n_planar = 0;
    const ullint ring_mask = (0x1 << ring_size);
    ullint* ring_ptr = tmp_ring_inclusion->data();
    for (int j = 0; j < j_pivot; j++) {
      ring_ptr[j_history[j]] |= ring_mask;
      tmp_ring_atoms->push_back(j_history[j]);
    }
    for (int k = k_pivot - 2; k >= 0; k--) {
      ring_ptr[k_history[k]] |= ring_mask;
      tmp_ring_atoms->push_back(k_history[k]);
    }
    tmp_ring_atom_bounds->push_back(tmp_ring_atoms->size());
    ring_count += 1;
  }
  
  // Mark the ring completion in each node
  linkdata[jpos].setRingCompletion(kbranch_in_j_atom);
  linkdata[kpos].setRingCompletion(jbranch_in_k_atom);
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::drawLewisStructures(const ValenceKit<double> &vk,
                                           const NonbondedKit<double> &nbk,
                                           const ChemicalDetailsKit &cdk) {

  // Loop over all molecules and find those which are similar
  const int3 init_midx = {-1, 0, 0};
  std::vector<int3> mol_index(cdk.nmol, init_midx);
  int nmi = 0;
  for (int i = 0; i < cdk.nmol; i++) {    
    if (mol_index[i].x >= 0) {
      continue;
    }
    mol_index[i].x = nmi;
    mol_index[i].y = 0;
    mol_index[i].z = 0;
    const int i_llim = cdk.mol_limits[i];
    const int i_hlim = cdk.mol_limits[i + 1];
    const int ni_atom = i_hlim - i_llim;
    for (int j = i + 1; j < cdk.nmol; j++) {
      const int j_llim = cdk.mol_limits[j];
      const int j_hlim = cdk.mol_limits[j + 1];
      if (j_hlim - j_llim != ni_atom) {
        continue;
      }

      // Check the Z numbers
      bool znums_differ = false;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        znums_differ = (znums_differ || cdk.z_numbers[moli_atom] != cdk.z_numbers[molj_atom]);
      }
      if (znums_differ) {
        continue;
      }
      
      // Check the charges
      bool charges_differ = false;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        charges_differ = (charges_differ ||
                          fabs(nbk.charge[moli_atom] - nbk.charge[molj_atom]) > 1.0e-4);
      }
      if (charges_differ) {
        continue;
      }

      // Check the bonds
      bool bonds_similar = true;
      for (int k = 0; k < ni_atom; k++) {
        const int moli_atom = cdk.mol_contents[i_llim + k];
        const int molj_atom = cdk.mol_contents[j_llim + k];
        bonds_similar = (bonds_similar &&
                        (vk.bond_asgn_bounds[moli_atom + 1] - vk.bond_asgn_bounds[moli_atom] == 
                         vk.bond_asgn_bounds[molj_atom + 1] - vk.bond_asgn_bounds[molj_atom]));
        if (bonds_similar) {
          const int offset = vk.bond_asgn_atoms[vk.bond_asgn_bounds[molj_atom]] -
                             vk.bond_asgn_atoms[vk.bond_asgn_bounds[moli_atom]];
          int jm = vk.bond_asgn_bounds[molj_atom];
          for (int im = vk.bond_asgn_bounds[moli_atom]; im < vk.bond_asgn_bounds[moli_atom + 1];
               im++) {
              bonds_similar = (bonds_similar &&
                               vk.bond_asgn_atoms[jm] - vk.bond_asgn_atoms[im] == offset);
            jm++;
          }
        }
      }
      if (bonds_similar == false) {
        continue;
      }
      
      // The molecules have been found to be the same.
      mol_index[j].x = nmi;
      mol_index[j].y = j_llim - i_llim;
      mol_index[j].z = vk.bond_asgn_terms[cdk.mol_contents[j_llim]] -
                       vk.bond_asgn_terms[cdk.mol_contents[i_llim]];
    }

    // Increment the number of molecules
    nmi++;
  }

  // Proceed through the list of molecules, using the pre-determined molecular matches.
  std::vector<bool> ls_covered(cdk.nmol, false);
  for (int i = 0; i < cdk.nmol; i++) {
    if (ls_covered[i]) {
      continue;
    }
    IndigoTable idg_tab(ag_pointer, i, temperature);
    const std::vector<CombineIDp> fc = idg_tab.getGroundStateFormalCharges();
    const std::vector<CombineIDp> bo = idg_tab.getGroundStateBondOrders();
    const std::vector<CombineIDp> fe = idg_tab.getGroundStateFreeElectrons();
    const int match_idx = mol_index[i].x;
    for (int j = i; j < cdk.nmol; j++) {
      if (mol_index[j].x == match_idx) {
        for (size_t k = 0; k < fc.size(); k++) {
          formal_charges.putHost(fc[k].y, fc[k].x + mol_index[j].y);
        }
        for (size_t k = 0; k < bo.size(); k++) {
          bond_orders.putHost(bo[k].y, bo[k].x + mol_index[j].z);
        }
        for (size_t k = 0; k < fe.size(); k++) {
          free_electrons.putHost(fe[k].y, fe[k].x + mol_index[j].y);
        }
        ls_covered[j] = true;
      }      
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findAromaticGroups(const ChemicalDetailsKit &cdk,
                                          const ValenceKit<double> &vk,
                                          const std::vector<int> &tmp_ring_atoms,
                                          const std::vector<int> &tmp_ring_atom_bounds,
                                          std::vector<int> *tmp_aromatic_group_bounds,
                                          std::vector<int> *tmp_aromatic_pi_electrons,
                                          std::vector<int> *tmp_aromatic_groups) {

  // Identify sp2- and possible sp2-centers within the ring atoms
  const int nring_atoms = tmp_ring_atoms.size();
  std::vector<bool> ring_sp2_character(nring_atoms, false);
  const double* bo_ptr = bond_orders.data();
  for (int pos = 0; pos < vk.nbond; pos++) {

    // Use a number just shy of 1 1/4 to accommodate double bonds distributed over up to four
    // atoms--it's hard to imagine such a case, but perhaps a sulfur atom making a tetrahedral
    // center between two aromatic rings at 90-degree angles to one another is physcially possible?
    // In the guanidino head group of arginine the double bond is distributed over three atoms.
    if (bo_ptr[pos] > 1.24999) {
      const int bi_atom = vk.bond_i_atoms[pos];
      const int bj_atom = vk.bond_j_atoms[pos];
      for (int i = 0; i < nring_atoms; i++) {
        ring_sp2_character[i] = (ring_sp2_character[i] || tmp_ring_atoms[i] == bi_atom ||
                                 tmp_ring_atoms[i] == bj_atom);
      }
    }
  }
  for (int i = 0; i < nring_atoms; i++) {
    const int atom_idx = tmp_ring_atoms[i];
    switch (cdk.z_numbers[atom_idx]) {
    case 1:
    case 9:
    case 17:
    case 35:

      // Hydrogen and each of the halogens are assumed to not participate in rings
      break;
    case 6:
    case 7:
    case 8:
    case 15:
    case 16:

      // Allow for a tiny amount of imprecision in whatever floating point arithmetic led to
      // the calculated free electron content.  Two free electrons can switch into the pi system
      // if this atom makes two or more bonds to other atoms, which will be the case for any atom
      // in a ring.  Even a nitrogen making two bonds to other ring atoms, with a hydrogen hanging
      // off of it (three bonds to other atoms in all) can donate a pair of electrons to the pi
      // system in the ring.  Allow one free electron to contribute, as a protonated histidine
      // ring will distribute a single positive formal charge and bonds of order 1.5 across both
      // nitrogens, leaving a lone free electron, on average, for each of them.  The individual
      // resonance states have either a full positive charge or a lone pair on either nitrogen,
      // and in either case the lone pair can contribute to an aromatic ring with a net positive
      // charge and six pi electrons.  Making inferences from the averages in this manner may
      // create some way for errors to creep in, but roll with it for now.
      ring_sp2_character[i] = (ring_sp2_character[i] ||
                               free_electrons.readHost(atom_idx) >= 0.999);
      break;
    default:
      break;
    }
  }
  
  // Determine fused rings, in preparation for aromaticity checks
  std::vector<int> ring_participation_bounds(atom_count + 1, 0);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i+1]; j++) {
      ring_participation_bounds[tmp_ring_atoms[j]] += 1;
    }
  }
  prefixSumInPlace<int>(&ring_participation_bounds, PrefixSumType::EXCLUSIVE, "ChemicalFeatures");
  std::vector<int> ring_participation(ring_participation_bounds[atom_count]);
  for (int i = 0; i < ring_count; i++) {
    for (int j = tmp_ring_atom_bounds[i]; j < tmp_ring_atom_bounds[i+1]; j++) {
      ring_participation[ring_participation_bounds[tmp_ring_atoms[j]]] = i;
      ring_participation_bounds[tmp_ring_atoms[j]] += 1;
    }
  }
  for (int i = atom_count; i > 0; i--) {
    ring_participation_bounds[i] = ring_participation_bounds[i - 1];
  }
  ring_participation_bounds[0] = 0;

  // Assess fused ring systems
  std::vector<int> fused_rings(ring_count);
  std::vector<bool> ring_covered(ring_count, false);
  std::vector<int> fused_ring_bounds(1, 0);
  int nr = 0;
  for (int i = 0; i < ring_count; i++) {
    if (ring_covered[i]) {
      continue;
    }
    ring_covered[i] = true;
    std::vector<int> new_rings(1, i);
    int current_ring_count;
    int previous_ring_count = 0;
    do {

      // Snapshot the current number of rings
      current_ring_count = new_rings.size();

      // Loop over all current rings and add more
      for (int rpos = previous_ring_count; rpos < current_ring_count; rpos++) {
        for (int j = tmp_ring_atom_bounds[new_rings[rpos]];
             j < tmp_ring_atom_bounds[new_rings[rpos] + 1]; j++) {
          const int tring_atom = tmp_ring_atoms[j];

          // Check for atoms that are also part of an additional ring
          for (int k = ring_participation_bounds[tring_atom] + 1;
               k < ring_participation_bounds[tring_atom + 1]; k++) {
            bool extra_ring_accounted = ring_covered[ring_participation[k]];
            const int total_ring_count = new_rings.size();
            for (int m = 0; m < total_ring_count; m++) {
              extra_ring_accounted = (extra_ring_accounted ||
                                      new_rings[m] == ring_participation[k]);
            }
            if (extra_ring_accounted == false) {
              new_rings.push_back(ring_participation[k]);
              ring_covered[ring_participation[k]] = true;
            }
          }
        }
      }

      // Update the number of rings previously known, in preparation for another pass
      previous_ring_count = current_ring_count;
    } while (new_rings.size() > current_ring_count);

    // Systems with only one ring are not fused rings
    if (new_rings.size() == 1) {
      continue;
    }

    // Contribute all contiguous rings to the fused ring list
    for (size_t j = 0; j < new_rings.size(); j++) {
      ring_covered[new_rings[j]] = true;
      fused_rings[nr] = new_rings[j];
      nr++;
    }
    fused_ring_bounds.push_back(nr);
  }

  // Record the number of fused ring systems, for convenience
  fused_ring_count = fused_ring_bounds.size() - 1;
  
  // Determine fused ring atom and bond content
  std::vector<int> fused_ring_atom_counts(fused_ring_count, 0);
  std::vector<int> fused_ring_bond_counts(fused_ring_count, 0);
  std::vector<int> fused_ring_atoms;
  std::vector<int> fused_ring_atom_bounds(1, 0);
  for (int i = 0; i < fused_ring_count; i++) {

    // Make a list of all atoms
    std::vector<int> atom_list;
    for (int j = fused_ring_bounds[i]; j < fused_ring_bounds[i + 1]; j++) {
      const int jring = fused_rings[j];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        bool found = false;
        const int atom_list_length = atom_list.size();
        const int tkr_atom = tmp_ring_atoms[k];
        for (int m = 0; m < atom_list_length; m++) {
          found = (found || atom_list[m] == tkr_atom);
        }
        if (found == false) {
          atom_list.push_back(tkr_atom);
        }
      }
    }
    fused_ring_atom_counts[i] = atom_list.size();
    fused_ring_atoms.insert(fused_ring_atoms.end(), atom_list.begin(), atom_list.end());
    fused_ring_atom_bounds.push_back(fused_ring_atoms.size());

    // Referencing each atom, make a list of all bonds
    for (int j = 0; j < fused_ring_atom_counts[i]; j++) {
      const int atom_idx = atom_list[j];
      for (int k = vk.bond_asgn_bounds[atom_idx]; k < vk.bond_asgn_bounds[atom_idx + 1]; k++) {
        bool found = false;
        const int partner_atom = vk.bond_asgn_atoms[k];
        for (int m = 0; m < fused_ring_atom_counts[i]; m++) {
          found = (found || atom_list[m] == partner_atom);
        }
        if (found) {
          fused_ring_bond_counts[i] += 1;
        }
      }
    }
  }
  
  // Mark rings as aromatic or not, based on the Lewis structure, rings, and fused_rings
  std::vector<int> pi_electron_count(ring_count);
  std::vector<int> lp_electron_count(ring_count);
  for (int i = 0; i < ring_count; i++) {
    const int rb_llim = tmp_ring_atom_bounds[i];
    const int rb_hlim = tmp_ring_atom_bounds[i + 1];
    double pi_electrons = 0.0;
    double lp_electrons = 0.0;
    for (int j = rb_llim; j < rb_hlim; j++) {

      // Compute the degree of pi bonding for bonds relating to this atom lying along the ring.
      const int atom_idx = tmp_ring_atoms[j];
      for (int k = vk.bond_asgn_bounds[atom_idx]; k < vk.bond_asgn_bounds[atom_idx + 1]; k++) {
        const int partner_idx = vk.bond_asgn_atoms[k];
        bool partner_in_ring = false;
        for (int m = rb_llim; m < rb_hlim; m++) {
          partner_in_ring = (partner_in_ring || tmp_ring_atoms[m] == partner_idx);
        }
        if (partner_in_ring) {
          const int bond_idx = vk.bond_asgn_terms[k];
          pi_electrons += 2.0 * (bond_orders.readHost(bond_idx) - 1.0);
        }
      }

      // Add the free electron content
      const double nfree_e = free_electrons.readHost(atom_idx);
      if (nfree_e >= 2.0) {
        lp_electrons += 2.0;
      }
      else if (nfree_e >= 0.0) {
        lp_electrons += nfree_e;
      }
    }
    pi_electron_count[i] = round(pi_electrons);
    lp_electron_count[i] = round(lp_electrons);
  }
  for (int i = 0; i < fused_ring_count; i++) {

    // Immediately eliminate fused rings that have less than five atoms, or that have one or more
    // atoms which cannot take on sp2 character.
    if (fused_ring_atom_counts[i] < 5) {
      continue;
    }
    const int frb_llim = fused_ring_bounds[i];
    const int frb_hlim = fused_ring_bounds[i + 1];
    bool ring_qualifies = true;
    for (int j = frb_llim; j < frb_hlim; j++) {
      const int jring = fused_rings[j];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        ring_qualifies = (ring_qualifies && ring_sp2_character[k]);
      }
    }
    if (ring_qualifies == false) {
      continue;
    }

    // The numbers of pi and lone pair electrons are the sum of each of the constituent rings,
    // without double-counting bonds or atoms shared by both rings.
    int pi_electrons = 0;
    int lp_electrons = 0;
    std::vector<int> ring_occupancy(fused_ring_atom_counts[i], 0);
    for (int j = frb_llim; j < frb_hlim; j++) {
      const int jring = fused_rings[j];      
      pi_electrons += pi_electron_count[jring];
      lp_electrons += lp_electron_count[jring];
      for (int k = tmp_ring_atom_bounds[jring]; k < tmp_ring_atom_bounds[jring + 1]; k++) {
        const int m_start = fused_ring_atom_bounds[i];
        for (int m = m_start; m < fused_ring_atom_bounds[i + 1]; m++) {
          ring_occupancy[m - m_start] += (fused_ring_atoms[m] == tmp_ring_atoms[k]);
        }
      }
    }
    for (int j = 0; j < fused_ring_atom_counts[i]; j++) {
      if (ring_occupancy[j] == 1) {
        continue;
      }
      const int n_extra = ring_occupancy[j] - 1;
      const int atomj_idx = fused_ring_atoms[fused_ring_atom_bounds[i] + j];
      const int nfree_e = static_cast<int>(round(free_electrons.readHost(atomj_idx)));
      if (nfree_e >= 2) {
        lp_electrons -= 2 * n_extra;
      }
      else if (nfree_e >= 0) {
        lp_electrons -= nfree_e * n_extra;
      }

      // Seek out bonds between atoms that are included in similar numbers of rings within
      // this fused system.
      for (int k = j + 1; k < fused_ring_atom_counts[i]; k++) {
        if (ring_occupancy[k] == ring_occupancy[j]) {
          const int atomk_idx = fused_ring_atoms[fused_ring_atom_bounds[i] + k];
          for (int m = vk.bond_asgn_bounds[atomk_idx]; m < vk.bond_asgn_bounds[atomk_idx + 1];
               m++) {
            if (vk.bond_asgn_atoms[m] == atomj_idx) {
              const double tbo = bond_orders.readHost(vk.bond_asgn_terms[m]);
              pi_electrons -= static_cast<int>(round(2.0 * (tbo - 1.0))) * n_extra;
            }
          }
          for (int m = vk.bond_asgn_bounds[atomj_idx]; m < vk.bond_asgn_bounds[atomj_idx + 1];
               m++) {
              const double tbo = bond_orders.readHost(vk.bond_asgn_terms[m]);
              pi_electrons -= static_cast<int>(round(2.0 * (tbo - 1.0))) * n_extra;
          }
        }
      }
    }

    // Is the ring itself aromatic?  Apply Huckel's rule to test whether there are 4n + 2 pi
    // electrons.
    if (((pi_electrons + lp_electrons) & 0x1) == 0 && (pi_electrons + lp_electrons) / 4 >= 1) {
      for (int j = fused_ring_atom_bounds[i]; j < fused_ring_atom_bounds[i + 1]; j++) {
        tmp_aromatic_groups->push_back(fused_ring_atoms[j]);
      }
      tmp_aromatic_group_bounds->push_back(tmp_aromatic_groups->size());
      tmp_aromatic_pi_electrons->push_back(pi_electrons + lp_electrons);
    }
  }
  for (int i = 0; i < ring_count; i++) {

    // Immediately eliminate rings that have less than five atoms, or that have one or more
    // atoms which cannot take on sp2 character.
    const int rb_llim = tmp_ring_atom_bounds[i];
    const int rb_hlim = tmp_ring_atom_bounds[i + 1];
    const int nring_atom = rb_hlim - rb_llim;
    if (nring_atom < 5) {
      continue;
    }
    bool ring_qualifies = true;
    for (int j = rb_llim; j < rb_hlim; j++) {
      ring_qualifies = (ring_qualifies && ring_sp2_character[j]);
    }
    if (ring_qualifies == false) {
      continue;
    }    
    const int pi_electrons = pi_electron_count[i];
    const int lp_electrons = lp_electron_count[i];
    
    // Is the ring itself aromatic?  Apply Huckel's rule to test whether there are 4n + 2 pi
    // electrons.
    if (((pi_electrons + lp_electrons) & 0x1) == 0 && (pi_electrons + lp_electrons) / 4 >= 1) {
      for (int j = rb_llim; j < rb_hlim; j++) {
        tmp_aromatic_groups->push_back(tmp_ring_atoms[j]);
      }
      tmp_aromatic_group_bounds->push_back(tmp_aromatic_groups->size());
    }
  }

  // Record the total number of groups, for convenience
  aromatic_group_count = tmp_aromatic_group_bounds->size() - 1;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::findChiralCenters(const NonbondedKit<double> &nbk,
                                                     const ValenceKit<double> &vk,
                                                     const ChemicalDetailsKit &cdk,
                                                     const CoordinateFrameReader &cfr) const {
  
  // Prepare to construct tree structures for each molecule, similar to what was done above.
  // However, the trees will not trace rings this time.  Rather, they will continue along the path
  // radiating outwards from any given atom which has four unique substituents, trying to determine
  // whether one is superior to another in terms of atomic numbers and bond orders.
  int max_branches = 0;
  for (int i = 0; i < atom_count; i++) {
    max_branches = std::max(nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i], max_branches);
  }
  const int max_molecule = ag_pointer->getLargestMoleculeSize();
  std::vector<std::vector<int>> all_branch_atoms(4, std::vector<int>(max_branches * max_molecule));
  std::vector<std::vector<BondedNode>> links(4, std::vector<BondedNode>(atom_count));
  std::vector<std::vector<int>> tree_positions(4, std::vector<int>(atom_count, -1));
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < max_molecule; j++) {
      links[i][j].setBranchPointer(&all_branch_atoms[i], j, max_branches);
    }
  }

  // Loop over all atoms, proceeding to explore all available bonds, until all atoms have been
  // either used as the start of chain / ring exploration or have been included in the exploration
  // initiated by some other atom.
  std::vector<std::vector<bool>> atom_touched(4, std::vector<bool>(atom_count, false));
  std::vector<int> result;
  for (int i = 0; i < atom_count; i++) {

    // Test whether this atom is a chiral candidate
    if (nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i] < 4) {
      continue;
    }
    bool candidate = true;
    int n_hydrogen = 0;
    int n_vs = 0;
    for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
      const int ex_znum = cdk.z_numbers[nbk.nb12x[j]];
      candidate = (candidate && ex_znum > 0);
      n_hydrogen += (ex_znum == 1);
      n_vs += (ex_znum == 0);
      candidate = (candidate && n_hydrogen < 2);
    }
    if (candidate == false || nbk.nb12_bounds[i + 1] - nbk.nb12_bounds[i] - n_vs != 4) {
      continue;
    }

    // Initiate the chains
    int n_real = 0;
    std::vector<int> layer_llim(4, 0);
    std::vector<int> layer_hlim(4, 1);
    std::vector<int> node_count(4, 1);
    std::vector<int> branch_ranks(4, 0);
    std::vector<int> branch_scores(4, 0);
    std::vector<int> touch_min(4);
    std::vector<int> touch_max(4);
    std::vector<bool> score_branch(4);
    int current_layer = 1;
    for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
      if (cdk.z_numbers[nbk.nb12x[j]] == 0) {
        continue;
      }
      links[n_real][0].addToTree(i, nbk.nb12x[j], 0, nbk, cdk);
      links[n_real][0].addBondOrder(vk, bond_orders);
      atom_touched[n_real][i] = true;
      atom_touched[n_real][nbk.nb12x[j]] = true;
      touch_min[n_real] = std::min(i, nbk.nb12x[j]);
      touch_max[n_real] = std::max(i, nbk.nb12x[j]);
      n_real++;
    }

    // Once a branch beats another, that dominance relationship must be maintained.  Keep a matrix
    // of the dominance relationships (i, j) : 0 = tie, +1 = i beats j, -1 = j beats i.  The matrix
    // will be antisymmetric.
    std::vector<int> chiral_dominance(16, 0);
    std::vector<int> parallel_growth(16, 0);
    bool advance_chains = scoreChiralBranches(links, layer_llim, layer_hlim, cdk,
                                              &chiral_dominance, &parallel_growth);
    while (advance_chains) {
      
      // Add more nodes
      for (int j = 0; j < 4; j++) {
        for (int k = layer_llim[j]; k < layer_hlim[j]; k++) {
          const int k_atom = links[j][k].getAtom();
          const int k_branch_count = links[j][k].getBranchCount();
          for (int m = 0; m < k_branch_count; m++) {
            const int m_atom = links[j][k].getBranchAtom(m);
            if (! atom_touched[j][m_atom]) {
              atom_touched[j][m_atom] = true;
              touch_min[j] = std::min(touch_min[j], m_atom);
              touch_max[j] = std::max(touch_max[j], m_atom);
              links[j][node_count[j]].addToTree(k_atom, m_atom, current_layer, nbk, cdk);
              links[j][node_count[j]].addBondOrder(vk, bond_orders);
              tree_positions[j][m_atom] = node_count[j];
              node_count[j] += 1;
            }
          }
        }
        layer_llim[j] = layer_hlim[j];
        layer_hlim[j] = node_count[j];
      }
      current_layer++;
      advance_chains = scoreChiralBranches(links, layer_llim, layer_hlim, cdk, &chiral_dominance,
                                           &parallel_growth);
    }

    // Check the dominance matrix: are the off-diagonal elements all nonzero?
    if (chiral_dominance[ 1] == 0 || chiral_dominance[ 2] == 0 || chiral_dominance[ 3] == 0 ||
        chiral_dominance[ 6] == 0 || chiral_dominance[ 7] == 0 || chiral_dominance[11] == 0) {
      continue;      
    }

    // Give the next chirality check a clean slate.  Use the touch_min and touch_max bounds
    // arrays on the range of atoms that became part of each branch's tree to keep us in O(N)
    // territory when doing things like the chirality of protein CA atoms.
    for (int j = 0; j < 4; j++) {
      for (int k = touch_min[j]; k < touch_max[j]; k++) {
        atom_touched[j][k] = false;
      }
    }

    // Check the dominance matrix for branch priorities.  Higher score = higher priority.
    std::vector<int2> priority(4, {0, 0});
    for (int j = 0; j < 4; j++) {
      priority[j].x += (chiral_dominance[j     ] == 1) + (chiral_dominance[j +  4] == 1) +
                       (chiral_dominance[j +  8] == 1) + (chiral_dominance[j + 12] == 1);
      priority[j].y = j;
    }
    std::sort(priority.begin(), priority.end(), [](int2 a, int2 b) { return a.x < b.x; });

    // Make the coordinate origin the atom center.  Align the atoms with respect to the axis
    // extending between the center and the root of the lowest priority branch.
    const int chiral_direction = getChiralOrientation(cfr, i, links[priority[0].y][0].getAtom(),
                                                      links[priority[3].y][0].getAtom(),
                                                      links[priority[2].y][0].getAtom(),
                                                      links[priority[1].y][0].getAtom());
    
    // Combine the chirality determination with the atom index (offset by one to prevent the
    // abiguity of zero and -zero), then contribute the result
    result.push_back(chiral_direction * (i + 1));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findRotatableBonds(const ValenceKit<double> &vk,
                                          const ChemicalDetailsKit &cdk,
                                          const NonbondedKit<double> &nbk,
                                          const std::vector<int> &ring_atoms,
                                          const std::vector<int> &ring_atom_bounds,
                                          std::vector<int> *tmp_rotatable_groups,
                                          std::vector<int> *tmp_rotatable_group_bounds) {

  // Prepare a table of atoms that are part of rings
  std::vector<bool> bond_in_ring(vk.nbond, false);
  for (int i = 0; i < ring_count; i++) {
    for (int j = ring_atom_bounds[i]; j < ring_atom_bounds[i + 1]; j++) {
      const int jatom = ring_atoms[j];
      for (int k = vk.bond_asgn_bounds[jatom]; k < vk.bond_asgn_bounds[jatom + 1]; k++) {
        const int katom = vk.bond_asgn_atoms[k];
        bool tb_in_ring = false;
        for (int m = ring_atom_bounds[i]; m < ring_atom_bounds[i + 1]; m++) {
          tb_in_ring = (tb_in_ring || katom == ring_atoms[m]);
        }
        bond_in_ring[vk.bond_asgn_terms[k]] = tb_in_ring;
      }
    }
  }
  
  // Scan over all bonds and accumulate a result
  Approx near_one(1.0, 0.21);
  std::vector<int2> rotators;
  std::vector<std::vector<int>> moving_lists;
  for (int pos = 0; pos < vk.nbond; pos++) {

    // Omit bonds within rings.  Those will be handled separately.  Otherwise, allow a liberal
    // definition of a single bond--a small amount of double-bond character might be permissible.
    if (bond_in_ring[pos] || bond_orders.readHost(pos) != near_one) {
      continue;
    }
    
    // Ensure that both ends have more branching from them, and that the branches are worth
    // rotating (more than just hydrogen atoms branching from them)
    const int atom_i = vk.bond_i_atoms[pos];
    const int atom_j = vk.bond_j_atoms[pos];
    int nbranch_i = 0;
    for (int i = nbk.nb12_bounds[atom_i]; i < nbk.nb12_bounds[atom_i + 1]; i++) {
      nbranch_i += (nbk.nb12x[i] != atom_j && cdk.z_numbers[nbk.nb12x[i]] > 1);
    }
    int nbranch_j = 0;
    for (int i = nbk.nb12_bounds[atom_j]; i < nbk.nb12_bounds[atom_j + 1]; i++) {
      nbranch_j += (nbk.nb12x[i] != atom_i && cdk.z_numbers[nbk.nb12x[i]] > 1);
    }    
    if (nbranch_i < 1 || nbranch_j < 1) {
      continue;
    }

    // Test the number of rotatable bonds with atom_i as the root and atom_j as the pivot (the
    // pivot atom is closest to the atoms that will move), then vice-versa.  List the atoms such
    // that the smallest number of atoms will move as a consequence of rotation about the bond.
    // Store the rotating atoms.
    const std::vector<int> ifirst = selectRotatingAtoms(ag_pointer, atom_i, atom_j);
    const std::vector<int> jfirst = selectRotatingAtoms(ag_pointer, atom_j, atom_i);
    if (ifirst.size() >= jfirst.size()) {
      rotators.push_back({atom_i, atom_j});
      moving_lists.push_back(jfirst);
    }
    else {
      rotators.push_back({atom_j, atom_i});
      moving_lists.push_back(ifirst);
    }
  }
  const size_t ngroup = rotators.size();
  tmp_rotatable_group_bounds->resize(ngroup + 1);
  int* gbounds_ptr = tmp_rotatable_group_bounds->data();
  int acc_size = 0;
  for (size_t i = 0; i < ngroup; i++) {
    gbounds_ptr[i] = acc_size;
    acc_size += 2 + moving_lists[i].size();
  }
  gbounds_ptr[ngroup] = acc_size;
  tmp_rotatable_groups->resize(acc_size);
  int* grp_ptr = tmp_rotatable_groups->data();
  size_t k = 0;
  for (size_t i = 0; i < ngroup; i++) {
    grp_ptr[k] = rotators[i].x;
    k++;
    grp_ptr[k] = rotators[i].y;
    k++;
    const size_t gsize = moving_lists[i].size();
    for (size_t j = 0; j < gsize; j++) {
      grp_ptr[k] = moving_lists[i][j];
      k++;
    }
  }

  // Signal that the rotatable groups have been mapped
  rotating_groups_mapped = true;
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::findHydrogenBondElements(const NonbondedKit<double> &nbk,
                                                const ChemicalDetailsKit &cdk,
                                                std::vector<int> *tmp_polar_h,
                                                std::vector<int> *tmp_hb_don,
                                                std::vector<int> *tmp_hb_acc) {
  std::vector<bool> polarh_covered(cdk.natom, false);
  std::vector<bool> donor_covered(cdk.natom, false);
  for (int i = 0; i < cdk.natom; i++) {
    if (cdk.z_numbers[i] == 1) {
      for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
        const int donor_atom = nbk.nb12x[j];
        const int donor_z = cdk.z_numbers[donor_atom];
        if (donor_z == 7 || donor_z == 8 || donor_z == 15 || donor_z == 16) {
          if (polarh_covered[i] == false) {
            tmp_polar_h->push_back(i);
            polarh_covered[i] = true;
          }
          if (donor_covered[donor_atom] == false) {
            tmp_hb_don->push_back(donor_atom);
            donor_covered[donor_atom] = true;
          }
        }
      }
    }
    else if ((cdk.z_numbers[i] ==  7 || cdk.z_numbers[i] ==  8 || cdk.z_numbers[i] == 15 |
              cdk.z_numbers[i] == 16) && free_electrons.readHost(i) > 1.25) {

      // Nitrogen, oxygen, phosphorus, and sulfur atoms with at least some lone pair occupancy
      // are considered hydrogen bond acceptors.
      tmp_hb_acc->push_back(i);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::repairPointers() {
  formal_charges.swapTarget(&double_data);
  free_electrons.swapTarget(&double_data);
  bond_orders.swapTarget(&double_data);
  planar_centers.swapTarget(&int_data);
  ring_atom_bounds.swapTarget(&int_data);
  ring_atoms.swapTarget(&int_data);
  aromatic_group_bounds.swapTarget(&int_data);
  aromatic_pi_electrons.swapTarget(&int_data);
  aromatic_groups.swapTarget(&int_data);
  chiral_centers.swapTarget(&int_data);
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void ChemicalFeatures::upload() {
  ring_inclusion.upload();
  int_data.upload();
  double_data.upload();
}
#endif

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getPlanarAtomCount() const {
  return planar_atom_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getRingCount() const {
  return ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getFusedRingCount() const {
  return fused_ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getMutableRingCount() const {
  return twistable_ring_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getAromaticGroupCount() const {
  return aromatic_group_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getPolarHydrogenCount() const {
  return polar_hydrogen_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getHydrogenBondDonorCount() const {
  return hbond_donor_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getHydrogenBondAcceptorCount() const {
  return hbond_acceptor_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getChiralCenterCount() const {
  return chiral_center_count;
}

//-------------------------------------------------------------------------------------------------
int ChemicalFeatures::getRotatableBondCount() const {
  return rotatable_bond_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getRingMask(const int min_ring_size,
                                                const int max_ring_size) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* rb_ptr = ring_atom_bounds.data();
  const int* ra_ptr = ring_atoms.data();
  for (int i = 0; i < ring_count; i++) {
    const int rsize = rb_ptr[i + 1] - rb_ptr[i];
    if (rsize >= min_ring_size && rsize <= max_ring_size) {
      for (int j = rb_ptr[i]; j < rb_ptr[i + 1]; j++) {
        const int mask_idx = ra_ptr[j] / n_bits;
        const int bit_idx  = ra_ptr[j] - (mask_idx * n_bits);
        result[mask_idx] |= (0x1 << bit_idx);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getAromaticMask(const int min_pi_electrons,
                                                    const int max_pi_electrons) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* amb_ptr = aromatic_group_bounds.data();
  const int* amp_ptr = aromatic_pi_electrons.data();
  const int* ama_ptr = aromatic_groups.data();
  for (int i = 0; i < ring_count; i++) {
    if (amp_ptr[i] >= min_pi_electrons && amp_ptr[i] <= max_pi_electrons) {
      for (int j = amb_ptr[i]; j < amb_ptr[i + 1]; j++) {
        const int mask_idx = ama_ptr[j] / n_bits;
        const int bit_idx  = ama_ptr[j] - (mask_idx * n_bits);
        result[mask_idx] |= (0x1 << bit_idx);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getPolarHydrogenList() const {
  return polar_hydrogens.readHost(0, polar_hydrogen_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getHydrogenBondDonorList() const {
  return hydrogen_bond_donors.readHost(0, hbond_donor_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::getHydrogenBondAcceptorList() const {
  return hydrogen_bond_acceptors.readHost(0, hbond_acceptor_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getPolarHydrogenMask() const {
  return numberSeriesToBitMask(polar_hydrogens, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getHydrogenBondDonorMask() const {
  return numberSeriesToBitMask(hydrogen_bond_donors, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getHydrogenBondAcceptorMask() const {
  return numberSeriesToBitMask(hydrogen_bond_acceptors, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ChemicalFeatures::listChiralCenters(const ChiralOrientation direction) const {
  std::vector<int> result;
  const int* chi_ptr = chiral_centers.data();
  int ncen;
  switch (direction) {
  case ChiralOrientation::RECTUS:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] < 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] < 0) {
        result[ncen] = 1 - chi_ptr[i];
        ncen++;
      }
    }
    break;
  case ChiralOrientation::SINISTER:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] > 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] > 0) {
        result[ncen] = chi_ptr[i] - 1;
        ncen++;
      }
    }
    break;
  case ChiralOrientation::NONE:
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      ncen += (chi_ptr[i] != 0);
    }
    result.resize(ncen);
    ncen = 0;
    for (int i = 0; i < chiral_center_count; i++) {
      if (chi_ptr[i] != 0) {
        result[ncen] = abs(chi_ptr[i]) - 1;
        ncen++;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> ChemicalFeatures::getChiralityMask(const ChiralOrientation direction) const {
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (atom_count + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);
  
  // Loop over all known rings in the system and select those meeting the size criterion
  const int* chi_ptr = chiral_centers.data();
  for (int i = 0; i < chiral_center_count; i++) {
    const int atom_index = abs(chi_ptr[i]) - 1;
    const int mask_idx = atom_index / n_bits;
    const int bit_idx  = atom_index - (mask_idx * n_bits);
    switch (direction) {
    case ChiralOrientation::RECTUS:
      if (chi_ptr[i] < 0) {
        result[mask_idx] |= (0x1 << bit_idx);
      }
      break;
    case ChiralOrientation::SINISTER:
      if (chi_ptr[i] > 0) {
        result[mask_idx] |= (0x1 << bit_idx);
      }
      break;
    case ChiralOrientation::NONE:      

      // Here, this indicates no preference for chirality.  The list still contains a distinct
      // subset of the aotms displaying chirality, not atoms that are achiral.
      result[mask_idx] |= (0x1 << bit_idx);
      break;      
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getFormalCharges() const {
  return formal_charges.readHost(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ChemicalFeatures::getBondOrders() const {
  return bond_orders.readHost(0, ag_pointer->getBondTermCount());
}

//-------------------------------------------------------------------------------------------------
std::vector<RotatorGroup> ChemicalFeatures::getRotatableBondGroups() const {
  std::vector<RotatorGroup> result(rotatable_bond_count);
  const int* rg_ptr = rotatable_groups.data();
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    result[i].root_atom  = rg_ptr[llim];
    result[i].pivot_atom = rg_ptr[hlim];
    result[i].rotatable_atoms.resize(hlim - llim - 2);
    int k = 0;
    for (int j = llim + 2; j < hlim; j++) {
      result[i].rotatable_atoms[k] = rg_ptr[j];
      k++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<RotatorGroup>
ChemicalFeatures::getRotatableBondGroups(const int cutoff, const int mol_index) const {

  // Collect all rotatable groups on a particular molecule larger than a stated cutoff size.
  // Order the results in descending order of the number of atoms that rotate as a consequence
  // of twisting the bond.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  if (mol_index < 0 || mol_index >= cdk.nmol) {
    rtErr("Molecule index " + std::to_string(mol_index) + " is invalid for a system with " +
          std::to_string(cdk.nmol) + " molecules.", "ChemicalFeatures", "getRotatableBondGroups");
  }
  int nrg = 0;
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    nrg += (cdk.mol_home[llim + 1] == mol_index && hlim - llim - 2 >= cutoff);
  }
  std::vector<RotatorGroup> result;
  result.reserve(nrg);
  const int* rg_ptr = rotatable_groups.data();
  nrg = 0;
  for (int i = 0; i < rotatable_bond_count; i++) {
    const int llim = rotatable_group_bounds.readHost(i);
    const int hlim = rotatable_group_bounds.readHost(i + 1);
    if (cdk.mol_home[llim + 1] == mol_index && hlim - llim - 2 >= cutoff) {
      RotatorGroup tg;
      tg.root_atom = rg_ptr[llim];
      tg.pivot_atom = rg_ptr[llim + 1];
      int k = 0;
      tg.rotatable_atoms.resize(hlim - llim - 2);
      for (int j = llim + 2; j < hlim; j++) {
        tg.rotatable_atoms[k] = rg_ptr[j];
        k++;
      }
      result.push_back(tg);
      nrg++;
    }
  }
  std::sort(result.begin(), result.end(),
            [](RotatorGroup a, RotatorGroup b) {
              return a.rotatable_atoms.size() > b.rotatable_atoms.size();
            });
  
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ChemicalFeatures::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
bool scoreChiralBranches(const std::vector<std::vector<BondedNode>> &links,
                         const std::vector<int> &layer_llim, const std::vector<int> &layer_hlim,
                         const ChemicalDetailsKit &cdk, std::vector<int> *chiral_dominance,
                         std::vector<int> *parallel_growth) {
  std::vector<bool> score_branch(4);
  int* cdom_ptr = chiral_dominance->data();
  score_branch[0] = (cdom_ptr[ 1] == 0 || cdom_ptr[ 2] == 0 || cdom_ptr[ 3] == 0);
  score_branch[1] = (cdom_ptr[ 1] == 0 || cdom_ptr[ 6] == 0 || cdom_ptr[ 7] == 0);
  score_branch[2] = (cdom_ptr[ 2] == 0 || cdom_ptr[ 6] == 0 || cdom_ptr[11] == 0);
  score_branch[3] = (cdom_ptr[ 3] == 0 || cdom_ptr[ 7] == 0 || cdom_ptr[11] == 0);

  // Determine if there is any work to do.  If there are no new layers to offer new scores,
  // return false to indicate that branches no longer need advancement.
  int max_leaves = 0;
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      max_leaves = std::max(layer_hlim[i] - layer_llim[i], max_leaves);
    }
  }
  if (max_leaves == 0) {
    return false;
  }

  // Check for parallel growth.  Two consecutive rounds of branches growing with the same atoms
  // implies that neither branch is superior to the other and thus the center is not chiral.
  int* pgrow_ptr = parallel_growth->data();
  for (int i = 0; i < 3; i++) {
    const int lyl_i = layer_llim[i];
    const int lyh_i = layer_hlim[i];
    for (int j = i + 1; j < 4; j++) {
      if (layer_llim[j] == lyl_i && layer_hlim[j] == lyh_i) {
        bool leaves_identical = true;
        for (int k = lyl_i; k < lyh_i; k++) {
          leaves_identical = (leaves_identical && links[i][k].getAtom() == links[j][k].getAtom());
        }
        if (leaves_identical) {
          pgrow_ptr[(4 * j) + i] += 1;
          pgrow_ptr[(4 * i) + j] += 1;
        }
        else {
          pgrow_ptr[(4 * j) + i] = 0;
          pgrow_ptr[(4 * i) + j] = 0;
        }
      }
    }
  }

  // Set up tables to score the various branches
  int max_znum = 0;
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      for (int j = layer_llim[i]; j < layer_hlim[i]; j++) {
        max_znum = std::max(cdk.z_numbers[links[i][j].getAtom()], max_znum);
      }
    }
  }
  std::vector<std::vector<double>> element_count(4, std::vector<double>(max_znum + 1, 0.0));
  for (int i = 0; i < 4; i++) {
    if (score_branch[i]) {
      for (int j = layer_llim[i]; j < layer_hlim[i]; j++) {
        element_count[i][cdk.z_numbers[links[i][j].getAtom()]] += links[i][j].getRootBondOrder();
      }
    }
  }
  for (int i = 0; i < 3; i++) {
    if (score_branch[i]) {
      for (int j = i + 1; j < 4; j++) {
        const int j4pi = (j * 4) + i;
        if (score_branch[j] && cdom_ptr[j4pi] == 0 && pgrow_ptr[j4pi] < 2) {
          const int i4pj = (i * 4) + j;
          for (int k = max_znum; k >= 1; k--) {
            if (cdom_ptr[j4pi] == 0) {
              if (element_count[i][k] > element_count[j][k] + constants::tiny) {
                cdom_ptr[j4pi] =  1;
                cdom_ptr[i4pj] = -1;
              }
              else if (element_count[i][k] < element_count[j][k] - constants::tiny) {
                cdom_ptr[j4pi] = -1;
                cdom_ptr[i4pj] =  1;
              }
            }
          }
        }
      }
    }
  }

  // If the dominance matrix still contains undecided elements, the search must keep going.
  return ((cdom_ptr[ 1] == 0 && pgrow_ptr[ 1] < 2) || (cdom_ptr[ 2] == 0 && pgrow_ptr[ 2] < 2) ||
          (cdom_ptr[ 3] == 0 && pgrow_ptr[ 3] < 2) || (cdom_ptr[ 6] == 0 && pgrow_ptr[ 6] < 2) ||
          (cdom_ptr[ 7] == 0 && pgrow_ptr[ 7] < 2) || (cdom_ptr[11] == 0 && pgrow_ptr[11] < 2));
}

//-------------------------------------------------------------------------------------------------
int getChiralOrientation(const CoordinateFrameReader &cfr, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {

  // Collect coordinates
  double ra[3], rb[3], rc[3], r_root[3], prja[3], prjb[3], prjc[3], acrb[3], bcrc[3];
  const double cax = cfr.xcrd[center_atom];
  const double cay = cfr.ycrd[center_atom];
  const double caz = cfr.zcrd[center_atom];
  ra[0] = cfr.xcrd[branch_a_atom] - cax;
  ra[1] = cfr.ycrd[branch_a_atom] - cay;
  ra[2] = cfr.zcrd[branch_a_atom] - caz;
  rb[0] = cfr.xcrd[branch_b_atom] - cax;
  rb[1] = cfr.ycrd[branch_b_atom] - cay;
  rb[2] = cfr.zcrd[branch_b_atom] - caz;
  rc[0] = cfr.xcrd[branch_c_atom] - cax;
  rc[1] = cfr.ycrd[branch_c_atom] - cay;
  rc[2] = cfr.zcrd[branch_c_atom] - caz;
  r_root[0] = cfr.xcrd[root_atom] - cax;
  r_root[1] = cfr.ycrd[root_atom] - cay;
  r_root[2] = cfr.zcrd[root_atom] - caz;

  // Remove the projections of each higher priority branch onto the lowest priority branch.
  project(ra, r_root, prja, 3);
  project(rb, r_root, prjb, 3);
  project(rc, r_root, prjc, 3);
  for (int i = 0; i < 3; i++) {
    ra[i] -= prja[i];
    rb[i] -= prjb[i];
    rc[i] -= prjc[i];
  }
  crossProduct(ra, rb, acrb);
  crossProduct(rb, rc, bcrc);

  // L- (S-) chirality occurs when the cross product of ra and rb, as well as the cross product of
  // rb and rc, following the right hand rule, point towards r_root.  However, L- chirality is
  // given a convention of + (as most amino acids, which have chiral centers, are exepcted to be
  // L-chiral).
  return ((2 * (dot(acrb, r_root, 3) < 0.0 && dot(bcrc, r_root, 3))) - 1);
}

} // namespace chemistry
} // namespace omni
