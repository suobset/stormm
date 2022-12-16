#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Constants/hpc_bounds.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Reporting/error_format.h"
#include "rmsd_plan.h"

namespace stormm {
namespace structure {

using card::HybridKind;
using chemistry::ChemicalFeatures;
using diskutil::getBaseName;
using math::roundUp;
  
//-------------------------------------------------------------------------------------------------
RMSDPlanReader::RMSDPlanReader(const int plan_count_in, const RMSDMethod strategy_in,
                               const double mass_fraction_in, const int* alignment_steps_in,
                               const int* core_atoms_in, const int* core_counts_in,
                               const int* core_starts_in, const int* symm_atoms_in,
                               const int4* symm_bounds_in, const int2* symm_ranges_in) :
    plan_count{plan_count_in}, strategy{strategy_in}, mass_fraction{mass_fraction_in},
    alignment_steps{alignment_steps_in}, core_atoms{core_atoms_in}, core_counts{core_counts_in},
    core_starts{core_starts_in}, symm_atoms{symm_atoms_in}, symm_bounds{symm_bounds_in},
    symm_ranges{symm_ranges_in}
{}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const RMSDMethod strategy_in, const double rmf_in) :
    plan_count{0}, general_strategy{strategy_in}, required_mass_fraction{rmf_in},
    alignment_steps{HybridKind::ARRAY, "rplan_align_steps"},
    asymmetric_core_atoms{HybridKind::ARRAY, "rplan_core_atoms"},
    asymmetric_core_counts{HybridKind::ARRAY, "rplan_core_counts"},
    asymmetric_core_starts{HybridKind::ARRAY, "rplan_core_starts"},
    symmetry_group_atoms{HybridKind::ARRAY, "rplan_symm_atoms"},
    symmetry_group_bounds{HybridKind::ARRAY, "rplan_symm_bounds"},
    symmetry_group_ranges{HybridKind::ARRAY, "rplan_symm_ranges"},
    ag_pointers{}
{}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                   const RMSDMethod strategy_in, const double rmf_in, const GpuDetails &gpu,
                   const int low_mol_idx, const int high_mol_idx) :
    RMSDPlan(strategy_in, rmf_in)
{
  // There is one plan that can be defined by the given information.
  plan_count = 1;

  // Get the formal charges, free electron content, chirality, and ring inclusions for the one
  // system.  Each system can, in fact, hold many molecules, but the constructor is set to 
  ChemicalFeatures chemfe(ag_in, cf_in);
  const AtomEquivalence eq_found(chemfe, nullptr, low_mol_idx, high_mol_idx);
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const AtomEquivalence &eq_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu) :
    RMSDPlan(strategy_in, rmf_in)
{
  // There is one plan defined by the given information
  plan_count = 1;

}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu, const int low_mol_idx,
                   const int high_mol_idx) :
    RMSDPlan(strategy_in, rmf_in)
{
  // There will be as many plans as the PhaseSpaceSynthesis has unique topologies.  It will be
  // assumed that a common range of molecules (again, defaulting to the first molecule) applies to
  // the desired RMSD calculations across all systems.
  ag_pointers = poly_ps_in.getUniqueTopologies();
  const std::vector<int> sys_top_indices = poly_ps_in.getUniqueTopologyIndices();
  plan_count = ag_pointers.size();
  std::vector<AtomEquivalence> eq_tables;
  eq_tables.reserve(plan_count);
  for (int i = 0; i < plan_count; i++) {
    const CoordinateFrame cf_example = poly_ps_in.exportCoordinates(sys_top_indices[i]);
    ChemicalFeatures chemfe(ag_pointers[i], cf_example);
    eq_tables.emplace_back(chemfe, nullptr, low_mol_idx, high_mol_idx);
  }
}

//-------------------------------------------------------------------------------------------------
RMSDPlan::RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in,
                   const std::vector<AtomEquivalence> &eq_list_in, const RMSDMethod strategy_in,
                   const double rmf_in, const GpuDetails &gpu) :
    RMSDPlan(strategy_in, rmf_in)
{
  // Check that there are enough plans to cover the various topologies in the synthesis, and the
  // topologies referenced by each plan should match those in the synthesis, at least in terms of
  // the atom counts.
  ag_pointers = poly_ps_in.getUniqueTopologies();
  if (ag_pointers.size() != eq_list_in.size()) {
    rtErr("One atom equivalence object must be provided for each unique topology in the "
          "synthesis.  " + std::to_string(ag_pointers.size()) + " unique topologies were "
          "detected, with " + std::to_string(eq_list_in.size()) + " atom equivalence table to "
          "support them.", "RMSDPlan");
  }
  plan_count = ag_pointers.size();
  for (int i = 0; i < plan_count; i++) {
    if (ag_pointers[i] != eq_list_in[i].getTopologyPointer() &&
        ag_pointers[i]->getAtomCount() != eq_list_in[i].getTopologyPointer()->getAtomCount()) {
      rtErr("The atom count for the unique topology index " + std::to_string(i) + " within the "
            "synthesis, (" + std::to_string(ag_pointers[i]->getAtomCount()) + " originating in "
            "file " + getBaseName(ag_pointers[i]->getFileName()) + "), does not match "
            "the atom count in the atom equivalence table (" +
            std::to_string(eq_list_in[i].getTopologyPointer()->getAtomCount()) + ", originating "
            "in file " + getBaseName(eq_list_in[i].getTopologyPointer()->getFileName()) + ").",
            "RMSDPlan");
    }
  }
}

//-------------------------------------------------------------------------------------------------
int RMSDPlan::getPlanCount() const {
  return plan_count;
}

//-------------------------------------------------------------------------------------------------
RMSDMethod RMSDPlan::getGeneralStrategy() const {
  return general_strategy;
}

//-------------------------------------------------------------------------------------------------
double RMSDPlan::getRequiredMassFraction() const {
  return required_mass_fraction;
}

//-------------------------------------------------------------------------------------------------
const RMSDPlanReader RMSDPlan::data(const HybridTargetLevel tier) const {
  return RMSDPlanReader(plan_count, general_strategy, required_mass_fraction,
                        alignment_steps.data(tier), asymmetric_core_atoms.data(tier),
                        asymmetric_core_counts.data(tier), asymmetric_core_starts.data(tier),
                        symmetry_group_atoms.data(tier), symmetry_group_bounds.data(tier),
                        symmetry_group_ranges.data(tier));
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::addSymmetryGroups(const std::vector<AtomEquivalence> &eq_tables) {

  // Compute the necessary sizes of each Hybrid object.
  const size_t ntab  = eq_tables.size();
  const size_t se_alignment_steps = alignment_steps.size();
  const size_t se_asymmetric_core_atoms = asymmetric_core_atoms.size();
  const size_t se_asymmetric_core_counts = asymmetric_core_counts.size();
  const size_t se_asymmetric_core_starts = asymmetric_core_starts.size();
  const size_t se_symmetry_group_atoms = symmetry_group_atoms.size();
  const size_t se_symmetry_group_bounds = symmetry_group_bounds.size();
  const size_t se_symmetry_group_ranges = symmetry_group_ranges.size();
  const size_t ne_alignment_steps = se_alignment_steps + ntab;
  size_t ne_asymmetric_core_atoms = se_asymmetric_core_atoms;
  const size_t ne_asymmetric_core_counts = se_asymmetric_core_counts + ntab;
  const size_t ne_asymmetric_core_starts = se_asymmetric_core_starts + ntab;
  size_t ne_symmetry_group_atoms = se_symmetry_group_atoms;
  size_t ne_symmetry_group_bounds = se_symmetry_group_bounds;
  const size_t ne_symmetry_group_ranges = se_symmetry_group_ranges + ntab;
  for (size_t i = 0; i < ntab; i++) {
    ne_asymmetric_core_atoms += roundUp(eq_tables[i].getAsymmetricAtoms().size(), warp_size_zu);
    const int ni_groups = eq_tables[i].getGroupCount();
    ne_symmetry_group_bounds += ni_groups;
    int igroup_natom = 0;    
    for (int j = 0; j < ni_groups; j++) {
      const int nij_gsize = eq_tables[i].getGroupSize(j) * eq_tables[i].getGroupOrder(j);
      igroup_natom += nij_gsize;
    }
    ne_symmetry_group_atoms += roundUp(static_cast<size_t>(igroup_natom), warp_size_zu);
  }

  // Resize all objects
  alignment_steps.resize(ne_alignment_steps);
  asymmetric_core_atoms.resize(ne_asymmetric_core_atoms);
  asymmetric_core_counts.resize(ne_asymmetric_core_counts);
  asymmetric_core_starts.resize(ne_asymmetric_core_starts);
  symmetry_group_atoms.resize(ne_symmetry_group_atoms);
  symmetry_group_bounds.resize(ne_symmetry_group_bounds);
  symmetry_group_ranges.resize(ne_symmetry_group_ranges);

  // Determine the alignment protocol, in the case that the system has symmetry-equivalent groups.
  for (size_t i = 0; i < ntab; i++) {
    if (eq_tables[i].getGroupCount() == 0) {
      alignment_steps.putHost(static_cast<int>(RMSDAlignmentProtocol::ALIGN_ALL),
                              se_alignment_steps + i);
    }
    else {
      for (int j = 0; j < eq_tables[i].getGroupCount(); j++) {
        alignment_steps.putHost(static_cast<int>(eq_tables[i].getGroupRule(j)),
                                se_alignment_steps + i);
      }
    }
  }
  
  // Fill out the core atoms array
  size_t asym_core_idx = se_asymmetric_core_atoms;
  for (size_t i = 0; i < ntab; i++) {
    const std::vector<int>& icore = eq_tables[i].getAsymmetricAtoms();
    const size_t isys_width = icore.size();
    for (size_t j = 0; j < isys_width; j++) {
      asymmetric_core_atoms.putHost(icore[j], asym_core_idx);
      asym_core_idx++;
    }
    asym_core_idx = roundUp(asym_core_idx, warp_size_zu);
  }
  
  // Fill out various bounds, offsets, and counts arrays
  size_t asym_counts_idx = se_asymmetric_core_counts;
  size_t asym_starts_idx = se_asymmetric_core_starts;
  int asym_core_starts_pos = asymmetric_core_starts.readHost(asym_starts_idx - 1);
  size_t symm_grp_bound_idx = se_symmetry_group_bounds;
  size_t symm_grp_range_idx = se_symmetry_group_ranges;
  const int4 tmp_sg = symmetry_group_bounds.readHost(se_symmetry_group_bounds - 1);
  int symmetry_group_fill = tmp_sg.x;
  for (size_t i = 0; i < ntab; i++) {
    asymmetric_core_counts.putHost(eq_tables[i].getAsymmetricAtoms().size(), asym_counts_idx);
    asym_core_starts_pos += asymmetric_core_counts.readHost(asym_counts_idx);
    asymmetric_core_starts.putHost(roundUp(asym_core_starts_pos, warp_size_int), asym_starts_idx);
    asym_counts_idx++;
    asym_starts_idx++;
    const int ni_groups = eq_tables[i].getGroupCount();
    for (int j = 0; j < ni_groups; j++) {
      const int nij_gsize  = eq_tables[i].getGroupSize(j);
      const int nij_gorder = eq_tables[i].getGroupOrder(j);
      const int4 sg = { symmetry_group_fill, symmetry_group_fill + (nij_gsize * nij_gorder),
                        nij_gsize, nij_gorder };
      symmetry_group_bounds.putHost(sg, symm_grp_bound_idx);
      symm_grp_bound_idx++;
    }
    symmetry_group_fill = roundUp(symmetry_group_fill, warp_size_int);
    int2 symmetry_group_range;
    symmetry_group_range.x = (i == 0) ? 0 :
                                        symmetry_group_ranges.readHost(symm_grp_range_idx - 1).x;
    symmetry_group_range.y = symmetry_group_range.x + ni_groups;
    symmetry_group_ranges.putHost(symmetry_group_range, symm_grp_range_idx);
    symm_grp_range_idx++;
  }
  
  // Add each symmetry group
  size_t sy_grp_idx = se_symmetry_group_atoms;
  for (size_t i = 0; i < ntab; i++) {

    // Determine the order in which to add each symmetry group.  The largest symmetry groups should
    // be added first, to ensure that cases where the core atoms need to be supplemented by
    // choosing the arrangements of one or more symmetry-related atom groups can be handled by
    // taking the first groups in the list.
    const int ni_group = eq_tables[i].getGroupCount();
    std::vector<int2> group_sizes(ni_group);
    for (int j = 0; j < ni_group; j++) {
      const int nij_gsize = eq_tables[i].getGroupSize(j) * eq_tables[i].getGroupOrder(j);
      group_sizes[j].x = j;
      group_sizes[j].y = nij_gsize;
    }
    std::sort(group_sizes.begin(), group_sizes.end(), [](int2 a, int2 b) { return a.x > b.x; });
    for (int j = 0; j < ni_group; j++) {
      const int grp_idx = group_sizes[j].x;
      const int nij_gsize  = eq_tables[i].getGroupSize(grp_idx);
      const int nij_gorder = eq_tables[i].getGroupOrder(grp_idx);
      for (int k = 0; k < nij_gorder; k++) {
        for (int m = 0; m < nij_gsize; m++) {
          symmetry_group_atoms.putHost(eq_tables[i].getSymmetryRelatedAtom(j, k, m), sy_grp_idx);
          sy_grp_idx++;
        }
      }
    }
    sy_grp_idx = roundUp(sy_grp_idx, warp_size_zu);
  }
}

//-------------------------------------------------------------------------------------------------
void RMSDPlan::writeSymmetryGroupCodes() {

  // Each atom may or may not belong to various symmetry groups, which will imply combinatorial
  // tests of the different symmetry-related partner atoms at various stages of computing the
  // best-fit RMSD.  The key is to make this information as compact as possible, and the way to
  // do that is to create bit-packed strings to encode the information about which symmetry groups
  // the atom is a part of.  Let the input coordinates be given { xyz0, xyz1, xyz2, ..., xyzN }
  // for N atoms and three Cartesian dimension x, y, and z.  View each atom in the molecule as a
  // virtual particle (not a virtual site in the sense of a massless particle, but a conceptual
  // entity) that could take the guise of one of a number of different coordinates from within the
  // given representation.  Atoms of topological indices 0, 1, and 2 might be three methyl
  // hydrogens constituting a symmetry group, with atom 3 being the carbon they are bound to.  This
  // methyl group itself may be symmetric with a second methyl group, with hydrogen atoms at
  // topological indices 7, 8, and 9 and the carbon atom at topological index 6.  First, the
  // decision must be made whether to swap the mythl groups, and then the decision can be made
  // in what order to take the hydrogen atoms.  In the final product, the virtual particles take
  // positions from the given coordinates that are then submitted to a standard RMSD calculation.
  // Virtual particle 3 could take coordinates from atoms 3 or 6.  Virtual particles 0, 1, and 2,
  // as well as 7, 8, and 9, could take coordinates from atoms 0, 1, and 2 in any order, without
  // replacement, or form atoms 7, 8, and 9, in any order without replacement.  The dependencies
  // of symmetry groups imply an order in which virtual particles take the coordinates from
  // particular atoms (processing symmetry groups which are not dependent on any other), then take
  // subsequent positions of other virtual particles (for dependent symmetry groups down the line
  // until there are no dependencies left to process).  In the present example, the groups of atoms
  // { 0, 1, 2, 3 } and { 7, 8, 9, 6 } take positions from atoms { 0, 1, 2, 3 } or { 7, 8, 9, 6 }
  // in the first step.  Next, virtual particles { 0, 1, 2 } choose whether to swap positions
  // among themselves, and virtual particles { 7, 8, 9 } go through the same process.  As a first
  // pass, the procedure is to find the best assignment among the highest level symmetry groups
  // (assigning coordinates from atoms in the given list to virtual particles), then filter down,
  // finding further improvements by attempting different assignments of the dependent groups.

}

} // namespace structure
} // namespace stormm
