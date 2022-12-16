#include <algorithm>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_abstracts.h"
#include "atom_equivalence.h"

namespace stormm {
namespace chemistry {

using math::sum;
using math::prefixSumInPlace;
using math::PrefixSumType;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
AtomRank::AtomRank(const AtomGraph *ag_in) :
    maximum_rank{-1}, ranks{}, rank_degeneracy_bounds{}, rank_instances{},
    ag_pointer{const_cast<AtomGraph*>(ag_in)}
{}
  
//-------------------------------------------------------------------------------------------------
AtomRank::AtomRank(const AtomGraph *ag_in, const std::vector<double> &formal_charges,
                   const std::vector<double> &free_electrons,
                   const std::vector<ullint> &ring_inclusion,
                   const std::vector<ChiralOrientation> &chiralities,
                   std::vector<int> *a_idx_tree, std::vector<int> *b_idx_tree,
                   std::vector<int> *a_zn_tree, std::vector<int> *b_zn_tree,
                   std::vector<double> *a_fc_tree, std::vector<double> *b_fc_tree,
                   std::vector<double> *a_fe_tree, std::vector<double> *b_fe_tree,
                   std::vector<ullint> *a_ri_tree, std::vector<ullint> *b_ri_tree,
                   std::vector<ChiralOrientation> *a_ch_tree,
                   std::vector<ChiralOrientation> *b_ch_tree,
                   std::vector<int> *a_coverage, std::vector<int> *b_coverage,
                   const int low_molecule_index, const int high_molecule_index) :
    AtomRank(ag_in)
{
  // Pull some critical abstracts from the topology and initialize storage vectors, beginning with
  // the result vector.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();
  ranks.resize(cdk.natom, -1);

  // Determine the range of molecules that will be studied
  const int actual_high_mol_idx = (high_molecule_index < 0) ? low_molecule_index + 1 :
                                                              high_molecule_index;

  // Loop over all molecules.  Separate copies of the same molecule will have different ranks for
  // each of their respective atoms.
  int rank_counter = 0;
  for (int mpos = low_molecule_index; mpos < actual_high_mol_idx; mpos++) {
  
    // Pre-compute the number of connections made by each atom to clean up the code and avoid
    // subsequent arithmetic.
    std::vector<int> nbxc(cdk.natom);
    for (int i = cdk.mol_limits[mpos]; i < cdk.mol_limits[mpos + 1]; i++) {
      const int atom_idx = cdk.mol_contents[i];
      nbxc[atom_idx] = nbk.nb12_bounds[atom_idx + 1] - nbk.nb12_bounds[atom_idx];
    }
  
    // Scanning through the topological order of atoms, find the first non-hydrogen atom and then
    // all atoms with the same view to the molecular bonding pattern.  These are the 0th ranked
    // atoms, and subsequent layers built from them will get ranks 1, 2, ..., M, M < N for a system
    // of N atoms.  If no heavy atoms is found, use the topology's first atom.
    std::vector<int> seed_atoms(1, cdk.mol_limits[mpos]);
    for (int i = cdk.mol_limits[mpos]; i < cdk.mol_limits[mpos + 1]; i++) {
      const int atom_idx = cdk.mol_contents[i];
      if (cdk.z_numbers[atom_idx] > 1) {
        seed_atoms[0] = atom_idx;
        break;
      }
    }
    const int seed_rank_zn    = cdk.z_numbers[seed_atoms[0]];
    const double seed_rank_fc = formal_charges[seed_atoms[0]];
    const double seed_rank_fe = free_electrons[seed_atoms[0]];
    const ullint seed_rank_ri = ring_inclusion[seed_atoms[0]];
    const ChiralOrientation seed_rank_ch = chiralities[seed_atoms[0]];
    const int seed_rank_nbxc = nbxc[seed_atoms[0]];
    for (int i = cdk.mol_limits[mpos]; i < cdk.mol_limits[mpos + 1]; i++) {
      const int atom_idx = cdk.mol_contents[i];

      // Perform this comparison ahead of the matchBondingPattern() call to avoid setup work that
      // it would perform to accomplish the same thing.
      if (atom_idx != seed_atoms[0] && cdk.z_numbers[atom_idx] == seed_rank_zn &&
          ring_inclusion[atom_idx] == seed_rank_ri && nbxc[atom_idx] == seed_rank_nbxc &&
          fabs(formal_charges[atom_idx] - seed_rank_fc) < 1.0e-6 &&
          fabs(free_electrons[atom_idx] - seed_rank_fe) < 1.0e-6 &&
          chiralities[atom_idx] == seed_rank_ch &&
          matchBondingPattern(ag_pointer, formal_charges, free_electrons, ring_inclusion,
                              chiralities, seed_atoms[0], atom_idx, a_idx_tree, b_idx_tree,
                              a_zn_tree, b_zn_tree, a_fc_tree, b_fc_tree, a_fe_tree, b_fe_tree,
                              a_ri_tree, b_ri_tree, a_ch_tree, b_ch_tree, a_coverage,
                              b_coverage)) {
        seed_atoms.push_back(atom_idx);
      }
    }

    // Expand outwards from the seed atom(s) until the entire molecule is covered.  Make a list of
    // atoms with each kind of property, including the rank(s) of the atom(s) they connect to in
    // the previous layer.  Given N atoms of a particular rank in the previous layer, up to N
    // connections to atoms with a particular set of properties in the new layer imply the same
    // rank among those connecting atoms.  If there are more than N such connections to atoms with
    // similar local properties, then all atoms in the present layer with those properties must be
    // tested for unique bonding patterns and thus distinct ranks.
    bool add_to_layers = true;
    std::vector<int> next_layer = seed_atoms;
    std::vector<int> next_layer_map(cdk.natom, -1);
    std::vector<std::vector<int>> origin_ranks(16);
    std::vector<int> rank_degeneracy(1, seed_atoms.size());
    std::vector<int> i_lookalikes(16);
    for (size_t i = 0; i < seed_atoms.size(); i++) {
      ranks[seed_atoms[i]] = rank_counter;
    }
    rank_counter++;
    while (next_layer.size() > 0) {
      const std::vector<int> current_layer = next_layer;
      const int n_current_layer = current_layer.size();
      next_layer.resize(0);
      origin_ranks.resize(0);
      for (int i = 0; i < n_current_layer; i++) {
        const int iatom = current_layer[i];
        const int irank = ranks[iatom];
        for (int j = nbk.nb12_bounds[iatom]; j < nbk.nb12_bounds[iatom + 1]; j++) {
          const int jatom = nbk.nb12x[j];

          // If the atom is not part of some previous layer and has not yet been added to the next
          // layer, add it now and start its list of originating ranks.  Mark the index of the next
          // layer at which this atom can be found, in case another atom in the current layer also
          // links to it.
          if (ranks[jatom] == -1) {
            if (next_layer_map[jatom] == -1) {
              next_layer_map[jatom] = next_layer.size();
              next_layer.push_back(jatom);
              origin_ranks.push_back(std::vector<int>(1, irank));
            }
            else {

              // Extend the list of ranks which link to this atom in the next layer.
              origin_ranks[next_layer_map[jatom]].push_back(irank);
            }
          }
        }
      }

      // Sort the origin ranks vectors for all atoms in the next layer to permit rapid comparisons.
      const int n_next_layer = next_layer.size();
      for (int i = 0; i < n_next_layer; i++) {
        switch (origin_ranks[i].size()) {
        case 0:
        case 1:
          break;
        case 2:
          if (origin_ranks[i][0] > origin_ranks[i][1]) {
            std::swap(origin_ranks[i][0], origin_ranks[i][1]);          
          }
          break;
        default:
          std::sort(origin_ranks[i].begin(), origin_ranks[i].end(),
                    []( int a, int b ) { return a < b; });
          break;
        }
      }

      // The next layer is now determined, with the origin ranks constituting a sort of additional
      // "property" for each atom, in addition to the formal charge, free electron content, ring
      // inclusion statistics, chirality, and overall number of atoms to which it bonds.  Use that
      // information, supplemented with explicit checks on the bonding patterns where necessary, to
      // determine ranks for each atom in the layer.
      std::vector<bool> nl_coverage(n_next_layer, false);
      for (int i = 0; i < n_next_layer; i++) {
        if (nl_coverage[i]) {
          continue;
        }
        nl_coverage[i] = true;
        const int iatom = next_layer[i];
        const int iatom_zn = cdk.z_numbers[iatom];
        const int iatom_nbxc = nbxc[iatom];
        const double iatom_fc = formal_charges[iatom];
        const double iatom_fe = free_electrons[iatom];
        const ullint iatom_ri = ring_inclusion[iatom];
        const ChiralOrientation iatom_ch = chiralities[iatom];
        int i_degeneracy = 1;
        const size_t iatom_norig = origin_ranks[i].size();
        i_lookalikes.resize(0);
        i_lookalikes.push_back(i);
        for (int j = i + 1; j < n_next_layer; j++) {
          const int jatom = next_layer[j];
          if (cdk.z_numbers[jatom] == iatom_zn && ring_inclusion[jatom] == iatom_ri &&
              nbxc[jatom] == iatom_nbxc && chiralities[jatom] == iatom_ch &&
              origin_ranks[j].size() == iatom_norig &&
              fabs(formal_charges[jatom] - iatom_fc) < 1.0e-6 &&
              fabs(free_electrons[jatom] - iatom_fe) < 1.0e-6) {
            bool origins_same = true;
            for (size_t k = 0; k < iatom_norig; k++) {
              origins_same = (origins_same &&
                              origin_ranks[i][k] == origin_ranks[j][k]);
            }
            if (origins_same) {
              i_degeneracy++;
              i_lookalikes.push_back(j);
            }
          }
        }
        switch (i_degeneracy) {
        case 1:

          // Only one such atom was found in the next layer.  Because no atoms in distinct layers
          // may share the same rank, this atom must be of a unique rank.
          ranks[iatom] = rank_counter;
          nl_coverage[i] = true;
          rank_counter++;
          break;
        default:

          // Check that there are no more of these atoms than the number of atoms occupying any of
          // the ranks that feed into them.  It may be possible to optimize this check further, but
          // such cases would be exceedingly rare, so use the most conservative condition.
          bool same_rank_safe = true;
          for (size_t j = 0; j < iatom_norig; j++) {
            same_rank_safe = (same_rank_safe &&
                              rank_degeneracy[origin_ranks[i][j]] >= i_degeneracy);
          }
          if (same_rank_safe) {
            for (int j = 0; j < i_degeneracy; j++) {
              const int jatom = next_layer[i_lookalikes[j]];
              ranks[jatom] = rank_counter;
              nl_coverage[i_lookalikes[j]] = true;
            }
            rank_counter++;
          }
          else {
            std::vector<bool> deg_coverage(i_degeneracy, false);
            for (int j = 0; j < i_degeneracy; j++) {
              const int jatom = next_layer[i_lookalikes[j]];
              if (deg_coverage[j]) {
                continue;
              }
              deg_coverage[j] = true;
              ranks[jatom] = rank_counter;
              for (int k = j + 1; k < i_degeneracy; k++) {
                const int katom = next_layer[i_lookalikes[k]];
                if (matchBondingPattern(ag_pointer, formal_charges, free_electrons, ring_inclusion,
                                        chiralities, jatom, katom, a_idx_tree, b_idx_tree,
                                        a_zn_tree, b_zn_tree, a_fc_tree, b_fc_tree, a_fe_tree,
                                        b_fe_tree, a_ri_tree, b_ri_tree, a_ch_tree, b_ch_tree,
                                        a_coverage, b_coverage)) {
                  ranks[katom] = rank_counter;
                  deg_coverage[k] = true;
                  nl_coverage[i_lookalikes[k]] = true;
                }
              }
              nl_coverage[i_lookalikes[j]] = true;
              rank_counter++;
            }
          }
          break;
        }
      }

      // With the ranks of all atoms in the next layer having been found, update the list of rank
      // degeneracies.  Also, clear the next layer's map to avoid pollution on the next iteration.
      rank_degeneracy.resize(rank_counter, 0);
      for (int i = 0; i < n_next_layer; i++) {
        rank_degeneracy[ranks[next_layer[i]]] += 1;
        next_layer_map[next_layer[i]] = -1;
      }
    }
  }
  
  // Compute the degeneracy of each atom rank throughout the entire system
  maximum_rank = rank_counter;
  rank_degeneracy_bounds.resize(maximum_rank + 1, 0);
  if (actual_high_mol_idx - low_molecule_index < cdk.nmol) {
    const size_t ilim = cdk.mol_limits[actual_high_mol_idx];
    for (size_t i = cdk.mol_limits[low_molecule_index]; i < ilim; i++) {
      const size_t atom_idx = cdk.mol_contents[i];
      rank_degeneracy_bounds[ranks[atom_idx]] += 1;
    }
  }
  else {
    for (int i = 0; i < cdk.natom; i++) {
      if (ranks[i] >= 0) {
        rank_degeneracy_bounds[ranks[i]] += 1;
      }
    }
  }
  prefixSumInPlace<int>(&rank_degeneracy_bounds, PrefixSumType::EXCLUSIVE);
  std::vector<int> rank_dg_counters = rank_degeneracy_bounds;
  rank_instances.resize(rank_degeneracy_bounds[maximum_rank]);
  if (actual_high_mol_idx - low_molecule_index < cdk.nmol) {
    const size_t ilim = cdk.mol_limits[actual_high_mol_idx];
    for (size_t i = cdk.mol_limits[low_molecule_index]; i < ilim; i++) {
      const size_t atom_idx = cdk.mol_contents[i];
      const int pos = rank_dg_counters[ranks[atom_idx]];
      rank_instances[pos] = atom_idx;
      rank_dg_counters[ranks[atom_idx]] = pos + 1;
    }
  }
  else {
    for (int i = 0; i < cdk.natom; i++) {
      if (ranks[i] >= 0) {
        const int pos = rank_dg_counters[ranks[i]];
        rank_instances[pos] = i;
        rank_dg_counters[ranks[i]] = pos + 1;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
AtomRank::AtomRank(const AtomGraph &ag, const std::vector<double> &formal_charges,
                   const std::vector<double> &free_electrons,
                   const std::vector<ullint> &ring_inclusion,
                   const std::vector<ChiralOrientation> &chiralities, const int low_molecule_index,
                   const int high_molecule_index) :
    AtomRank()
{
  // Allocate space for the matchBondingPattern() function
  const int natom = ag.getAtomCount();
  std::vector<int> a_idx_tree(natom), b_idx_tree(natom), a_zn_tree(natom), b_zn_tree(natom);
  std::vector<double> a_fc_tree(natom), b_fc_tree(natom), a_fe_tree(natom), b_fe_tree(natom);
  std::vector<ullint> a_ri_tree(natom), b_ri_tree(natom);
  std::vector<ChiralOrientation> a_ch_tree(natom), b_ch_tree(natom);
  std::vector<int> a_coverage(natom), b_coverage(natom);
  *this = AtomRank(ag.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                   chiralities, &a_idx_tree, &b_idx_tree, &a_zn_tree, &b_zn_tree, &a_fc_tree,
                   &b_fc_tree, &a_fe_tree, &b_fe_tree, &a_ri_tree, &b_ri_tree, &a_ch_tree,
                   &b_ch_tree, &a_coverage, &b_coverage, low_molecule_index, high_molecule_index);
}

//-------------------------------------------------------------------------------------------------
AtomRank::AtomRank(const AtomGraph *ag, const ChemicalFeatures &chemfe,
                   const int low_molecule_index, const int high_molecule_index) :
    AtomRank(*ag, chemfe.getFormalCharges(), chemfe.getFreeElectrons(), chemfe.getRingInclusion(),
             chemfe.getAtomChirality(), low_molecule_index, high_molecule_index)
{}

//-------------------------------------------------------------------------------------------------
AtomRank::AtomRank(const AtomGraph &ag, const ChemicalFeatures &chemfe,
                   const int low_molecule_index, const int high_molecule_index) :
    AtomRank(ag, chemfe.getFormalCharges(), chemfe.getFreeElectrons(), chemfe.getRingInclusion(),
             chemfe.getAtomChirality(), low_molecule_index, high_molecule_index)
{}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& AtomRank::getRanks() const {
  return ranks;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomRank::getRanks(const int low_index, const int high_index) const {
  if (low_index < 0 || high_index >= ag_pointer->getAtomCount() || low_index >= high_index) {
    rtErr("The range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          "is invalid for a list of ranks for " + std::to_string(ag_pointer->getAtomCount()) +
          "atoms.", "AtomRank", "getRank");
  }
  std::vector<int> result(high_index - low_index);
  for (int i = low_index; i < high_index; i++) {
    result[i - low_index] = ranks[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomRank::getRank(const int atom_index) const {
  return ranks[atom_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomRank::getAtomsWithRank(const int rank_value) const {
  if (rank_value < 0 || rank_value >= maximum_rank) {
    rtErr("Rank " + std::to_string(rank_value) + " is invalid for a list of " +
          std::to_string(maximum_rank) + " ranks.", "AtomRank", "getAtomsWithRank");
  }
  const int imin = rank_degeneracy_bounds[rank_value];
  const int imax = rank_degeneracy_bounds[rank_value + 1];
  std::vector<int> result(imax - imin);
  for (int i = imin; i < imax; i++) {
    result[i - imin] = rank_instances[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomRank::getRankPartners(const int atom_index) const {
  if (atom_index < 0 || atom_index >= ag_pointer->getAtomCount()) {
    rtErr("Atom index " + std::to_string(atom_index) + " is out of range for a system with " +
          std::to_string(ag_pointer->getAtomCount()) + " atoms.", "AtomRank", "getRankPartners");
  }
  if (ranks[atom_index] < 0) {
    rtErr("No rank has yet been computed for atom index " + std::to_string(atom_index) + ".",
          "AtomRank", "getRankPartners");
  }
  return getAtomsWithRank(ranks[atom_index]);
}

//-------------------------------------------------------------------------------------------------
AtomEquivalence::AtomEquivalence() :
    group_count{0}, group_atoms{}, group_sizes{}, group_orders{}, group_bounds{}, dependencies{},
    dependency_bounds{}, asymmetric_atoms{}, group_rules{}, ag_pointer{nullptr}
{}

//-------------------------------------------------------------------------------------------------
AtomEquivalence::AtomEquivalence(const AtomGraph *ag_in, const std::vector<double> &formal_charges,
                                 const std::vector<double> &free_electrons,
                                 const std::vector<ullint> &ring_inclusion,
                                 const std::vector<ChiralOrientation> &chiralities,
                                 StopWatch *timer, const int low_molecule_index,
                                 const int high_molecule_index) :
    AtomEquivalence()
{
  // Assign any time up until this point to miscellaneous
  int aeq_timings, dep_timings;
  if (timer != nullptr) {
    aeq_timings = timer->addCategory("Atom equivalence calculation");
    dep_timings = timer->addCategory("Atom equivalence dependencies");
    timer->assignTime(0);
  }
  
  // Set the topology pointer
  ag_pointer = const_cast<AtomGraph*>(ag_in);

  // Determine the upper limit of molecules to assess for symmetry
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const int actual_high_mol_idx = (high_molecule_index == -1) ? low_molecule_index + 1 :
                                                                high_molecule_index;
  if (low_molecule_index < 0 || actual_high_mol_idx > cdk.nmol ||
      low_molecule_index >= actual_high_mol_idx) {
    rtErr("Molecule range " + std::to_string(low_molecule_index) + " to " +
          std::to_string(actual_high_mol_idx) + " is invalid for a system with " +
          std::to_string(cdk.nmol) + " molecules.", "AtomEquivalence");
  }

  // Pre-allocate various arrays and pass them down so that they do not need to be repeatedly
  // allocated by recursive function calls in and under findEquivalentAtoms().
  std::vector<int> a_idx_tree(cdk.natom), b_idx_tree(cdk.natom), a_zn_tree(cdk.natom);
  std::vector<int> b_zn_tree(cdk.natom), a_coverage(cdk.natom), b_coverage(cdk.natom);
  std::vector<int> domain_coverage(cdk.natom), allowed_atoms(cdk.natom);
  std::vector<int> candidate_hopper(cdk.natom), domain_assignment(cdk.natom);
  std::vector<int> jumbled_groups(cdk.natom), aligned_groups(cdk.natom), layer_bounds(cdk.natom);
  std::vector<double> a_fc_tree(cdk.natom), b_fc_tree(cdk.natom), a_fe_tree(cdk.natom);
  std::vector<double> b_fe_tree(cdk.natom);
  std::vector<ullint> a_ri_tree(cdk.natom), b_ri_tree(cdk.natom);
  std::vector<ChiralOrientation> a_ch_tree(cdk.natom), b_ch_tree(cdk.natom);

  // Compute atom ranks
  const AtomRank arnks(ag_pointer, formal_charges, free_electrons, ring_inclusion, chiralities,
                       &a_idx_tree, &b_idx_tree, &a_zn_tree, &b_zn_tree, &a_fc_tree, &b_fc_tree,
                       &a_fe_tree, &b_fe_tree, &a_ri_tree, &b_ri_tree, &a_ch_tree, &b_ch_tree,
                       &a_coverage, &b_coverage, low_molecule_index, actual_high_mol_idx);
  
  // Find top-level symmetry groups.  Once these have been identified, the search will continue for
  // sub-symmetries within each group.
  group_count = 0;
  group_bounds.push_back(0);
  std::vector<int> map_to_subset(cdk.natom, -1);
  for (int hmol = low_molecule_index; hmol < actual_high_mol_idx; hmol++) {
    std::vector<int> hmol_atoms(cdk.mol_limits[hmol + 1] - cdk.mol_limits[hmol]);
    const int hmol_llim = cdk.mol_limits[hmol];
    const int hmol_hlim = cdk.mol_limits[hmol + 1];
    const int hmol_natom = hmol_hlim - hmol_llim;
    for (int i = hmol_llim; i < hmol_hlim; i++) {
      hmol_atoms[i - hmol_llim] = cdk.mol_contents[i];
      map_to_subset[cdk.mol_contents[i]] = i - hmol_llim;
    }
    findEquivalentAtoms(hmol_atoms, &map_to_subset, formal_charges, free_electrons, ring_inclusion,
                        chiralities, arnks, &domain_coverage, &allowed_atoms, &candidate_hopper,
                        &domain_assignment, &jumbled_groups, &aligned_groups, &layer_bounds);
  }
  if (timer != nullptr) timer->assignTime(aeq_timings);
  
  // Determine the dependencies of symmetry-equivalent groups on one another.  Having the atom
  // indices pre-sorted is again helpful in this respect.
  findDependencies();

  // Determine the list of asymmetric atoms (those which are not in any symmetry group)
  std::vector<bool> in_symmetry_group(cdk.natom, false);
  const int glim = group_bounds[group_count];
  for (int i = 0; i < glim; i++) {
    in_symmetry_group[group_atoms[i]] = true;
  }
  int nasym = 0;
  for (int i = 0; i < cdk.natom; i++) {
    nasym += (! in_symmetry_group[i]);
  }
  asymmetric_atoms.reserve(nasym);
  for (int i = 0; i < cdk.natom; i++) {
    if (! in_symmetry_group[i]) {
      asymmetric_atoms.push_back(i);
    }
  }
  if (timer != nullptr) timer->assignTime(dep_timings);
}

//-------------------------------------------------------------------------------------------------
AtomEquivalence::AtomEquivalence(const AtomGraph &ag_in, const std::vector<double> &formal_charges,
                                 const std::vector<double> &free_electrons,
                                 const std::vector<ullint> &ring_inclusion,
                                 const std::vector<ChiralOrientation> &chiralities,
                                 StopWatch *timer, const int low_molecule_index,
                                 const int high_molecule_index) :
    AtomEquivalence(ag_in.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                    chiralities, timer, low_molecule_index, high_molecule_index)
{}

//-------------------------------------------------------------------------------------------------
AtomEquivalence::AtomEquivalence(const ChemicalFeatures &chemfe_in, StopWatch *timer,
                                 const int low_molecule_index, const int high_molecule_index) :
    AtomEquivalence(chemfe_in.getTopologyPointer(), chemfe_in.getFormalCharges(),
                    chemfe_in.getFreeElectrons(), chemfe_in.getRingInclusion(),
                    chemfe_in.getAtomChirality(), timer, low_molecule_index, high_molecule_index)
{}

//-------------------------------------------------------------------------------------------------
int AtomEquivalence::getGroupCount() const {
  return group_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomEquivalence::getGroup(const int group_index) const {
  const size_t llim = group_bounds[group_index];
  const size_t hlim = group_bounds[group_index + 1];
  std::vector<int> result(hlim - llim);
  for (size_t i = llim; i < hlim; i++) {
    result[i - llim] = group_atoms[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int AtomEquivalence::getSymmetryRelatedAtom(const int group_index, const int domain_index,
                                            const int atom_index) const {
  const size_t lookup = group_bounds[group_index] + (domain_index * group_sizes[group_index]) +
                        atom_index;
  return group_atoms[lookup];
}

//-------------------------------------------------------------------------------------------------
const int* AtomEquivalence::getGroupPointer(const int group_index) const {
  return &group_atoms.data()[group_bounds[group_index]];
}

//-------------------------------------------------------------------------------------------------
int AtomEquivalence::getGroupSize(const int group_index) const {
  return group_sizes[group_index];
}

//-------------------------------------------------------------------------------------------------
int AtomEquivalence::getGroupOrder(const int group_index) const {
  return group_orders[group_index];
}

//-------------------------------------------------------------------------------------------------
EquivalenceSwap AtomEquivalence::getGroupRule(int group_index) const {
  return group_rules[group_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomEquivalence::getGroupDependencies(int group_index) const {
  const size_t llim = dependency_bounds[group_index];
  const size_t hlim = dependency_bounds[group_index + 1];
  std::vector<int> result(hlim - llim);
  for (size_t i = llim; i < hlim; i++) {
    result[i - llim] = dependencies[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& AtomEquivalence::getAsymmetricAtoms() const {
  return asymmetric_atoms;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* AtomEquivalence::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::findEquivalentAtoms(const std::vector<int> &subset_atoms,
                                          std::vector<int> *map_to_subset,
                                          const std::vector<double> &formal_charges,
                                          const std::vector<double> &free_electrons,
                                          const std::vector<ullint> &ring_inclusion,
                                          const std::vector<ChiralOrientation> &chiralities,
                                          const AtomRank &arnks, std::vector<int> *domain_coverage,
                                          std::vector<int> *allowed_atoms,
                                          std::vector<int> *candidate_hopper,
                                          std::vector<int> *domain_assignment,
                                          std::vector<int> *jumbled_groups,
                                          std::vector<int> *aligned_groups,
                                          std::vector<int> *layer_bounds) {

  // Track whether atoms in the subset have been included in symmetry groups during this recursive
  // function call.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const int nsubset_atom = subset_atoms.size();
  std::vector<bool> in_symmetry_group(nsubset_atom, false);
  int* map_to_subset_ptr = map_to_subset->data();
  const AtomGraph& ag = *ag_pointer;
  const int* atom_rank_ptr = arnks.getRanks().data();
  for (int i = 0; i < nsubset_atom; i++) {
    const int atom_i = subset_atoms[i];
    const int zn_i = cdk.z_numbers[atom_i];
    const double fc_i = formal_charges[atom_i];
    const double fe_i = free_electrons[atom_i];
    const ullint ri_i = ring_inclusion[atom_i];
    const ChiralOrientation ch_i = chiralities[atom_i];
    if (in_symmetry_group[i]) {
      continue;
    }
    const std::vector<int> tmp_eq_partners = arnks.getRankPartners(atom_i);
    const size_t n_tmp_partners = tmp_eq_partners.size();
    std::vector<int> eq_partners;
    for (size_t j = 0; j < n_tmp_partners; j++) {
      for (int k = 0; k < nsubset_atom; k++) {

        // Messy to figure out which atoms are in the subset and correlate those to the known
        // equivalence partners.  But, this method is a brute-force way to do it.
        if (subset_atoms[k] == tmp_eq_partners[j] && in_symmetry_group[k] == false) {
          eq_partners.push_back(tmp_eq_partners[j]);
        }
      }
    }
    if (eq_partners.size() >= 2) {

      // Find the extent of the equivalent groups that the equivalent atoms may represent, as
      // defined by the groups of atoms that extend outward from each atom to the edges of the
      // molecule and backwards toward the rest of the molecule, up to the point that the groups
      // bump into each other.  This will work whether the atoms came from an entire molecule or
      // a subset of it.  The pointer arrays fed into drawEquivalentGroups() are chosen based on
      // the data type, are for scratch work only, and may emerge modified.  Nothing that happens
      // to the arrays will have any consequence for recursive calls to this function.
      const int last_group_count = group_count;
      drawEquivalentGroups(eq_partners, subset_atoms, formal_charges, free_electrons,
                           ring_inclusion, chiralities, arnks, domain_coverage, allowed_atoms,
                           candidate_hopper, domain_assignment, jumbled_groups, aligned_groups,
                           layer_bounds);

      // The drawEquivalentGroups() function will have constructed one or more new groups and
      // incremented the counter.  To log the atoms that may have been included in the new
      // group(s), scan over the group_atoms array based on the appropriate bounds.
      for (int j = group_bounds[last_group_count]; j < group_bounds[group_count]; j++) {
        in_symmetry_group[map_to_subset_ptr[group_atoms[j]]] = true;
      }

      // Recursive function calls for each symmetry-related subset of atoms
      for (int j = last_group_count; j < group_count; j++) {
        std::vector<int> inner_subset_atoms(group_sizes[j]);
        if (group_sizes[j] < 2) {
          continue;
        }
        std::vector<int2> outer_map_to_subset(group_sizes[j]);
        for (int k = 0; k < group_orders[j]; k++) {

          // Record the current mapping of atoms to the subset, so that it may be put back
          // after returning from the recursive function call and thus avoid the need to
          // allocate a new array of the total number of topological atoms.  The x member of
          // each tuple stores the topological index of the atom whose mapping was modified, and
          // the y member stores the value that it should be returned t after exiting the
          // recursive function call.
          for (int m = 0; m < group_sizes[j]; m++) {
            const int atom_m = group_atoms[group_bounds[j] + (k * group_sizes[j]) + m];
            inner_subset_atoms[m] = atom_m;
            outer_map_to_subset[m].x = atom_m;
            outer_map_to_subset[m].y = map_to_subset_ptr[atom_m];
            map_to_subset_ptr[atom_m] = m;
          }

          // Recursively call this function.
          findEquivalentAtoms(inner_subset_atoms, map_to_subset, formal_charges, free_electrons,
                              ring_inclusion, chiralities, arnks, domain_coverage, allowed_atoms,
                              candidate_hopper, domain_assignment, jumbled_groups, aligned_groups,
                              layer_bounds);

          // Restore the previous subset mapping.
          for (int m = 0; m < group_sizes[j]; m++) {
            map_to_subset_ptr[outer_map_to_subset[m].x] = outer_map_to_subset[m].y;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::drawEquivalentGroups(const std::vector<int> &partners,
                                           const std::vector<int> &subset_atoms,
                                           const std::vector<double> &formal_charges,
                                           const std::vector<double> &free_electrons,
                                           const std::vector<ullint> &ring_inclusion,
                                           const std::vector<ChiralOrientation> &chiralities,
                                           const AtomRank &arnks,
                                           std::vector<int> *domain_coverage,
                                           std::vector<int> *allowed_atoms,
                                           std::vector<int> *candidate_hopper,
                                           std::vector<int> *domain_assignment,
                                           std::vector<int> *jumbled_groups,
                                           std::vector<int> *aligned_groups,
                                           std::vector<int> *layer_bounds) {

  // Get the abstracts locally to avoid de-referencing to the greatest extent possible.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const NonbondedKit<double> nbk = ag_pointer->getDoublePrecisionNonbondedKit();

  // Clear the coverage for the entire subset of atoms.  All equivalence partners will be within
  // the subset.  At the same time, mark a separate array, itself capable of holding a mask for
  // all atoms in the topology, with special values to indicate that atoms are within the subset
  // at all.  At the end of this routine, this scratch array will be cleared to ensure that the
  // space of allowed atoms is not polluted for future calls to this routine.
  int* domain_coverage_ptr = domain_coverage->data();
  int* allowed_atoms_ptr = allowed_atoms->data();
  int* domain_assignment_ptr = domain_assignment->data();
  const int* atom_rank_ptr = arnks.getRanks().data();
  const int n_subset = subset_atoms.size();
  const int subset_code = -4000;
  for (int i = 0; i < n_subset; i++) {
    domain_coverage_ptr[subset_atoms[i]] = 0;

    // The setting of -4000 indicates that this atom is one of the allowed settings.  No other
    // procedure using the array passed in for allowed_atoms will generate such a value.
    allowed_atoms_ptr[subset_atoms[i]] = subset_code;
    domain_assignment_ptr[subset_atoms[i]] = -1;
  }

  // Set pointers to pre-allocated data arrays
  int* jumbled_group_ptr = jumbled_groups->data();
  int* aligned_group_ptr = aligned_groups->data();
  int* candidate_hopper_ptr = candidate_hopper->data();

  // Start with each partner atom and expand outwards until hitting the edges of the molecule or
  // the coverage of another partner.  Once no additional atoms can be added, the groups are
  // complete.  Each partner should arrive at the extent of its symmetry group simultaneously.
  const int n_partners = partners.size();
  for (int i = 0; i < n_partners; i++) {
    domain_coverage_ptr[partners[i]] = 1;
    domain_assignment_ptr[partners[i]] = i;
    jumbled_group_ptr[i] = partners[i];
  }
  bool add_to_domains = true;
  int atoms_per_domain = 1;
  int last_layer_start = 0;
  int last_layer_end = 1;
  if (layer_bounds->size() < 8) {
    layer_bounds->resize(8);
  }
  int* layer_bounds_ptr = layer_bounds->data();
  layer_bounds_ptr[0] = 0;
  layer_bounds_ptr[1] = 1;
  int n_layers = 1;
  int candidate_space = candidate_hopper->size();
  while (add_to_domains) {

    // Loop over all atoms in this group and propose atoms to add to each group.  The number of
    // candidates will be the same for each partner as the layers of their respective domains grow
    // in step with one another.  The number of candidates found will be used as an index into the
    // hopper array, then preserved as the total number of candidates for each partner for use in
    // the subsequent loop.
    int n_candidates;
    for (int i = 0; i < n_partners; i++) {
      n_candidates = 0;
      for (int j = last_layer_start; j < last_layer_end; j++) {
        const int jatom = jumbled_group_ptr[(j * n_partners) + i];
        for (int k = nbk.nb12_bounds[jatom]; k < nbk.nb12_bounds[jatom + 1]; k++) {

          // Only make each atom a candidate if some previous pass has not already covered it.
          if (domain_coverage_ptr[nbk.nb12x[k]] == 0 &&
              allowed_atoms_ptr[nbk.nb12x[k]] == subset_code) {

            // Check to ensure that the hopper array does not overflow.  If that is a danger,
            // resize the array and reset the data pointer.
            if ((n_candidates + 1) * n_partners > candidate_space) {
              candidate_space = (n_candidates + 8) * n_partners;
              candidate_hopper->resize(candidate_space);
              candidate_hopper_ptr = candidate_hopper->data();
            }
            candidate_hopper_ptr[(n_candidates * n_partners) + i] = nbk.nb12x[k];
            n_candidates++;
          }
        }
      }
    }

    // Loop back over all proposed candidates and assess the coverage.  If two or more partners'
    // domains propose to incorporate the same candidate, the coverage will be incremented more
    // than once and the atom should go into none of the domains.  The coverage array will remain
    // in the super-subscribed (e.g. +2 for two competing domains, +3 for three competing domains)
    // state so that future passes don't try to incorporate the atom at all.
    for (int i = 0; i < n_partners; i++) {
      for (int j = 0; j < n_candidates; j++) {
        const int candidate_atom = candidate_hopper_ptr[(j * n_partners) + i];
        domain_coverage_ptr[candidate_atom] += 1;
      }
    }

    // Commit atoms that fit neatly within a single equivalence partner's domain.  Atoms that are
    // claimed by two competing domains will get dropped.  As before, use atoms_per_domain to
    // increment the index counter for adding new atoms, resetting it for each equivalence partner,
    // then hold the value from the last group as the new total size of each equivalence domain.
    for (int i = 0; i < n_partners; i++) {
      atoms_per_domain = last_layer_end;
      for (int j = 0; j < n_candidates; j++) {
        const int candidate_atom = candidate_hopper_ptr[(j * n_partners) + i];
        if (domain_assignment_ptr[candidate_atom] < 0) {
          if (domain_coverage_ptr[candidate_atom] == 1) {
            jumbled_group_ptr[(atoms_per_domain * n_partners) + i] = candidate_atom;
            domain_assignment_ptr[candidate_atom] = i;
            atoms_per_domain++;
          }
          else {

            // Search the neighbors of any atom that was covered more than once--if all coverage
            // came from the same domain, it can be included in that domain as well.
            int belongs_to_domain = -1;
            bool competing_domains = false;
            for (int k = nbk.nb12_bounds[candidate_atom]; k < nbk.nb12_bounds[candidate_atom + 1];
                 k++) {
              if (allowed_atoms_ptr[nbk.nb12x[k]] == subset_code &&
                  domain_assignment_ptr[nbk.nb12x[k]] >= 0) {
                if (belongs_to_domain == -1) {
                  belongs_to_domain = domain_assignment_ptr[nbk.nb12x[k]];
                }
                else if (belongs_to_domain != domain_assignment_ptr[nbk.nb12x[k]]) {
                  competing_domains = true;
                }
              }
            }
            if (competing_domains == false) {
              jumbled_group_ptr[(atoms_per_domain * n_partners) + i] = candidate_atom;
              domain_assignment_ptr[candidate_atom] = i;
              atoms_per_domain++;
            }
          }
        }
      }
    }

    // Record the bounds of this new layer in each partner's domain.  Continue adding layers if
    // the domains are still growing.
    if (atoms_per_domain > last_layer_end) {
      n_layers++;
      if (n_layers >= layer_bounds->size()) {
        layer_bounds->resize(n_layers + 16);
        layer_bounds_ptr = layer_bounds->data();
      }
      layer_bounds_ptr[n_layers] = atoms_per_domain;
      last_layer_start = last_layer_end;
      last_layer_end = atoms_per_domain;
    }
    else {
      add_to_domains = false;
    }
  }

  // The jumbled groups were formed in strides of n_partners.  Rearrange the aligned groups to
  // put each domain into a contiguous space, then allocate scratch space for out-of-place
  // re-ordering of each layer.
  for (int i = 0; i < atoms_per_domain; i++) {
    aligned_group_ptr[i] = jumbled_group_ptr[(i * n_partners)];
  }
  for (int i = 1; i < n_partners; i++) {
    aligned_group_ptr[(i * atoms_per_domain)] = jumbled_group_ptr[i]; 
  }
  int largest_layer = 0;
  for (int i = 0; i < n_layers; i++) {
    largest_layer = std::max(largest_layer, layer_bounds_ptr[i + 1] - layer_bounds_ptr[i]);
  }
  std::vector<bool> links_found(largest_layer);
  std::vector<int> correspondence_candidates(largest_layer);

  // Covert the stored pointer into a reference to the topology
  const AtomGraph& ag = *ag_pointer;
  
  // Sort the domains, layer by layer, to arrange symmetric atoms at each position.  In effect, the
  // prior pass only determined the extend of each domain, even though the atom indices that each
  // domain will ultimately contain needed to be recorded in order to accomplish that.  The atoms
  // in the jumbled_groups array are valid but not in trustworthy order, however.  Retrace the
  // tree of domain atoms for the first equivalent atom partner, noting how many atoms each
  // branches to in the subsequent layer (if two atoms of layer K connect to an atom in layer K+1,
  // the first atom of layer K making the connection will be credited).  This tree will then be
  // used as a template for retracing the same order in other symmetry-equivalent domains, a
  // process that is, in principle, still potentially explosive in its complexity but likely to
  // remain well bounded in practice.
  for (int i = 0; i < n_layers - 1; i++) {

    // Assume that the correspondence for this layer has been worked out and that, for all atoms in
    // it, corresponding atoms in each domain are placed in the same order. (This is implicitly
    // true for the first layer.) Loop over all of the non-bonded connections for the atoms in the
    // ith layer of the first domain, verify that they are present in the next layer, find the
    // corresponding atoms in each of the other domains, and then order them as appropriate.
    const int llim = layer_bounds_ptr[i];
    const int hlim = layer_bounds_ptr[i + 1];
    const int xlim = layer_bounds_ptr[i + 2];
    for (int j = hlim; j < xlim; j++) {
      links_found[j - hlim] = false;
    }
    for (int j = llim; j < hlim; j++) {

      // Loop over all atoms to which the one in the present layer bonds.  Verify that they are in
      // the next layer, then seek out the corresponding atoms in other domains.  Mark that each
      // atom has been scouted, to ensure that, in the case of multimple atoms in layer i bonding
      // to the same atom in layer i + 1, effort is not duplicated.
      const size_t atomj_zero = aligned_group_ptr[j];
      for (int k = nbk.nb12_bounds[atomj_zero]; k < nbk.nb12_bounds[atomj_zero + 1]; k++) {
        const int atomk = nbk.nb12x[k];
        int next_layer_kloc = -1;
        for (int m = hlim; m < xlim; m++) {
          if (aligned_group_ptr[m] == atomk) {
            next_layer_kloc = m;
          }
        }
        if (next_layer_kloc >= 0 && links_found[next_layer_kloc - hlim] == false) {

          // Determine the properties of this atom, then go looking for it in the corresponding
          // layer of each of the other domains.
          const int atomk_zn = cdk.z_numbers[atomk];
          const int atomk_links = nbk.nb12_bounds[atomk + 1] - nbk.nb12_bounds[atomk];
          const double atomk_fc = formal_charges[atomk];
          const double atomk_fe = free_electrons[atomk];
          const ullint atomk_ri = ring_inclusion[atomk];
          for (int dom_idx = 1; dom_idx < n_partners; dom_idx++) {
            const int atomj_next = aligned_group_ptr[(dom_idx * atoms_per_domain) + j];
            int ncorr = 0;
            for (int m = nbk.nb12_bounds[atomj_next]; m < nbk.nb12_bounds[atomj_next + 1]; m++) {
              const int atomm = nbk.nb12x[m];
              bool atomm_in_next_layer = false;
              for (int n = hlim; n < xlim; n++) {
                atomm_in_next_layer = (atomm_in_next_layer ||
                                       jumbled_group_ptr[(n * n_partners) + dom_idx] == atomm);
              }
              if (atomm_in_next_layer) {
                const int atomm_zn = cdk.z_numbers[atomm];
                if (atomm_zn == atomk_zn) {
                  const int atomm_links = nbk.nb12_bounds[atomm + 1] - nbk.nb12_bounds[atomm];
                  const double atomm_fc = formal_charges[atomm];
                  const double atomm_fe = free_electrons[atomm];
                  const ullint atomm_ri = ring_inclusion[atomm];
                  if (atomk_links == atomm_links && atomk_ri == atomm_ri &&
                      fabs(atomk_fc - atomm_fc) < 1.0e-6 && fabs(atomk_fe - atomm_fe) < 1.0e-6) {
                    correspondence_candidates[ncorr] = atomm;
                    ncorr++;
                  }
                }
              }
            }
            if (ncorr == 1) {

              // There is only one possible solution.  Mark that as the corresponding atom, then
              // erase the atom from the "jumbled" list as it is no longer available.
              const int atomm = correspondence_candidates[0];
              aligned_group_ptr[(dom_idx * atoms_per_domain) + next_layer_kloc] = atomm;
              for (int n = hlim; n < xlim; n++) {
                if (jumbled_group_ptr[(n * n_partners) + dom_idx] == atomm) {
                  jumbled_group_ptr[(n * n_partners) + dom_idx] = -1;
                }
              }
            }
            else {
              for (int m = 0; m < ncorr; m++) {
                const int atomm = correspondence_candidates[m];
                if (atom_rank_ptr[atomk] == atom_rank_ptr[atomm]) {

                  // The properties of the atom linked in the first domain and the atom linked
                  // in domain dom_idx are identical, as are their bonding patterns.  Transfer the
                  // atom from the "jumbled" list to the "aligned" list.
                  aligned_group_ptr[(dom_idx * atoms_per_domain) + next_layer_kloc] = atomm;
                  for (int n = hlim; n < xlim; n++) {
                    if (jumbled_group_ptr[(n * n_partners) + dom_idx] == atomm) {
                      jumbled_group_ptr[(n * n_partners) + dom_idx] = -1;
                    }
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Rearrange the aligned groups: put each domain's atom indices in increasing topological order.
  // This will help later, when making comparisons between the atoms in proposed symmetry groups
  // and those which are already known.
  std::vector<int> buffer(atoms_per_domain);
  for (int i = 0; i < n_partners; i++) {
    for (int j = 0; j < atoms_per_domain; j++) {
      buffer[j] = aligned_group_ptr[(i * atoms_per_domain) + j];
    }
    std::sort(buffer.begin(), buffer.end(), [](int a, int b) { return a < b; }); 
    for (int j = 0; j < atoms_per_domain; j++) {
      aligned_group_ptr[(i * atoms_per_domain) + j] = buffer[j];
    }
  }
  
  // Check for borders between groups
  std::vector<bool> touch_table(n_partners * n_partners, false);
  for (int i = 0; i < n_partners; i++) {
    for (int j = 0; j < atoms_per_domain; j++) {
      const int atom_ij = aligned_group_ptr[(i * atoms_per_domain) + j];
      for (int k = nbk.nb12_bounds[atom_ij]; k < nbk.nb12_bounds[atom_ij + 1]; k++) {
        const int atom_k = nbk.nb12x[k];
        if (domain_assignment_ptr[atom_k] < 0) {
          for (int m = nbk.nb12_bounds[atom_k]; m < nbk.nb12_bounds[atom_k + 1]; m++) {
            const int atom_m = nbk.nb12x[m];
            const int other_domain = domain_assignment_ptr[atom_m];
            if (other_domain >= 0) {
              touch_table[(other_domain * n_partners) + i] = true;
              touch_table[(i * n_partners) + other_domain] = true;
            }
          }
        }
      }
    }
  }

  // Record the equivalent groups
  std::vector<bool> group_assigned(n_partners, false);
  std::vector<int> swappable_domains;
  swappable_domains.reserve(n_partners);
  for (int i = 0; i < n_partners; i++) {
    if (group_assigned[i]) {
      continue;
    }
    group_assigned[i] = true;
    swappable_domains.resize(1);
    swappable_domains[0] = i;
    int n_swappable = 1;
    for (int j = i + 1; j < n_partners; j++) {
      if (touch_table[(i * n_partners) + j]) {
        group_assigned[j] = true;
        swappable_domains.push_back(j);
        n_swappable++;
      }
    }
    if (n_swappable == 1) {
      rtWarn("A group of atoms was determined to be symmetric only with itself.",
             "AtomEquivalence", "drawEquivalentGroups");
      freeForAllSymmetryGroup(*aligned_groups, atoms_per_domain, swappable_domains);
    }
    else if (n_swappable == n_partners) {
      freeForAllSymmetryGroup(*aligned_groups, atoms_per_domain, swappable_domains);
    }
    else {

      // Determine whether this is a free-for-all group or some sort of rotary arrangement
      bool ffa = true;
      for (int j = 0; j < n_swappable; j++) {
        const int jdom = swappable_domains[j];
        for (int k = j + 1; k < n_swappable; k++) {
          const int kdom = swappable_domains[k];
          ffa = (ffa && touch_table[(jdom * n_partners) + kdom]);
        }
      }
      if (ffa) {
        freeForAllSymmetryGroup(*aligned_groups, atoms_per_domain, swappable_domains);
      }
      else {
        rotarySymmetryGroup(*aligned_groups, atoms_per_domain, swappable_domains, touch_table);
      }
    }
  }
  
  // Reset the allowed atoms array, erasing any entries containing the special code.
  for (int i = 0; i < n_subset; i++) {
    allowed_atoms_ptr[subset_atoms[i]] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::addSymmetryGroup(const std::vector<int> &all_domains, const int domain_size,
                                       const std::vector<int> &selected_domains,
                                       const EquivalenceSwap plan) {
  const int n_domains = selected_domains.size();

  // Check that this group has not already been defined
  for (int i = 0; i < group_count; i++) {
    if (group_sizes[i] == domain_size && group_orders[i] == n_domains) {
      std::vector<bool> taken(n_domains, false);
      std::vector<bool> matched(n_domains, false);
      for (int j = 0; j < n_domains; j++) {
        if (taken[j]) {
          continue;
        }
        const int known_domain_start = group_bounds[i] + (domain_size * j);
        for (int k = 0; k < n_domains; k++) {
          const int tba_domain_start = domain_size * selected_domains[k];
          bool match = true;
          for (int m = 0; m < domain_size; m++) {
            match = (match &&
                     group_atoms[known_domain_start + m] == all_domains[tba_domain_start + m]);
          }
          if (match) {
            taken[j] = true;
            matched[k] = true;
          }
        }
      }
      bool all_matched = true;
      for (int j = 0; j < n_domains; j++) {
        all_matched = (all_matched && matched[j]);
      }
      if (all_matched) {
        return;
      }
    }
  }

  // Add this symmetry group
  const int atom_limit = group_bounds[group_count] + (n_domains * domain_size);
  group_atoms.resize(atom_limit);
  int gpos = group_bounds[group_count];
  std::vector<int> buffer(domain_size);
  for (int i = 0; i < n_domains; i++) {
    const int adpos = selected_domains[i] * domain_size;
    for (int j = 0; j < domain_size; j++) {
      group_atoms[gpos + j] = all_domains[adpos + j];
    }
    gpos += domain_size;
  }
  group_sizes.push_back(domain_size);
  group_orders.push_back(n_domains);
  group_bounds.push_back(atom_limit);
  group_rules.push_back(plan);
  group_count += 1;
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::freeForAllSymmetryGroup(const std::vector<int> &all_domains,
                                              const int domain_size,
                                              const std::vector<int> &selected_domains) {
  addSymmetryGroup(all_domains, domain_size, selected_domains, EquivalenceSwap::FREE_FOR_ALL);
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::rotarySymmetryGroup(const std::vector<int> &all_domains,
                                          const int domain_size,
                                          const std::vector<int> &selected_domains,
                                          const std::vector<bool> &touch_table) {

  // Rotary symmetry is recognized for four or more swappable groups A, B, C, D, ... such that
  // A <-> B and B <-> C, but A <-!-> C.  The chain must close, but this requirement is implicitly
  // satisfied or the groups would not be symmetry-equivalent in the first place.  However, if the
  // chain-link connections from one group to the next are not found, use the approximation that
  // all groups are interchangeable.  This is a fallback onto the notion that correlations among
  // particles will be self-enforcing, especially when there are already contiguous groups of
  // particles to compare.
  const int n_domains = selected_domains.size();
  const int total_domains = all_domains.size() / domain_size;
  std::vector<bool> coverage(total_domains, false);
  int idom = selected_domains[0];
  for (int i = 0; i < n_domains; i++) {
    coverage[idom] = true;

    // Find the number of connections
    int n_connect = 0;
    for (int j = 0; j < total_domains; j++) {
      n_connect += (j != i && touch_table[(idom * total_domains) + j]);
    }
    if (n_connect != 2) {
      addSymmetryGroup(all_domains, domain_size, selected_domains, EquivalenceSwap::FREE_FOR_ALL);
      return;
    }
    
    // Find the next link
    for (int j = 0; j < total_domains; j++) {
      if (j != i && coverage[j] == false) {
        coverage[j] = true;
        idom = j;
        break;
      }
    }
  }
  bool all_links_found = true;
  for (int i = 0; i < n_domains; i++) {
    all_links_found = (all_links_found && coverage[selected_domains[i]]);
  }
  if (all_links_found) {
    addSymmetryGroup(all_domains, domain_size, selected_domains, EquivalenceSwap::ROTARY);
  }
  else {
    addSymmetryGroup(all_domains, domain_size, selected_domains, EquivalenceSwap::FREE_FOR_ALL);
  }
}

//-------------------------------------------------------------------------------------------------
void AtomEquivalence::findDependencies() {
  std::vector<int2> minmax(group_count, { ag_pointer->getAtomCount(), -1 });
  for (int i = 0; i < group_count; i++) {
    const int jlim = group_bounds[i + 1];
    for (int j = group_bounds[i]; j < jlim; j++) {
      minmax[i].x = std::min(minmax[i].x, group_atoms[j]);
      minmax[i].y = std::max(minmax[i].y, group_atoms[j]);
    }
  }
  dependency_bounds.resize(1);
  dependency_bounds[0] = 0;
  for (int i = 0; i < group_count; i++) {
    const int imin = minmax[i].x;
    const int imax = minmax[i].y;
    const int icnt = group_sizes[i];
    const int iord = group_orders[i];
    for (int j = 0; j < group_count; j++) {
      if (group_sizes[j] < icnt && minmax[j].x >= imin && minmax[j].y <= imax) {
        const int jcnt = group_sizes[j];
        const int jord = group_orders[j];

        // Attempt to find each domain of the jth symmetry group within a domain of the ith group
        bool all_enclosed = true;
        int k = 0;
        while (k < jord && all_enclosed) {
          const int jgroup_llim = group_bounds[j] + (k * jcnt);
          const int jgroup_hlim = jgroup_llim + jcnt;
          bool domain_enclosed = false;
          int m = 0;
          while (m < iord && domain_enclosed == false) {
            const int igroup_llim = group_bounds[i] + (m * icnt);
            const int igroup_hlim = igroup_llim + icnt;
            int itrack = igroup_llim;
            bool all_atoms_enclosed = true;
            for (int n = jgroup_llim; n < jgroup_hlim; n++) {
              while (group_atoms[itrack] < group_atoms[n] && itrack < igroup_hlim) {
                itrack++;
              }
              if (itrack == igroup_hlim || group_atoms[itrack] != group_atoms[n]) {
                all_atoms_enclosed = false;
                n = jgroup_hlim;
                break;
              }
            }
            domain_enclosed = (domain_enclosed || all_atoms_enclosed);
            m++;
          }
          all_enclosed = (all_enclosed && domain_enclosed);
          k++;
        }
        if (all_enclosed) {
          dependencies.push_back(j);
        }
      }
    }
    dependency_bounds.push_back(dependencies.size());
  }
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph *ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b, std::vector<int> *a_idx_tree,
                         std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree,
                         std::vector<int> *b_zn_tree, std::vector<double> *a_fc_tree,
                         std::vector<double> *b_fc_tree, std::vector<double> *a_fe_tree,
                         std::vector<double> *b_fe_tree, std::vector<ullint> *a_ri_tree,
                         std::vector<ullint> *b_ri_tree,
                         std::vector<ChiralOrientation> *a_ch_tree,
                         std::vector<ChiralOrientation> *b_ch_tree,
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage) {

  // If the two atoms have the same index, return true.
  if (atom_a == atom_b) {
    return true;
  }
  
  // Check that both atoms are part of the same molecule
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const int mol_ab = cdk.mol_home[atom_a];
  if (mol_ab != cdk.mol_home[atom_b]) {
    rtErr("Atom indices " + std::to_string(atom_a) + " and " + std::to_string(atom_b) +
          " reside in molecules " + std::to_string(cdk.mol_home[atom_a]) + " and " +
          std::to_string(cdk.mol_home[atom_b]) + ", respectively.  Bonding patterns can only "
          "match for atoms in the same molecule.", "matchBondingPatterns");
  }
  
  // Get the non-bonded abstract to track connectivity in the topology
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();

  // Set pointers for various storage arrays.  They are passed in to avoid needing to re-allocate
  // large amounts of memory for finding equivalencies in smaller and smaller subdivisions of the
  // system.
  int* a_idx_tptr = a_idx_tree->data();
  int* b_idx_tptr = b_idx_tree->data();
  int* a_zn_tptr = a_zn_tree->data();
  int* b_zn_tptr = b_zn_tree->data();
  double* a_fc_tptr = a_fc_tree->data();
  double* b_fc_tptr = b_fc_tree->data();
  double* a_fe_tptr = a_fe_tree->data();
  double* b_fe_tptr = b_fe_tree->data();
  ullint* a_ri_tptr = a_ri_tree->data();
  ullint* b_ri_tptr = b_ri_tree->data();
  ChiralOrientation* a_ch_tptr = a_ch_tree->data();
  ChiralOrientation* b_ch_tptr = b_ch_tree->data();
  int* a_coverage_ptr = a_coverage->data();
  int* b_coverage_ptr = b_coverage->data();
  
  // Initialize trees for each atom
  const int mol_natom = cdk.mol_limits[mol_ab + 1] - cdk.mol_limits[mol_ab];
  for (int k = cdk.mol_limits[mol_ab]; k < cdk.mol_limits[mol_ab + 1]; k++) {
    a_coverage_ptr[cdk.mol_contents[k]] = 0;
    b_coverage_ptr[cdk.mol_contents[k]] = 0;
  }
  a_coverage_ptr[atom_a] = 1;
  b_coverage_ptr[atom_b] = 1;
  a_idx_tptr[0] = atom_a;
  b_idx_tptr[0] = atom_b;
  bool add_to_trees = true;
  bool previous_layer_same = false;
  int last_layer_start = 0;
  int last_layer_end = 1;
  int n_atree_atoms = 1;
  int n_btree_atoms = 1;
  while (add_to_trees) {
    for (int i = last_layer_start; i < last_layer_end; i++) {
      const int a_tree_atom = a_idx_tptr[i];
      for (int j = nbk.nb12_bounds[a_tree_atom]; j < nbk.nb12_bounds[a_tree_atom + 1]; j++) {
        if (a_coverage_ptr[nbk.nb12x[j]] == 0) {
          const size_t next_atom = nbk.nb12x[j];
          a_idx_tptr[n_atree_atoms] = next_atom;
          a_zn_tptr[n_atree_atoms] = cdk.z_numbers[next_atom];
          a_fc_tptr[n_atree_atoms] = formal_charges[next_atom];
          a_fe_tptr[n_atree_atoms] = free_electrons[next_atom];
          a_ri_tptr[n_atree_atoms] = ring_inclusion[next_atom];
          a_ch_tptr[n_atree_atoms] = chiralities[next_atom];
          a_coverage_ptr[next_atom] = 1;
          n_atree_atoms++;
        }
      }
      const int b_tree_atom = b_idx_tptr[i];
      for (int j = nbk.nb12_bounds[b_tree_atom]; j < nbk.nb12_bounds[b_tree_atom + 1]; j++) {
        if (b_coverage_ptr[nbk.nb12x[j]] == 0) {
          const size_t next_atom = nbk.nb12x[j];
          b_idx_tptr[n_btree_atoms] = next_atom;
          b_zn_tptr[n_btree_atoms] = cdk.z_numbers[next_atom];
          b_fc_tptr[n_btree_atoms] = formal_charges[next_atom];
          b_fe_tptr[n_btree_atoms] = free_electrons[next_atom];
          b_ri_tptr[n_btree_atoms] = ring_inclusion[next_atom];
          b_ch_tptr[n_btree_atoms] = chiralities[next_atom];
          b_coverage_ptr[next_atom] = 1;
          n_btree_atoms++;
        }
      }
    }
    if (n_atree_atoms != n_btree_atoms) {
      return false;
    }
    
    // If the previous check passed, both trees are still the same size.  Update the layer bounds,
    // then check sums of the various Z numbers, free electrons, ring inclusions, and formal
    // charges.  Even if there are atoms in 60+ membered rings, the sums of unsigned long long ints
    // should still overflow in the same manner.
    last_layer_start = last_layer_end;
    last_layer_end = n_atree_atoms;
    int zn_sum_a = 0;
    int zn_sum_b = 0;
    int ch_sum_a = 0;
    int ch_sum_b = 0;
    double fc_sum_a = 0.0;
    double fc_sum_b = 0.0;
    double fe_sum_a = 0.0;
    double fe_sum_b = 0.0;
    for (int i = last_layer_start; i < last_layer_end; i++) {
      zn_sum_a += a_zn_tptr[i];
      zn_sum_b += b_zn_tptr[i];
    }
    if (zn_sum_a != zn_sum_b) {
      return false; 
    }
    for (int i = last_layer_start; i < last_layer_end; i++) {
      fc_sum_a += a_fc_tptr[i];
      fc_sum_b += b_fc_tptr[i];
      fe_sum_a += a_fe_tptr[i];
      fe_sum_b += b_fe_tptr[i];
      ch_sum_a += static_cast<int>(a_ch_tptr[i]);
      ch_sum_b += static_cast<int>(b_ch_tptr[i]);
    }
    if (fabs(fc_sum_a - fc_sum_b) > 1.0e-6 || fabs(fe_sum_a - fe_sum_b) > 1.0e-6 ||
        ch_sum_a != ch_sum_b) {
      return false; 
    }

    // Check to see whether the trees are adding atoms in the same order, and can thus be expected
    // to do so forever more.  While the sum of integer atom indices in a very, very large system
    // could, in principle, overflow, the results can still be expected to be equal for two
    // identical series of numbers, and int addition is faster than int => long long int conversion
    // with long long int addition.
    if (sum<int>(&a_idx_tptr[last_layer_start], last_layer_end - last_layer_start) ==
        sum<int>(&b_idx_tptr[last_layer_start], last_layer_end - last_layer_start)) {

      // Because all entries in each array will be unique, just check that each entry of one is
      // present in the other.
      bool all_covered = true;
      int search_start = last_layer_start;
      for (int i = last_layer_start; i < last_layer_end; i++) {
        const int a_tree_atom = a_idx_tptr[i];
        int j = search_start;
        bool found = false;
        while (j < last_layer_end && (! found)) {
          found = (a_tree_atom == b_idx_tptr[j]);
          search_start += (found && j == search_start);
          j++;
        }
        all_covered = (all_covered && found);
      }

      // If two successive layers are identical, then it can be inferred that subsequent layers
      // will only continue to add the same atoms with the same properties and connections.
      if (previous_layer_same) {
        if (all_covered) {
          return true;
        }
        else {
          previous_layer_same = false;
        }
      }
      else {
        previous_layer_same = all_covered;
      }
    }

    // Did the trees grow?
    add_to_trees = (last_layer_end > last_layer_start);
  }

  // If the trees stop growing, the entirety of the molecule is covered.  If nothing has indicated
  // a mismatch by that point, the bonding patterns must be identical.
  return true;
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b, std::vector<int> *a_idx_tree,
                         std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree,
                         std::vector<int> *b_zn_tree, std::vector<double> *a_fc_tree,
                         std::vector<double> *b_fc_tree, std::vector<double> *a_fe_tree,
                         std::vector<double> *b_fe_tree, std::vector<ullint> *a_ri_tree,
                         std::vector<ullint> *b_ri_tree,
                         std::vector<ChiralOrientation> *a_ch_tree,
                         std::vector<ChiralOrientation> *b_ch_tree,
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage) {
  return matchBondingPattern(ag.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                             chiralities, atom_a, atom_b, a_idx_tree, b_idx_tree, a_zn_tree,
                             b_zn_tree, a_fc_tree, b_fc_tree, a_fe_tree, b_fe_tree, a_ri_tree,
                             b_ri_tree, a_ch_tree, b_ch_tree, a_coverage, b_coverage);
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b) {
  const int natom = ag.getAtomCount();
  std::vector<int> a_idx_tree(natom), b_idx_tree(natom), a_zn_tree(natom), b_zn_tree(natom);
  std::vector<double> a_fc_tree(natom), b_fc_tree(natom), a_fe_tree(natom), b_fe_tree(natom);
  std::vector<ullint> a_ri_tree(natom), b_ri_tree(natom);
  std::vector<ChiralOrientation> a_ch_tree(natom), b_ch_tree(natom);
  std::vector<int> a_coverage(natom), b_coverage(natom);
  return matchBondingPattern(ag.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                             chiralities, atom_a, atom_b, &a_idx_tree, &b_idx_tree, &a_zn_tree,
                             &b_zn_tree, &a_fc_tree, &b_fc_tree, &a_fe_tree, &b_fe_tree,
                             &a_ri_tree, &b_ri_tree, &a_ch_tree, &b_ch_tree, &a_coverage,
                             &b_coverage);
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph &ag, const ChemicalFeatures &chemfe, const int atom_i,
                         const int atom_j) {
  return matchBondingPattern(ag, chemfe.getFormalCharges(), chemfe.getFreeElectrons(),
                             chemfe.getRingInclusion(), chemfe.getAtomChirality(), atom_i, atom_j);
}

} // namespace chemistry
} // namespace stormm
