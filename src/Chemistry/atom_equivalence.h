// -*-c++-*-
#ifndef STORMM_ATOM_EQUIVALENCE_H
#define STORMM_ATOM_EQUIVALENCE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "UnitTesting/stopwatch.h"
#include "chemistry_enumerators.h"

namespace stormm {
namespace chemistry {

using topology::AtomGraph;
using testing::StopWatch;
  
class AtomEquivalence {
public:

  /// \brief The constructor relies on a topology pointer plus arrays of formal charges, bond
  ///        orders, and ring inclusion as will have been derived by a ChemicalFeatures object.
  ///        Due to this object possibly feeding back into the ChemicalFeatures, however, its
  ///        constructor only takes the distilled data.
  ///
  /// \param ag_in           Topology for the system of interest (a pointer will be retained by
  ///                        the object)
  /// \param formal_charges  Array of formal charges for all atoms in the entire topology
  /// \param free_electrons  Array of free electron content for all atoms in the entire topology
  /// \param ring_inclusion  Array of ring inclusion marks for all atoms in the entire topology
  /// \param chiralities     Array of chiral statuses for all atoms in the entire topology
  /// \param timer           Wall time tracker to monitor time spent in various stages of the
  ///                        calculation
  /// \{
  AtomEquivalence();
  
  AtomEquivalence(const AtomGraph &ag_in, const std::vector<double> &formal_charges,
                  const std::vector<double> &free_electrons,
                  const std::vector<ullint> &ring_inclusion,
                  const std::vector<int> &chiral_centers, StopWatch *timer = nullptr);
  /// \}

  /// The default copy and move constructors, as well as copy and move assignment operators, will
  /// be adequate for this object which contains only Standard Template Library components and no
  /// const member variables.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object to which the current one shall be assigned
  /// \{
  AtomEquivalence(const AtomEquivalence &original) = default;
  AtomEquivalence(AtomEquivalence &&original) = default;
  AtomEquivalence& operator=(const AtomEquivalence &original) = default;
  AtomEquivalence& operator=(AtomEquivalence &&original) = default;
  /// \}
  
  /// \brief Get the number of atom equivalence groups.
  int getGroupCount() const;

  /// \brief Get the atoms of one group as an independent Standard Template Library vector.
  ///
  /// \param group_index  The group of interest
  std::vector<int> getGroup(int group_index) const;

  /// \brief Get a pointer to the atom indices of a particular group.
  ///
  /// \param group_index  The group of interest  
  const int* getGroupPointer(int group_index) const;

  /// \brief Get the size of a particular group.
  ///
  /// \param group_index  The group of interest  
  int getGroupSize(int group_index) const;

  /// \brief Get the order of a particular group
  ///
  /// \param group_index  The group of interest  
  int getGroupOrder(int group_index) const;
  
  /// \brief Get a list of a group's dependencies (other groups that the group in question contains
  ///        completely).  When evaluating dependent groups, the procedure must be to work from the
  ///        outside in, swapping 
  ///
  /// \param group_index  The group of interest  
  std::vector<int> getGroupDependencies(int group_index) const;

  /// \brief Get a const reference to the list of atoms outside any symmetry group.
  const std::vector<int>& getAsymmetricAtoms() const;

  /// \brief Get a const pointer to the topology for which these equivalencies apply.
  const AtomGraph* getTopologyPointer() const;

private:
  int group_count;                     ///< The number of equivalent atom groups
  std::vector<int> group_atoms;        ///< Indices of atoms in each equivalent group.  For a group
                                       ///<   with K atoms and N-fold symmetry, the equivalent
                                       ///<   atoms of the first group instance are listed
                                       ///<   [ A(1,i), A(1,i+1), ..., A(1,k) ] followed by those
                                       ///<   of the second instance [ A(2,i), ..., A(2,k) ] and
                                       ///<   so on up to [ A(N,i), ..., A(N,k) ].
  std::vector<int> group_sizes;        ///< The number of atoms in each equivalent atom group
  std::vector<int> group_orders;       ///< The number of symmetric instances of each group
  std::vector<int> group_bounds;       ///< Bounds array for group_atoms above
  std::vector<int> dependencies;       ///< Indices of dependent groups for each group in the
                                       ///<   molecule.  If group 4 is completely subsumed within
                                       ///<   group 3, then the dependencies of group 3 will
                                       ///<   include group 4.
  std::vector<int> dependency_bounds;  ///< Bounds array for dependencies above, group_count + 1 in
                                       ///<   length
  std::vector<int> asymmetric_atoms;   ///< List of atoms that participate in no symmetry group

  /// The manner by which each group can exchange the coordinates of its symmetric subunits
  std::vector<EquivalenceSwap> group_rules;

  /// Pointer to the topology this is describing
  AtomGraph *ag_pointer;

  /// \brief Find equivalent atoms from within a subset of the topology, then call the
  ///        drawEquivalentGroups() function below.
  ///
  /// \param ag                System topology
  /// \param subset_atoms      Topological indices for the subset of atoms that this call will look
  ///                          at, i.e. { 4, 5, 6, 7, 21, 22, 23, 24 }.
  /// \param map_to_subset     Mapping of each topological atom its place in the subset, i.e. for
  ///                          the example above, { x, x, x, x, 0, 1, 2, 3, x, x, x, ..., 4, 5,
  ///                          6, 7 }.  Not all indices of this array can be trusted at all times.
  ///                          If the atom is not actually in the subset, then the map will not
  ///                          contain a valid index into the subset_atoms array.
  /// \param formal_charges    Array of formal charges for all atoms in the entire topology
  /// \param free_electrons    Array of free electron content for all atoms in the entire topology
  /// \param ring_inclusion    Array of ring inclusion marks for all atoms in the entire topology
  /// \param chiralities       Array of chiral statuses for all atoms in the entire topology
  /// \param domain_coverage   Array (with space for all topological atom indices) for recording
  ///                          the extent to which various symmetry-related groups expanding
  ///                          outwards from each partner have covered the available atoms
  /// \param allowed_atoms     Array (also with space for all topological atom indices) for
  ///                          recording which atoms of the topology are to get consideration for
  ///                          building the symmetry-equivalent partner domains that will form
  ///                          the basis for symmetry groups
  /// \param candidate_hopper  Space reserved for atoms that might be found to contribute to the
  ///                          next layer of each symmetry-related group
  /// \param jumbled_groups    List of non-overlapping atoms making up each symmetry-equivalent
  ///                          group.  This list is filtered by passing through the
  ///                          candidate_hopper array, eliminating atoms that result in
  ///                          double-coverage.
  /// \param aligned_groups    Re-arrangement of jumbled groups that places atoms from each
  ///                          symmetry-equivalent domain, layer by layer working outward from the
  ///                          symmetric partner atom, into a contiguous stretch of memory
  /// \param a_idx_tree        Tree of topological atom indices branching out from the first (A)
  ///                          atom of a bond pattern comparison.  This and subsequent arrays are
  ///                          pre-allocated for the call to findEquivalentAtoms() so that they
  ///                          can be passed down through recursive calls, onward to
  ///                          drawEquivalentGroups(), and finally used in calls to the
  ///                          matchBondingPattern() function.  Re-allocation for each recursive
  ///                          call is then not needed.
  /// \param b_idx_tree        Tree of topological atom indices branching out from the second (B)
  ///                          atom of a bond pattern comparison
  /// \param a_zn_tree         Tree of atomic numbers found in a bond pattern comparison when
  ///                          branching from the first atom
  /// \param b_zn_tree         Tree of atomic numbers found in a bond pattern comparison when
  ///                          branching from the second atom
  /// \param a_fc_tree         Tree of formal charges branching from the first atom in a bond
  ///                          pattern comparison
  /// \param b_fc_tree         Tree of formal charges branching from the second atom in a bond
  ///                          pattern comparison
  /// \param a_fe_tree         Tree of free electron contents branching from the first atom in a
  ///                          bond pattern comparison
  /// \param b_fe_tree         Tree of free electron contents branching from the second atom in a
  ///                          bond pattern comparison
  /// \param a_ri_tree         Tree of ring inclusions branching from the first atom in a bond
  ///                          pattern comparison
  /// \param b_ri_tree         Tree of ring inclusions branching from the second atom in a bond
  ///                          pattern comparison
  /// \param a_ch_tree         Tree of chiral designations branching from the first atom in a bond
  ///                          pattern comparison
  /// \param b_ch_tree         Tree of chiral designations branching from the second atom in a bond
  ///                          pattern comparison
  /// \param a_coverage        Coverage of the first atom's tree in a bond pattern comparison
  /// \param b_coverage        Coverage of the second atom's tree in a bond pattern comparison
  void findEquivalentAtoms(const AtomGraph &ag, const std::vector<int> &subset_atoms,
                           std::vector<int> *map_to_subset,
                           const std::vector<double> &formal_charges,
                           const std::vector<double> &free_electrons,
                           const std::vector<ullint> &ring_inclusion,
                           const std::vector<ChiralOrientation> &chiralities,
                           std::vector<int> *domain_coverage, std::vector<int> *allowed_atoms,
                           std::vector<int> *candidate_hopper, std::vector<int> *domain_assignment,
                           std::vector<int> *jumbled_groups, std::vector<int> *aligned_groups,
                           std::vector<int> *layer_bounds, std::vector<int> *a_idx_tree,
                           std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree,
                           std::vector<int> *b_zn_tree, std::vector<double> *a_fc_tree,
                           std::vector<double> *b_fc_tree, std::vector<double> *a_fe_tree,
                           std::vector<double> *b_fe_tree, std::vector<ullint> *a_ri_tree,
                           std::vector<ullint> *b_ri_tree,
                           std::vector<ChiralOrientation> *a_ch_tree,
                           std::vector<ChiralOrientation> *b_ch_tree,
                           std::vector<int> *a_coverage, std::vector<int> *b_coverage);

  /// \brief Draw the equivalent groups for a molecule, based on a trusted list of partners each
  ///        representing a distinct collection of symmetry-related atoms.  The strategy is to
  ///        begin with a list of atoms that have been found to have equivalent bonding patterns:
  ///        starting from each atom and exploring outwards throughout the molecule, the elements,
  ///        bonds, and electronic structures of each atom are identical.  The search expands
  ///        outwards from these atoms, accumulating symmetry-equivalent "domains" for each partner
  ///        atom, until further expansion would cause the domains to overlap.  If two partners'
  ///        domains both attempt to include the same atom at one time, it will become part of
  ///        neither domain.  Each call to this function may make more than one group of equivalent
  ///        atoms given the set of partners: only those domains which come into contact with one
  ///        another will be considered interchangeable.  If there are symmetric partner atoms A,
  ///        B, C, and D, with A and B at one end of the molecule and C and D at the other, such
  ///        the domains of A and B or those of C and D might collide, and be able to expand no
  ///        further, well before the domains of A and C or those of B and D ever come into
  ///        contact.  In such a case, the domains of A and B would be considered interchangeable,
  ///        and those of C and D as well, but A and B would not be allowed to exchange with
  ///        C or D.
  ///
  /// Argument descriptions follow from findEquivalentAtoms() above (many of the arrays are meant
  /// to be passed on to calls there), with the additions of:
  ///
  /// \param partners          A list of equivalent atoms with one entry spanning each equivalent
  ///                          group
  /// \param subset_atoms      A list of atoms within which the equivalence groups can be drawn
  /// \param layer_bounds      Working array for the bounds of each layer forming the domains of
  ///                          equivalent atom groups
  void drawEquivalentGroups(const std::vector<int> &partners,
                            const std::vector<int> &subset_atoms,
                            const std::vector<double> &formal_charges,
                            const std::vector<double> &free_electrons,
                            const std::vector<ullint> &ring_inclusion,
                            const std::vector<ChiralOrientation> &chiralities,
                            std::vector<int> *domain_coverage, std::vector<int> *allowed_atoms,
                            std::vector<int> *candidate_hopper,
                            std::vector<int> *domain_assignment, std::vector<int> *jumbled_groups,
                            std::vector<int> *aligned_groups, std::vector<int> *layer_bounds,
                            std::vector<int> *a_idx_tree, std::vector<int> *b_idx_tree,
                            std::vector<int> *a_zn_tree, std::vector<int> *b_zn_tree,
                            std::vector<double> *a_fc_tree, std::vector<double> *b_fc_tree,
                            std::vector<double> *a_fe_tree, std::vector<double> *b_fe_tree,
                            std::vector<ullint> *a_ri_tree, std::vector<ullint> *b_ri_tree,
                            std::vector<ChiralOrientation> *a_ch_tree,
                            std::vector<ChiralOrientation> *b_ch_tree,
                            std::vector<int> *a_coverage, std::vector<int> *b_coverage);

  /// \brief Add a selection of the symmetry-equivalent domains as an interchangeable collection
  ///        of atom groups, with a plan for how to do the swapping.
  ///
  /// \param all_domains       The collection of atom topological indices for all symmetry-related
  ///                          domains of interest.  This may include more domains than just the
  ///                          selected few which can be swapped.  If a dumbell-shaped molecule has
  ///                          three symmtry-equivalent moieties on either end, the all_domains
  ///                          array will contain the atoms of all six moieties, but only three
  ///                          will be listed in selected_domains.
  /// \param domain_size       The size of each domain.  Each domain is presented as a contiguous
  ///                          stream of atm indices in all_domains.
  /// \param selected_domains  The domains which are free to swap
  /// \param plan              Prescribed method for interchanging domains
  void addSymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                        const std::vector<int> &selected_domains, EquivalenceSwap plan);
  
  /// \brief Construct a symmetry group in which free interchanges of all of its subunits are
  ///        permissible.  Descriptions of parameters follow from addSymmetryGroup() above.
  void freeForAllSymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                               const std::vector<int> &selected_domains);

  /// \brief Construct a symmetry group in which interchanges of the subunits ("domains)" are
  ///        permitted in the manner K -> K+1, K+1 -> K+2, K -> K-1, etc., with N total
  ///        permutations for a ring of N subunits.  This function will first verify that there is
  ///        a ring structure defined by the connectivity between subunits, and if not the group
  ///        will be approximated as a "free for all" symmetry group.  Descriptions of parameters
  ///        follow from addSymmetryGroup(), with the addition of:
  ///
  /// \param touch_table  Table indicating whether any pair of domains / subunits touch one another
  void rotarySymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                           const std::vector<int> &selected_domains,
                           const std::vector<bool> &touch_table);

  /// \brief Analyze the detected symmetry groups and determine which ones are subsumed within
  ///        others.  When making plans to create 
  void findDependencies();
};

/// \brief Beginning with two distinct atoms in a molecule within a topology, proceed throughout
///        the molecular bonding pattern verifying that the atoms one encounters when stepping
///        outward from each atom have identical properties.
///
/// Overloaded:
///   - Accept a collection of pre-allocated vectors to use a scratch space
///   - Allocate temporary storage space for all work
///
/// \param ag              System topology
/// \param formal_charges  Array of formal charges for all atoms in the entire topology
/// \param bond_orders     Array of bond orders for all bonds in the topology
/// \param free_electrons  Array of ree electron content for all atoms in the entire topology
/// \param ring_inclusion  Array of ring inclusion marks for all atoms in the entire topology
/// \param chiralities     Array of chiral statuses for all atoms in the entire topology
/// \param atom_a          The first atom to compare   
/// \param atom_b          The second atom to compare
/// \param a_idx_tree      Tree of topological atom indices branching out from the first (A) atom
/// \param b_idx_tree      Tree of topological atom indices branching out from the second (B) atom
/// \param a_zn_tree       Tree of atomic numbers found when branching from the first atom
/// \param b_zn_tree       Tree of atomic numbers found when branching from the second atom
/// \param a_fc_tree       Tree of formal charges branching from the first atom
/// \param b_fc_tree       Tree of formal charges branching from the second atom
/// \param a_fe_tree       Tree of free electron contents branching from the first atom
/// \param b_fe_tree       Tree of free electron contents branching from the second atom
/// \param a_ri_tree       Tree of ring inclusions branching from the first atom
/// \param b_ri_tree       Tree of ring inclusions branching from the second atom
/// \param a_ch_tree       Tree of chiral designations branching from the first atom
/// \param b_ch_tree       Tree of chiral designations branching from the second atom
/// \param a_coverage      Coverage of the first atom's tree
/// \param b_coverage      Coverage of the second atom's tree
/// \{
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
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage);

bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, int atom_i,
                         int atom_j);
/// \}

} // namespace chemistry
} // namespace stormm

#endif
