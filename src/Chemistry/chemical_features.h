// -*-c++-*-
#ifndef OMNI_CHEMICAL_FEATURES_H
#define OMNI_CHEMICAL_FEATURES_H

#include "Accelerator/hybrid.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "chemistry_enumerators.h"

namespace omni {
namespace chemistry {

using card::Hybrid;
using card::HybridTargetLevel;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;

/// \brief An unguarded struct for delivering information about a rotatable group of atoms.
struct RotatorGroup {
  int root_atom;                     ///< The root atom helping to define the rotatable bond axis
  int pivot_atom;                    ///< Pivot atom completing the rotatable bond axis (more
                                     ///<   proximal to atoms that will turn than the root atom)
  std::vector<int> rotatable_atoms;  ///< List of all atoms that turn as a consequence of twisting
                                     ///<   about the rotatable bond axis
};
  
/// \brief A struct to serve as a tracker of progress through a molecule in the search for loops.
///        Proceeding forward in the search, every link will have come from one and only one
///        previous link, but could go in multiple directions thereafter.
class BondedNode {
public:
  
  /// \brief Basic constructor simply initializes the members to blank values
  BondedNode();

  /// \brief Take a pre-existing std::vector<int> in order to set the object's own next_links
  ///        pointer.  This is in the interest of not repeatedly allocating and freeing many tiny
  ///        std::vector<int> objects as BondedNodes are created and destroyed to grow a chain of
  ///        atoms.
  ///
  /// \param vi            Vector of integers in which this BondedNode will store some of its data
  /// \param pos           Indicator of the space in the storage vector allocate to this object
  /// \param max_branches  The maximum number of branches that this chain link shall support
  void setBranchPointer(std::vector<int> *vi, size_t pos, size_t max_branches);

  /// \brief Set this as the next link in the chain, or as the chain initiator.
  ///
  /// \param previous_in    Index of the previous BondedNode leading to this one, or -1 if this is
  ///                       intended to be the chain initiator
  /// \param current_atom   Topological index of the current atom, which this BondedNode describes
  /// \param current_layer  The current layer in the growing tree
  /// \param nbk            Nonbonded interactions abstract for the overarching topology (needed
  ///                       for list of 1:2 exclusions, which indicate bonds)
  /// \param cdk            Atom, residue, and molecule details taken from the original topology
  void addToTree(int previous_in, int current_atom, int current_layer,
                 const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk);

  /// \brief Add the bond order between this atom and its previous atom
  ///
  /// \param vk           Valence interactions abstract from the original topology
  /// \param bond_orders  Orders of bonds determined by Lewis structure drawing
  void addBondOrder(const ValenceKit<double> &vk, const Hybrid<double> &bond_orders);
  
  /// \brief Get the previous atom index
  int getPreviousAtom() const;

  /// \brief Get the central atom for this BondedNode.
  int getAtom() const;

  /// \brief Get the layer of which this BondedNode is a part.
  int getLayer() const;

  /// \brief Get the number of branches associated with this BondedNode's central atom (this is the
  ///        total number of other atoms that it bonds to, less the previous atom in the chain)
  int getBranchCount() const;

  /// \brief Get the topological atom index of a specific branch that can come off of this
  ///        BondedNode.
  ///
  /// \param index  The branch index of interest (could be any topological atom)
  int getBranchAtom(int index) const;

  /// \brief Get the branch index of a specific atom in a BondedNode.
  ///
  /// \param search_index  The topological atom index to search for within the node's branches
  int findBranchAtom(int search_index) const;

  /// \brief Get the ring completion status for a particular node
  ///
  /// \param branch_index  The branch index of the partner atom with which this node forms a ring
  ///                      (this index is not a topological index, and must be ascertained in
  ///                      advance by findBranchAtom() above)
  uint getRingCompletion(int branch_index) const;

  /// \brief Get the bond order of the bond between this atom and the previous atom in the tree.
  double getRootBondOrder() const;

  /// \brief Set the ring completion for this object (this relies on no atom having more bonds than
  ///        the number of bits in an unsigned int)
  ///
  /// \param partner_index  Topological index of the atom to which a ring has been formed
  void setRingCompletion(int partner_index);

  /// \brief Wipe clean the ring completion for this object
  void wipeRingCompletion();

private:
  int previous_atom_index;   ///< Topological index of the previous atom
  int atom_index;            ///< Topological index of this atom
  int layer_index;           ///< Index of the layer to which this BondedNode belongs.  The first
                             ///<   atom in the tree is layer zero.
  double root_bond_order;    ///< Order by which this atom bonds to the previous atom (initalized
                             ///<   to zero by default, filled with addBondOrder())
  int branch_count;          ///< The number of branching particles coming off of this atom,
                             ///<   excluding the previous atom
  int* branch_atoms;         ///< Pointer to data array containing indices of all branch atoms
  uint rings_completed;      ///< Bitmask indicating, in the nth bit, whether a ring has been
                             ///<   completed by following the nth branch
};
  
/// \brief An object to store information about chemical motifs: participation in rings, planarity,
///        chirality, aromaticity, conjugation, planarity, and bonds with different rotatability.
class ChemicalFeatures {
public:

  /// \brief The constructor requires a topology and some coordinate set.
  ///
  /// Overloaded:
  ///   - Create a blank object
  ///   - Create with a light CoordinateFrame object
  ///   - Create with a heavy-duty PhaseSpace object holding the coordinates (this will just be
  ///     stripped down to a CoordinateFrame object)
  ///
  /// \param ag_in           Pointer to the system topology.  This topology will not be modified by
  ///                        submitting it to this constructor, but it is needed as a constant
  ///                        pointer so that the object itself can store a valid pointer to the
  ///                        original topology (passing by const reference would not create a valid
  ///                        pointer).
  /// \param cfr             Coordinates of the system
  /// \param ps              Coordinates of the system (a CoordinateFrameReader will be extracted)
  /// \param map_group_in    Indicator of whether to map rotatable groups (this is an O(N^2)
  ///                        algorithm in memory as well as computation, as larger structures will
  ///                        have more rotatable bonds and the entirety of the structure must be
  ///                        traced in relation to each of them)
  /// \param temperature_in  Temperature at which to take Boltzmann weights of different resonance
  ///                        states
  /// \{
  ChemicalFeatures();

  ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrameReader &cfr,
                   MapRotatableGroups map_group_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0);

  ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrame &cf,
                   MapRotatableGroups map_group_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0);

  ChemicalFeatures(const AtomGraph *ag_in, const PhaseSpace &ps,
                   MapRotatableGroups map_group_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0);
  /// \}

  /// \brief Copy and move constructors
  ///
  /// \param original  The ChemicalFeatures object to copy
  /// \{
  ChemicalFeatures(const ChemicalFeatures &original);
  ChemicalFeatures(ChemicalFeatures &&original);
  /// \}
  
  /// \brief Copy assignment and move assignment operators
  ///
  /// \param other     The ChemicalFeatures object to copy (a different name for a better semantic
  ///                  fit in the context of the = sign)
  /// \{
  ChemicalFeatures& operator=(const ChemicalFeatures &other);
  ChemicalFeatures& operator=(ChemicalFeatures &&other);
  /// \}
  
#ifdef OMNI_USE_HPC
  /// \brief Upload data to the GPU.
  void upload();
#endif

  /// \brief Get the number of planar atoms.
  int getPlanarAtomCount() const;

  /// \brief Get the number of rings in the system.
  int getRingCount() const;

  /// \brief Get the number of fused rings in the system.
  int getFusedRingCount() const;

  /// \brief Get the number of malleable, non-aromatic twistable ringsin the system.
  int getMutableRingCount() const;

  /// \brief Get the number of aromatic groups in the system.
  int getAromaticGroupCount() const;

  /// \brief Get the number of polar hydrogens in the system.
  int getPolarHydrogenCount() const;

  /// \brief Get the number of hydrogen bond donors in the system.
  int getHydrogenBondDonorCount() const;

  /// \brief Get the number of hydrogen bond acceptors in the system.
  int getHydrogenBondAcceptorCount() const;

  /// \brief Get the number of chiral centers in the system.
  int getChiralCenterCount() const;

  /// \brief Get the number of rotatable bonds in the system.
  int getRotatableBondCount() const;
  
  /// \brief Return a mask of rings within a given size range for this system.
  ///
  /// \param min_ring_size  The minimum number of atoms in the rings that will be reported
  /// \param max_ring_size  The maximum number of atoms in the rings that will be reported
  std::vector<uint> getRingMask(int min_ring_size, int max_ring_size) const;
  
  /// \brief Return a mask of atomatic atoms in the system.
  ///
  /// \param min_pi_electrons  Minimum number of electrons in the aromatic ring system to report
  /// \param max_pi_electrons  Maximum number of electrons in the aromatic ring system to report
  std::vector<uint> getAromaticMask(int min_pi_electrons, int max_pi_electrons) const;

  /// \brief Get a list of all polar hydrogen atoms.
  std::vector<int> getPolarHydrogenList() const;

  /// \brief Get a list of all hydrogen bond donor atoms.
  std::vector<int> getHydrogenBondDonorList() const;

  /// \brief Get a list of all hydrogen bond acceptor atoms.
  std::vector<int> getHydrogenBondAcceptorList() const;

  /// \brief Get a bit mask of all polar hydrogen atoms in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getPolarHydrogenMask() const;

  /// \brief Get a bit mask of all hydrogen bond donors in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getHydrogenBondDonorMask() const;

  /// \brief Get a bit mask of all hydrogen bond acceptors in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getHydrogenBondAcceptorMask() const;

  /// \brief List the chiral centers in a system, using topological indices.
  ///
  /// \param direction  Preferred chiral orientation of the centers to return (D-, L-, or both)
  std::vector<int> listChiralCenters(ChiralOrientation direction = ChiralOrientation::NONE) const;
  
  /// \brief Return a mask of chiral centers in the system.
  ///
  /// \param direction  Allows one to select R- (D-), S- (L-), or both chiralities for the maks
  std::vector<uint> getChiralityMask(ChiralOrientation direction) const;

  /// \brief Return a vector containing the formal charges on all particles in the system (this
  ///        includes virtual sites, which will have formal charges of zero since they are not real
  ///        atoms).
  std::vector<double> getFormalCharges() const;

  /// \brief Return a vector containing the orders of all bonds in the system (this includes bonds
  ///        to virtual sites only if they are defined in the topology's connectivity, and such
  ///        bonds will have order zero).
  std::vector<double> getBondOrders() const;

  /// \brief Get the atom endpoints of a rotatable bond.  The bond root atom is returned in the x
  ///        member of the tuple, the pivot atom (the second atom, closest to atoms that will turn)
  ///        is listed in the y member.
  ///
  /// Overloaded:
  ///   - Get all rotatable bonds without re-ordering the list.
  ///   - Get a list of all rotatable bonds upon which a minimum number of atoms turn.
  ///   - Get a list of rotatable bonds for which the pivot is within a specific cutoff of the
  ///     current conformation's center of mass.
  ///   - Get a list of rotatable bonds for which the pivot is within a specific cutoff of the
  ///     center of mass of some atom mask (the mask must be supplied as a raw bitmask of the
  ///     entire system, a std::vector of unsigned ints, to prevent making a circular dependency
  ///     whereby the AtomMask object depends on ChemicalFeatures and vice-versa).
  ///   - Get a list of the N largest rotatable groups of atoms.
  ///
  /// \param index   The index of the rotatable bond to obtain.  Such a request is naive as to what
  ///                the list actually contains, but a bounds check will be applied to ensure that
  ///                the request is valid.
  /// \param cutoff  The threshold at which to accept rotatable bond groups.  The meaning depends
  ///                the value of the choice enumeration (see below).  If the choice is
  ///                COM_PRXIMITY, then cutoff is a distance with units of Angstroms.  If the
  ///                choice is GROUP_SIZE, then cutoff is a minimium number of rotating atoms.
  /// \param choice     Specification of the means for discriminating different rotatable bond
  ///                   groups.  See the description for cutoff, above.
  /// \param mol_index  The molecule of interest (the system may have multiple molecules).
  /// \{
  std::vector<RotatorGroup> getRotatableBondGroups() const;
  std::vector<RotatorGroup> getRotatableBondGroups(int cutoff, int mol_index = 0) const;
  /// \}

  /// \brief Get a pointer to the AtomGraph which built this object.
  const AtomGraph* getTopologyPointer() const;
  
private:
  int atom_count;              ///< Number of atoms in the system
  int planar_atom_count;       ///< Number of atoms at the centers of improper dihedrals
  int ring_count;              ///< Number of unique ring systems found in the topology
  int fused_ring_count;        ///< Number of fused rings in the system
  int twistable_ring_count;    ///< Number of rings which can undergo conformational changes, i.e.
                               ///<   boat / chair.  The conformations are not called boat or chair
                               ///<   per se, but provided in tabulated lists based on the size of
                               ///<   the ring.  This indicates the number of such rings, which
                               ///<   will be a subset of all rings (excluding those which are
                               ///<   aromatic).
  int conjugated_group_count;  ///< Number of conjugated systems found in the topology
  int aromatic_group_count;    ///< Number of aromatic groups found in the topology
  int polar_hydrogen_count;    ///< Number of polar hydrogens in the system, across all molecules
  int hbond_donor_count;       ///< Number of hydrogen bond donors, across all molecules
  int hbond_acceptor_count;    ///< Number of hydrogen bond acceptors, across all molecules
  int chiral_center_count;     ///< Number of chiral centers in the topology
  int rotatable_bond_count;    ///< Number of fully rotatable single bonds
  int double_bond_count;       ///< Number of detected double bonds
  int triple_bond_count;       ///< Number of detected triple bonds
  int max_ring_size;           ///< The maximum size of a ring, based on the size of a long long
                               ///<   unsigned integer (stored for reference)
  double temperature;          ///< The temperature at which these chemical features were
                               ///<   determined (this can influence the Lewis structure)
  bool rotating_groups_mapped; ///< Flag to indicate that rotating groups have been mapped
  
  /// List of atoms which constitute planar centers (needs no bounds array, simply a list of
  /// unique atoms at the centers of improper dihedrals)
  Hybrid<int> planar_centers;
  
  /// Indicates, for each atom i in the nth bit, whether the atom is part of an n-membered ring.
  /// Rings of up to 60 atoms can be tracked in this manner.  The highest four bits are reserved
  /// to indicate how many unique rings an atom takes part in (up to 15, which should be chemically
  /// impossible except in the case of a bonded crystal).  For atoms in chemically bonded crystal
  /// lattices, the ring numbers could technically explode (every atom being considered part of a
  /// multiple higher-order rings.  For this reason, the algorithm for detecting rings will be
  /// recursive but always refer back to the array itself for the current status of the ring
  /// occupancy.  It will cut off if all 60 slots have been filled, and declare the number of
  /// simultaneously occupied rings to be 15 (a sort of inf representation).
  Hybrid<ullint> ring_inclusion;

  /// Bounds for the rings that have been found
  Hybrid<int> ring_atom_bounds;
  
  /// A list of all rings in the system.  The starting point of the kth ring is given by
  /// ring_atom_bounds[k] and the size of the ring can be computed by referencing
  /// ring_atom_bonds[k + 1].
  Hybrid<int> ring_atoms;

  /// Bounds array for the aromatic atoms list below
  Hybrid<int> aromatic_group_bounds;

  /// Counts of pi electrons in each detected aromatic group
  Hybrid<int> aromatic_pi_electrons;
  
  /// List of atoms indicating groups of aromaticity in the system
  Hybrid<int> aromatic_groups;

  /// List of polar hydrogen atoms
  Hybrid<int> polar_hydrogens;
  
  /// List of hydrogen bond donor atoms (while the order of donors and actual polar hydrogens
  /// may be similar, the fact that any one donor could have more than one polar hydrogen attached
  /// to it prevents there being a 1:1 mapping between the arrays)
  Hybrid<int> hydrogen_bond_donors;

  /// List of hydrogen bond acceptor atoms
  Hybrid<int> hydrogen_bond_acceptors;
  
  /// List of chiral centers (needs no bounds array, but each center is a five-membered tuple
  /// consisting of the topological index of the center itself, followed by the four substituent
  /// atoms, in order of decreasing IUPAC score, which define its chirality.  Positive values in
  /// this list indicate that an atom is L-chiral, while negative values indicate that an atom
  /// with the index equal to the absolute value is D-chiral.
  Hybrid<int> chiral_centers;

  /// List of rotating atoms, and the endpoints of bonds defining the axes of rotation.  The first
  /// atom in each group is the "root," distal to the atoms that move from the "pivot," the second
  /// atom listed in each group.  Neither the root atom or the pivot actually moves, but subsequent
  /// atoms in each group do.  Use the rotatable_group_bounds array to determine the extent of each
  /// group.
  Hybrid<int> rotatable_groups;

  /// Bounds array for rotatable groups.  Each group will have at least three members: the root
  /// atom, the pivot atom, and one or more atoms beyond the pivot that rotate as a consequence of
  /// twisting about the rotatable bond.
  Hybrid<int> rotatable_group_bounds;

  /// List of atoms that move (by a C2 symmetry operation) to invert a chiral center.  These groups
  /// are chosen to minimize the number of moving atoms involved in flipping a chiral center.  For
  /// some molecules like proteins, the list can get very long, so a flag is provided for disabling
  /// the search if it is not needed.
  Hybrid<int> invertible_groups;
  
  /// Bounds array for chiral inversion groups.  Each group will have at least two members,
  /// covering two branches of the chiral center that will undergo a C2 symmetry operation in order
  /// to invert the center.
  Hybrid<int> invertible_group_bounds;

  /// First atoms of the heaviest chain branching from each chiral center
  Hybrid<int> anchor_a_branches;
  
  /// First atoms of the second-heaviest chain branching from each chiral center
  Hybrid<int> anchor_b_branches;
  
  /// Formal charges determined for all atoms in the system, based on a Lewis structure drawing
  /// with the Indigo method
  Hybrid<double> formal_charges;

  /// Bond orders determined for all bonds in the system, based on a Lewis structure drawing with
  /// the Indigo method
  Hybrid<double> bond_orders;

  /// Free electron content of each atom in the molecule, averaged over resonance states
  Hybrid<double> free_electrons;

  /// Integer data (the Hybrid<int> arrays above are POINTER-kind and will targer this storage)
  Hybrid<int> int_data;  

  /// Double-precision data, targeted by the POINTER-kind Hybrid<double> objects above
  Hybrid<double> double_data;
  
  /// Pointer to the topology that this object describes (needed, if not just for reference, to
  /// obtain the actual atom indices by the various bonds enumerated in some of the lists above)
  const AtomGraph *ag_pointer;

  /// \brief Find all planar atoms in the system based on improper dihedral terms.  This prepares
  ///        data arrays that will later be loaded into the ChemicalFeatures object, but does not
  ///        directly modify the object's members.
  ///
  /// \param vk  Abstract of valence terms from the topology
  std::vector<int> findPlanarAtoms(const ValenceKit<double> &vk) const;

  /// \brief 

  /// \brief Trace all rings in the system based on a tree structure linked list.
  ///
  /// \param nbk                   Nonbonded abstract from the original topology
  /// \param cdk                   Atom, residue, and molecule details taken from the original
  ///                              topology
  /// \param tmp_ring_inclusions   Developing list of ring inlusions.  The nth bit of the kth
  ///                              array element indicates whether the kth atom participates in a
  ///                              ring of size n.
  /// \param tmp_ring_atoms        Growing list of atoms in all rings found thus far
  /// \param tmp_ring_atom_bounds  Bounds array for tmp_ring_atoms.  Must be initialized with a
  ///                              leading zero so that the push_back method can add the boundaries
  ///                              of each successive group and arrive at a capped, exclusive
  ///                              prefix sum of all ring sizes.
  void traceTopologicalRings(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                             std::vector<ullint> *tmp_ring_inclusion,
                             std::vector<int> *tmp_ring_atoms,
                             std::vector<int> *tmp_ring_atom_bounds);
  
  /// \brief Draw a ring containing two atoms, based on the histories of each atom's branch.
  ///
  /// \param j_atom                The first atom that is definitely in the ring, links to k_atom
  /// \param k_atom                The second atom that is definitely in the ring, links to j_atom
  /// \param tree_positions        Positions of various atoms in the linked list of BondedNode
  ///                              objects
  /// \param node_count            The number of nodes in the linked list object
  /// \param links                 The linked list object storing information on what atoms bond
  /// \param tmp_ring_inclusion    Array recording whether each atom takes part in a ring of size
  ///                              n, as shown by the nth bit of each ullint 
  /// \param tmp_ring_atoms        Growing list of known ring-forming atoms
  /// \param tmp_ring_atom_bounds  Bounds array for tmp_ring_atoms
  void markRingAtoms(int j_atom, int k_atom, const std::vector<int> &tree_positions,
                     int node_count, std::vector<BondedNode> *links,
                     std::vector<ullint> *tmp_ring_inclusion, std::vector<int> *tmp_ring_atoms,
                     std::vector<int> *tmp_ring_atom_bounds,
                     const ChemicalDetailsKit &cdk);

  /// \brief Draw Lewis structures over the entire topology using the Indigo method.  Lewis
  ///        structures will be drawn for each unique molecule and copied otherwise, but that
  ///        requires a list of all unique molecules.
  ///
  /// \param vk   Information on valence interactions, taken from the original topology
  /// \param nbk  Information on non-bonded interactions, taken from the original topology
  /// \param cdk  Atom, residue, and molecule details taken from the original topology
  void drawLewisStructures(const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                           const ChemicalDetailsKit &cdk);

  /// \brief Find groups of aromatic atoms.  This will also enumerate the number of fused rings in
  ///        the system as such systems must be considered for aromaticity even if their component
  ///        rings fall short.
  ///
  /// \param cdk                        Chemical details abstract from the original topology, with
  ///                                   atomic numbers of all atoms
  /// \param vk                         Valence term abstract from the original topology
  /// \param tmp_ring_atoms             List of atoms involved in rings.  The associated bounds
  ///                                   array is not needed in this case, as the function just
  ///                                   steps over all atoms that are in some ring.
  /// \param tmp_ring_atom_bounds       Bounds array for each ring in the topology
  /// \param tmp_aromatic_groups        Growing list of aromatic groups detected in the topology
  /// \param tmp_aromatic_group_bounds  Developing bounds array for aromatic groups.  Must be
  ///                                   pre-intialized with a leading zero to allow use of the
  ///                                   push_back method to grow the array.
  void findAromaticGroups(const ChemicalDetailsKit &cdk, const ValenceKit<double> &vk,
                          const std::vector<int> &tmp_ring_atoms,
                          const std::vector<int> &tmp_ring_atom_bounds,
                          std::vector<int> *tmp_aromatic_group_bounds,
                          std::vector<int> *tmp_aromatic_pi_electrons,
                          std::vector<int> *tmp_aromatic_groups);

  /// \brief Find chiral centers in the system.  This will use methods similar to the ring
  ///        detection system.  Apply IUPAC rules based on substituent atomic numbers and bond
  ///        orders.  Isotopes will not be applied.  This returns a list of detected chiral
  ///        centers for later incorporation into the actual object.
  ///
  /// \param nbk  Nonbonded system details, abstract taken from the original topology
  /// \param vk   Valence term abstract from the original topology
  /// \param cdk  Chemical details abstract from the original topology
  /// \param cf   Coordinates for the system, needed to eventually determine R- or S-chirality
  std::vector<int> findChiralCenters(const NonbondedKit<double> &nbk, const ValenceKit<double> &vk,
                                     const ChemicalDetailsKit &cdk,
                                     const CoordinateFrameReader &cfr) const;

  /// \brief Find rotatable bonds in the system, those with bond order of 1.0 and nontrivial groups
  ///        sprouting from either end, and return a vector of the atom indices at either end.
  ///
  /// \param vk                          Valence term abstract from the original topology
  /// \param cdk                         Chemical details of the system (for atomic numbers)
  /// \param nbk                         Nonbonded system details, abstract taken from the
  ///                                    original topology
  /// \param ring_atoms                  List of atoms in rings, indexing the original topology
  /// \param ring_atom_bounds            Bounds array for the ring_atoms array
  /// \param tmp_rotatable_groups        List of all atoms involved in rotation about a rotatable
  ///                                    bond, including the endpoints of the bond itself in the
  ///                                    first two slots.  This vector is assembled and returned.
  /// \param tmp_rotatable_group_bounds  Bounds array for tmp_rotatable_groups, assembled and
  ///                                    returned.
  void findRotatableBonds(const ValenceKit<double> &vk, const ChemicalDetailsKit &cdk,
                          const NonbondedKit<double> &nbk, const std::vector<int> &ring_atoms,
                          const std::vector<int> &ring_atom_bounds,
                          std::vector<int> *tmp_rotatable_groups,
                          std::vector<int> *tmp_rotatable_group_bounds);

  /// \brief Find invertible groups in the system, those comprising two branches of a chiral center
  ///        that have the fewest possible atoms.  The number of invertible groups is the number of
  ///        chiral centers.
  ///
  /// \param tmp_chiral_centers           Array of pre-determined chiral centers
  /// \param tmp_anchor_a_branches        Roots of the largest branch on each chiral center (in
  ///                                     terms of the number of atoms), not to be inverted
  /// \param tmp_anchor_b_branches        Roots of the second largest branch on each chiral center,
  ///                                     not to be inverted
  /// \param tmp_invertible_groups        Array of atoms in invertible groups (assembled and
  ///                                     returned)
  /// \param tmp_invertible_group_bounds  Bounds array for tmp_invertible_groups (assembled and
  ///                                     returned)
  void findInvertibleGroups(const std::vector<int> tmp_chiral_centers,
                            std::vector<int> *tmp_anchor_a_branches,
                            std::vector<int> *tmp_anchor_b_branches,
                            std::vector<int> *tmp_invertible_groups,
                            std::vector<int> *tmp_invertible_group_bounds);
  
  /// \brief Find polar heavy atoms that can act as hydrogen bond donors, and label polar hydrogens
  ///        in the process.  Unlike findRotatableBonds above, this is a fast evaluation and will
  ///        be performed automatically with the ChemicalFeatures object construction.
  ///
  /// \param nbk          Non-bonded interaction abstract from the original topology
  /// \param cdk          Chemical details of the system (for atomic numbers)
  /// \param tmp_polar_h  Array for listing polar hydrogen atoms (modified and returned)
  /// \param tmp_hb_don   Array for listing hydrogen bond donor atoms (modified and returned)
  /// \param tmp_hb_acc   Array for listing hydrogen bond acceptor atoms (modified and returned)
  void findHydrogenBondElements(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                std::vector<int> *tmp_polar_h, std::vector<int> *tmp_hb_don,
                                std::vector<int> *tmp_hb_acc);

  /// \brief Reset POINTER-kind Hybrid objects to target the appropriate ARRAY-kind object in
  ///        the copy constructor and copy assignment operator.
  void repairPointers();
};

/// \brief Score the four branches of a chiral molecule.  This is called by the findChiralCenters()
///        member function of the ChemicalFeatures object, but written as a free function as it
///        would not benefit from any of the object's member variables.  Returns true if scoring
///        indicates that more branch exploration is needed and might successfully discriminate
///        among the four branches, false otherwise.
///
/// \param links             Quartet of trees describing branches out of the putative chiral
///                          center.  Also contains information on bond orders.
/// \param layer_llim        Lower bounds for each branch's atoms added since the previous pass
/// \param layer_hlim        Upper bounds for each branch's atoms added since the previous pass
/// \param cdk               Chemical details of the system, an abstract taken from the original
///                          topology.  Contains atomic numbers for all atoms.
/// \param chiral_dominance  Dominance matrix describing whether one branch beats another.
///                          Modified and returned.
/// \param parallel_growth   Matrix indicating whether each branch grows in parallel with another
bool scoreChiralBranches(const std::vector<std::vector<BondedNode>> &links,
                         const std::vector<int> &layer_llim, const std::vector<int> &layer_hlim,
                         const ChemicalDetailsKit &cdk, std::vector<int> *chiral_dominance,
                         std::vector<int> *parallel_growth);

/// \brief Determine the orientation of a chiral center given the priorities of its various
///        branches.  Only the first atom of each branch is critical at this point.  Return +11 for
///        L-chiral centers (S-) or -1 for D-chiral (R-) centers.  This is then used as a
///        multiplier for the chiral center atom index in subsequent masks or atom retrievals.
///
/// \param cfr            Coordinates of the entire system
/// \param center_atom    Index of the center atom
/// \param root_atom      Index of the "root" atom--the very lowest priority branch
/// \param branch_a_atom  Index of the highest priority branch
/// \param branch_b_atom  Index of the next highest priority branch
/// \param branch_c_atom  Index of the lowest priority branch, aside from the root branch
int getChiralOrientation(const CoordinateFrameReader &cfr, int center_atom, int root_atom,
                         int branch_a_atom, int branch_c_atom, int branch_d_atom);

} // namespace chemistry
} // namespace omni

#endif
