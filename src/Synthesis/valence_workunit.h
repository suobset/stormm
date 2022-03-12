// -*-c++-*-
#ifndef OMNI_VALENCE_WORKUNIT_H
#define OMNI_VALENCE_WORKUNIT_H

#include <string>
#include <vector>
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "atomgraph_synthesis.h"

namespace omni {
namespace synthesis {

using restraints::RestraintApparatus;
using restraints::RestraintApparatusDpReader;
using topology::AtomGraph;
using topology::ConstraintKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
  
/// \brief Object to track how different valence terms in a topology are delegated.  Valence work
///        units may evaluate a valence term without being responsible for moving both atoms, or
///        even for moving any of the atoms at all.  Each valence term is only fully delegated
///        once valence work units that evaluate it are responsible for moving all atoms that the
///        valence term contains.  In order for a work unit to move any aotm, it must evaluate all
///        valence terms that include that atom.
class ValenceDelegator {
public:

  /// \brief The object is constructed based on a single topology and oversees the construction of
  ///        an array of valence work units.
  ///
  /// Overloaded:
  ///   - Constructor for an empty object
  ///   - Take the original data structures by const pointer or const reference
  ///
  /// \param ag  The topology containing valence terms to delegate among work units
  /// \param ra  The complete collection of restraints applicable to this topology
  /// \{
  ValenceDelegator();
  ValenceDelegator(const AtomGraph &ag, const RestraintApparatus &ra);
  ValenceDelegator(const AtomGraph *ag, const RestraintApparatus *ra);
  /// \}
  
private:
  int atom_count;                 ///< The number of atoms in the system overall (taken from the
                                  ///<   topology)

  // The following are all readable as column-format matrices, the value of element (y, k)
  // indicating the index of the kth work unit (W.U.) containing atom I, J, ..., M in bond / angle /
  // ... / CMAP term y.  Most atoms will be present in only one work unit.  The number of rows in
  // each matrix is the number of bond / angle / ... / CMAP terms.
  std::vector<int> bond_i_presence;  ///< Indices of W.U.'s containing atom I of each bond
  std::vector<int> bond_j_presence;  ///< Indices of W.U.'s containing atom J of each bond
  std::vector<int> angl_i_presence;  ///< Indices of W.U.'s containing atom I of each bond angle
  std::vector<int> angl_j_presence;  ///< Indices of W.U.'s containing atom J of each bond angle
  std::vector<int> angl_k_presence;  ///< Indices of W.U.'s containing atom K of each bond angle
  std::vector<int> dihe_i_presence;  ///< Indices of W.U.'s containing atom I of each dihedral
  std::vector<int> dihe_j_presence;  ///< Indices of W.U.'s containing atom J of each dihedral
  std::vector<int> dihe_k_presence;  ///< Indices of W.U.'s containing atom K of each dihedral
  std::vector<int> dihe_l_presence;  ///< Indices of W.U.'s containing atom L of each dihedral
  std::vector<int> ubrd_i_presence;  ///< Indices of W.U.'s containing atom I of each Urey-Bradley
  std::vector<int> ubrd_k_presence;  ///< Indices of W.U.'s containing atom K of each Urey-Bradley
  std::vector<int> cimp_i_presence;  ///< Indices of W.U.'s containing atom I of each
                                     ///<   CHARMM improper dihedral
  std::vector<int> cimp_j_presence;  ///< Indices of W.U.'s containing atom J of each
                                     ///<   CHARMM improper dihedral
  std::vector<int> cimp_k_presence;  ///< Indices of W.U.'s containing atom K of each
                                     ///<   CHARMM improper dihedral
  std::vector<int> cimp_l_presence;  ///< Indices of W.U.'s containing atom L of each
                                     ///<   CHARMM improper dihedral
  std::vector<int> cmap_i_presence;  ///< Indices of W.U.'s containing atom I of each CMAP term
  std::vector<int> cmap_j_presence;  ///< Indices of W.U.'s containing atom J of each CMAP term
  std::vector<int> cmap_k_presence;  ///< Indices of W.U.'s containing atom K of each CMAP term
  std::vector<int> cmap_l_presence;  ///< Indices of W.U.'s containing atom L of each CMAP term
  std::vector<int> cmap_m_presence;  ///< Indices of W.U.'s containing atom M of each CMAP term

  // The work unit noted in element (y, k) of the following array indicates that the yth virtual
  // site or one of its frame atoms is present the work unit (W.U.) with the stated index.
  std::vector<int> vs_presence;      ///< Indices of W.U.'s containing each virtual site
  std::vector<int> vsf1_presence;    ///< Indices of W.U.'s containing the first frame atom of each
                                     ///<   virtual site (the "parent" atom)
  std::vector<int> vsf2_presence;    ///< Indices of W.U.'s containing the second frame atom of each
                                     ///<   virtual site
  std::vector<int> vsf3_presence;    ///< Indices of W.U.'s containing the third frame atom of each
                                     ///<   virtual site
  std::vector<int> vsf4_presence;    ///< Indices of W.U.'s containing the fourth frame atom of each
                                     ///<   virtual site

  // The following arrays are again to be interpreted as matrices, with a number of rows equal to
  // the number of atoms in constrained groups or the total number of fast waters.  The work unit
  // noted in element (y, k) indicates that the yth constrained atom is present the work unit (W.U.)
  // with the stated index.  The group constraints bounds array will be needed to determine whether
  // a particular constrained group is completely present in one work unit or another.
  std::vector<int> cnst_n_presence;  ///< Indices of W.U.'s containing the nth constrained atom from
                                     ///<   the list of all atoms in constrained groups
  std::vector<int> sett_ox_presence; ///< Indices of W.U.'s containing the oxygen atom of a
                                     ///<   particular SETTLE group
  std::vector<int> sett_h1_presence; ///< Indices of W.U.'s containing the oxygen atom of a
                                     ///<   particular SETTLE group
  std::vector<int> sett_h2_presence; ///< Indices of W.U.'s containing the oxygen atom of a
                                     ///<   particular SETTLE group

  // The following arrays are matrices, with a number of rows equal to the number of each type of
  // restraint.  The work unit noted in element (y, k) indicates that the atom in some restraint of
  // the member name-implied type is present the the work unit (W.U.) with the stated index.
  std::vector<int> rposn_i_presence; ///< Indices of W.U.'s containing specific atoms subject to
                                     ///<   positional restraints
  std::vector<int> rbond_i_presence; ///< Indices of W.U.'s containing the first atom of each
                                     ///<   distance restraint
  std::vector<int> rbond_j_presence; ///< Indices of W.U.'s containing the second atom of each
                                     ///<   distance restraint
  std::vector<int> rangl_i_presence; ///< Indices of W.U.'s containing the first atom of each
                                     ///<   three-point angle restraint
  std::vector<int> rangl_j_presence; ///< Indices of W.U.'s containing the second atom of each
                                     ///<   three-point angle restraint
  std::vector<int> rangl_k_presence; ///< Indices of W.U.'s containing the third atom of each
                                     ///<   three-point angle restraint
  std::vector<int> rdihe_i_presence; ///< Indices of W.U.'s containing the first atom of each
                                     ///<   four-point dihedral angle restraint
  std::vector<int> rdihe_j_presence; ///< Indices of W.U.'s containing the second atom of each
                                     ///<   four-point dihedral angle restraint
  std::vector<int> rdihe_k_presence; ///< Indices of W.U.'s containing the third atom of each
                                     ///<   four-point dihedral angle restraint
  std::vector<int> rdihe_l_presence; ///< Indices of W.U.'s containing the fourth atom of each
                                     ///<   four-point dihedral angle restraint
  
  // The following lists are quick to construct and contain a double-, triple-, or even a
  // quintuple-counting of all information in the topology's own "atom assignment" arrays for the
  // same valence terms.  Is there a better way?
  std::vector<int> bond_affector_list;    ///< List of all harmonic bonds affecting a given atom
  std::vector<int> bond_affector_bounds;  ///< Bounds array for bond_affector_list
  std::vector<int> angl_affector_list;    ///< List of all bond angles affecting a given atom
  std::vector<int> angl_affector_bounds;  ///< Bounds array for angl_affector_list
  std::vector<int> dihe_affector_list;    ///< List of all cosine dihedrals affecting a given atom
  std::vector<int> dihe_affector_bounds;  ///< Bounds array for dihe_affector_list
  std::vector<int> ubrd_affector_list;    ///< List of Urey-Bradley terms affecting a given atom
  std::vector<int> ubrd_affector_bounds;  ///< Bounds array for ubrd_affector_list
  std::vector<int> cimp_affector_list;    ///< List of CHARMM impropers affecting a given atom
  std::vector<int> cimp_affector_bounds;  ///< Bounds array for cimp_affector_list
  std::vector<int> cmap_affector_list;    ///< List of all CMAP terms affecting a given atom
  std::vector<int> cmap_affector_bounds;  ///< Bounds array for cmap_affector_list

  // The virtual site affector arrays will be pre-allcoated to hold virtual sites with up to four
  // frame atoms apiece.
  std::vector<int> vste_affector_list;    ///< List of all virtual sites that name a given atom as
                                          ///<   one of their frame atoms
  std::vector<int> vste_affector_bounds;  ///< Bounds array for vste_affector_list

  // SHAKE and RATTLE groups likewise must have all atoms present in a work group in order to
  // evaluate bond constraints.
  std::vector<int> cnst_affector_list;    ///< List of all constraint groups affecting any atom
                                          ///<   (these must be detected through a special search
                                          ///<   so that only a work unit with all affected atoms
                                          ///<   present will be tasked with evaluating the
                                          ///<   constraints).
  std::vector<int> cnst_affector_bounds;  ///< Bounds array for cnst_affector_list

  // SETTLE groups are special constraint groups for a fast rigid water implementation, kept
  // distinct because of the different treatment they get in the dynamics loop.
  std::vector<int> sett_affector_list;    ///< List of all SETTLE groups affecting any atom
                                          ///<   (only one work group will be tasked to evaluate
                                          ///<   each SETTLE fast water)
  std::vector<int> sett_affector_bounds;  ///< Bounds array for sett_affector_list

  // Positional restraints apply to a single atom, but multiple restraints could still, in theory,
  // restrain a single atom to a ring or other confined region of space.  A bounds array is still
  // necessary.  Other restraints affect multiple atoms and map like their valence term analogs.
  std::vector<int> rposn_affector_list;   ///< List of all positional restraints affecting any atom
  std::vector<int> rposn_affector_bounds; ///< Bounds array for rposn_affector_list
  std::vector<int> rbond_affector_list;   ///< List of all bond restraints affecting any atom
  std::vector<int> rbond_affector_bounds; ///< Bounds array for rbond_affector_list
  std::vector<int> rangl_affector_list;   ///< List of all three-point restraints affecting any atom
  std::vector<int> rangl_affector_bounds; ///< Bounds array for rangl_affector_list
  std::vector<int> rdihe_affector_list;   ///< List of all four-point restraints affecting any atom
  std::vector<int> rdihe_affector_bounds; ///< Bounds array for rdihe_affector_list
  
  // Individual atoms must leave a record of their whereabouts in the valence work units for rapid
  // retrieval of their locations
  std::vector<int> work_unit_assignment_count;  ///< The numbers of work units in which each atom
                                                ///<   can be found
  std::vector<int> work_unit_presence;          ///< Lists of the work units in which each atom is
                                                ///<   found.  This is a column-format matrix with
                                                ///<   atom_count rows and a number of columns
                                                ///<   expanded as needed to accommodate the largest
                                                ///<   entry in work_unit_assignments.

  /// \brief Allocate the necessary space for this work unit
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void allocate(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                const ConstraintKit &cnk, const RestraintApparatusDpReader &rar);

  /// \brief Fill the arrays describing how different atoms are affected by each potential term
  ///        (including restraints), each virtual site, and each constraint.
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void fillAffectorArrays(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                          const ConstraintKit &cnk, const RestraintApparatusDpReader &rar);
};
  
/// \brief An object to collect the components of a valence work unit (which will also track frozen
///        atoms to implement coordinate updates, velocity updates, and constraints).  While the
///        work unit is encoded in the AtomGraphSynthesis object, the assembly is best done by a
///        dedicated object with plenty of its own methods.
class ValenceWorkUnit {
public:

  /// \brief The constructor takes a specific input topology (multiple systems using the same
  ///        topology in an AtomGraphSynthesis can thereby take the same valence work unit and
  ///        translate the atom indices as appropriate rather than regenerating the work unit
  ///        for many equivalent systems).  Mapping starts from a specific atom and proceeds until
  ///        a maximum number of atoms has been accumulated in order to process as many related
  ///        valence terms as possible.
  ///
  /// \param ag         Topology to work from
  /// \param seed_atom  The first atom to incorporate into the work unit.  Subsequent atoms will
  ///                   be either bonded in some chain to the seed, retracing previous topological
  ///                   indices whereby previous work units left some atoms behind, or jumping
  ///                   forward to the next new molecule.
  /// \param max_atoms  The maximum number of atoms to accumulate in the work unit
  ValenceWorkUnit(const AtomGraph &ag, ValenceDelegator *vdel, int seed_atom,
                  int max_atoms = 768);
  
  /// \brief Add a new atom to a work unit
  
  
private:
  int atom_count;
  int max_atom_index;
  int min_atom_index;
  std::vector<int> atom_import_list;
  
};

} // namespace topology
} // namespace omni

#endif
