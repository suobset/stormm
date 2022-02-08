// -*-c++-*-
#ifndef OMNI_VALENCE_WORKUNIT_H
#define OMNI_VALENCE_WORKUNIT_H

#include <string>
#include <vector>
#include "atomgraph.h"
#include "atomgraph_synthesis.h"

namespace omni {
namespace topology {

/// \brief Object to track how different valence terms in a topology are delegated.  Valence work
///        units may evaluate a valence term without being responsible for moving both atoms, or
///        even for moving any of the atoms at all.  Each valence term is only fully delegated
///        once valence work units that evaluate it are responsible for moving all atoms that the
///        valence term contains.  In order for a work unit to move any aotm, it must evaluate all
///        valence terms that include that atom.
struct ValenceDelegation {

  /// \brief The object is constructed based on a single topology and oversees the construction of
  ///        an array of valence work units.
  ///
  /// \param ag  The topology containing valence terms to delegate among work units
  ValenceDelegation(const AtomGraph &ag);
  
private:
  int atom_count;                 ///< The number of atoms in the system overall (taken from the
                                  ///<   topology)

  // The yth element of each of the following arrays indicates that atom I, J, ..., M in bond /
  // angle / ... / CMAP term y is moved by the work unit noted in the yth array element.
  std::vector<int> bond_i_atoms;  ///< Indices of work units that move atom I of each bond
  std::vector<int> bond_j_atoms;  ///< Indices of work units that move atom J of each bond
  std::vector<int> angl_i_atoms;  ///< Indices of work units that move atom I of each bond angle
  std::vector<int> angl_j_atoms;  ///< Indices of work units that move atom J of each bond angle
  std::vector<int> angl_k_atoms;  ///< Indices of work units that move atom K of each bond angle
  std::vector<int> dihe_i_atoms;  ///< Indices of work units that move atom I of each dihedral
  std::vector<int> dihe_j_atoms;  ///< Indices of work units that move atom J of each dihedral
  std::vector<int> dihe_k_atoms;  ///< Indices of work units that move atom K of each dihedral
  std::vector<int> dihe_l_atoms;  ///< Indices of work units that move atom L of each dihedral
  std::vector<int> ubrd_i_atoms;  ///< Indices of work units moving atom I of each Urey-Bradley
  std::vector<int> ubrd_k_atoms;  ///< Indices of work units moving atom K of each Urey-Bradley
  std::vector<int> cimp_i_atoms;  ///< Indices of work units moving atom I of each CHARMM improper
  std::vector<int> cimp_j_atoms;  ///< Indices of work units moving atom J of each CHARMM improper
  std::vector<int> cimp_k_atoms;  ///< Indices of work units moving atom K of each CHARMM improper
  std::vector<int> cimp_l_atoms;  ///< Indices of work units moving atom L of each CHARMM improper
  std::vector<int> cmap_i_atoms;  ///< Indices of work units moving atom I of each CMAP term
  std::vector<int> cmap_j_atoms;  ///< Indices of work units moving atom J of each CMAP term
  std::vector<int> cmap_k_atoms;  ///< Indices of work units moving atom K of each CMAP term
  std::vector<int> cmap_l_atoms;  ///< Indices of work units moving atom L of each CMAP term
  std::vector<int> cmap_m_atoms;  ///< Indices of work units moving atom M of each CMAP term

  /// The work unit noted in the yth element of the following array places the yth virtual site.
  std::vector<int> virtual_site_placement;  
  
  // The following lists are quick to construct and contain a double-counting of all information
  // in the topology's own "atom assignment" arrays for the same valence terms.  Is there a better
  // way?
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
};
  
/// \brief An object to collect the components of a valence work unit (which will also track frozen
///        atoms to implement coordinate updates, velocity updates, and constraints).  While the
///        work unit is encoded in the AtomGraphSynthesis object, the assembly is best done by a
///        dedicated object with plenty of its own methods.
struct ValenceWorkUnit {

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
  ValenceWorkUnit(const AtomGraph &ag, ValenceDelegation *vdel, int seed_atom,
                  int max_atoms = 768);

  /// \brief Test whether this object contains a particular atom
  
private:
  int atom_count;
  std::vector<int> atom_import_list;
  
};

} // namespace topology
} // namespace omni

#endif
