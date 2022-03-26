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

/// \brief The minimium value for the atom limit in a valence work unit--there should be no need
///        to force the size below 128, which was the maximum size in pmemd.cuda.
constexpr int minimum_valence_work_unit_atoms = 128;

/// \brief The maximum value for the atom limit in a valence work unit--any higher and the amount
///        of __shared__ memory in a block of 1024 threads might need to be stretched.
constexpr int maximum_valence_work_unit_atoms = 768;
  
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
  /// \param ag_in  Pointer to the topology containing valence terms to delegate among work units
  /// \param ra_in  Pointer to the complete collection of restraints applicable to the system
  ValenceDelegator(const AtomGraph *ag_in, const RestraintApparatus *ra_in = nullptr);

  /// \brief Get the number of work units to which a particular atom is (currently) assigned.
  ///
  /// \param atom_index  The atom of interest, as indexed in the original topology
  int getAtomAssignmentCount(int atom_index) const;

  /// \brief Get the index of the first unassigned atom in the topology.
  int getFirstUnassignedAtom() const;

  /// \brief Get a list of bonds affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getBondAffectors(int atom_index) const;
  
  /// \brief Get a list of harmonic angle terms affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getAngleAffectors(int atom_index) const;

  /// \brief Get a list of cosine-based dihedrals affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getDihedralAffectors(int atom_index) const;
  
  /// \brief Get a list of Urey-Bradley harmonic angles affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getUreyBradleyAffectors(int atom_index) const;
  
  /// \brief Get a list of CHARMM improper dihedral terms affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getCharmmImproperAffectors(int atom_index) const;
  
  /// \brief Get a list of CMAP terms affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getCmapAffectors(int atom_index) const;
  
  /// \brief Get a list of positional restraints affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getPositionalRestraintAffectors(int atom_index) const;
  
  /// \brief Get a list of distance restraints affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getDistanceRestraintAffectors(int atom_index) const;

  /// \brief Get a list of three-point angle restraints affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getAngleRestraintAffectors(int atom_index) const;

  /// \brief Get a list of four-point dihedral restraints affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getDihedralRestraintAffectors(int atom_index) const;

  /// \brief Get a list of SETTLE constraint groups affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getSettleGroupAffectors(int atom_index) const;

  /// \brief Get a list of SHAKE or RATTLE constraint groups affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getConstraintGroupAffectors(int atom_index) const;

  /// \brief Get a list of virtual sites affecting a given atom.
  ///
  /// param atom_index  Topological index of the atom of interest
  std::vector<int> getVirtualSiteAffectors(int atom_index) const;

  /// \brief Mark the addition of an atom to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index   Index of the ValenceWorkUnit receiving the new atom
  /// \param atom_index  Index of the atom to add, referencing the original topology
  void markAtomAddition(int vwu_index, int atom_index);

  /// \brief Mark the addition of a harmonic bond term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param bond_term_index  Index of the bond term to add, referencing the original topology
  void markBondAddition(int vwu_index, int bond_term_index);

  /// \brief Mark the addition of a harmonic angle term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param angl_term_index  Index of the angle term to add, referencing the original topology
  void markAngleAddition(int vwu_index, int angl_term_index);

  /// \brief Mark the addition of a cosine-based dihedral term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param dihe_term_index  Index of the dihedral term to add, referencing the original topology
  void markDihedralAddition(int vwu_index, int dihe_term_index);

  /// \brief Mark the addition of a Urey-Bradley term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param ubrd_term_index  Index of the Urey-Bradley term to add, referencing the original
  ///                         topology
  void markUreyBradleyAddition(int vwu_index, int ubrd_term_index);

  /// \brief Mark the addition of a CHARMM improper dihedral term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param ubrd_term_index  Index of the CHARMM improper dihedral term to add, referencing the
  ///                         original topology
  void markCharmmImproperAddition(int vwu_index, int ubrd_term_index);
  
  /// \brief Mark the addition of a CMAP splined surface term to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param cmap_term_index  Index of the CMAP term to add, referencing the original topology
  void markCmapAddition(int vwu_index, int cmap_term_index);

  /// \brief Mark the addition of a positional restraint to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param posn_rstr_index  Index of the positional restraint to add, referencing the original
  ///                         restraint apparatus
  void markPositionalRestraintAddition(int vwu_index, int posn_rstr_index);

  /// \brief Mark the addition of a distance restraint to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param dist_rstr_index  Index of the distance restraint to add, referencing the original
  ///                         restraint apparatus
  void markDistanceRestraintAddition(int vwu_index, int dist_rstr_index);

  /// \brief Mark the addition of a three-point angle restraint to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param angl_rstr_index  Index of the angle restraint to add, referencing the original
  ///                         restraint apparatus
  void markAngleRestraintAddition(int vwu_index, int angl_rstr_index);

  /// \brief Mark the addition of a four-point dihedral restraint to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index        Index of the ValenceWorkUnit receiving the new atom
  /// \param dihe_rstr_index  Index of the dihedral restraint to add, referencing the original
  ///                         restraint apparatus
  void markDihedralRestraintAddition(int vwu_index, int dihe_rstr_index);

  /// \brief Mark the addition of a constraint group to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index         Index of the ValenceWorkUnit receiving the new atom
  /// \param cnst_group_index  Index of the constraint group to add, referencing the original
  ///                          topology
  void markConstraintGroupAddition(int vwu_index, int cnst_group_index);
  
  /// \brief Mark the addition of a SETTLE group to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index         Index of the ValenceWorkUnit receiving the new atom
  /// \param sett_group_index  Index of the SETTLE group to add, referencing the original
  ///                          topology
  void markSettleGroupAddition(int vwu_index, int sett_group_index);

  /// \brief Mark the addition of a virtual site frame to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index    Index of the ValenceWorkUnit receiving the new atom
  /// \param vsite_index  Index of the virtual site frame to add, referencing the original
  ///                     topology's list of virtual sites (not its list of atoms)
  void markVirtualSiteAddition(int vwu_index, int vsite_index);

private:
  int atom_count;                 ///< The number of atoms in the system overall (taken from the
                                  ///<   topology)
  int first_unassigned_atom;      ///< The first atom index with no current work unit assignments
  int max_presence_allocation;    ///< The allocation size for the maximum number of work units
                                  ///<   to which any atom can be assigned.

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
  std::vector<int> rposn_affector_list;    ///< List of positional restraints affecting any atom
  std::vector<int> rposn_affector_bounds;  ///< Bounds array for rposn_affector_list
  std::vector<int> rbond_affector_list;    ///< List of bond restraints affecting any atom
  std::vector<int> rbond_affector_bounds;  ///< Bounds array for rbond_affector_list
  std::vector<int> rangl_affector_list;    ///< List of three-point restraints affecting any atom
  std::vector<int> rangl_affector_bounds;  ///< Bounds array for rangl_affector_list
  std::vector<int> rdihe_affector_list;    ///< List of four-point restraints affecting any atom
  std::vector<int> rdihe_affector_bounds;  ///< Bounds array for rdihe_affector_list

  // Individual atoms must leave a record of their whereabouts in the valence work units for rapid
  // retrieval of their locations
  std::vector<int> work_unit_assignment_count;  ///< The numbers of work units in which each atom
                                                ///<   can be found
  std::vector<int> work_unit_presence;          ///< Lists of the work units in which each atom is
                                                ///<   found.  This is a column-format matrix with
                                                ///<   atom_count columns and a number of rows
                                                ///<   expanded as needed to accommodate the
                                                ///<   largest entry in work_unit_assignments.
  std::vector<int> assigned_update_work_units;  ///< List indices for the work units assigned to
                                                ///<   update (and log) each atom's position and
                                                ///<   velocity in the master coordinate set

  /// Pointers to the original topology and restraint apparatus that formed this object
  const AtomGraph *ag_pointer;
  const RestraintApparatus *ra_pointer;

  /// \brief Make whatever needed space for the arrays indicating which work units any particular
  ///        atom is present in.  This is not a compact array with a bounds list, but rather a
  ///        padded array with space for each atom, due to the frequency with which it might be
  ///        updated.
  ///
  /// \param n_units  The maximum number of work units that any atom might be a part of.  Allocate
  ///                 an array for this many inclusions.
  void resizePresenceArrays(int n_units);
  
  /// \brief Allocate the necessary space for this work unit
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void allocate();

  /// \brief Fill the arrays describing how different atoms are affected by each potential term
  ///        (including restraints), each virtual site, and each constraint.
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void fillAffectorArrays(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                          const ConstraintKit<double> &cnk, const RestraintApparatusDpReader &rar);
};
  
/// \brief An object to collect the components of a valence work unit (which will also track frozen
///        atoms to implement coordinate updates, velocity updates, and constraints).  While the
///        work unit is encoded in the AtomGraphSynthesis object, the assembly is best done by a
///        dedicated object with plenty of its own methods operating on a single topology
///        (AtomGraph).  All systems in the AtomGraphSynthesis are designed to function
///        independently of one another--the only difference is that they have consensus tables of
///        most parameters and differen atom indexing.  Translating a valence work unit into a
///        list of instructions within an AtomGraphSynthesis is therefore a critical member
///        function of this class.
class ValenceWorkUnit {
public:

  /// \brief The constructor takes a specific input topology (multiple systems using the same
  ///        topology in an AtomGraphSynthesis can thereby take the same valence work unit and
  ///        translate the atom indices as appropriate rather than regenerating the work unit
  ///        for many equivalent systems).  Mapping starts from a specific atom and proceeds until
  ///        a maximum number of atoms has been accumulated in order to process as many related
  ///        valence terms as possible.
  ///
  /// \param ag_in          Pointer to the topology to work from
  /// \param ra_in          Pointer to the collection of restraints applicable to the topology
  /// \param vdel_in        Valence delegator managing the creation of this valence work unit
  /// \param list_index_in  Index of this unit in a larger list (the unit should remember its own
  ///                       index number, for the purposes of coordinating with other work units)
  /// \param seed_atom_in   The first atom to incorporate into the work unit.  Subsequent atoms
  ///                       will be either bonded in some chain to the seed, retracing previous
  ///                       topological indices whereby previous work units left some atoms
  ///                       behind, or jumping forward to the next new molecule.
  /// \param max_atoms_in   The maximum number of atoms to accumulate in the work unit
  ValenceWorkUnit(const AtomGraph *ag_in, const RestraintApparatus *ra_in,
                  ValenceDelegator *vdel_in, int list_index_in, int seed_atom_in,
                  int max_atoms_in = 768);

  /// \brief Get the number of atoms currently involved in this work unit.
  int getAtomCount() const;

  /// \brief Get the list index of this work unit.
  int getListIndex() const;

  /// \brief Get the minimum topological atom index of any used by this work unit.
  int getMinAtomIndex() const;

  /// \brief Get the maximum topological atom index of any used by this work unit.
  int getMaxAtomIndex() const;
  
  /// \brief Get the maximum atom count that this work unit can hold.
  int getMaxAtoms() const;

  /// \brief Get the pointer to the ValenceDelegator managing the creation of this object.
  ValenceDelegator* getDelegatorPointer();

  /// \brief Get a pointer to the topology for which this work unit applies.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a pointer to the restraint collection for which this work unit applies.
  const RestraintApparatus* getRestraintApparatusPointer() const;

  /// \brief Get the topological index of one of the work unit's imported atoms.
  ///
  /// \param index  The local index of the imported atom of interest
  int getImportedAtom(int index) const;
  
  /// \brief Set the list index of this work unit, in the event that the list of work units for
  ///        a particular topology needs to be re-ordered.
  ///
  /// \param list_index_in  The new list index for the work unit
  void setListIndex(int list_index_in);

  /// \brief Set the atom limit for a valence work unit.  This can be useful in situations where
  ///        it is desirable to form several work units out of a single molecule, despite there
  ///        being enough room in just one to hold all atoms of the molecule.
  ///
  /// \param new_limit  The new limit on the number of atoms.  This cannot be lower than the
  ///                   number of atoms already in the work unit.
  void setAtomLimit(int new_limit);
  
  /// \brief Add a new atom to a work unit.  This will update the associated ValenceDelegator and
  ///        all assignments therein.
  ///
  /// \param atom_index  Index of the atom of interest
  void addNewAtom(int atom_index);

  /// \brief Add a new harmonic bond term to the work unit.  Record the fact in the associated
  ///        delegator.
  ///
  /// \param bond_term_index  Index of the bond term to add, based on the original topology
  void addNewBondTerm(int bond_term_index);

  /// \brief Add a new harmonic bond angle term to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param angl_term_index  Index of the angle term to add, based on the original topology
  void addNewAngleTerm(int angl_term_index);

  /// \brief Add a new cosine-based dihedral term to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param dihe_term_index  Index of the dihedral term to add, based on the original topology
  void addNewDihedralTerm(int dihe_term_index);

  /// \brief Add a new Urey-Bradley harmonic angle term to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param ubrd_term_index  Index of the Urey-Bradley term to add, based on the original topology
  void addNewUreyBradleyTerm(int ubrd_term_index);

  /// \brief Add a new CHARMM harmonic improper dihedral term to the work unit.  Record the fact
  ///        in the associated delegator.
  ///
  /// \param cimp_term_index  Index of the CHARMM improper to add, based on the original topology
  void addNewCharmmImproperTerm(int cimp_term_index);

  /// \brief Add a new CMAP bicubic spline surface term to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param cmap_term_index  Index of the CMAP to add, based on the original topology
  void addNewCmapTerm(int cmap_term_index);

  /// \brief Add a new positional restraint to the work unit.  Record the fact in the associated
  ///        delegator.
  ///
  /// \param posn_rstr_index  The index of the positional restraint to add
  void addNewPositionalRestraint(int posn_rstr_index);
  
  /// \brief Add a new distance restraint to the work unit.  Record the fact in the associated
  ///        delegator.
  ///
  /// \param dist_rstr_index  The index of the distance restraint to add
  void addNewDistanceRestraint(int dist_rstr_index);

  /// \brief Add a new three-point angle restraint to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param angl_rstr_index  The index of the angle restraint to add
  void addNewAngleRestraint(int angl_rstr_index);

  /// \brief Add a new four-point dihedral angle restraint to the work unit.  Record the fact in
  ///        the associated delegator.
  ///
  /// \param dihe_rstr_index  The index of the dihedral angle restraint to add
  void addNewDihedralRestraint(int dihe_rstr_index);

  /// \brief Add a new constraint group to the work unit.  Record the fact in the associated
  ///        delegator.
  ///
  /// \param cnst_group_index  The index of the constraint group to add
  void addNewConstraintGroup(int cnst_group_index);

  /// \brief Add a new SETTLE rigid water group to the work unit.  Record the fact in the
  ///        associated delegator.
  ///
  /// \param sett_group_index  The index of the constraint group to add
  void addNewSettleGroup(int sett_group_index);

  /// \brief Add a new virtual site, referencing the virtual site index (as opposed to the actual
  ///        atom index of the virtual site particle) in the original topology.  Record the fact
  ///        in the associated delegator.
  ///
  /// \param vsite_index  The index of the virtual site to add
  void addNewVirtualSite(int vsite_index);

  /// \brief Assign atom update responsibilities to a work unit, subject to the stipulations that
  ///        no other work unit has been assigned to update the atom's position and velocity, and
  ///        that the work unit has all of the atoms needed to compute the forces that will affect
  ///        the atom itself.  If the atom is part of a constraint group, the work unit must also
  ///        have all atoms needed to compute forces on other atoms in the constraint group, and
  ///        if the atom is a virtual site the work unit will also need to have all frame atoms
  ///        plus whatever atoms are needed to compute forces on them.  This function will scan
  ///        over all atoms in the work unit as it is found, checking that the criteria are met.
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  /// \param cnk  Constraint group abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  void assignUpdateTasks(const ValenceKit<double> &vk, const RestraintApparatusDpReader &rar,
                         const ConstraintKit<double> &cnk, const VirtualSiteKit<double> &vsk);
  
private:
  int atom_count;                     ///< Number of atoms in the work unit
  int bond_term_count;                ///< Number of bond terms in the work unit
  int angl_term_count;                ///< Number of angle terms in the work unit
  int dihe_term_count;                ///< Number of cosine-based dihedral terms in the work unit
  int ubrd_term_count;                ///< Number of Urey-Bradley terms in the work unit
  int cimp_term_count;                ///< Number of CHARMM harmonic improper dihedral terms in
                                      ///<   the work unit
  int cmap_term_count;                ///< Number of CMAP terms in the work unit
  int rposn_term_count;               ///< Number of positional restraints handled by this work
                                      ///<   unit
  int rbond_term_count;               ///< Number of distance restraints handled by this work unit
  int rangl_term_count;               ///< Number of angle restraints handled by this work unit
  int rdihe_term_count;               ///< Number of dihedral restraints handled by this work unit
  int cnst_group_count;               ///< Number of SHAKE or RATTLE groups managed by this work
                                      ///<   unit (excludes SETTLE-constrained waters)
  int sett_group_count;               ///< Number of SETTLE-constrained rigid waters managed by
                                      ///<   this work unit
  int vste_count;                     ///< Number of virtual sites managed by this work unit
  int list_index;                     ///< Index of the work unit in a larger list of similar
                                      ///<   objects coordinating to cover an entire topology
  int min_atom_index;                 ///< Lowest topological index of any imported atom
  int max_atom_index;                 ///< Highest topological index of any imported atom
  int atom_limit;                     ///< Largest number of atoms that this work unit can hold

  /// The list of imported atoms, indicating indices into the original topology.  The position of
  /// each atom in this list indicates its local index within the work unit, as referenced by
  /// subsequent ????_(i,j,k,...)_atoms arrays.
  std::vector<int> atom_import_list;

  // Valence terms for the work unit, typical force field elements
  std::vector<int> bond_term_list;    ///< List of harmonic bonds for which this work unit is
                                      ///<   responsible (more than one work unit may be tasked
                                      ///<   with computing any of the relevant energy terms.
                                      ///<   One and only one work unit will be tasked with moving
                                      ///<   each (mobile) atom.  If an atom is not mobile, no
                                      ///<   work unit will be tasked with moving it and terms
                                      ///<   pertaining to it may or may not be computed.
  std::vector<int> angl_term_list;    ///< List of harmonic angle terms to be computed by this
                                      ///<   work unit, indexed into the original topology
  std::vector<int> dihe_term_list;    ///< List of cosine-based dihedral terms to be computed by
                                      ///<   this work unit, indexed into the original topology
  std::vector<int> ubrd_term_list;    ///< List of Urey-Bradley harmonic angle terms to be computed
                                      ///<   by this work unit, indexed into the original topology
  std::vector<int> cimp_term_list;    ///< List of CHARMM harmonic improper dihedral terms to be
                                      ///<   computed by this work unit
  std::vector<int> cmap_term_list;    ///< List of CMAP terms to be computed by this work unit
  std::vector<int> bond_i_atoms;      ///< List of I atoms in each harmonic bond computed by this
                                      ///<   work unit, tracking the order in bond_term_list but
                                      ///<   indicating the local indices of imported atoms
  std::vector<int> bond_j_atoms;      ///< List of J atoms in each harmonic bond computed by this
                                      ///<   work unit, indicating local atom indices
  std::vector<int> angl_i_atoms;      ///< List of local indices for I atoms in each angle term
  std::vector<int> angl_j_atoms;      ///< List of local indices for J atoms in each angle term
  std::vector<int> angl_k_atoms;      ///< List of local indices for K atoms in each angle term
  std::vector<int> dihe_i_atoms;      ///< Local indices for I atoms in each dihedral term
  std::vector<int> dihe_j_atoms;      ///< Local indices for J atoms in each dihedral term
  std::vector<int> dihe_k_atoms;      ///< Local indices for K atoms in each dihedral term
  std::vector<int> dihe_l_atoms;      ///< Local indices for L atoms in each dihedral term
  std::vector<int> ubrd_i_atoms;      ///< Local indices for I atoms in each Urey-Bradley term
  std::vector<int> ubrd_k_atoms;      ///< Local indices for K atoms in each Urey-Bradley term
  std::vector<int> cimp_i_atoms;      ///< Local indices for I atoms in each CHARMM improper term
  std::vector<int> cimp_j_atoms;      ///< Local indices for J atoms in each CHARMM improper term
  std::vector<int> cimp_k_atoms;      ///< Local indices for K atoms in each CHARMM improper term
  std::vector<int> cimp_l_atoms;      ///< Local indices for L atoms in each CHARMM improper term
  std::vector<int> cmap_i_atoms;      ///< Local indices for I atoms in each CMAP term
  std::vector<int> cmap_j_atoms;      ///< Local indices for J atoms in each CMAP term
  std::vector<int> cmap_k_atoms;      ///< Local indices for K atoms in each CMAP term
  std::vector<int> cmap_l_atoms;      ///< Local indices for L atoms in each CMAP term
  std::vector<int> cmap_m_atoms;      ///< Local indices for M atoms in each CMAP term

  // Restraint terms for this work unit
  std::vector<int> rposn_term_list;  ///< Positional restraint terms, indexing into the original
                                     ///<   restraint apparatus
  std::vector<int> rbond_term_list;  ///< Distance restraint terms, indexing into the original
                                     ///<   restraint apparatus
  std::vector<int> rangl_term_list;  ///< Three-point angle restraint terms, indexing into the
                                     ///<   original restraint apparatus
  std::vector<int> rdihe_term_list;  ///< Four-point dihedral restraint terms, indexing into the
                                     ///<   original restraint apparatus
  std::vector<int> rposn_atoms;      ///< Local indices for atoms subject to each positional
                                     ///<   restraint term in this work unit
  std::vector<int> rbond_i_atoms;    ///< Local indices for I atoms subject to distance restraints
  std::vector<int> rbond_j_atoms;    ///< Local indices for J atoms subject to distance restraints
  std::vector<int> rangl_i_atoms;    ///< Local indices for I atoms subject to angle restraints
  std::vector<int> rangl_j_atoms;    ///< Local indices for J atoms subject to angle restraints
  std::vector<int> rangl_k_atoms;    ///< Local indices for K atoms subject to angle restraints
  std::vector<int> rdihe_i_atoms;    ///< Local indices for I atoms subject to dihedral restraints
  std::vector<int> rdihe_j_atoms;    ///< Local indices for J atoms subject to dihedral restraints
  std::vector<int> rdihe_k_atoms;    ///< Local indices for K atoms subject to dihedral restraints
  std::vector<int> rdihe_l_atoms;    ///< Local indices for L atoms subject to dihedral restraints

  // Constraint groups for this work unit
  std::vector<int> cnst_group_list;    ///< List of constraint groups, indexing into the group
                                       ///<   enumerated in the original topology, that this
                                       ///<   work unit is responsible for enforcing (the work
                                       ///<   unit must be responsible for moving all atoms in
                                       ///<   any such constraint group)
  std::vector<int> sett_group_list;    ///< List of fast rigid water SETTLE groups, indexing
                                       ///<   into the groups enumerated in the original
                                       ///<   topology, assigned to this work unit
  std::vector<int> cnst_group_atoms;   ///< Local indices of atoms in all constrained groups.
                                       ///<   The bounds of this list, delineating separate groups,
                                       ///<  are found in cnst_group_bounds.
  std::vector<int> cnst_group_bounds;  ///< Bounds array for cnst_group_atoms
  std::vector<int> sett_ox_atoms;      ///< Local indices of oxygen atoms in each of this work
                                       ///<   unit's SETTLE groups
  std::vector<int> sett_h1_atoms;      ///< Local indices of the first hydrogen atoms in each of
                                       ///<   this work unit's SETTLE groups
  std::vector<int> sett_h2_atoms;      ///< Local indices of the second hydrogen atoms in each of
                                       ///<   this work unit's SETTLE groups

  // Virtual sites in this work unit
  std::vector<int> virtual_site_list;   ///< List of virtual sites, indexing into the original
                                        ///<   topology's list of virtual sites (i.e. the 5th
                                        ///<   virtual site, which could be atom index 19 in a box
                                        ///<   of TIP4P-Eq water)
  std::vector<int> vsite_atoms;         ///< Local indices for each of the respective virtual
                                        ///<   particles in this work unit's above list
  std::vector<int> vsite_parent_atoms;  ///< Local indices of this work unit's virtual site parent
                                        ///<   atoms
  std::vector<int> vsite_frame2_atoms;  ///< Local indices of this work unit's virtual site second
                                        ///<   frame atoms
  std::vector<int> vsite_frame3_atoms;  ///< Local indices of this work unit's virtual site third
                                        ///<   frame atoms
  std::vector<int> vsite_frame4_atoms;  ///< Local indices of this work unit's virtual site fourth
                                        ///<   frame atoms
  
  ValenceDelegator *vdel_pointer;        ///< Pointer to a delegator managing this object's
                                         ///<   creation
  const AtomGraph *ag_pointer;           ///< Pointer to the topology to which this object pertains
  const RestraintApparatus *ra_pointer;  ///< Pointer to the restraint apparatus to which this
                                         ///<   object pertains
};

/// \brief Build a series of valence work units to cover a topology.
///
/// \param ag  The topology of interest
/// \param ra  Restraints linked to the topology
std::vector<ValenceWorkUnit>
buildValenceWorkUnits(const AtomGraph *ag, const RestraintApparatus *ra,
                      int max_atoms_per_vwu = maximum_valence_work_unit_atoms);

} // namespace topology
} // namespace omni

#endif
