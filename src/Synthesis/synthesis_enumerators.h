// -*-c++-*-
#ifndef STORMM_SYNTHESIS_ENUMERATORS_H
#define STORMM_SYNTHESIS_ENUMERATORS_H

namespace stormm {
namespace synthesis {
  
/// \brief Enumerate the different kinds of non-bonded work.  All-to-all calculations for systems
///        in isolated boundary conditions have two options, based on how large the system or
///        systems are.  Neighbor list-based calculations for systems in periodic boundary
///        conditions have only one option.
enum class NbwuKind {
  TILE_GROUPS,  ///< Work units process small lists of all-to-all non-bonded tiles
  SUPERTILES,   ///< Work units process large, regular arrangements of all-to-all non-bonded tiles
  HONEYCOMB,    ///< Neighbor list-based spatial decomposition, with each work unit being a
                ///<   locally clustered patch of atoms
  UNKNOWN       ///< Placeholder for the non-bonded work type until it is properly set
};

/// \brief Enumerate the different tasks that this work unit can perform.
enum class VwuTask {
  BOND = 0,  ///< Harmonic bond interactions
  ANGL,      ///< Harmonic angle interactions
  DIHE,      ///< Cosine-based dihedral interactions
  UBRD,      ///< Urey-Bradley harmonic stretching angle terms
  CBND,      ///< Composite bond term count, including harmonic bonds and Urey-Bradley angle
             ///<   stretching terms
  CIMP,      ///< CHARMM improper dihedrals
  CDHE,      ///< Composite dihedral term count (includes multiple dihedrals that share the same
             ///<   atoms as well as CHARMM impropers in one list)
  CMAP,      ///< CMAP spline-based surface potentials
  INFR14,    ///< Inferred 1:4 attenuated non-bonded interactions
  RPOSN,     ///< Positional restraints
  RBOND,     ///< Distance-based restraints
  RANGL,     ///< Three-point angle restraints
  RDIHE,     ///< Four-point dihedral restraints
  VSITE,     ///< Virtual sites
  SETTLE,    ///< SETTLE analytic (fast, rigid water) constraints
  CGROUP,    ///< Constraint groups (hub-and-spoke, Hydrogen-on-heavy atom constraint clusters)
  ALL_TASKS  ///< Total number of distinct tasks that the ValenceWorkUnit can perform (for array
             ///<   sizing purposes)
};

/// \brief Define the space to be allocated to any ValenceWorkUnit abstract (see the enumerator
///        below).  This number must grow in step with the enumerator.
constexpr int vwu_abstract_length = 32;
  
/// \brief Create a map of the contents of the ValenceWorkUnit instruction set abstract, a list of
///        which will appear in any AtomGraphSynthesis object (vwu_instruction_sets member
///        variable) for use by HPC applications.  The abstract is just a series of limits stored
///        in an array of integers, one abstract per stride of tuples, based on this enumerator.
enum class VwuAbstractMap {
  IMPORT = 0,  ///< Atom import indexing limits (these are the bounds of a list of atom
               ///<   indices to import, not the low end of a consecutive series of topological
               ///<   atoms to import)
  MANIPULATE,  ///< Limits for the imported atom movement mask (this will be about 1/32nd the size
               ///<   of the import list)
  CBND,        ///< Composite bond instruction index limits
  ANGL,        ///< Harmonic angle instruction limits
  CDHE,        ///< Composite dihedral instruction limits
  CMAP,        ///< CMAP instruction limits
  INFR14,      ///< Inferred 1:4 attenuated non-bonded instruction limits
  RPOSN,       ///< Positional restraint instruction limits
  RBOND,       ///< Distance restraint instruction limits
  RANGL,       ///< Three-point angle restraint instruction limits
  RDIHE,       ///< Four-point dihedral angle restraint instruction limits
  VSITE,       ///< Virtual site placement and force transfer instruction limits
  SETTLE,      ///< SETTLE constraint group instruction limits
  CGROUP,      ///< Hub-and-spoke constraint group instruction limits
  CBND_NRG,    ///< Limits of the composite bond energy accumulation bit mask
  ANGL_NRG,    ///< Limits of the harmonic angle energy accumulation bit mask
  CDHE_NRG,    ///< Limits of the composite dihedral energy accumulation bit mask
  CMAP_NRG,    ///< Limits of the CMAP energy accumulation bit mask
  INFR14_NRG,  ///< Limits of the 1:4 attenuated interaction energy accumulation bit mask
  RPOSN_NRG,   ///< Limits of the positional restraint energy accumulation bit mask
  RBOND_NRG,   ///< Limits of the distance restraint energy accumulation bit mask
  RANGL_NRG,   ///< Limits of the three-point angle restraint energy accumulation bit mask
  RDIHE_NRG,   ///< Limits of the four-point dihedral restraint energy accumulation bit mask
  SYSTEM_ID,   ///< System to which the valence work unit pertains, referencing the list of
               ///<   systems within the AtomGraphSynthesis.  Only the x member of this element
               ///<   has meaning.
};

/// \brief There are two modes in which an evaluation of ValenceWorkUnit instructions can proceed.
enum class VwuGoal {
  ACCUMULATE,     ///< The objective is a global ccumulation of forces.  This case, the bit
                  ///<   strings that dictact whether a ValenceWorkUnit logs the energy due to a
                  ///<   particular interaction will dictate whether the ValenceWorkUnit
                  ///<   evaluates a term at all.  The local force accumulations will be added
                  ///<   back to a global array.  Particles cannot be moved under this mode.
  MOVE_PARTICLES  ///< The objective to is to move particles.  In this mode, global energy
                  ///<   accumulators will still be valid as the work units continue to follow
                  ///<   their bit strings dictating which interactions they are responsible for
                  ///<   logging, but forces accumulated locally may contain redundant effects
                  ///<   needed to carry out a local atom move and therefore cannot be pooled at
                  ///<   the end of the calculation.  Atoms will move and their global positions
                  ///<   will be updated according to separate bit strings dictating which work
                  ///<   units are responsible for each update.
};

} // namespace synthesis
} // namespace stormm

#endif
