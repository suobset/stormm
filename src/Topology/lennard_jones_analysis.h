// -*-c++-*-
#ifndef LENNARD_JONES_ANALYSIS_H
#define LENNARD_JONES_ANALYSIS_H

#include <map>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"

namespace stormm {
namespace topology {

using energy::VdwCombiningRule;
using synthesis::AtomGraphSynthesis;

/// \brief A struct to encode two atom types and the Lennard-Jones parameters by which they
///        interact.
struct PairLJInteraction {
public:

  /// \brief The default constructor will initialize atom types and parameters, but they are almost
  ///        certain to register as invalid in later searches if left unspecified.
  /// \{
  PairLJInteraction(const char4 type_a_in = { ' ', ' ', ' ', ' ' },
                    const char4 type_b_in = { ' ', ' ', ' ', ' ' }, double lja_in = 0.0,
                    double ljb_in = 0.0);

  /// \brief The default copy and move constructions as well as assignment operators are valid,
  ///        making this object easy to manipulate.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of an assignment statement
  /// \{
  PairLJInteraction(const PairLJInteraction &original) = default;
  PairLJInteraction(PairLJInteraction &&original) = default;
  PairLJInteraction& operator=(const PairLJInteraction &original) = default;
  PairLJInteraction& operator=(PairLJInteraction &&original) = default;
  /// \}

  // All member variables are public.
  char4 type_a;  ///< Atom type of the first atom in the pair
  char4 type_b;  ///< Atom type of the second atom in the pair
  double lja;    ///< Lennard-Jones A coefficient for the interaction, as in U = A/r^12 - B/r^6
  double ljb;    ///< Lennard-Jones B coefficient for the interaction
};
  
/// \brief A class to sort and hold details of a topology's Lennard-Jones interactions.  The
///        Lennard-Jones table is parsed for a combining rule, tables of sigma and epsilon
///        parameters for self-interactions are kept, and a list of pair-specific sigma and
///        epsilon parameters are stored according to the type indices they modify.
class LennardJonesAnalysis {
public:

  /// \brief The constructor accepts a topology.  Supplying additional topologies will expand the
  ///        analysis to cover multiple systems and create a consensus of their respective
  ///        Lennard-Jones parameters with arrays to index the original topologies' atoms into the
  ///        unified tables.
  /// \{
  LennardJonesAnalysis(const NonbondedKit<double> &nbk,
                       const std::vector<std::vector<char4>> &atom_type_aliases_in);
  /// \}

  /// Composed of Standard Template Library objects with no pointers to repair, the default copy
  /// and move constructors are valid for this object, as are the default copy and move assignment
  /// operators.
  ///
  /// \param original  A pre-existing object to copy or move
  /// \param other     An object to place on the right hand side of the assignment statement
  /// \{
  LennardJonesAnalysis(const LennardJonesAnalysis &original) = default;
  LennardJonesAnalysis(LennardJonesAnalysis &&original) = default;
  LennardJonesAnalysis& operator=(const LennardJonesAnalysis &original) = default;
  LennardJonesAnalysis& operator=(LennardJonesAnalysis &&original) = default;
  /// \}

  /// \brief Get the number of Lennard-Jones indices, the number of distinct Lennard Jones
  ///        parameter sets which are needed to describe all interactions.
  int getLJTypeCount() const;

  /// \brief Get the number of atom types, including all types with equivalent Lennard-Jones
  ///        parameters.
  int getAtomTypeCount() const;

  /// \brief Get the most prevalent Lennard-Jones rule, as judged by the formula that fits the
  ///        most combinations of distinct sigma and epsilon parameters.
  VdwCombiningRule getMostPrevalentCombiningRule() const;

  /// \brief Get the atom type names associated with a particular Lennard-Jones parameter set,
  ///        type index in one of the topologies, or aliases of a particular type name.
  ///
  /// Overloaded:
  ///   - Indicate the index of the Lennard-Jones type within the consensus parameter tables
  ///   - Indicate the index of the topology and the Lennard-Jones type index within that topology
  ///   - Indicate the topology by pointer and the Lennard-Jones type index within that topology
  ///   - Indicate the sigma and epsilon self-interaction parameters, to within 1.0e-4 precision
  ///   - Indicate the name of one of the atom types so that all aliases can be found
  ///
  /// \param consensus_index  Index of the Lennard-Jones interaction type with the consensus table
  /// \param ag_query_index   Index of the original AtomGraph within the list of referenced
  ///                         topologies (ag_pointers)
  /// \param lj_type_index    Lennard-Jones interaction index from within the indicated topology
  /// \param ag_query         The topology to search for within the list of referenced topologies
  ///                         (failing to find such a topology produces a runtime error)
  /// \param sigma_query      The self-interaction sigma parameter of interest
  /// \param epsilon_query    The self-interaction epsilon parameter of interest
  /// \param tolerance        The tolerance for finding a Lennard-Jones parameter match, if
  ///                         supplying sigma and epsilon explicitly
  /// \param atom_type_query  Name of one of the atom type aliases, provided so that all others
  ///                         can be found
  /// \{
  const std::vector<char4>& getLJAliases(int consensus_index) const;

  const std::vector<char4>& getLJAliases(double sigma_query, double epsilon_query,
                                         double tolerance = 1.0e-4) const;

  const std::vector<char4>& getLJAliases(const char4 atom_type_query) const;
  /// \}  

  /// \brief Get the Lennard-Jones self interaction sigma for a particular interaction index from
  ///        within the consensus tables.
  ///
  /// \param consensus_index  The interaction type index of interest.  It may not be obvious what
  ///                         this is if the analysis includes more than one topology, but it can
  ///                         still be useful to have an accessor that accepts this sort of input.
  double getLJSigma(int consensus_index) const;
  
  /// \brief Get the Lennard-Jones self interaction epsilon for a particular interaction index from
  ///        within the consensus tables.
  ///
  /// \param consensus_index  The interaction type index of interest
  double getLJEpsilon(int consensus_index) const;

  /// \brief Get the Lennard-Jones sigma and epsilon parameters for the self interaction of a
  ///        specific interaction index within the consensus tables or atom type name.  The sigma
  ///        parameter appears in the "x" member of the result and the epsilon parameter in the
  ///        "y" member.
  ///
  /// Overloaded:
  ///   - Specify the index from within the consensus tables
  ///   - Specify one of the names of an atom type alias bearing the parameters of interest
  ///
  /// \param consensus_index
  /// \{
  double2 getLJParameters(int consensus_index) const;
  double2 getLJParameters(const char4 atom_type_query) const;
  /// \}

  /// \brief Add a new set of Lennard-Jones interactions to the list and expand the internal tables
  ///        with all unique parameters it may contain, as well as new atom type names for existing
  ///        parameters.  The new topology will be checked to ensure that it does not redefine the
  ///        parameters associated with known atom type names, or contradict the prevailing
  ///        combining rules of topologies already incorporated into the tables.
  ///
  /// \param nbk                The non-bonded abstract of the topology responsible for the new
  ///                           Lennard-Jones parameter set.  While this library is included in
  ///                           the AtomGraph library itself and therefore cannot reference the
  ///                           AtomGraph class in any member functions, the abstracts can be
  ///                           passed in much the same way that other information reduced to its
  ///                           basic form can be handled by this library.
  /// \param othr_sigma         Sigma parameters determined for the new Lennard-Jones parameters
  /// \param othr_epsilon       Epsilon parameters determined for the new Lennard-Jones parameters
  /// \param othr_type_aliases  Type aliases for the new set of Lennard-Jones parameters
  void addSet(const NonbondedKit<double> &nbk,
              const std::vector<std::vector<char4>> &othr_type_aliases);

private:

  /// Total number of unique Lennard-Jones indices in the consensus set
  int lj_type_count;

  /// Total number of unique atom type names in the consensus set.  Multiple atom type names may
  /// take the same parameters, for the self interaction as well as in combination with other atom
  /// types, making atom_type_count >= lj_type_count.
  int atom_type_count;

  /// Set as the most common combining rule (goemetric or Lorentz-Berthelot) that is apparent in
  /// the first topology used to contruct the object.  Even if the Lennard-Jones parameters of the
  /// topology technically fall under the NBFIX enumeration, there will likely be a rule which
  /// the majority of the off-diagonal elements adhere to.  NBFix details will be identified from
  /// with the original topology and any subsequent topologies added to the object, then added to
  /// the list of edits (see below).
  VdwCombiningRule prevalent_rule;

  /// Set as the absolute rule for the combinations of Lennard-Jones parameters.  This will only
  /// register Lorentz-Berthelot or geometric rules if all parameters in the matrix obey them.
  /// Otherwise, the rule will be declared NBFIX.
  VdwCombiningRule absolute_rule;
  
  /// Sigma parameters for the self-interactions of each type found within any of the topologies
  /// referenced in ag_pointers
  std::vector<double> sigma;

  /// Epsilon parameters for the self-interactions of each type found within any of the topologies
  /// referenced in ag_pointers
  std::vector<double> epsilon;

  /// The mixed parameters for consensus tables, taken directly from the referenced topologies if
  /// possible or constructed using the prevailing mixing rule.
  /// \{
  std::vector<double> lja_coeff;
  std::vector<double> ljb_coeff;
  /// \}

  /// List of all NBFix pair-specifc interactions that break the standard combining rule found in
  /// prevalent_rule above
  std::vector<PairLJInteraction> edits;

  /// Bounds array for ag_index_origins above
  std::vector<int> ag_index_origins_bounds;

  /// Arrays of aliases for each atom type name, i.e. if N and NA have the same Lennard-Jones
  /// parameters and no NBFix detail applies to one of them in particular, they are aliases.
  std::vector<std::vector<char4>> atom_type_aliases;

  /// A map to indicate the Lennard-Jones type index, within the consensus tables, for any
  /// particular atom type name.  The names are converted to uint prior to becoming map keys, to
  /// expedite comparisons and thus lookup.
  std::map<uint, int> atom_type_map;
  
  /// Maps between the Lennard-Jones interaction indices in each referenced topology and the
  /// consensus parameter tables.  This table should be used when the Lennard-Jones interaction
  /// type index in a specific topology is known (the outer dimension is the number of referenced
  /// topologies), and the goal is to find the index of that type as represented in the consensus
  /// tables.
  std::vector<std::vector<int>> set_to_consensus_map;

  /// Maps between the Lennard-Jones interaction types in the consensus tables and their instances
  /// in each of the referenced topologies.  The "x" member of each tuple indicates the index of
  /// some referenced topology and the "y" member indicates the Lennard-Jones interaction index in
  /// that topology.  The kth list in the array pertains to all instances of the kth consensus
  /// Lennard-Jones atom type.  Keeping this as a vector of vectors eliminates the need for a
  /// bounds array and thus makes it more easily extensible with new topologies.
  std::vector<std::vector<int2>> consensus_to_set_map;
};

/// \brief Determine the Lennard-Jones combining rule in effect.  The function will accept either
///        single- or double-precision data, but internally it uses double-precision calculations.
///
/// Overloaded:
///   - Provide the A and B coefficient arrays by C-style arrays with trusted length
///   - Provide the A and B coefficient arrays as Standard Template Library vectors
///   - Provide the A and B coefficient arrays as Hybrid objects
///
/// \param lja             Array of Lennard-Jones A coefficients
/// \param ljb             Array of Lennard-Jones B coefficients
/// \param lj_type_count   The number of Lennard-Jones types, indicating the trusted lengths of
///                        lja and ljb by the number's square
/// \param policy          Protocol in the event that there are only one (fewer) atom types
/// \param seek_prevalent  A rare boolean input that, if set to TRUE, will avoid declaring a
///                        parameter set "NBFIX" and instead declare it to be "GEOMETRIC" or
///                        "LORENTZ_BERTHELOT" depending on whether more cases of either combining
///                        rule can be found.  A parameter set might still be decalred "NBFIX" if
///                        no instances of the other combining rules can describe any off-diagonal
///                        interactions in the matrix.
/// \{
template <typename T>
VdwCombiningRule inferCombiningRule(const T* lja, const T* ljb, int lj_type_count,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);

template <typename T>
VdwCombiningRule inferCombiningRule(const std::vector<T> &lja, const std::vector<T> &ljb,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);

template <typename T>
VdwCombiningRule inferCombiningRule(const Hybrid<T> &lja, const Hybrid<T> &ljb,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);
/// \}

/// \brief Locate the Lennard-Jones types where pair-specific combining rules apply.  This will
///        return a Standard Template Library vector of booleans reading TRUE if the Lennard-Jones
///        type of each index has pair-specific interactions.  Overloading and descriptions of
///        input parameters follow from inferCombiningRule(), above.
/// \{
template <typename T>
std::vector<bool> findPairSpecificParticipation(const T* lja, const T* ljb, int lj_type_count);

template <typename T>
VdwCombiningRule findPairSpecificParticipation(const std::vector<T> &lja,
                                               const std::vector<T> &ljb);

template <typename T>
VdwCombiningRule findPairSpecificParticipation(const Hybrid<T> &lja, const Hybrid<T> &ljb);
/// \}
  
} // namespace topology
} // namespace stormm

#include "lennard_jones_analysis.tpp"

#endif
