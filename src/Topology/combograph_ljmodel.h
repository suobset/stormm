// -*-c++-*-
#ifndef STORMM_COMBOGRAPH_LJMODEL_H
#define STORMM_COMBOGRAPH_LJMODEL_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "atomgraph.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_analysis.h"

namespace stormm {
namespace topology {

using energy::getEnumerationName;
using parse::operator==;
using synthesis::AtomGraphSynthesis;
using topology::AtomGraph;
  
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
  //                                                                                               
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
  
/// \brief An object to manage the combination of Lennard-Jones parameters from two topologies.
///        This can be used to insert "NBFix" pair-specific (off-rule) parameters where neither
///        topology by itself might contain the requisite information.  It can also be used to
///        build parameter tables for a union of the two underlying topologies and produce arrays
///        indices into that table for a new topology.
class ComboGraphLJModel {
public:

  /// \brief The constructor accepts two topologies or two Lennard-Jones parameter sets.  Edits
  ///        to the Lennard-Jones parameters can be provided to the constructor or entered later.
  ///
  /// \param poly_ag_othr  
  /// \{
  ComboGraphLJModel(const AtomGraph *base_ag_in, const AtomGraph *ag_othr_in = nullptr,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &base_ag_in, const AtomGraph &ag_othr_in,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &base_ag_in, const AtomGraphSynthesis &poly_ag_othr,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});
  
  ComboGraphLJModel(const AtomGraph *base_ag_in, const std::vector<double> &lj_a_in,
                    const std::vector<double> &lj_b_in,
                    const std::vector<char4> &lj_type_names = {},
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &base_ag_in, const std::vector<double> &lj_a_in,
                    const std::vector<double> &lj_b_in,
                    const std::vector<char4> &lj_type_names = {},
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});
  /// \}

  /// \brief Get the number of Lennard-Jones parameter matrix sets.
  int getSetCount() const;

  /// \brief Get the default combining rule for mixing parameters of different topologies.
  VdwCombiningRule getCombiningRule() const;

  /// \brief Get the size of a particular combination matrix.  The matrix has as many rows as the
  ///        base topology has atom types and as many columns as the new (other) topology has
  ///        atom types.
  ///
  /// \param index  The combination of interest
  int2 getMatrixSize(int index) const;

  /// \brief Get the matrix of Lennard-Jones A coefficients between the base topology and another
  ///        topology.  The base atom types control rows of the matrix, while the other topology's
  ///        atom types control columns.  The matrix has column format, like others in STORMM.
  ///
  /// \param index  The combination of interest
  const std::vector<double>& getACoefficients(int index) const;
  
  /// \brief Get the matrix of Lennard-Jones B coefficients between the base topology and another
  ///        topology.  The base atom types control rows of the matrix, while the other topology's
  ///        atom types control columns.  The matrix has column format, like others in STORMM.
  ///
  /// \param index  The combination of interest
  const std::vector<double>& getBCoefficients(int index) const;

  /// \brief Add an interaction matrix to the list of combinations.
  ///
  /// Overloaded:
  ///   - Provide a combination of vectors or C-style arrays for the new parameter set
  ///   - Provide the new topology by const pointer or const reference
  ///
  /// \param lja_in         Array of Lennard-Jones A coefficients, i.e. U = A/r^12 + B/r^6, to be
  ///                       combined per the object's default rule
  /// \param ljb_in         Array of Lennard-Jones B coefficients, to be combined per the object's
  ///                       default rule
  /// \param lj_type_count  The number of Lennard-Jones types, implying a trusted length for lja_in
  ///                       and ljb_in, if provided as C-style arrays
  /// \param ag_new         The topology to incorporate.  A null pointer will abort the inclusion.
  /// \param edits          A list of pair-specific Lennard-Jones parameters to apply
  /// \{
  template <typename T>
  void addCombination(const T* lja_in, const T* ljb_in, int lj_type_count,
                      const char4* lj_type_names = nullptr,
                      const std::vector<PairLJInteraction> &edits = {});

  template <typename T>
  void addCombination(const std::vector<T> &lja_in, const std::vector<T> &ljb_in,
                      const std::vector<char4> &lj_type_names = {},
                      const std::vector<PairLJInteraction> &edits = {});

  void addCombination(const AtomGraph *ag_othr, const std::vector<PairLJInteraction> &edits = {});

  void addCombination(const AtomGraph &ag_othr, const std::vector<PairLJInteraction> &edits = {});
  /// \}
  
private:

  /// The number of Lennard-Jones parameter sets being tracked by this object, and the length of
  /// arrays set_sizes, set_lja, set_ljb, set_edits, and new_ag_pointers
  int set_count;                             

  /// The default combining rule to be applied to interactions in each topology.
  VdwCombiningRule default_rule;
  
  /// The rule governing Lennard-Jones parameter mixing in the original topology.  This will be
  /// inferred from the parameter matrices contained therein.
  VdwCombiningRule base_topology_rule;

  /// The number of Lennard-Jones types in the base topology
  int base_atom_type_count;
  
  /// The numbers of Lennard-Jones atom types in each additional topology.
  std::vector<int> othr_atom_type_counts;

  /// Rules inferred for Lennard-Jones parameter mixing in all other topologies.  The ith index of
  /// this array corresponds to the secondary topology in the ith combination.
  std::vector<VdwCombiningRule> othr_topology_rules;
  
  /// The rules by which each combination can be taken to occur.
  std::vector<VdwCombiningRule> set_rules;
  
  /// Lennard-Jones A coefficients for each set, given as a column-format matrix with a number of
  /// rows given in the 'x" member of the set_sizes array and a number of columns given in the "y"
  /// member.
  std::vector<std::vector<double>> set_lja;

  /// Lennard-Jones B coefficients for each parameter matrix
  std::vector<std::vector<double>> set_ljb;

  /// The edits that occurred in each set.  The presence of any edits (a nonzero array length) will
  /// change the combining rule for that set to "NBFIX."
  std::vector<std::vector<PairLJInteraction>> set_edits;

  /// Pointer to the original topology, whose atom types (or types derived thereof) feed into each
  /// row of the resulting parameter matrices.
  AtomGraph *base_ag_pointer;

  /// Pointers to new topologies, whose atom types (or types derived thereof) appear in the columns
  /// of the resulting parameter matrices.
  std::vector<AtomGraph*> othr_ag_pointers;

  /// \brief Validate the index of some requested matrix set.
  ///
  /// \param index  The index of interest
  void validateSetIndex(int index) const;
};
  
} // namespace topology
} // namespace stormm

#include "combograph_ljmodel.tpp"

#endif
