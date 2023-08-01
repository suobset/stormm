#include "copyright.h"
#include "combograph_ljmodel.h"

namespace stormm {
namespace topology {

using topology::NonbondedKit;
  
//-------------------------------------------------------------------------------------------------
PairLJInteraction::PairLJInteraction(const char4 type_a_in, const char4 type_b_in,
                                     const double lja_in, const double ljb_in) :
    type_a{type_a_in}, type_b{type_b_in}, lja{lja_in}, ljb{ljb_in}
{}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph *base_ag_in, const AtomGraph *ag_othr_in,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    set_count{0},
    default_rule{default_rule_in},
    base_topology_rule{VdwCombiningRule::LORENTZ_BERTHELOT},
    base_atom_type_count{0},
    othr_atom_type_counts{},
    othr_topology_rules{},
    set_lja{}, set_ljb{}, set_edits{},
    base_ag_pointer{const_cast<AtomGraph*>(base_ag_in)},
    othr_ag_pointers{}
{
  if (base_ag_pointer == nullptr) {
    rtErr("A base topology must be specified in order to create parameter combination matrices.",
          "ComboGraphLJModel");
  }
  if (ag_othr_in == nullptr) {
    return;
  }

  // The default mixing rule cannot be NBFix (this would imply that there is, in effect, no rule).
  switch (default_rule) {
  case VdwCombiningRule::LORENTZ_BERTHELOT:
  case VdwCombiningRule::GEOMETRIC:
    break;    
  case VdwCombiningRule::NBFIX:
    rtErr("A default combining rule of " + getEnumerationName(default_rule) + " would imply that "
          "there is no combining rule and is therefore invalid.", "ComboGraphLJModel");
  }
  
  // Infer the Lennard-Jone combination rules present in the original topology.
  const NonbondedKit<double> base_nbk = base_ag_pointer->getDoublePrecisionNonbondedKit();
  base_topology_rule = inferCombiningRule(base_nbk.lja_coeff, base_nbk.ljb_coeff,
                                          base_nbk.n_lj_types);
  base_atom_type_count = base_nbk.n_lj_types;
  addCombination(ag_othr_in, edits);
}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &base_ag_in, const AtomGraph &ag_othr_in,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(base_ag_in.getSelfPointer(), ag_othr_in.getSelfPointer(), default_rule_in,
                      edits)
{}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &base_ag_in,
                                     const AtomGraphSynthesis &poly_ag_othr,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(base_ag_in.getSelfPointer(), nullptr, default_rule_in, edits)
{}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph *base_ag_in,
                                     const std::vector<double> &lj_a_in,
                                     const std::vector<double> &lj_b_in,
                                     const std::vector<char4> &lj_type_names,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(base_ag_in, nullptr, default_rule_in, edits)
{
  addCombination(lj_a_in, lj_b_in, lj_type_names, edits);
}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &base_ag_in,
                                     const std::vector<double> &lj_a_in,
                                     const std::vector<double> &lj_b_in,
                                     const std::vector<char4> &lj_type_names,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
  ComboGraphLJModel(base_ag_in.getSelfPointer(), lj_a_in, lj_b_in, lj_type_names, default_rule_in,
                    edits)
{}

//-------------------------------------------------------------------------------------------------
int ComboGraphLJModel::getSetCount() const {
  return set_count;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule ComboGraphLJModel::getCombiningRule() const {
  return default_rule;
}

//-------------------------------------------------------------------------------------------------
int2 ComboGraphLJModel::getMatrixSize(const int index) const {
  validateSetIndex(index);
  return {base_atom_type_count, othr_atom_type_counts[index] };
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ComboGraphLJModel::getACoefficients(const int index) const {
  validateSetIndex(index);
  return set_lja[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ComboGraphLJModel::getBCoefficients(const int index) const {
  validateSetIndex(index);
  return set_ljb[index];
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::validateSetIndex(const int index) const {
  if (index < 0 || index >= set_count) {
    rtErr("Combination index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(set_count) + " topology combinations.", "ComboGraphLJModel",
          "validateSetIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::addCombination(const AtomGraph *ag_othr,
                                       const std::vector<PairLJInteraction> &edits) {
  if (ag_othr == nullptr) {
    return;
  }
  const NonbondedKit<double> othr_nbk = ag_othr->getDoublePrecisionNonbondedKit();
  addCombination(othr_nbk.lja_coeff, othr_nbk.ljb_coeff, othr_nbk.n_lj_types);

  // Replace the nullptr appended by the topology-free overload of this function.
  othr_ag_pointers.back() = const_cast<AtomGraph*>(ag_othr);
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::addCombination(const AtomGraph &ag_othr,
                                       const std::vector<PairLJInteraction> &edits) {
  addCombination(ag_othr.getSelfPointer(), edits);
}

} // namespace topology
} // namespace stormm
