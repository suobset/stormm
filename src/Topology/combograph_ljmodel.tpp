// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
void ComboGraphLJModel::addCombination(const T* lja_in, const T* ljb_in, const int lj_type_count,
                                       const char4* lj_type_names,
                                       const std::vector<PairLJInteraction> &edits) {
  othr_topology_rules.push_back(inferCombiningRule(lja_in, ljb_in, lj_type_count));

  // Compute the combination of the original (base) topology's Lennard-Jones parameters with the
  // new topology's Lennard-Jones parameters, according to the default rule.  The default rule for
  // mixing the two topologies could, in principle, be different than the combining rules apparent
  // in either topology itself.  Each atom type of the new topology controls a column of the
  // resulting matrix, each atom type of the base topology controls a row.
  const NonbondedKit<double> base_nbk = base_ag_pointer->getDoublePrecisionNonbondedKit();
  std::vector<double> combi_lja(base_nbk.n_lj_types * lj_type_count);
  std::vector<double> combi_ljb(base_nbk.n_lj_types * lj_type_count);
  const std::vector<double> base_sig = base_ag_pointer->getLennardJonesSigma<double>();
  const std::vector<double> base_eps = base_ag_pointer->getLennardJonesEpsilon<double>();
  for (int j = 0; j < lj_type_count; j++) {
    const double othr_sig  = sqrt(cbrt(lja_in[j] / ljb_in[j]));
    const double nsig_three = othr_sig * othr_sig * othr_sig;
    const double othr_eps = 0.25 * ljb_in[j] / (nsig_three * nsig_three);
    for (int i = 0; i < base_nbk.n_lj_types; i++) {
      const double ij_eps = sqrt(othr_eps * base_eps[i]);
      double ij_sig;
      switch (default_rule) {
      case VdwCombiningRule::LORENTZ_BERTHELOT:
        ij_sig = 0.5 * (othr_sig + base_sig[i]);
        break;
      case VdwCombiningRule::GEOMETRIC:
        ij_sig = sqrt(othr_sig * base_sig[i]);        
        break;
      case VdwCombiningRule::NBFIX:

        // This invalid case is trapped in the constructor.
        break;
      }
      double ij_sig_six = (ij_sig * ij_sig * ij_sig);
      ij_sig_six *= ij_sig_six;
      const size_t ij_idx = (j * base_nbk.n_lj_types) + i;
      combi_ljb[ij_idx] = 4.0 * ij_eps * ij_sig_six;
      combi_lja[ij_idx] = combi_ljb[ij_idx] * ij_sig_six;
    }      
  }

  // Apply the Lennard-Jones edits
  VdwCombiningRule combi_rule = default_rule;
  if (edits.size() > 0) {
    const std::vector<char4> base_type_names = base_ag_pointer->getAtomTypeNameTable();
    const size_t nedit = edits.size();
    for (size_t i = 0; i < nedit; i++) {

      // The atom types are not required to be found.
      for (int j = 0; j < base_nbk.n_lj_types; j++) {
        for (int k = 0; k < lj_type_count; k++) {
          if ((edits[i].type_a == base_type_names[j] && edits[i].type_b == lj_type_names[k]) ||
              (edits[i].type_b == base_type_names[j] && edits[i].type_a == lj_type_names[k])) {

            // The types to be edited are matched.  Assign the new parameters and note that the
            // combining rule must now be assumed to be "NBFIX".
            const size_t jk_idx = (k * base_nbk.n_lj_types) + j;
            combi_lja[jk_idx] = edits[i].lja;
            combi_ljb[jk_idx] = edits[i].ljb;
            combi_rule = VdwCombiningRule::NBFIX;
          }
        }
      }
    }
  }
  
  // Push the new set sizes, parameter tables, and edits to the object's growing arrays.
  othr_atom_type_counts.push_back(lj_type_count);
  set_rules.push_back(combi_rule);
  set_lja.push_back(combi_lja);
  set_ljb.push_back(combi_ljb);
  set_edits.push_back(edits);

  // Set the source of the new tables as the null pointer.  If there was an actual topology
  // providing the input Lennard-Jones parameters, a pointer to it will be supplied by one of the
  // calling overloads and replace the nullptr.
  othr_ag_pointers.push_back(nullptr);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ComboGraphLJModel::addCombination(const std::vector<T> &lja_in, const std::vector<T> &ljb_in,
                                       const std::vector<char4> &lj_type_names,
                                       const std::vector<PairLJInteraction> &edits) {
  const int lj_type_count = inferLennardJonesTypeCount(lja_in.size(), ljb_in.size(),
                                                       "addCombination");

  // Trap cases where edits are provided but the names of atom types are not.
  if (edits.size() > 0 && static_cast<int>(lj_type_names.size()) != lj_type_count) {
    rtErr("Indicating Lennard-Jones edits requires that the names of atom types be provided "
          "for matching purposes.", "ComboGraphLJModel", "AddCombination");
  }
  addCombination(lja_in.data(), ljb_in.data(), lj_type_count, lj_type_names.data(), edits);
}

} // namespace topology
} // namespace stormm
