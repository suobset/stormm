#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_analysis.h"
#include "lennard_jones_analysis.h"
#include "topology_util.h"

namespace stormm {
namespace topology {

using parse::char4ToString;
using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
PairLJInteraction::PairLJInteraction(const char4 type_a_in, const char4 type_b_in,
                                     const double lja_in, const double ljb_in) :
    type_a{type_a_in}, type_b{type_b_in}, lja{lja_in}, ljb{ljb_in}
{}

//-------------------------------------------------------------------------------------------------
LennardJonesAnalysis::
LennardJonesAnalysis(const NonbondedKit<double> &nbk,
                     const std::vector<std::vector<char4>> &atom_type_aliases_in) :
    lj_type_count{nbk.n_lj_types},
    atom_type_count{0}, prevalent_rule{}, absolute_rule{}, sigma{},
    epsilon{}, lja_coeff{}, ljb_coeff{}, edits{},
    atom_type_aliases{atom_type_aliases_in},
    atom_type_map{}, set_to_consensus_map{}, consensus_to_set_map{}
{
  // Factor out the sigma and epsilon parameters
  sigma.resize(nbk.n_lj_types);
  epsilon.resize(nbk.n_lj_types);
  for (int i = 0; i < lj_type_count; i++) {
    const size_t diag_idx = (i * nbk.n_lj_types) + i;
    sigma[i] = sqrt(cbrt(nbk.lja_coeff[diag_idx] / nbk.ljb_coeff[diag_idx]));
    epsilon[i] = 0.25 * nbk.ljb_coeff[diag_idx] / pow(sigma[i], 6.0);
  }

  // Group the atom types involved in each Lennard-Jones interaction
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const int n_type_names = static_cast<int>(atom_type_aliases[i].size());
    atom_type_count += n_type_names;
    for (int j = 0; j < n_type_names; j++) {
      std::map<uint, int>::iterator it = atom_type_map.find(char4ToUint(atom_type_aliases[i][j]));
      if (it != atom_type_map.end()) {
        rtErr("Atom type " + char4ToString(atom_type_aliases[i][j]) + " controls multiple "
              "Lennard-Jones parameter sets in the consensus tables.", "LennardJonesAnalysis");
      }
      atom_type_map[char4ToUint(atom_type_aliases[i][j])] = i;
    }
  }

  // Extract the Lennard-Jones A and B coefficients directly from the first topology.
  prevalent_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types,
                                              ExceptionResponse::DIE, true);
  absolute_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types,
                                              ExceptionResponse::DIE, false);
  const size_t nlj_squared = nbk.n_lj_types * nbk.n_lj_types;
  lja_coeff.resize(nlj_squared);
  ljb_coeff.resize(nlj_squared);
  for (size_t i = 0; i < nlj_squared; i++) {
    lja_coeff[i] = nbk.lja_coeff[i];
    ljb_coeff[i] = nbk.ljb_coeff[i];
  }

  // Identify any edits that would be characteristic of NBFix.
  for (int i = 1; i < nbk.n_lj_types; i++) {
    for (int j = 0; j < i; j++) {
      const size_t ji_idx = (nbk.n_lj_types * i) + j;
      double ij_eps = epsilon[i] * epsilon[j];
      ij_eps = (ij_eps > constants::tiny) ? sqrt(ij_eps) : 0.0;
      double ij_sig;
      switch (prevalent_rule) {
      case VdwCombiningRule::GEOMETRIC:
        ij_sig = sigma[i] * sigma[j];
        ij_sig = (ij_sig > constants::tiny) ? sqrt(ij_sig) : 0.0;
        break;
      case VdwCombiningRule::LORENTZ_BERTHELOT:
        ij_sig = 0.5 * (sigma[i] + sigma[j]);
        break;
      case VdwCombiningRule::NBFIX:
        edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                          nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        continue;
      }
      double ij_sig_six = ij_sig * ij_sig * ij_sig;
      ij_sig_six *= ij_sig_six;
      const double pred_ljb = 4.0 * ij_eps * ij_sig_six;
      const double pred_lja = pred_ljb * ij_sig_six;
      if (pred_ljb < constants::small) {
        if (fabs(pred_lja - nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs(pred_ljb - nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        }
      }
      else {
        if (fabs((pred_lja - nbk.lja_coeff[ji_idx]) / nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs((pred_ljb - nbk.ljb_coeff[ji_idx]) / nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        }
      }
    }
  }

  // Set the origins of each Lennard-Jones type index in the original topology.
  set_to_consensus_map.resize(1);
  set_to_consensus_map[0].resize(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    set_to_consensus_map[0][i] = i;
  }
  
  // Initialize the consensus map, an explanation of which Lennard-Jones type index in the
  // consensus tables each Lennard-Jones index in one of the included topologies references,
  // and a reciprocal lookup for which Lennard-Jones indices in one exemplary topology 
  consensus_to_set_map.resize(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    consensus_to_set_map[i].resize(1);
    consensus_to_set_map[i][0] = { 0, i };
  }
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getLJTypeCount() const {
  return lj_type_count;
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getAtomTypeCount() const {
  return atom_type_count;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule LennardJonesAnalysis::getMostPrevalentCombiningRule() const {
  return prevalent_rule;
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const int consensus_index) const {
  if (consensus_index < 0 || consensus_index >= lj_type_count) {
    rtErr("Index " + std::to_string(consensus_index) + " is invalid for a collection of " +
          std::to_string(lj_type_count) + " unique Lennard-Jones interactions.",
          "LennardJonesAnalysis", "getLJAliases");
  }
  return atom_type_aliases[consensus_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const double sigma_query,
                                                             const double epsilon_query,
                                                             const double tolerance) const {

  // Find a match and return immediately, if possible.
  int best_match = -1;
  double best_error = 4.0 * tolerance * tolerance;
  for (int i = 0; i < lj_type_count; i++) {
    if (fabs(sigma[i] - sigma_query) < tolerance && fabs(epsilon[i] - epsilon_query) < tolerance) {
      const double dsig = sigma[i] - sigma_query;
      const double deps = epsilon[i] - epsilon_query;
      const double derr = (dsig * dsig) + (deps * deps);
      if (derr < best_error) {
        best_match = i;
        best_error = derr;
      }
    }
  }
  if (best_match == -1) {
    rtErr("No sigma and epsilon parameter combination " +
          realToString(sigma_query, 7, 4, NumberFormat::STANDARD_REAL) + " (sigma) and " +
          realToString(epsilon_query, 7, 4, NumberFormat::STANDARD_REAL) + " (epsilon) was "
          "found in the consensus tables.");
  }
  return atom_type_aliases[best_match];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const char4 atom_type_query) const {
  for (int i = 0; i < lj_type_count; i++) {
    const size_t n_names = atom_type_aliases[i].size();
    for (size_t j = 0; j < n_names; j++) {
      if (atom_type_aliases[i][j] == atom_type_query) {
        return atom_type_aliases[i];
      }
    }
  }
  rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
        "consensus tables.", "LennardJonesAnalysis", "getLJAliases");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJSigma(const int consensus_index) const {
  return sigma[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJEpsilon(const int consensus_index) const {
  return epsilon[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const int consensus_index) const {
  return { sigma[consensus_index], epsilon[consensus_index] };
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const char4 atom_type_query) const {
  for (int i = 0; i < lj_type_count; i++) {
    const size_t n_names = atom_type_aliases[i].size();
    for (size_t j = 0; j < n_names; j++) {
      if (atom_type_aliases[i][j] == atom_type_query) {
        return { sigma[i], epsilon[i] };
      }
    }
  }
  rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
        "consensus tables.", "LennardJonesAnalysis", "getLJSParameters");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addSet(const NonbondedKit<double> &nbk,
                                  const std::vector<std::vector<char4>> &othr_type_aliases) {

  // Obtain sigma and epsilon parameters for the new set.
  std::vector<double> othr_sigma(nbk.n_lj_types);
  std::vector<double> othr_epsilon(nbk.n_lj_types);
  for (int i = 0; i < lj_type_count; i++) {
    const size_t diag_idx = (i * nbk.n_lj_types) + i;
    othr_sigma[i] = sqrt(cbrt(nbk.lja_coeff[diag_idx] / nbk.ljb_coeff[diag_idx]));
    othr_epsilon[i] = 0.25 * nbk.ljb_coeff[diag_idx] / pow(othr_sigma[i], 6.0);
  }

  // Obtain a table of which atoms in the current table have pair-specific combining rules.
  const std::vector<bool> base_has_nbfix_terms =
    findPairSpecificParticipation(lja_coeff.data(), ljb_coeff.data(), lj_type_count, absolute_rule,
                                  prevalent_rule);
  const std::vector<bool> othr_has_nbfix_terms =
    findPairSpecificParticipation(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types, absolute_rule,
                                  prevalent_rule);
  
  // Check for overlap in the atom types.  Make a list of all atom type names that are found in
  // both topologies, look up their indices, and check how they interact in each case.
  std::vector<int> existing_ljt, othr_ljt;
  std::vector<int> ljt_correspondence(nbk.n_lj_types, -1);
  std::vector<bool> unique_types(nbk.n_lj_types, true);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const int jlim = othr_type_aliases.size();
    int lj_type_footprint = -1;
    for (int j = 0; j < jlim; j++) {
      std::map<uint, int>::iterator it = atom_type_map.find(char4ToUint(othr_type_aliases[i][j]));
      if (it != atom_type_map.end()) {
        const int guess_lj_type = it->second;
        unique_types[i] = false;

        // An atom type representing the new topology's ith Lennard-Jones type is represented in
        // the existing tables.  But does it have the same Lennard-Jones properties in the
        // existing tables?
        const int current_table_ljidx = it->second;
        othr_ljt.push_back(i);
        existing_ljt.push_back(current_table_ljidx);

        // Ensure that overlapping atom types are confined to specific Lennard-Jones types.
        if (lj_type_footprint == -1) {
          lj_type_footprint = current_table_ljidx;
        }
        else if (current_table_ljidx != lj_type_footprint) {
          std::string i_alias_list, cfoot_alias_list, ccurr_alias_list;
          for (size_t k = 0; k < othr_type_aliases[i].size(); k++) {
            i_alias_list += char4ToString(othr_type_aliases[i][k]) + " ";
          }
          for (size_t k = 0; k < atom_type_aliases[lj_type_footprint].size(); k++) {
            cfoot_alias_list += char4ToString(atom_type_aliases[lj_type_footprint][k]) + " ";
          }
          for (size_t k = 0; k < atom_type_aliases[current_table_ljidx].size(); k++) {
            ccurr_alias_list += char4ToString(atom_type_aliases[current_table_ljidx][k]) + " ";
          }
          rtErr(char4ToString(othr_type_aliases[i][j]) + " should map to the same consensus "
                "Lennard-Jones type as other atom types in its group [ " + i_alias_list +
                "].  Instead it maps to [ " + cfoot_alias_list + "] as well as [ " +
                ccurr_alias_list + "], and possibly the Lennard-Jones interactions of other "
                "groups in the consensus tables.", "LennardJonesAnalysis", "addTopology");
        }
      }
    }
    if (unique_types[i]) {

      // If the atom type was not found, check that there are no pair-specific combinations of the
      // incoming Lennard-Jones type within its own topology.  Such a finding would enforce type
      // uniqueness.  If not, then check that the Lennard-Jones parameters are indeed unique and
      // that no existing pair-specific parameters might confound the equivalence.
      if (othr_has_nbfix_terms[i] == false) {
        for (int j = 0; j < lj_type_count; j++) {
          if (fabs(sigma[j] - othr_sigma[i]) <= 1.0e-6 &&
              fabs(epsilon[j] - othr_epsilon[i]) <= 1.0e-6 && base_has_nbfix_terms[j]) {
            unique_types[i] = false;
            ljt_correspondence[i] = j;
          }
        }
      }
    }
  }
  
  // Check the Lennard-Jones A and B coefficients for all interaction types that are found in
  // common between the two systems.  If there are discrepancies, this is an error.
  const int n_common = othr_ljt.size();
  std::vector<double> common_lja(n_common * n_common), common_ljb(n_common * n_common);
  for (int i = 0; i < n_common; i++) {
    const int exst_ljt_i = existing_ljt[i];
    const int othr_ljt_i = othr_ljt[i];
    for (int j = 0; j < n_common; j++) {
      const int exst_ljt_j = existing_ljt[j];
      const int othr_ljt_j = othr_ljt[j];
      const double lja_exst = lja_coeff[(exst_ljt_j * lj_type_count) + exst_ljt_i];
      const double ljb_exst = ljb_coeff[(exst_ljt_j * lj_type_count) + exst_ljt_i];
      const double lja_othr = nbk.lja_coeff[(othr_ljt_j * nbk.n_lj_types) + othr_ljt_i];
      const double ljb_othr = nbk.ljb_coeff[(othr_ljt_j * nbk.n_lj_types) + othr_ljt_i];
      if (lja_exst < constants::tiny) {
        if (fabs(lja_exst - lja_othr) >= 1.0e-3) {
          std::string exst_type_list_i("[ "), exst_type_list_j("[ ");
          for (size_t k = 0; k < atom_type_aliases[exst_ljt_i].size(); k++) {
            const char4 exst_atyp = atom_type_aliases[exst_ljt_i][k];
            bool found = false;
            for (size_t m = 0; m < othr_type_aliases[othr_ljt_i].size(); m++) {
              found = (found || exst_atyp == othr_type_aliases[othr_ljt_i][m]);
            }
            if (found == false) {
              exst_type_list_i += char4ToString(exst_atyp) + " ";
            }
          }
          for (size_t k = 0; k < atom_type_aliases[exst_ljt_j].size(); k++) {
            const char4 exst_atyp = atom_type_aliases[exst_ljt_j][k];
            bool found = false;
            for (size_t m = 0; m < othr_type_aliases[othr_ljt_j].size(); m++) {
              found = (found || exst_atyp == othr_type_aliases[othr_ljt_j][m]);
            }
            if (found == false) {
              exst_type_list_j += char4ToString(exst_atyp) + " ";
            }
          }
          exst_type_list_i += "]";
          exst_type_list_j += "]";
          rtErr("Atom types " + exst_type_list_i + " and " + exst_type_list_j + "interact with "
                "different properties in the consensus tables than they do in the new set.",
                "LennardJonesAnalysis", "addTopology");
        }
      }
    }
  }

  // Having survived the type comparisons, the combination procedure should now isolate the
  // new types and create new tables.

  // CHECK
  printf("Lennard Jones Types to atom types:\n");
  for (int i = 0; i < nbk.n_lj_types; i++) {
    printf("  %2d : ", i);
    for (size_t j = 0; j < othr_type_aliases[i].size(); j++) {
      const char4 tmp_atyp = othr_type_aliases[i][j];
      printf(" %c%c%c%c", tmp_atyp.x, tmp_atyp.y, tmp_atyp.z, tmp_atyp.w);
    }
    printf("\n");
  }
  // END CHECK

  // CHECK
  printf("Uniqueness of added types: ");
  for (int i = 0; i < nbk.n_lj_types; i++) {
    printf(" %c", unique_types[i] ? 'T' : 'F');
  }
  printf("\n");
  // END CHECK

}

} // namespace topology
} // namespace stormm
