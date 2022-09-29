#include <string>
#include "copyright.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "topology_util.h"

namespace stormm {
namespace topology {

using math::findBin;
using math::prefixSumInPlace;
using math::PrefixSumType;
using parse::char4ToString;
  
//-------------------------------------------------------------------------------------------------
std::string writeAtomList(const std::vector<int> &atom_list, const ChemicalDetailsKit &cdk) {
  std::string result;
  const size_t natom = atom_list.size(); 
  for (size_t i = 0; i < natom; i++) {
    const int atom_idx = atom_list[i];
    const int res_idx = findBin(cdk.res_limits, atom_idx, cdk.nres);
    result += std::to_string(cdk.atom_numbers[atom_idx]) + " " +
              char4ToString(cdk.atom_names[atom_idx]) + " " +
              char4ToString(cdk.res_names[res_idx]) + " " +
              std::to_string(cdk.res_numbers[atom_idx]);
    if (i < natom - 1LLU) {
      result += ", ";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const NonbondedKit<double> &nbk, const int atom_i, const int atom_j) {
  const int hlim = nbk.nb12_bounds[atom_i + 1];
  bool result = false;
  for (int i = nbk.nb12_bounds[atom_i]; i < hlim; i++) {
    result = (result || (atom_j == nbk.nb12x[i]));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const AtomGraph &ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  return isBonded(nbk, atom_i, atom_j);
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const AtomGraph *ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  return isBonded(nbk, atom_i, atom_j);
}

} // namespace topology
} // namespace stormm
