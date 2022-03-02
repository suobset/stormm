#include "DataTypes/omni_vector_types.h"
#include "set_ops.h"
#include "vector_ops.h"

namespace omni {
namespace math {

//-------------------------------------------------------------------------------------------------
std::vector<int> getSubsetIndexPattern(const std::vector<int> &x_subset, int *n_subset_indices) {

  // Obtain the maximum index
  const size_t nidx = x_subset.size();
  const int max_index = maxValue(x_subset);
  std::vector<bool> index_found(max_index + 1, false);
  for (size_t i = 0; i < nidx; i++) {
    index_found[x_subset[i]] = true;
  }
  int new_index_counter = 0;
  std::vector<int> result(max_index + 1, -1);
  for (size_t i = 0; i <= max_index; i++) {
    if (index_found[i]) {
      result[i] = new_index_counter;     
      new_index_counter++;
    }
  }
  if (n_subset_indices != nullptr) {
    *n_subset_indices = new_index_counter;
  }
  return result;
}

} // namespace math
} // namespace omni
