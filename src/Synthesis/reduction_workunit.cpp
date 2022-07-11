#include "Constants/hpc_bounds.h"
#include "Reporting/error_format.h"
#include "reduction_workunit.h"

namespace omni {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
ReductionWorkUnit::ReductionWorkUnit(const int atom_start_in, const int atom_end_in,
                                     const int result_index_in, const int dependency_start_in,
                                     const int dependency_end_in) :
  atom_start{atom_start_in}, atom_end{atom_end_in}, result_index{result_index_in},
  dependency_start{dependency_start_in}, dependency_end{dependency_end_in}
{}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getAtomStart() const {
  return atom_start;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getAtomEnd() const {
  return atom_end;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getResultIndex() const {
  return result_index;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getDependencyStart() const {
  return dependency_start;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getDependencyEnd() const {
  return dependency_end;
}

//-------------------------------------------------------------------------------------------------
std::vector<ReductionWorkUnits> buildReductionWorkUnits(const std::vector<int> &atom_starts,
                                                        const std::vector<int> &atom_counts,
                                                        const GpuDetails &gpu,
                                                        GpuLaunch *launcher) {
  if (atom_starts.size() != atom_counts.size()) {
    rtErr("Starting indices were provided for the atoms of " + std::to_string(atom_starts.size()) +
          " systems, but atom counts for " + std::to_string(atom_counts.size()) + " systems were "
          "provided.", "buildReductionWorkUnits");
  }
  const int nsys = atom_starts.size();
  
  // Compare the number of systems to the number of streaming multiprocessors on the GPU, when
  // blocks of various sizes are used.
  for (int block_size = tiny_block_size; block_size < large_block_size; block_size *= 2) {

  }

  // Commit this wisdom to the kernel launch guide
  launcher->setReductionGridDims(best_block_size, best_grid_size);
}

}
}
