#ifndef OMNI_NONBONDED_WORKUNIT_H
#define OMNI_NONBONDED_WORKUNIT_H

#include "Topology/atomgraph.h"
#include "Potential/static_exclusionmask.h"

namespace omni {
namespace synthesis {

/// \brief Management for the non-bonded delegation of one or more topologies.
class NonbondedDelegator {
public:

private:
};
  
/// \brief Collect a series of tiles for non-bonded computations as well as the required atom
///        imports to carry them out.  This will accomplish the task of planning the non-bonded
///        computation, given a single topology or a synthesis of topologies, to optimize the
///        thread utilization on a GPU.
class NonbondedWorkUnit {
public:

  /// \brief The constructor assumes a 16x16 tile and accepts a target number of tiles for each
  ///        work unit to perform.  There will be a few "magic" numbers: eight tiles, the minimum
  ///        for occupying 256 threads on NVIDIA or AMD GPUs, sixteen tiles, the most that can be
  ///        accomplished with 64 x 64 imported atoms, sixty-four, for a much larger tile that
  ///        approaches the maximum __shared__ memory usage for imported atoms and Generalized Born
  ///        force computations.
  ///
  /// \param 
  NonbondedWorkUnit(int tile_count_in);

private:
  int abscissa_atom_count;
  int ordinate_atom_count;
  int tile_count;
  std::vector<uint2> tile_instructions;
  const AtomGraph *ag_pointer;
  const StaticExclusionMask *se_pointer;
};

/// \brief Create a list of work units of a consistent size so as to optimize the load balancing
///        on a given resource (GPU, on or more CPUs, etc.).
///
/// Overloaded:
///   - Take a single topology or a list of topologies
///   - Take one or more static exclusion masks corresponding to the input topologies
///   - Accept a target number of work units
///   - Accept a GPU specifications from which to infer a target number of work units
///
/// \param ag_in                The input topology
/// \param ag_list_in           List of topologies
/// \param se_in                Static exclusion mask corresponding to the input topology
/// \param se_list_in           List of static exclusion masks corresponding to each input
///                             topology
/// \param topology_indices_in  Indicies of the topologies referenced by each system
/// \param target_nbwu_count    The target number of work units to produce.  The goal will be to
///                             get as close to this number without going over, or if the number
///                             must be exceeded due to absolute limits on the size of any one
///                             work unit, to produce a number of work units that comes as close
///                             to a multiple of this number as possible, without going over.
/// \param gpu                  GPU specifications (HPC compilation only)
/// \{
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const <AtomGraph*> &ag_in, const StaticExclusionMask *se_in,
                        const std::vector<int> &topology_indices_in = { 0 },
                        int target_nbwu_count = static_cast<int>(mega));

std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const AtomGraph *ag_in, const StaticExclusionMask *se_in,
                        const std::vector<int> &topology_indices_in = { 0 },
                        int target_nbwu_count = static_cast<int>(mega));
/// \}

} // namespace synthesis
} // namespace omni

#endif
