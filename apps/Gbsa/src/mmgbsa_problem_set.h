// -*-c++-*-
#ifndef MMGBSA_PROBLEM_SET_H
#define MMGBSA_PROBLEM_SET_H

#include <string>
#include <vector>
#include "../../../src/copyright.h"
#include "../../../src/Accelerator/hybrid.h"
#include "../../../src/Accelerator/gpu_details.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Constants/fixed_precision.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Trajectory/phasespace.h"
#include "mmgbsa_enumerators.h"

namespace mmgbsa {

using stormm::card::default_hpc_format;
using stormm::card::GpuDetails;
using stormm::card::HybridFormat;
using stormm::constants::ExceptionResponse;
#ifndef STORMM_USE_HPC
using stormm::data_types::int2;
#endif
using stormm::numerics::default_force_scale_bits;
using stormm::numerics::default_globalpos_scale_bits;
using stormm::numerics::default_velocity_scale_bits;
using stormm::synthesis::AtomGraphSynthesis;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SystemCache;
using stormm::topology::AtomGraph;
using stormm::trajectory::PhaseSpace;

//-------------------------------------------------------------------------------------------------
// A class to manage the various systems and complexes involved in the MM-GBSA calculation.  The
// problem set begins with an instance of the uliganded (apo) receptor, followed by a series of
// ligand poses isolated by themselves, and finishes with a series of complexes of each ligand and
// the receptor, matching the order in which the individual ligands are presented.  While many of
// the ligand poses may result in the same structure after energy minimization, their small size
// means that these problems will result in relatively small costs to the overall calculation.
//-------------------------------------------------------------------------------------------------
class MMGBSAProblemSet {
public:

  // The constructor requires an existing SystemCache object, containing the data from
  // user-supplied files.
  MMGBSAProblemSet(const SystemCache *sysche_pointer_in);

  // With no pointers to repair (and only Standard Template Library components for arrays), the
  // default copy and move constructors as well as compy and move assignment operators all apply.
  //
  // Arguments:
  //   original:  The original object to copy or move
  //   other:     Another object placed on the right hand side of the assignment statement
  MMGBSAProblemSet(const MMGBSAProblemSet &original) = default;
  MMGBSAProblemSet(MMGBSAProblemSet &&original) = default;
  MMGBSAProblemSet& operator=(const MMGBSAProblemSet &original) = default;
  MMGBSAProblemSet& operator=(MMGBSAProblemSet &&original) = default;
  
  // Get the number of problems involved in the calculation.
  int getProblemCount() const;

  // Get the number of states of the receptor.
  int getReceptorStates() const;

  // Get the number of complexes involved in the calculation.
  int getComplexCount() const;

  // Get the cache index which is found to be the receptor.
  int getReceptorTopologyCacheIndex() const;

  // Obtain the ligand topology index in the system cache used in each problem.  The result may be
  // -1 if there is no ligand, i.e. in the first instance for the _apo_ receptor.
  //
  // Overloaded:
  //   - Get the topology index of a specific problem
  //   - Get the topology indices for the ligands in all problems
  //
  // Arguments:
  //   problem_idx:  The index of the problem of interest.  A check will be applied to ensure that
  //                 the index is valid.
  int getLigandTopologyCacheIndex(int problem_idx) const;
  const std::vector<int>& getLigandTopologyCacheIndex() const;

  // Obtain the (initial) ligand coordinate index used to create each problem.  Overloading and
  // descriptions of input variables follow from getLigandTopologyIndex(), above.
  int getLigandCoordinateCacheIndex(int problem_idx) const;
  const std::vector<int>& getLigandCoordinateCacheIndex() const;

  // Get a pointer to the topology used by one of the problems.
  //
  // Arguments:
  //   problem_idx:  The index of the problem of interest.  A check will be applied to ensure that
  //                 the index is valid.
  const AtomGraph* getProblemTopologyPointer(int problem_idx) const;

  // Get the array of all unique topology pointers, an indication of what is involved in any
  // synthesis that emerges from the object.
  const std::vector<AtomGraph*>& getUniqueTopologies() const;
  
  // Get a const reference to one of the initial coordinate sets for a specific problem.
  // Descriptions of input parameters follow from getProblemTopologyPointer(), above.
  const PhaseSpace& getProblemInitialCoordinates(int problem_idx) const;

  // Get the unique topology index for a particular system.
  //
  // Arguments:
  //   problem_idx:  The index of the problem of interest.  A check will be applied to ensure that
  //                 the index is valid.
  int getMasterTopologyIndex(int problem_idx) const;
  
  // Get a list of the topology indices for each system.
  const std::vector<int>& getMasterTopologyIndices() const;

  // Create a topology synthesis for the set of MM-GBSA problems at hand.
  AtomGraphSynthesis
  exportTopologySynthesis(const GpuDetails &gpu = null_gpu,
                          ExceptionResponse policy = ExceptionResponse::DIE) const;

  // Create a coordinate synthesis for the set of MM-GBSA problems at hand.  The synthesis has
  // systems ordered as follows: receptor states, ligand structures for problems 0, 1, 2, ...,
  // complex structures for problems 0, 1, 2, ...
  PhaseSpaceSynthesis
  exportCoordinateSynthesis(const GpuDetails &gpu = null_gpu,
                            HybridFormat fmt = default_hpc_format,
                            int gpos_scale_bits = default_globalpos_scale_bits,
                            int vel_scale_bits = default_velocity_scale_bits,
                            int frc_scale_bits = default_force_scale_bits) const;

  // Get the system cache pointer referenced for all topologies and coordinates used to create the
  // MM-GBSA problem set.
  const SystemCache* getCachePointer() const;
  
private:
  int problem_count;                  // The total number of molecular systems (receptor, free
                                      //   ligands, or complexes) involved in the problem set
  int receptor_states;                // The number of independent states of the receptor.  Its
                                      //   total energy is averaged over all such states.
  int receptor_topology_cache_index;  // Index of the receptor's topology in the list maintained
                                      //   by the system cache.  This topology will be combined
                                      //   with various ligand topologies as demanded by each
                                      //   problem.
  int unique_ligand_count;            // The number of unique ligands (many poses of the same
                                      //   ligand count as one for the purposes of this tally)

  // The ligand topology index in the cache's condensed list for each problem.  This array will
  // read -1 if there is no ligand, in the first element where the _apo_ receptor is handled.
  std::vector<int> ligand_topology_index;

  // The ligand coordinate index in the cache's list for each problem.  As above, this array will
  // read -1 if there is no ligand.
  std::vector<int> ligand_coordinate_index;

  // Array to indicate whether the receptor is present in system (e.g. it is the _apo_ receptor
  // in the first problem, or it is a complex of any sort)
  std::vector<bool> receptor_is_present;

  // Array of unique receptor / ligand combinations, each a topology that is not present in the
  // system cache and must be created in the program.  The "x" member of the tuple refers to the
  // cache index of the receptor topology and will be the value set forth in
  // receptor_topology_cache_index.  The "y" member of the tuple holds the cache index of the
  // ligand topology.
  std::vector<int2> complex_topology_index;

  // Array of unique complex topologies, each entry in the list corresponding to the tuple in
  // complex_topology_index, above.
  std::vector<AtomGraph> custom_topologies;

  // The array of all receptor, ligand pose, and complex systems, ready for translation into a
  // synthesis
  std::vector<PhaseSpace> all_systems;

  // The master array of topology pointers for each system.  These can map into the original
  // system cache, or into the array of custom, complex topologies maintained by this object.
  // There is one entry in this array for each unique topology, so that it can be used in making
  // a topology synthesis without triggering warnings about unused topologies.  Use the array
  // all_master_topology_indices below to find a specific problem's topology by an index into this
  // array.
  std::vector<AtomGraph*> master_topology_list;

  // The master array of topological indices, referencing elements of master_topology_array,
  // needed to describe each entry in the member variable all_systems above.
  std::vector<int> all_master_topology_indices;

  // Pointer to the cache of user-supplied systems referenced in laying out the object
  SystemCache *sysche_pointer;

  // Validate a request for a problem index within the object.
  //
  // Arguments:
  //   problem_idx:  Index of the problem of interest
  void validateProblemIndex(int problem_idx) const;
};

} // namespace mmgbsa

#endif
