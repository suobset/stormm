#include <algorithm>
#include "../../../src/copyright.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/Math/vector_ops.h"
#include "mmgbsa_problem_set.h"

namespace mmgbsa {

using stormm::numerics::default_localpos_scale_bits;
using stormm::stmath::incrementingSeries;
using stormm::stmath::reduceUniqueValues;
using stormm::errors::rtErr;
using stormm::errors::rtWarn;
  
//-------------------------------------------------------------------------------------------------
MMGBSAProblemSet::MMGBSAProblemSet(const SystemCache *sysche_pointer_in) :
    problem_count{0}, receptor_states{0}, receptor_topology_cache_index{0}, unique_ligand_count{0},
    ligand_topology_index{}, ligand_coordinate_index{}, receptor_is_present{},
    complex_topology_index{}, custom_topologies{}, all_systems{}, master_topology_list{},
    all_master_topology_indices{}, sysche_pointer{const_cast<SystemCache*>(sysche_pointer_in)}
{
  // Look for the ligands and the receptor
  const std::vector<int> rec_cache_idx = sysche_pointer->getMatchingSystemIndices("receptor");
  receptor_topology_cache_index = sysche_pointer->getSystemTopologyIndex(rec_cache_idx[0]);
  const std::vector<int> lig_cache_idx = sysche_pointer->getMatchingSystemIndices("ligand");

  // The number of receptors here referes ot the total number of receptor snapshots, even though
  // there can only be one definition of the receptor (one receptor topology).
  receptor_states = rec_cache_idx.size();
  
  // The number of ligands here refers to the total number of ligand poses, and is greater than or
  // at least equal to the number of unique ligand topologies.
  const int n_ligand = lig_cache_idx.size();
  problem_count = (2 * n_ligand) + receptor_states;
  ligand_topology_index.resize(problem_count);
  ligand_coordinate_index.resize(problem_count);
  receptor_is_present.resize(problem_count);
  for (int i = 0; i < receptor_states; i++) {
    ligand_topology_index[i] = -1;
    ligand_coordinate_index[i] = -1;
  }
  receptor_is_present[0] = true;
  for (int i = 0; i < n_ligand; i++) {
    const int lig_idx = i + receptor_states;
    const int cmp_idx = lig_idx + n_ligand;
    ligand_topology_index[lig_idx] = sysche_pointer->getSystemTopologyIndex(lig_cache_idx[i]);
    ligand_topology_index[cmp_idx] = ligand_topology_index[lig_idx];
    ligand_coordinate_index[lig_idx] = lig_cache_idx[i];
    ligand_coordinate_index[cmp_idx] = lig_cache_idx[i];
    receptor_is_present[lig_idx] = false;
    receptor_is_present[cmp_idx] = true;
  }

  // Search the ligand topology indices for unique entries, and when they are found in conjunction
  // with receptor presence, note that a new topology, the fusion of both ligand and receptor, must
  // be created.
  std::vector<int> unique_ligand_topidx = reduceUniqueValues(ligand_topology_index);
  unique_ligand_topidx.erase(std::find(unique_ligand_topidx.begin(),
                                       unique_ligand_topidx.end(), -1));
  unique_ligand_count = unique_ligand_topidx.size();
  std::vector<int> unique_ligtop_map(problem_count);
  for (int i = 0; i < receptor_states; i++) {
    unique_ligtop_map[i] = -1;
  }
  for (int i = receptor_states; i < n_ligand + receptor_states; i++) {
    if (i > receptor_states &&
        ligand_topology_index[i] == ligand_topology_index[i - receptor_states]) {
      unique_ligtop_map[i] = unique_ligtop_map[i - receptor_states];
      unique_ligtop_map[i + n_ligand] = unique_ligtop_map[i + n_ligand - receptor_states];
    }
    else {
      for (int j = 0; j < unique_ligand_count; j++) {
        if (ligand_topology_index[i] == unique_ligand_topidx[j]) {
          unique_ligtop_map[i] = j;
          unique_ligtop_map[i + n_ligand] = j;
        }
      }
    }
  }
  
  // Create fused topologies of each unique ligand topology and the receptor.
  complex_topology_index.reserve(unique_ligand_count);
  custom_topologies.reserve(unique_ligand_count);
  const AtomGraph& rec_top_ref = sysche_pointer->getTopology(receptor_topology_cache_index);
  for (int i = 0; i < unique_ligand_count; i++) {
    complex_topology_index[i] = { receptor_topology_cache_index, unique_ligand_topidx[i] };
    const AtomGraph& lig_top_ref = sysche_pointer->getTopology(unique_ligand_topidx[i]);
    custom_topologies.emplace_back(rec_top_ref, lig_top_ref);
  }

  // Create the complete list of initial coordinates
  all_systems.reserve(problem_count);
  for (int i = 0; i < receptor_states; i++) {
    all_systems.push_back(sysche_pointer->getCoordinates(rec_cache_idx[i]));
  }
  for (int i = 0; i < n_ligand; i++) {
    all_systems.push_back(sysche_pointer->getCoordinates(lig_cache_idx[i]));
  }
  for (int i = 0; i < n_ligand; i++) {
    all_systems.emplace_back(all_systems[0], rec_top_ref,
                             sysche_pointer->getCoordinates(lig_cache_idx[i]),
                             sysche_pointer->getSystemTopology(lig_cache_idx[i]));
  }
  
  // Create the complete list of topology pointers for each system.
  master_topology_list.reserve((2 * unique_ligand_count) + 1);
  master_topology_list.push_back(const_cast<AtomGraph*>(rec_top_ref.getSelfPointer()));
  for (int i = 0; i < unique_ligand_count; i++) {
    const AtomGraph* ilig_ag = sysche_pointer->getTopologyPointer(unique_ligand_topidx[i]);
    master_topology_list.push_back(const_cast<AtomGraph*>(ilig_ag));
  }
  for (int i = 0; i < unique_ligand_count; i++) {
    master_topology_list.push_back(&custom_topologies[i]);
  }
  all_master_topology_indices.resize(problem_count);
  for (int i = 0; i < receptor_states; i++) {
    all_master_topology_indices[i] = 0;
  }
  for (int i = 0; i < n_ligand; i++) {
    all_master_topology_indices[i + receptor_states] = unique_ligtop_map[i + receptor_states] + 1;
    all_master_topology_indices[i + receptor_states + n_ligand] =
      unique_ligtop_map[i + receptor_states + n_ligand] + 1 + unique_ligand_count;
  }
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getProblemCount() const {
  return problem_count;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getReceptorStates() const {
  return receptor_states;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getComplexCount() const {
  return (problem_count - receptor_states) / 2;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getReceptorTopologyCacheIndex() const {
  return receptor_topology_cache_index;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getLigandTopologyCacheIndex(const int problem_idx) const {
  validateProblemIndex(problem_idx);
  return ligand_topology_index[problem_idx];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MMGBSAProblemSet::getLigandTopologyCacheIndex() const {
  return ligand_topology_index;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getLigandCoordinateCacheIndex(const int problem_idx) const {
  validateProblemIndex(problem_idx);
  return ligand_coordinate_index[problem_idx];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MMGBSAProblemSet::getLigandCoordinateCacheIndex() const {
  return ligand_coordinate_index;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* MMGBSAProblemSet::getProblemTopologyPointer(const int problem_idx) const {
  validateProblemIndex(problem_idx);
  return master_topology_list[all_master_topology_indices[problem_idx]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*>& MMGBSAProblemSet::getUniqueTopologies() const {
  return master_topology_list;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace& MMGBSAProblemSet::getProblemInitialCoordinates(const int problem_idx) const {
  validateProblemIndex(problem_idx);
  return all_systems[problem_idx];
}

//-------------------------------------------------------------------------------------------------
int MMGBSAProblemSet::getMasterTopologyIndex(const int problem_idx) const {
  validateProblemIndex(problem_idx);
  return all_master_topology_indices[problem_idx];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MMGBSAProblemSet::getMasterTopologyIndices() const {
  return all_master_topology_indices;
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis
MMGBSAProblemSet::exportTopologySynthesis(const GpuDetails &gpu,
                                          const ExceptionResponse policy) const {
  return AtomGraphSynthesis(master_topology_list, all_master_topology_indices, policy, gpu);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
MMGBSAProblemSet::exportCoordinateSynthesis(const GpuDetails &gpu, const HybridFormat fmt,
                                            const int gpos_scale_bits, const int vel_scale_bits,
                                            const int frc_scale_bits) const {
  return PhaseSpaceSynthesis(all_systems, incrementingSeries<int>(0, problem_count),
                             master_topology_list, all_master_topology_indices,
                             gpos_scale_bits, default_localpos_scale_bits, vel_scale_bits,
                             frc_scale_bits, fmt, gpu);
}

//-------------------------------------------------------------------------------------------------
const SystemCache* MMGBSAProblemSet::getCachePointer() const {
  return sysche_pointer;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAProblemSet::validateProblemIndex(const int problem_idx) const {
  if (problem_idx < 0 || problem_idx >= problem_count) {
    rtErr("Problem index " + std::to_string(problem_idx) + " is invalid for a collection of " +
          std::to_string(problem_count) + " problems (" + std::to_string((problem_count - 1) / 2) +
          " poses of " + std::to_string(unique_ligand_count) + " unique ligands against a single "
          "receptor).", "MMGBSAProblemSet", "validateProblemIndex");
  }
}

} // namespace mmgbsa
