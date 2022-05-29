#include "nonbonded_workunit.h"

namespace omni {
namespace synthesis {

using energy::supertile_length;
using energy::tile_length;
using energy::tile_lengths_per_supertile;

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask *se_in,
                                     const std::vector<int2> &tile_list) {
    atom_limit{se_in->atom_count},
    imports{std::vector<int>(16)},
    tile_count{static_cast<int>(tile_list.size())},
    tile_instructions{std::vector<uint2>(tile_count, 0)},
    score{0},
{
  const int n_supert = se_in->supertile_stride_count;
  std::vector<int> tmp_imports(2 * tile_count);
  for (int i = 0; i < tile_count; i++) {
    if (tile_list[i].x < 0 || tile_list[i].y < 0 || tile_list[i].x >= se_in->atom_count ||
        tile_list[i].y >= se_in->atom_count) {
      rtErr("A tile with lower left corner (" + std::to_string(tile_list[i].x) + ", " +
            std::to_string(tile_list[i].y) + ") cannot be part of a system with " +
            std::to_string(se_in->atom_count) + " atoms.", "NonbondedWorkUnit");
    }
    const int sti = tile_list[i].x / supertile_length;
    const int stj = tile_list[i].y / supertile_length;
    const int stij_map_index = se_in->supertile_map_idx[(stj * n_supert) + sti];
    const int ti = (tile_list[i].x - (sti * supertile_length)) / tile_length;
    const int tj = (tile_list[i].y - (stj * supertile_length)) / tile_length;
    const int tij_map_index = se_in->tile_map_idx[stij_map_index +
                                                  (tj * tile_lengths_per_supertile) + ti];
    tile_instructions[j].y = tij_map_index;

    // Divide the starting atom indices by the tile length to help exploit the compressed sorting
    // feature in reduceUniqueValues().  The numbers will be re-inflated later.
    tmp_imports[(2 * i)    ] = tile_list[i].x / tile_length;
    tmp_imports[(2 * i) + 1] = tile_list[i].y / tile_length;
  }
  reduceUniqueValues(tmp_imports);
  const int nbatch = tmp_imports.size();
  for (int i = 0; i < nbatch; i++) {
    tmp_imports[i] *= tile_length;
    imports[i] = tmp_imports[i];
  }
  for (int i = nbatch; i < 16; i++) {
    imports[i] = -1;
  }

  // Compute the workunit effort score, an estimate of how long one work unit will take relative
  // to others.  Work units with the highest effort will be placed first, to backfill idle
  // processes with shorter work units.
  score = nbatch + (tile_count * 8);
  
  // The limited number of imports, and the manipulations that need to happen to assemble the
  // work unit instructions, likely makes it preferable to just search them all rather than
  // calling a binary search function.
  for (int i = 0; i < nbatch; i++) {
    const int import_idx = imports[i];
    const uint absc_mask = 16 * i;
    const uint ordi_mask = (absc_mask << 16);
    for (int j = 0; j < tile_count; j++) {
      if (tile_list[j].x == import_idx) {
        tile_instructions[j].x |= absc_mask;
      }
      else if (tile_list[j].y == import_idx) {
        tile_instructions[j].x |= ordi_mask;        
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask *se_in,
                                     const int abscissa_start, const int ordinate_start) :
    atom_limit{se_in->atom_count},
    imports{},
    tile_count{0},
    tile_instructions{},
    score{0}
{
  imports.resize(2);
  imports[0] = abscissa_start;
  imports[1] = ordinate_start;
  
  // Compute the tile count and the work unit score
  if (abscissa_start == ordinate_start) {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * (tile_lengths_per_supertile + 1) / 2;
      score = tile_lengths_per_supertile;
    }
    else {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;  
      tile_count = ntile_rem * (ntile_rem + 1) / 2;
      score = ntile_rem;
    }
  }
  else {
    if (abscissa_start < atom_limit - supertile_length &&
        ordinate_start < atom_limit - supertile_length) {
      tile_count = tile_lengths_per_supertile * tile_lengths_per_supertile;
      score = tile_lengths_per_supertile * 2;
    }
    else if (abscissa_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - ordinate_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      score = tile_lengths_per_supertile + ntile_rem;
    }
    else if (ordinate_start < atom_limit - supertile_length) {
      const int ntile_rem = (atom_limit - abscissa_start + tile_length - 1) / tile_length;
      tile_count = tile_lengths_per_supertile * ntile_rem;
      score = tile_lengths_per_supertile + ntile_rem;
    }
    else {
      rtErr("A trapezoidal tile should not exist in a supertile-based work unit.",
            "NonbondedWorkUnit");
    }
  }
  score += time_count * 8;
}
 
//-------------------------------------------------------------------------------------------------
size_t estimateNonbondedWorkUnitCount(const std::vector<int> &atom_counts, const int nbwu_tile) {
  const size_t nsys = atom_counts.size();
  size_t result = 0LLU;
  const int big_nbwu = tile_lengths_per_supertile * tile_lengths_per_supertile;
  for (size_t i = 0; i < atom_counts; i++) {
    const size_t ntile_side = (nbwu_tile == big_nbwu) ? atom_counts[i] / supertile_length :
                                                        atom_counts[i] / tile_length;
    result += ntile_side * (ntile_side + 1LLU) / 2LLU;
  }
  if (nbwu_tile != 8 && nbwu_tile != 16 && nbwu_tile != 32 && nbwu_tile != 64) {
    const size_t nbwu_tile_zu = static_cast<size_t>(nbwu_tile);
    return (result + nbwu_tile_zu - 1LLU) / nbwu_tile_zu;
  }
  else if (nbwu_tile == big_nbwu) {
    return result;
  }
  else {
    rtErr("Only 8, 16, 32, 64, or " + std::to_string(big_nbwu) + " tiles may be the target size "
          "of any one NonbondedWorkUnit.", "estimateNonbondedWorkUnitCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
enumerateNonbondedWorkUnits(const std::vector<AtomGraph*> &ag_list,
                            const std::vector<StaticExclusionMask*> se_list,
                            const std::vector<int> &topology_indices,
                            const int target_tile_count) {
  std::vector<NonbondedWorkUnit> result;
  result.reserve(nbwu_count);
  const int nsys = topology_indices.size();
  const int2 zero2 = { 0, 0 };
  const int full_tile_absc_count = (target_tile_count == 8) ? 2 :
    (target_tile_count == 
  
  for (int sysid = 0; sysid < nsys; sysid++) {
    const int natom = ag_list[sysid]->getAtomCount();
    
    // Take care of extremely small systems
    if (natom <= tile_length) {
      result.emplace_back(se_list[sysid], std::vector<int2>(1, zero2));
      continue;
    }
    const int total_tile_lengths = (natom + tile_length - 1) / tile_length;
    const int last_tile_start = tile_length * (total_tile_lengths - 1);
    std::vector<bool> coverage(total_tile_lengths * total_tile_lengths, false);
    for (int i = 0; i < total_tile_lengths; i++) {
      for (int j = i + 1; j < total_tile_lengths; j++) {
        coverage[(j * total_tile_lengths) + i] = true;
      }
    }

    // The strategy is to list tiles in the typical order, keeping track of import requirements.
    // Get as close as possible to the target work unit size, keeping to multiples of eight if at
    // all possible, then commit work units to the list as they reach the target size or it
    // becomes impossible to grow the work unit without exceeding the import limits.
    const int max_imports = small_block_max_atoms / tile_length;
    int absc_import_min = 0;
    int absc_import_max = 0;
    int ordi_import_min = 0;
    int ordi_import_max = 0;
    std::vector<int2> tile_list(target_tile_count);
    int current_tile_count = 0;
    int i = 0;
    int j = 0;
    const int ordi_stride = (target_tile_count ==  8) ? 2 :
                            (target_tile_count == 16) ? 4 :
                            (target_tile_count == 32) ? 6 :
                            (target_tile_count == 64) ? 8;
    int ordi_llim = 0;
    int ordi_hlim = std::min(ordi_stride, total_tile_lengths);
    int absc_llim = 0;
    int absc_hlim = 1;
    int i = 0;
    int j = 0;
    while (ordi_lim <= total_tile_lengths) {
      for (i = 0; i < 
    }

    for (int i = 0; i < total_tile_lengths; i++) {
      for (j = 0; j <= i; j++) {
        
      }
      
      tile_list[current_tile_count].x = tile_length * i;
      tile_list[current_tile_count].y = tile_length * i;
      current_tile_count++;
      if (current_tile_count == target_tile_count) {
        current_tile_count = 0;
        result.emplace_back(se_list[sysid], tile_list);
      }
    }
    for (int i = total_tile_lengths - 1; i >= 0; i--) {
      tile_list[current_tile_count].x = last_tile_start;
      if (current_tile_count == target_tile_count) {
        current_tile_count = 0;
        result.emplace_back(se_list[sysid], tile_list);
      }
    }
    
    
  }
}

//-------------------------------------------------------------------------------------------------

     
//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const std::vector<AtomGraph*> &ag_list,
                        const std::vector<StaticExclusionMask*> se_list,
                        const std::vector<int> &topology_indices, const int target_nbwu_count) {

  // Determine the optimal overall size for work units.  Given that this process is guided by
  // static (as opposed to forward) exclusion masks, this is a matter of how many atoms are
  // present in all topologies and the number of tiles it would take to cover them all.
  const size_t nsys = topology_indices.size();
  size_t pair_acc = 0LLU;
  std::vector<int> atom_counts(nsys);
  for (size_t i = 0; i < nsys; i++) {
    atom_counts[i] = ag_list[topology_indices[i]]->getAtomCount();
  }
  for (
    const size_t natom_zu = atom_counts[i];
    pair_acc += natom_zu * natom_zu;
  }

  // Try work units of eight tiles
  const int huge_nbwu_size = tile_lengths_per_supertile * tile_lengths_per_supertile;
  const size_t tiny_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, 8);
  const size_t small_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, 16);
  const size_t medium_wu_count = estimateNonbondedWorkUnitCount(atom_counts, 32);
  const size_t large_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, 64);
  const size_t huge_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, huge_nbwu_size);
  if (tiny_wu_count < mega) {
    
  }
  else if (small_wu_count < mega) {
    
  }
  else if (medium_wu_count < mega) {

  }
  else if (large_wu_count < mega) {

  }
  else {

    // Use the huge tiles (this will imply the large block size)
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const AtomGraph* ag, const StaticExclusionMask* se,
                        const int target_nbwu_count) {
  const std::vector<AtomGraph*> agv(1, const_cast<AtomGraph*>(ag));
  const std::vector<StaticExclusionMask*> sev(1, const_cast<StaticExclusionMask*>(se));
  return buildNonbondedWorkUnits(agv, sev, std::vector<int>(1, 0), target_nbwu_count);
}

} // namespace synthesis
} // namespace omni
