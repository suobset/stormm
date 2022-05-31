#include "Constants/scaling.h"
#include "Math/vector_ops.h"
#include "nonbonded_workunit.h"

namespace omni {
namespace synthesis {

using constants::kilo;
using constants::mega;
using energy::StaticExclusionMaskReader;
using energy::supertile_length;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
using math::accumulateBitmask;
using math::reduceUniqueValues;
  
//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask &se,
                                     const std::vector<int3> &tile_list) :
    tile_count{static_cast<int>(tile_list.size())},
    score{0},
    imports{std::vector<int>(small_block_max_imports)},
    import_size_keys{std::vector<int>(small_block_max_imports / 4, 0)},
    tile_instructions{std::vector<uint2>(tile_count)}
{
  // In this form of the constructor, one system is present with a single StaticExclusionMask.
  // Get the total number of atoms and inflate the tile starting positions to obtain the atom
  // import limits.
  const StaticExclusionMaskReader ser = se.data();
  const int n_supert = se.getSuperTileStrideCount();
  std::vector<int> tmp_imports(2 * tile_count);
  for (int i = 0; i < tile_count; i++) {
    if (tile_list[i].x < 0 || tile_list[i].y < 0 || tile_list[i].x * tile_length >= ser.natom ||
        tile_list[i].y * tile_length >= ser.natom) {
      rtErr("A tile with lower left corner (" + std::to_string(tile_list[i].x * tile_length) +
            ", " + std::to_string(tile_list[i].y * tile_length) + ") cannot be part of a system "
            "with " + std::to_string(ser.natom) + " atoms.", "NonbondedWorkUnit");
    }
    const int sti = tile_list[i].x / tile_lengths_per_supertile;
    const int stj = tile_list[i].y / tile_lengths_per_supertile;
    const int stij_map_index = ser.supertile_map_idx[(stj * n_supert) + sti];
    const int ti = tile_list[i].x - (sti * tile_lengths_per_supertile);
    const int tj = tile_list[i].y - (stj * tile_lengths_per_supertile);
    const int tij_map_index = ser.tile_map_idx[stij_map_index +
                                               (tj * tile_lengths_per_supertile) + ti];
    tile_instructions[i].y = tij_map_index;
    tmp_imports[(2 * i)    ] = tile_list[i].x;
    tmp_imports[(2 * i) + 1] = tile_list[i].y;
  }
  reduceUniqueValues(&tmp_imports);
  const int n_imports = tmp_imports.size();
  for (int i = 0; i < n_imports; i++) {
    tmp_imports[i] *= tile_length;
    imports[i] = tmp_imports[i];
  }
  for (int i = n_imports; i < small_block_max_imports; i++) {
    imports[i] = -1;
  }
  
  // Compute the workunit effort score, an estimate of how long one work unit will take relative
  // to others.  Work units with the highest effort will be placed first, to backfill idle
  // processes with shorter work units.
  score = n_imports + (tile_count * 8);
  
  // The limited number of imports, and the manipulations that need to happen to assemble the
  // work unit instructions, likely makes it preferable to just search them all rather than
  // calling a binary search function.  Each instruction's y member has already been filled out.
  for (int i = 0; i < n_imports; i++) {
    const int import_idx = imports[i] / tile_length;
    const uint absc_mask = 16 * i;
    const uint ordi_mask = (absc_mask << 16);
    for (int j = 0; j < tile_count; j++) {
      if (tile_list[j].x == import_idx) {
        tile_instructions[j].x |= absc_mask;
      }
      if (tile_list[j].y == import_idx) {
        tile_instructions[j].x |= ordi_mask;
      }
    }
  }

  // Fill out the size keys for each import, a series of bit-packed integers with one import
  // batch size (up to tile_length atoms, i.e. 16) per eight bits.
  if (constants::int_bit_count_int < 32) {
    rtErr("Descriptors for non-bonded work units require that signed integers be of size >= 32 "
          "bits.", "NonbondedWorkUnit");
  }
  for (int i = 0; i < n_imports; i++) {
    const int key_idx = i / 4;
    const int key_pos = i - (key_idx * 4);
    int import_batch_size;
    if (imports[i] >= 0) {
      import_batch_size = std::min(ser.natom - imports[i], tile_length);
    }
    else {
      import_batch_size = 0;
    }
    import_size_keys[key_idx] |= (import_batch_size << (8 * key_pos));
  }
}

//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const StaticExclusionMask &se,
                                     const int abscissa_start, const int ordinate_start) :
    tile_count{0},
    score{0},
    imports{},
    import_size_keys{},
    tile_instructions{}
{
  imports.resize(2);
  imports[0] = abscissa_start;
  imports[1] = ordinate_start;
  
  // Compute the tile count and the work unit score
  const int atom_limit = se.getAtomCount();
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
  score += tile_count * 8;
}

//-------------------------------------------------------------------------------------------------
int NonbondedWorkUnit::getTileCount() const {
  return tile_count;
}

//-------------------------------------------------------------------------------------------------
int4 NonbondedWorkUnit::getTileLimits(const int index) const {
  const int absc_idx = (tile_instructions[index].x & 0xffff) / tile_length;
  const int ordi_idx = ((tile_instructions[index].x >> 16) & 0xffff) / tile_length;
  const int absc_key = absc_idx / 4;
  const int absc_pos = absc_idx - (absc_key * 4);
  const int ordi_key = ordi_idx / 4;
  const int ordi_pos = ordi_idx - (ordi_key * 4);
  const int absc_length = ((import_size_keys[absc_key] >> (8 * absc_pos)) & 0xff);
  const int ordi_length = ((import_size_keys[ordi_key] >> (8 * ordi_pos)) & 0xff);
  return { imports[absc_idx], imports[absc_idx] + absc_length,
           imports[ordi_idx], imports[ordi_idx] + ordi_length };
}

//-------------------------------------------------------------------------------------------------
size_t estimateNonbondedWorkUnitCount(const std::vector<int> &atom_counts, const int nbwu_tile) {
  const size_t nsys = atom_counts.size();
  size_t result = 0LLU;
  for (size_t i = 0; i < nsys; i++) {
    const size_t ntile_side = (nbwu_tile == huge_nbwu_tiles) ? atom_counts[i] / supertile_length :
                                                               atom_counts[i] / tile_length;
    result += ntile_side * (ntile_side + 1LLU) / 2LLU;
  }
  if (nbwu_tile == tiny_nbwu_tiles || nbwu_tile == small_nbwu_tiles ||
      nbwu_tile == medium_nbwu_tiles || nbwu_tile == large_nbwu_tiles) {
    const size_t nbwu_tile_zu = static_cast<size_t>(nbwu_tile);
    return (result + nbwu_tile_zu - 1LLU) / nbwu_tile_zu;
  }
  else if (nbwu_tile == huge_nbwu_tiles) {
    return result;
  }
  else {
    rtErr("Only 8, 16, 32, 64, or " + std::to_string(huge_nbwu_tiles) + " tiles may be the target "
          "size of any one NonbondedWorkUnit.", "estimateNonbondedWorkUnitCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool addTileToWorkUnitList(int3* tile_list, int* import_coverage, int *import_count,
                           int *current_tile_count, const int ti, const int tj, const int sysid) {
  const int c_ic  = *import_count;
  if (c_ic + 2 - import_coverage[ti] - import_coverage[tj] > small_block_max_imports) {
    return false;
  }
  else {
    const size_t c_tile = *current_tile_count;
    tile_list[c_tile].x = ti;
    tile_list[c_tile].y = tj;
    tile_list[c_tile].z = sysid;
    *import_count = c_ic + 2 - import_coverage[ti] - import_coverage[tj];
    import_coverage[ti] = 1;
    import_coverage[tj] = 1;
    *current_tile_count += 1;
    return true;
  }
  __builtin_unreachable();
}
     
//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se) {

  // Determine the optimal overall size for work units.  Given that this process is guided by
  // static (as opposed to forward) exclusion masks, this is a matter of how many atoms are
  // present in all topologies and the number of tiles it would take to cover them all.
  const int nsys = poly_se.getSystemCount();
  std::vector<int> atom_counts(nsys);
  for (int i = 0; i < nsys; i++) {
    atom_counts[i] = poly_se.getAtomCount(i);
  }

  // Try work units of eight tiles
  const size_t tiny_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, tiny_nbwu_tiles);
  const size_t small_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, small_nbwu_tiles);
  const size_t medium_wu_count = estimateNonbondedWorkUnitCount(atom_counts, medium_nbwu_tiles);
  const size_t large_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, large_nbwu_tiles);
  const size_t huge_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, huge_nbwu_tiles);
  if (tiny_wu_count < 16 * kilo) {

    // Implement tiny work units
  }
  else if (small_wu_count < 32 * kilo) {
    
  }
  else if (medium_wu_count < 64 * kilo) {

  }
  else if (large_wu_count < 128 * kilo) {

  }
  else {

    // Use the huge tiles (this will imply the large block size)
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMask &se) {
  const std::vector<int> atom_counts = { se.getAtomCount() };
  const size_t tiny_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, tiny_nbwu_tiles);
  const size_t small_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, small_nbwu_tiles);
  const size_t medium_wu_count = estimateNonbondedWorkUnitCount(atom_counts, medium_nbwu_tiles);
  const size_t large_wu_count  = estimateNonbondedWorkUnitCount(atom_counts, large_nbwu_tiles);
  const size_t huge_wu_count   = estimateNonbondedWorkUnitCount(atom_counts, huge_nbwu_tiles);
  if (tiny_wu_count < 16 * kilo) {
    return enumerateNonbondedWorkUnits(se, tiny_nbwu_tiles, tiny_wu_count, atom_counts);
  }
  else if (small_wu_count < 32 * kilo) {
    return enumerateNonbondedWorkUnits(se, small_nbwu_tiles, small_wu_count, atom_counts);
  }
  else if (medium_wu_count < 64 * kilo) {
    return enumerateNonbondedWorkUnits(se, medium_nbwu_tiles, medium_wu_count, atom_counts);
  }
  else if (large_wu_count < 128 * kilo) {
    return enumerateNonbondedWorkUnits(se, large_nbwu_tiles, large_wu_count, atom_counts);
  }
  else {
    return enumerateNonbondedWorkUnits(se, huge_nbwu_tiles, huge_wu_count, atom_counts);
  }
  __builtin_unreachable();
}

} // namespace synthesis
} // namespace omni
