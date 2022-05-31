#include "static_mask_synthesis.h"

namespace omni {
namespace synthesis {

using card::HybridKind;
using energy::StaticExclusionMaskReader;
using energy::tile_length;
using energy::tiles_per_supertile;

//-------------------------------------------------------------------------------------------------
SeMaskSynthesisReader::SeMaskSynthesisReader(const int nsys_in, const int* atom_counts_in,
                                             const int nsupertile_in, const int ntile_in,
                                             const int* supertile_map_idx_in,
                                             const int* supertile_map_bounds_in,
                                             const int* tile_map_idx_in,
                                             const uint* mask_data_in) :
    nsys{nsys_in}, atom_counts{atom_counts_in}, nsupertile{nsupertile_in}, ntile{ntile_in},
    supertile_map_idx{supertile_map_idx_in}, supertile_map_bounds{supertile_map_bounds_in},
    tile_map_idx{tile_map_idx_in}, mask_data{mask_data_in}
{}

//-------------------------------------------------------------------------------------------------
StaticExclusionMaskSynthesis::
StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask*> &base_masks,
                             const std::vector<int> &topology_indices) :
    atom_counts{},
    unique_supertile_count{0},
    unique_tile_count{0},
    supertile_map_indices{HybridKind::ARRAY, "sesyn_supertile_idx"},
    supertile_map_bounds{HybridKind::ARRAY, "sesyn_supertile_bnd"},
    tile_map_indices{HybridKind::ARRAY, "sesyn_tile_idx"},
    all_masks{HybridKind::ARRAY, "sesyn_mask_data"}
{
  build(base_masks, topology_indices);
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMaskSynthesis::
StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask> &base_masks,
                             const std::vector<int> &topology_indices) :
    atom_counts{},
    unique_supertile_count{0},
    unique_tile_count{0},
    supertile_map_indices{HybridKind::ARRAY, "sesyn_supertile_idx"},
    supertile_map_bounds{HybridKind::ARRAY, "sesyn_supertile_bnd"},
    tile_map_indices{HybridKind::ARRAY, "sesyn_tile_idx"},
    all_masks{HybridKind::ARRAY, "sesyn_mask_data"}
{
  const size_t ntop = base_masks.size();
  std::vector<StaticExclusionMask*> base_mask_ptrs(ntop);
  for (size_t topid = 0; topid < ntop; topid++) {
    base_mask_ptrs[topid] = const_cast<StaticExclusionMask*>(&base_masks[topid]);
  }
  build(base_mask_ptrs, topology_indices);
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMaskSynthesis::
StaticExclusionMaskSynthesis(const std::vector<AtomGraph*> &base_topologies,
                             const std::vector<int> &topology_indices) :
    atom_counts{},
    unique_supertile_count{0},
    unique_tile_count{0},
    supertile_map_indices{HybridKind::ARRAY, "sesyn_supertile_idx"},
    supertile_map_bounds{HybridKind::ARRAY, "sesyn_supertile_bnd"},
    tile_map_indices{HybridKind::ARRAY, "sesyn_tile_idx"},
    all_masks{HybridKind::ARRAY, "sesyn_mask_data"}
{
  const size_t ntop = base_topologies.size();
  std::vector<StaticExclusionMask> base_masks;
  std::vector<StaticExclusionMask*> base_mask_ptrs;
  base_masks.reserve(ntop);
  base_mask_ptrs.resize(ntop);
  for (size_t topid = 0; topid < ntop; topid++) {
    base_masks.emplace_back(base_topologies[topid]);
    base_mask_ptrs[topid] = &base_masks[topid];
  }
  build(base_mask_ptrs, topology_indices);
}

//-------------------------------------------------------------------------------------------------
void StaticExclusionMaskSynthesis::build(const std::vector<StaticExclusionMask*> &base_masks,
                                         const std::vector<int> &topology_indices) {

  // Compute the numbers of unique tiles and other properties of the array of unique masks
  unique_supertile_count = 1;
  unique_tile_count = 1;
  const int ntop = base_masks.size();
  std::vector<int> supertile_index_offsets(ntop, 0);
  std::vector<int> tile_index_offsets(ntop, 0);
  for (int topid = 0; topid < ntop; topid++) {
    supertile_index_offsets[topid] = unique_supertile_count * tiles_per_supertile;
    tile_index_offsets[topid] = unique_tile_count * 2 * tile_length;
    unique_supertile_count += base_masks[topid]->getUniqueSuperTileCount() - 1;
    unique_tile_count += base_masks[topid]->getUniqueTileCount() - 1;
  }
  tile_map_indices.resize(tiles_per_supertile * unique_supertile_count);
  all_masks.resize(unique_tile_count * 2 * tile_length);
  
  // Lay out bounds arrays and basic descriptors spanning all systems
  const int nsys = topology_indices.size();
  atom_counts.resize(nsys);
  supertile_map_bounds.resize(nsys + 1);
  int supertile_acc = 0;
  supertile_map_bounds.putHost(supertile_acc, 0);
  for (int sysid = 0; sysid < nsys; sysid++) {
    atom_counts.putHost(base_masks[topology_indices[sysid]]->getAtomCount(), sysid);
    const int stcount = base_masks[topology_indices[sysid]]->getSuperTileStrideCount();
    supertile_acc += stcount * stcount;
    supertile_map_bounds.putHost(supertile_acc, sysid + 1);
  }
  supertile_map_indices.resize(supertile_acc);

  // Compute the supertile map indices based on the forthcoming consensus tables
  int* stile_map_idx_ptr = supertile_map_indices.data();
  for (int sysid = 0; sysid < nsys; sysid++) {
    const StaticExclusionMaskReader ser = base_masks[topology_indices[sysid]]->data();
    const int stcount = ser.supertile_stride_count;
    const int stm_bound = supertile_map_bounds.readHost(sysid);
    const int sttop_offset = supertile_index_offsets[topology_indices[sysid]];
    for (int i = 0; i < stcount * stcount; i++) {

      // Supertile 0 in any topology's static exclusion mask is supertile 0 in the synthesis
      // exclusion mask, but other supertiles are offset based on prior topologies incorporated
      // into the synthesis.
      if (ser.supertile_map_idx[i] == 0) {
        stile_map_idx_ptr[stm_bound + i] = 0;
      }
      else {
        stile_map_idx_ptr[stm_bound + i] = ser.supertile_map_idx[i] + sttop_offset -
                                           tiles_per_supertile;
      }
    }
  }
    
  // Transfer supertile information to the synthesis
  int* tile_map_idx_ptr = tile_map_indices.data();
  uint* all_masks_ptr = all_masks.data();
  for (int i = 0; i < tiles_per_supertile; i++) {
    tile_map_idx_ptr[i] = 0;
  }
  for (int i = 0; i < 2 * tile_length; i++) {
    all_masks_ptr[i] = 0U;
  }
  int tdata_idx = 2 * tile_length;
  int tbase = tiles_per_supertile;
  for (int topid = 0; topid < ntop; topid++) {
    const StaticExclusionMaskReader ser = base_masks[topid]->data();

    // Copy the unique tile data
    const int nuni_t = (base_masks[topid]->getUniqueTileCount() - 1) * 2 * tile_length;
    for (int i = 0; i < nuni_t; i++) {
      all_masks_ptr[tdata_idx + i] = ser.mask_data[i + (2 * tile_length)];
    }
    tdata_idx += nuni_t;
    
    // Skip the first supertile, which corresponds to zero interactions.  That is already in the
    // list of unique consensus supertiles from work above.
    const int nuni_st = base_masks[topid]->getUniqueSuperTileCount() * tiles_per_supertile;
    const int ttop_offset = tile_index_offsets[topid];
    for (int i = tiles_per_supertile; i < nuni_st; i++) {
      if (ser.tile_map_idx[i] == 0) {
        tile_map_idx_ptr[i + tbase - tiles_per_supertile] = 0;
      }
      else {
        tile_map_idx_ptr[i + tbase - tiles_per_supertile] = ser.tile_map_idx[i] + ttop_offset -
                                                            (2 * tile_length);
      }
    }
    tbase += (base_masks[topid]->getUniqueSuperTileCount() - 1) * tiles_per_supertile;
  }
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMaskSynthesis::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int StaticExclusionMaskSynthesis::getAtomCount(const int index) const {
  return atom_counts.readHost(index);
}
  
//-------------------------------------------------------------------------------------------------
SeMaskSynthesisReader StaticExclusionMaskSynthesis::data(const HybridTargetLevel tier) const {
  return SeMaskSynthesisReader(system_count, atom_counts.data(tier), unique_supertile_count,
                               unique_tile_count, supertile_map_indices.data(tier),
                               supertile_map_bounds.data(tier), tile_map_indices.data(tier),
                               all_masks.data(tier));
}

} // namespace synthesis
} // namespace omni
