#include "nonbonded_workunit.h"

namespace omni {
namespace synthesis {

using energy::tile_length;
  
//-------------------------------------------------------------------------------------------------
NonbondedWorkUnit::NonbondedWorkUnit(const AtomGraph *ag_in, const StaticExclusionMask *se_in) {
    abscissa_atom_count{0}, ordinate_atom_count{0}, tile_count{0}, tile_instructions{},
    ag_pointer{ag_in}, se_pointer{se_in}
{}

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
