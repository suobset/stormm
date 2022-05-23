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
int getAllSystemsWorkUnitCount(const std::vector<int> &atom_counts, const int nbwu_tile) {

  // There are pre-determined tile counts for the non-bonded work units, each with their own
  // width and height in terms of the tile layout.
  int absc_tile_dim, ordi_tile_dim;
  switch (nbwu_tile) {
  case 8:
    absc_tile_dim = 4;
    ordi_tile_dim = 2;
  case 16:
    absc_tile_dim = 4;
    ordi_tile_dim = 4;
  case 32:
    absc_tile_dim = 4;
    ordi_tile_dim = 8;
  case 64:
    absc_tile_dim = 8;
    ordi_tile_dim = 8;
  default:
    rtErr("Only 8, 16, 32, or 64 tiles may be the target size of any one NonbondedWorkUnit.",
          "getAllSystemsWorkUnitCount");
  }
  
  const size_t nsys = atom_counts.size();
  int result = 0;
  for (size_t i = 0LLU; i < nsys; i++) {
    const int ntile_side  = atom_counts[i] / tile_length;
    const int nabsc_whole = tile_length / absc_tile_dim;
    const int nordi_whole = tile_length / ordi_tile_dim;
    const int nabsc_rem   = tile_length - (nabsc_whole * absc_tile_dim);

    // If there is a partial tile along the side, find out how many work units the trim can
    // decompose into.  The abscissa and ordinate atom arrays must add to at most 256, which
    // means sixteen tile lengths in both directions, combined.  Each work unit along the trim
    // must have no more than the number of tiles in any other work unit.
    int absc_idx = nabsc_whole * absc_tile_dim;
    int ordi_idx = 0;
    while (ordi_idx < ntile_side) {

    }
  }
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
  size_t tile_pair = 
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
