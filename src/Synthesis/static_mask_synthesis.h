// -*-c++-*-
#ifndef OMNI_STATIC_MASK_SYNTHESIS_H
#define OMNI_STATIC_MASK_SYNTHESIS_H

#include "Accelerator/hybrid.h"
#include "Potential/static_exclusionmask.h"
#include "Topology/atomgraph.h"

namespace omni {
namespace synthesis {

using card::Hybrid;
using energy::StaticExclusionMask;
using topology::AtomGraph;
  
class StaticExclusionMaskSynthesis {
public:

  /// \brief The constructor requires a list of pre-computed StaticExclusionMask objects, each of
  ///        which will contain pointers to the original topologies, and a list of indices stating
  ///        which masks to apply in a given order for the synthesis of systems.
  ///
  /// Overloaded:
  ///   - Accept a list of static exclusion mask object pointers
  ///   - Accept a list of static exclusion mask objects
  ///   - Accept a list of topology pointers, from which exclusion masks will be constructed
  ///
  /// \param base_masks        List of exclusion masks for all unique topologies in the synthesis
  /// \param base_topologies   List of topologies from which to construct the exclusion masks
  ///                          (providing this removes the need to generate the list of exclusion
  ///                          masks outside this function)
  /// \param topology_indices  List of indices into base_masks indicating how to compile the
  ///                          synthesis of systems
  /// \{
  StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask*> &base_masks,
                               const std::vector<int> &topology_indices);

  StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask> &base_masks,
                               const std::vector<int> &topology_indices);

  StaticExclusionMaskSynthesis(const std::vector<AtomGraph*> &base_toplogies,
                               const std::vector<int> &topology_indices);
  /// \}

private:
  Hybrid<int> atom_counts;            ///< Atom counts for all systems
  int unique_supertile_count;         ///< The number of unique supertiles, essentially a sum of
                                      ///<   unique supertiles in all unique masks but omitting the
                                      ///<   zero supertile for all but the first mask.
  int unique_tile_count;              ///< Number of unique tiles compiled from all unique masks
  Hybrid<int> supertile_map_indices;  ///< Indicies into tile_map_indices where each supertile's
                                      ///<   list of tile map indices can be found.  These indices
                                      ///<   are pre-inflated by 256 (the number of tiles in a
                                      ///<   supertile).
  Hybrid<int> supertile_map_bounds;   ///< Bounds array for the supertile_map_indices array above.
                                      ///<   This is the substantive change from single-topology
                                      ///<   static exclusion masks: the bounds array provides
                                      ///<   offsets for each system, avoiding a O(N^2) explosion
                                      ///<   in the supertile indexing, even though such an
                                      ///<   explosion would be somewhat muted by the large size
                                      ///<   of supertiles.  The presence of this array also
                                      ///<   negates the need for a hypothetical array of starting
                                      ///<   indices for each system's atoms in some concatenated
                                      ///<   list, which would not match up with the concatenation
                                      ///<   in a PhaseSpaceSynthesis anyway.
  Hybrid<int> tile_map_indices;       ///< Indices into the all_masks array where each tile's
                                      ///<   exclusions are to be found.  These indices are
                                      ///<   pre-inflated by 32 (see above).
  Hybrid<uint> all_masks;             ///< The actual mask data.  Each tile's exclusion map is
                                      ///<   a stretch of 32 unsigned integers in this array (16
                                      ///<   for the tile's sending atoms and 16 for the tile's
                                      ///<   receiving atoms).  

  /// \brief Encapsulate the contents of the constructor to permit multiple pathways for creating
  ///        the object.
  ///
  /// \param base_masks        List of exclusion masks for all unique topologies in the synthesis
  /// \param topology_indices  List of indices into base_masks indicating how to compile the
  ///                          synthesis of systems
  void build(const std::vector<StaticExclusionMask*> &base_masks,
             const std::vector<int> &topology_indices);
};
  
} // namespace synthesis
} // namespace omni

#endif
