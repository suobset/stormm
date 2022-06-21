// -*-c++-*-
#ifndef OMNI_NONBONDED_WORKUNIT_H
#define OMNI_NONBONDED_WORKUNIT_H

#include "Topology/atomgraph.h"
#include "Potential/static_exclusionmask.h"
#include "static_mask_synthesis.h"
#include "synthesis_enumerators.h"

namespace omni {
namespace synthesis {

using energy::StaticExclusionMask;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
  
/// \brief The maximum number of imported atoms in a "small" nonbonded block (256 threads)
constexpr int small_block_max_atoms = 320;

/// \brief The maximum number of tile sides' worth of atoms that can be imported
constexpr int small_block_max_imports = small_block_max_atoms / tile_length;

/// \brief The number of tiles that a small block will be designed to process simultaneously
constexpr int small_block_tile_width = 8;

/// \brief Tile counts for various sizes of non-bonded work units that use the small block size
/// \{
constexpr int tiny_nbwu_tiles   =     small_block_tile_width;
constexpr int small_nbwu_tiles  = 2 * small_block_tile_width;
constexpr int medium_nbwu_tiles = 4 * small_block_tile_width;
constexpr int large_nbwu_tiles  = 8 * small_block_tile_width;
/// \}

/// \brief Tile count for non-bonded work units using the large block size
constexpr int huge_nbwu_tiles   = tile_lengths_per_supertile * tile_lengths_per_supertile;

/// \brief Length of the abstract for a non-bonded work unit based on all-to-all tile groups
constexpr int tile_groups_wu_abstract_length = 64;

/// \brief Length of the abstract for a non-bonded work unit based on all-to-all supertiles
constexpr int supertile_wu_abstract_length = 8;

/// \brief Collect a series of tiles for non-bonded computations as well as the required atom
///        imports to carry them out.  This will accomplish the task of planning the non-bonded
///        computation, given a single topology or a synthesis of topologies, to optimize the
///        thread utilization on a GPU.
class NonbondedWorkUnit {
public:

  /// \brief The constructor accepts an exclusion mask and a list of nonbonded interaction tiles
  ///        to compute.  Non-bonded interaction tiles go according to "abscissa atoms" and
  ///        "ordinate atoms," although both abscissa and ordinate indices refer to the same list
  ///        of imported atoms.  There are three cases to consider:
  ///
  ///        Static exclusion mask, tiny up to large work units:
  ///        This will be the majority of the cases with implicit-solvent systems.  In practice,
  ///        the non-bonded work unit becomes a list of 28 integers representing the total number
  ///        of tile imports (up to 20) followed by the starting locations of each tile's atoms to
  ///        import (and convert to floating point representations) from within what may be a
  ///        concatenated list of atoms for many systems.  Five integers following these import
  ///        starting positions indicate the numbers of atoms to import in each tile (usually the
  ///        tile length, but fewer if the tile is incomplete at the end of the system's list of
  ///        atoms).  The number of atoms after each starting point is a matter of the tile_length
  ///        constant, and to some degree this couples it to the block size used for smaller
  ///        non-bonded work units.  There are then low and high limits on the tile computations to
  ///        perform, indexing into a master list of uint2 tuples, each listing the abscissa atom
  ///        start and ordinate atom starts in bits 1-16 and 17-32 of the x member (referring to
  ///        an index of local atom imports) and the tile exclusion mask index (referring to an
  ///        index of a separate master list out in main memory) in the y member.  These cases
  ///        will run with a 256-thread block size.
  ///
  ///        Static exclusion mask, huge work units:
  ///        This will handle cases of implicit-solvent systems with sizes so large that millions
  ///        of smaller work units would be required to cover everything.  In practice, the
  ///        non-bonded work unit is reduced to a list of 4 integers, now representing the lower
  ///        limits of the abscissa and ordinate atoms to import, and the numbers of atoms to
  ///        import along each axis, in a supertile for which the work unit is to compute all
  ///        interactions.  There is no meaningful list of all interactions in this case, as it
  ///        might be prohibitive even to store such a thing.  Instead, the work unit will proceed
  ///        over all tiles in the supertile after computing whether it lies along the diagonal.
  ///        This will require a larger block size (512 threads minimum, up to 768 depending on
  ///        the architecture).
  ///
  ///        Forward exclusion mask:
  ///        This will handle cases of neighbor list-based non-bonded work units.  The work unit
  ///        assumes a honeycomb-packed image of all atoms in or about the primary unit cell (some
  ///        honeycomb pencils will straddle the unit cell boundary but their positions will be
  ///        known as part of the decomposition).  The work unit will consist of thirty integers:
  ///        seven starting locations of atom imports, seven bit-packed integers detailing the
  ///        lengths of each stretch of atoms (first 20 bits) and the obligatory y- and z- imaging
  ///        moves to make with such atoms (last 12 bits), seven integers detailing segments of
  ///        each stretch of imported atoms to replicate in +x (and where in the list of imported
  ///        atoms to put them), and finally seven integers detailing segments of each stretch of
  ///        imported atoms to replicate in -x (and where to put the replicas).  The final two
  ///        integers state the starting and ending indices of a list of tile instructions to
  ///        process.  The tile instructions for the neighbor list-based non-bonded work units are
  ///        much more complex than those for non-bonded work based on a static exclusion mask.
  ///
  /// \param se              Static exclusion mask, or synthesis thereof, for one or more systems
  ///                        in isolated boundary conditions
  /// \param tile_list       Paired with a static exclusion mask and non-huge tiles, the specific
  ///                        list of tiles to include for this work unit in the x- and y-members,
  ///                        plus the system index number in the z member (if more than one system
  ///                        is present in a synthesis).  Among them, the tiles must not call for
  ///                        importing more than small_block_max_atoms atoms.
  /// \param abscissa_start  The abscissa axis start of the supertile to process.  Paired with a
  ///                        static exclusion mask in extreme circumstances of very large implicit
  ///                        solvent systems.  This will call for importing 512 atoms
  ///                        (2 x supertile_length) in the most general case and will require
  ///                        larger thread blocks.  If computing for a synthesis of static
  ///                        exclusion masks, the abscissa starting point is given as a relative
  ///                        index within the local system to which the supertile belongs.
  /// \param ordinate_start  The ordinate axis start of the supertile to process.  If computing for
  ///                        a synthesis of static exclusion masks, the ordinate starting point is
  ///                        given as a relative index within the local system to which the
  ///                        supertile belongs.
  /// \{
  NonbondedWorkUnit(const StaticExclusionMask &se, const std::vector<int3> &tile_list);
  
  NonbondedWorkUnit(const StaticExclusionMask &se, int abscissa_start, int ordinate_start);

  NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se,
                    const std::vector<int3> &tile_list);

  NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se, int abscissa_start, int ordinate_start,
                    int system_index);
  /// \}

  /// \brief Get the tile count of this work units.
  int getTileCount() const;

  /// \brief Get the abscissa and ordinate atom limits for a tile from within this work unit.
  ///        The abscissa limits are returned in the x and y members, the ordinate limits in the
  ///        z and w members.
  ///
  /// \param index  Index of the tile from within the work unit
  int4 getTileLimits(int index) const;

  /// \brief Get the list of tile instructions for this work unit.
  const std::vector<uint2>& getTileInstructions() const;
  
  /// \brief Get an abstract for this work unit, layed out to work within an AtomGraphSynthesis.
  ///
  /// \param instruction_start  The starting point of instructions for this group of tiles, if
  ///                           the work unit will fit on a small thread block, having less than
  ///                           huge_nbwu_tiles tiles.  Otherwise no instruction start is needed.
  std::vector<int> getAbstract(int instruction_start = 0) const;
  
private:
  int tile_count;                          ///< Number of tiles to be processed by this work unit
  NbwuKind kind;                           ///< The type of non-bonded work unit.  For isolated
                                           ///<   boundary conditions, there is a choice between
                                           ///<   TILE_GROUPS and SUPERTILES.
  int score;                               ///< The estimated effort score of this work unit
  std::vector<int> imports;                ///< List of starting positions for atoms that must be
                                           ///<   cached in order to process all tiles in this
                                           ///<   work unit
  std::vector<int> import_system_indices;  ///< System indices of each imported tile.  All atoms
                                           ///<   in any one tile pertain to the same system, but
                                           ///<   each imported tile within the work unit may be
                                           ///<   part of a different system.
  std::vector<int> import_size_keys;       ///< Bit-packed integers with the number of atoms to
                                           ///<   import after each starting position
  std::vector<uint2> tile_instructions;    ///< Instructions for processing each tile based on an
                                           ///<   appropriate exclusion mask

  /// \brief Find the system to which a particular import belongs.  This function is only used in
  ///        the event that the non-bonded work unit pertains to synthesis of systems, otherwise
  ///        the system index is obviously zero.
  ///
  /// \param ser         Reader abstract for the exclusion mask synthesis
  /// \param atom_index  Index of the first imported atom in the tile of interest
  int getImportSystemIndex(const SeMaskSynthesisReader &ser, int atom_index);
};

/// \brief Add one tile to the growing list that will eventually define a work unit, if it is
///        possible to do so.  This encapsualtes some work at a rather high cost in terms of
///        pointers for multiple returned values.  Return true if the tile addition was successful
///        or false if not.
///
/// \param tile_list           The growing list of tiles, pre-allocated to be able to hold any
///                            additions that this function might make.  Appended and returned.
/// \param import_coverage     Array of 0's and 1's (an int, rather than boolean, array for
///                            versatility in passing to and from this function and for direct
///                            summation with other integers).  Edited and returned.
/// \param import_count        The total number of imported tiles' worth of atoms.  Updated and
///                            returned.
/// \param current_tile_count  The current number of tiles, updarted and returned.
/// \param ti                  Tile index of the abscissa atoms (starting atom index divided by
///                            tile_length)
/// \param tj                  Tile index of the ordinate atoms (starting atom index divided by
///                            tile_length)
/// \param sysid               System index from which the tiles are derived
bool addTileToWorkUnitList(int3* tile_list, int* import_coverage, int *import_count,
                           int *current_tile_count, const int ti, const int tj, const int sysid);

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
buildNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se);

std::vector<NonbondedWorkUnit> buildNonbondedWorkUnits(const StaticExclusionMask &se);
/// \}

} // namespace synthesis
} // namespace omni

#include "nonbonded_workunit.tpp"

#endif
