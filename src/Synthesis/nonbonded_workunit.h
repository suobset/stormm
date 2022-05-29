// -*-c++-*-
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

  /// \brief The constructor accepts an exclusion mask and a list of nonbonded interaction tiles
  ///        to compute.  Non-bonded interaction tiles go according to "abscissa atoms" and
  ///        "ordinate atoms," although both abscissa and ordinate indices refer to the same list
  ///        of imported atoms.  There are three cases to consider:
  ///
  ///        Static exclusion mask, tiny up to large work units:
  ///        This will be the majority of the cases with implicit-solvent systems.  In practice,
  ///        the non-bonded work unit becomes a list of 19 integers representing the upper limit
  ///        of atoms in the particular system and up to 16 starting locations of tile atoms to
  ///        import (and convert to floating point representations) from within what may be a
  ///        concatenated list of atoms for many systems.  The number of atoms after each starting
  ///        point is a matter of the tile_length constant, and to some degree this couples it to
  ///        the block size used for smaller non-bonded work units.  There are then low and high
  ///        limits on the tile computations to perform, indexing into a master list of uint2
  ///        tuples, each listing the abscissa atom start and ordinate atom starts in bits 1-16 and
  ///        17-32 of the x member (referring to an index of local atom imports) and the tile
  ///        exclusion mask index (referring to an index of a separate master list out in main
  ///        memory) in the y member.  These cases will run with a 256-thread block size.
  ///
  ///        Static exclusion mask, huge work units:
  ///        This will handle cases of implicit-solvent systems with sizes so large that millions
  ///        of smaller work units would be required to cover everything.  In practice, the
  ///        non-bonded work unit is reduced to a list of 4 integers, now representing the upper
  ///        limit of atoms in the system and the lower limits of the abscissa and ordinate atoms
  ///        to import in a supertile for which the work unit is to compute all interactions.
  ///        There is no meaningful list of all interactions in this case, as it might be
  ///        prohibitive even to store such a thing.  Instead, the work unit will proceed over all
  ///        tiles in the supertile after computing whether it lies along the diagonal.  This will
  ///        require a larger block size (512 threads minimum, up to 768 depending on the
  ///        architecture).
  ///
  ///        Forward exclusion mask:
  ///        This will handle cases of neighborlist-based non-bonded work units.  The work unit
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
  /// \param se_in        Static exclusion mask for a system in isolated boundary conditions
  /// \param tile_list    Paired with a static exclusion mask and non-huge tiles, the specific list
  ///                     of lower left corners for tiles to include in this work unit.  Among
  ///                     them, the tiles must not call for importing more than 256 atoms.
  /// \param tile_corner  The lower limits of the supertile to process.  Paired with a static
  ///                     exclusion mask in extreme circumstances of very large implicit solvent
  ///                     systems.  This will call for importing 512 atoms (2 x supertile_length)
  ///                     in the most general case and will require larger thread blocks.
  NonbondedWorkUnit(const StaticExclusionMask *se_in, const std::vector<int2> &tile_list);

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
