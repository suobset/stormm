// -*-c++-*-
#ifndef OMNI_SYSTEM_PREP_H
#define OMNI_SYSTEM_PREP_H

#include <vector>
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace omni {
namespace synthesis {

/// \brief Produce definitive lists of unique starting coordinates and the topologies and restraint
///        vectors that govern their movement.  This unrolls collections of frames, replicas, and
///        interpolations of various topologies.
///
/// \param file_io_input  Object built from a &files namelist and possible command-line edits
/// \param top_cache      Cache of topologies to be produced by unrolling the &files input
///                       (modified and returned)
/// \param crd_cache      Cache of (initial) coordinates to be produced by unrolling the &files
///                       input (modified and returned)
/// \param top_indices    List of topology indices linking initial coordinates to the topologies
///                       that govern them (modified and returned)
void unrollSystemLists(const FilesControls &fcon, std::vector<AtomGraph> *top_cache,
                       std::vector<PhaseSpace> *crd_cache, std::vector<int> *top_indices);
  
} // namespace synthesis
} // namespace omni

#endif
