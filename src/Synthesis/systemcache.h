// -*-c++-*-
#ifndef OMNI_SYSTEM_PREP_H
#define OMNI_SYSTEM_PREP_H

#include <vector>
#include "Constants/behavior.h"
#include "Namelists/nml_files.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"

namespace omni {
namespace synthesis {

using constants::ExceptionResponse;
using namelist::FilesControls;
using topology::AtomGraph;
using trajectory::PhaseSpace;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::detectCoordinateFileKind;

/// \brief Simple but central struct to store the primitive form of the collection of systems that
///        most OMNI applications will work with.  There are topologies, coordinate sets, and a
///        list of indices to indicate how they mix and match.  This object exists on the CPU and
///        in host memory only, but the topology and coordinate contents have CUDA-capable memory
///        components.
struct SystemCache {

  /// \brief The constructor takes the output of &files namelist reading, perhaps with some
  ///        additional command-line edits.  An empty object overload is available because it is
  ///        trivial and works well for delegation.
  ///
  /// \param file_io_input  Object built from a &files namelist and possible command-line edits
  /// \param policy         Response to bad user input, i.e. files of the wrong type
  /// \{
  SystemCache();
  SystemCache(const FilesControls &fcon, ExceptionResponse policy);
  /// \}

  /// \brief Get a pointer to a topology in the cache.
  ///
  /// Overloaded:
  ///   - Get a pointer to a specific topology (this will apply a bounds check)
  ///   - Get a pointer to the array for all topologies (this will not apply a bounds check)
  ///   - Const and non-const versions of each
  ///
  /// \param int index  Index of the requested topology
  /// \{
  const AtomGraph* getTopologyPointer(int index) const;
  AtomGraph* getTopologyPointer(int index);
  const AtomGraph* getTopologyPointer() const;
  AtomGraph* getTopologyPointer();
  /// \}

  /// \brief Get a pointer to a set of coordinates, velocities, and forces in the cache.
  ///
  /// Overloaded:
  ///   - Get a pointer to a specific PhaseSpace object (this will apply a bounds check)
  ///   - Get a pointer to the array for all PhaseSpace objects (this will not check bounds)
  ///   - Const and non-const versions of each
  ///
  /// \param int index  Index of the requested PhaseSpace object
  /// \{
  const PhaseSpace* getPhaseSpacePointer(int index) const;
  PhaseSpace* getPhaseSpacePointer(int index);
  const PhaseSpace* getPhaseSpacePointer() const;
  PhaseSpace* getPhaseSpacePointer();
  /// \}

private:
  /// An array of all topologies to be read by the system: all free topologies and all topologies
  /// read as part of a MoleculeSystem.
  std::vector<AtomGraph> topology_cache;

  /// An array of all coordinate sets to be read by the system: all free coordinate sets and all
  /// coordinates read as part of a MoleculeSystem.
  std::vector<PhaseSpace> coordinates_cache;

  /// The vector of all topology indices guiding each simulation.  This may contain repeats, if
  /// the various MoleculeSystem objects contain the same topology, but the list will be reduced
  /// when composing the synthesis objects.
  std::vector<int> topology_indices;
};
  
} // namespace synthesis
} // namespace omni

#endif
