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
class SystemCache {
public:

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

  /// \brief Get the number of systems
  int getSystemCount() const;

  /// \brief Get the number of topologies
  int getTopologyCount() const;

  /// \brief Get the topology index of one of the coordinate sets contained in this cache.  This
  ///        will apply a bounds check to the coordinate index query.  This function should be used
  ///        to access topologies in the output of a getTopologyReference() or getTopologyPointer()
  ///        call when the full array is returned without function input arguments.  Calling either
  ///        getTopologyReference() or getTopologyPointer() with an integer argument implicitly
  ///        performs the operation of getting the appropriate topology list index, as would be
  ///        found by this function, and thus returns a reference or pointer to the appropriate
  ///        topology.
  ///
  /// \param int coord_index  Index of the PhaseSpace entry object of interest
  int getSystemTopologyIndex(int coord_index) const;

  /// \brief Get the index of a coordinate set which provides an example of the system that one of
  ///        the topologies in the cache describes.
  ///
  /// \param int topology_index  Index of the topology of interest
  int getCoordinateExample(int topology_index) const;

  /// \brief Get pointers to one or all of the unique topologies from within the compact topology
  ///        cache.
  ///
  /// Overloaded:
  ///   - Get a pointer to one topology
  ///   - Get a vector of pointers to all unique topologies in the SystemCache.
  ///
  /// \param int topology_index  Index of the topology of interest
  /// \{
  const AtomGraph* getTopologyPointer(int topology_index) const;
  std::vector<const AtomGraph*> getTopologyPointer() const;
  /// \}

  /// \brief Get a reference to one of the unique topologies from within the topology cache, based
  ///        on an index into the cache itself.
  ///
  /// \param int topology_index  Index of the topology of interes
  const AtomGraph& getTopologyReference(int topology_index) const;  
  
  /// \brief Get a pointer to a topology in the cache associated with a particular coordinate set.
  ///
  /// Overloaded:
  ///   - Get a pointer to the topology for a specific coordinate set (this will apply a bounds
  ///     check)
  ///   - Get pointers to topologies for all coordinate sets (this will not apply a bounds check)
  ///   - Const and non-const versions of each
  ///
  /// \param int index  Index of the PhaseSpace entry of interest
  /// \{
  const AtomGraph* getSystemTopologyPointer(int index) const;
  AtomGraph* getSystemTopologyPointer(int index);
  std::vector<const AtomGraph*> getSystemTopologyPointer() const;
  std::vector<AtomGraph*> getSystemTopologyPointer();
  /// \}

  /// \brief Get a const std::vector of non-const (const-casted) topology pointers.  This is for
  ///        the specific case of feeding other functions that require std::vectors of topology
  ///        pointers in this format.
  const std::vector<AtomGraph*> getSystemTopologyPointerCC() const;

  /// \brief Return a reference to the topology that describes a particular set of coordinates
  ///        within the SystemCache.  All overloads apply a bounds check.
  ///
  /// Overloaded:
  ///   - Return a const reference to the relevant topology in a const SystemCache
  ///   - Return a non-const reference to the relevant topology in a non-const SystemCache
  ///   - Return a reference to the entire array of AtomGraph objects (this is a raw take on the
  ///     underlying topologies and will need to be accessed according to the indexing present in
  ///     the topology_indices array)
  ///
  /// \param int index  Index of the PhaseSpace entry of interest
  /// \{
  const AtomGraph& getSystemTopologyReference(const int index) const;
  AtomGraph& getSystemTopologyReference(const int index);
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
  const PhaseSpace* getCoordinatePointer(int index) const;
  PhaseSpace* getCoordinatePointer(int index);
  std::vector<const PhaseSpace*> getCoordinatePointer() const;
  std::vector<PhaseSpace*> getCoordinatePointer();
  /// \}

  /// \brief Return a reference to a set of coordinates, velocities, and forces in the cache.
  ///
  /// Overloaded:
  ///   - Return a const reference to an object in a const SystemCache
  ///   - Return a non-const reference to an object in a non-const SystemCache
  ///   - Return a reference to the entire array of PhaseSpace objects
  ///
  /// \param int index  Index of the PhaseSpace object of interest
  /// \{
  const PhaseSpace& getCoordinateReference(const int index) const;
  PhaseSpace& getCoordinateReference(const int index);
  const std::vector<PhaseSpace>& getCoordinateReference() const;
  std::vector<PhaseSpace>& getCoordinateReference();
  /// \}

  /// \brief Get the number of systems described by a topology of the given index in this cache.
  ///
  /// \param topology_index  The index of the topology of interest
  int getTopologyCaseCount(int topology_index) const;

  /// \brief Get a list of indices for all systems in the cache described by the topology of
  ///        interest.
  ///
  /// \param topology_index  The index of the topology of interest
  std::vector<int> getTopologicalCases(int topology_index) const;

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

  /// Each topology describes one or more systems in the cache.  This array will store the index
  /// of a coordinate set which serves as an example of each topology.
  std::vector<int> example_indices;

  /// Locations of all systems using a given topology.  This list is accessed using the bounds
  /// array in the following member variable: the content between indices given by elements [k]
  /// and [k + 1] of the bounds array shows all systems using the kth topology in the cache's
  /// list.
  std::vector<int> topology_cases;

  /// Bounds array for topology_cases above
  std::vector<int> topology_case_bounds;
};
  
} // namespace synthesis
} // namespace omni

#endif
