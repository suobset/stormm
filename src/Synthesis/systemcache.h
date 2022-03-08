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
  ///        to access topologies in the output of a getTopologyReference() call when the full
  ///        array is returned without function input arguments, but not to access the output og
  ///        getTopologyPointer() as that array is returned indexing the topology pointers with
  ///        the coordinate systems (just use coord_index as seen in input to this function to
  ///        access that array directly).
  ///
  /// \param int coord_index  Index of the PhaseSpace entry object of interest
  int getTopologyIndex(int coord_index) const;

  /// \brief Get the index of a coordinate set which provides an example of the system that one of
  ///        the topologies in the cache describes.
  ///
  /// \param int topology_index  Index of the topology of interest
  int getCoordinateExample(int topology_index) const;
  
  /// \brief Get a pointer to a topology in the cache associated with a particular coordinate set.
  ///
  /// Overloaded:
  ///   - Get a pointer to the topology for a specific coordinate set (this will apply a bounds
  ///     check)
  ///   - Get pointers to topologies for all coordinate sets (this will not apply a bounds check)
  ///   - Const and non-const versions of each
  ///
  /// \param int index  Index of the PhaseSpace entry object of interest
  /// \{
  const AtomGraph* getTopologyPointer(int index) const;
  AtomGraph* getTopologyPointer(int index);
  std::vector<const AtomGraph*> getTopologyPointer() const;
  std::vector<AtomGraph*> getTopologyPointer();
  /// \}

  /// \brief Get a const std::vector of non-const (const-casted) topology pointers.  This is for
  ///        the specific case of feeding other functions that require std::vectors of topology
  ///        pointers in this format.
  const std::vector<AtomGraph*> getTopologyPointerCC() const;

  /// \brief Return a reference to the topology that describes a particular set of coordinates
  ///        within the SystemCache.  All overloads apply a bounds check.
  ///
  /// Overloaded:
  ///   - Return a const reference to a const SystemCache
  ///   - Return a non-const reference to a non-const SystemCache
  ///   - Return a reference to the entire array of AtomGraph objects (this is a raw take on the
  ///     underlying topologies and will need to be accessed according to the indexing present in
  ///     the topology_indices array)
  ///
  /// \param int index  Index of the PhaseSpace entry object of interest
  /// \{
  const AtomGraph& getTopologyReference(const int index) const;
  AtomGraph& getTopologyReference(const int index);
  const std::vector<AtomGraph>& getTopologyReference() const;
  std::vector<AtomGraph>& getTopologyReference();
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
  /// \param int index  Index of the PhaseSpace entry object of interest
  /// \{
  const PhaseSpace& getCoordinateReference(const int index) const;
  PhaseSpace& getCoordinateReference(const int index);
  const std::vector<PhaseSpace>& getCoordinateReference() const;
  std::vector<PhaseSpace>& getCoordinateReference();
  /// \}

  /// \brief Get the indexing of a given 
  
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
};
  
} // namespace synthesis
} // namespace omni

#endif
