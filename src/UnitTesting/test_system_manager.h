// -*-c++-*-
#ifndef STORMM_TEST_SYSTEM_MANAGER_H
#define STORMM_TEST_SYSTEM_MANAGER_H

#include "copyright.h"
#include <string>
#include <vector>
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace testing {

using constants::ExceptionResponse;
using topology::AtomGraph;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;
  
/// \brief Provide the means to read a series of topology and coordinate files, then organize them
///        into corresponding lists of AtomGraph and PhaseSpace objects.  The object can also
///        export coordinates in other formats.  Error checking is provided to abort the (test)
///        program or issue a warning, and the object can signal whether dependent tests should
///        be run or aborted.
class TestSystemManager {
public:

  /// \brief The constructor takes lists of topologies, then lists of coordinates, with a
  ///        one-to-one correspondence of files in each list.
  ///
  /// Overloaded:
  ///   - Create a blank object (so that this object can be easily containerized)
  ///   - Accept straight lists of file names for both topologies and the coordinates
  ///   - Accept a common base name and extension for the topologies and the coordinates
  ///   - Accept a common base name for the topologies and the coordinates
  /// \{
  TestSystemManager();

  TestSystemManager(const std::string &topology_base_in, const std::string &topology_extn_in,
                    const std::vector<std::string> &topology_names_in,
                    const std::string &coordinate_base_in, const std::string &coordinate_extn_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);

  TestSystemManager(const std::string &topology_base_in,
                    const std::vector<std::string> &topology_names_in,
                    const std::string &coordinate_base_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);

  TestSystemManager(const std::vector<std::string> &topology_names_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);
  /// \}

  /// \brief Get the number of systems in the object.
  int getSystemCount() const;
  
  /// \brief Return the base path for topologies
  std::string getTopologyBasePath() const;
  
  /// \brief Return the base path for coordinates
  std::string getCoordinateBasePath() const;

  /// \brief Return the extension for topologies
  std::string getTopologyExtension() const;

  /// \brief Return the extension for coordinates
  std::string getCoordinateExtension() const;
  
  /// \brief Return the full name of a topology file according to some index in the list.
  std::string getTopologyFile(int index) const;

  /// \brief Return the full name of a coordinate file according to some index in the list.
  std::string getCoordinateFile(int index) const;

  /// \brief Get the planned course of action for subsequent tests in the event that a file is
  ///        non-existent, or for any reason unreadable.
  TestPriority getTestingStatus() const;

  /// \brief Get a copy of the coordinates for one system as a PhaseSpace object.
  ///
  /// \param index  The system of interest
  PhaseSpace exportPhaseSpace(int index) const;
  
  /// \brief Get a copy of the coordinates for one system as a CoordinateFrame object.
  ///
  /// \param index  The system of interest
  CoordinateFrame exportCoordinateFrame(int index) const;

  /// \brief Get a const reference to the coordinates for one system.
  ///
  /// \param index  The system of interest
  const PhaseSpace& viewCoordinates(int index) const;

  /// \brief Get a const reference to the topology for one or more systems.
  ///
  /// Overloaded:
  ///   - Get a reference to a single topology
  ///   - Get a reference to the vector of all topologies in order
  ///
  /// \param index  The system of interest
  /// \{
  const AtomGraph& getTopologyReference(int index) const;
  const std::vector<AtomGraph>& getTopologyReference() const;
  /// \}
  
  /// \brief Get a pointer to one or more topologies
  ///
  /// Overloaded:
  ///   - Get a pointer to a single topology
  ///   - Get a const vector of pointers to all topologies in order
  ///
  /// \param index  Identifier of a specific topology of interest
  /// \{
  const AtomGraph* getTopologyPointer(int index) const;
  const std::vector<AtomGraph*> getTopologyPointer();
  /// \}
  
private:
  int system_count;              ///< The number of systems managed by this object
  std::string topology_base;     ///< Common base name for all topology files (a directory
                                 ///<   character will be placed between this and any subsequent
                                 ///<   topology file names)
  std::string topology_extn;     ///< Common extension for all topology files (a dot '.' will come
                                 ///<   between this and any preceding coordinate file names)
  std::string coordinate_base;   ///< Common base name for all coordinate files
  std::string coordinate_extn;   ///< Common extension for all coordinate files
  TestPriority fault_response;   ///< Action to take if files from either list are missing
  TestPriority all_go_response;  ///< Action to take if all files are present and parseable
  bool fault_found;              ///< Set to TRUE if any files are unreadable for any reason (if
                                 ///<   any element of the subsequent <bool>-type vectors are TRUE)

  // Vectors of names and flags describing the file parsing
  std::vector<std::string> topology_names;    ///< Names for all topology files
  std::vector<std::string> coordinate_names;  ///< Names for all coordinate files
  std::vector<bool> topology_exist;           ///< Indications that topology files were found
  std::vector<bool> topology_success;         ///< Indications of successful topology parsing
  std::vector<bool> coordinate_exist;         ///< Indications that coordinate files were found
  std::vector<bool> coordinate_success;       ///< Indications of successful coordinate parsing

  /// All topologies (this array contains empty topologies if the corresponding files were
  /// unreadable)
  std::vector<AtomGraph> all_topologies;

  /// All coordinates (this array contains empty coordinates objects if the corresponding files
  /// were unreadable)
  std::vector<PhaseSpace> all_coordinates;

  /// \brief Test that an index does not exceed the valid range for this object.
  ///
  /// \param index   The index of a topology or coordinate file to access
  /// \param caller  Name of the calling function, for error tracing purposes
  void checkIndexing(int index, const char* caller) const;
};

} // namespace testing
} // namespace stormm

#endif
