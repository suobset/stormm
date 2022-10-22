// -*-c++-*-
#ifndef STORMM_SYSTEM_PREP_H
#define STORMM_SYSTEM_PREP_H

#include <vector>
#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Constants/behavior.h"
#include "FileManagement/file_enumerators.h"
#include "MoleculeFormat/mdlmol.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_restraint.h"
#include "Potential/forward_exclusionmask.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"
#include "UnitTesting/stopwatch.h"

namespace stormm {
namespace synthesis {

using chemistry::ChemicalFeatures;
using chemistry::MapRotatableGroups;
using constants::ExceptionResponse;
using diskutil::PrintSituation;
using energy::ForwardExclusionMask;
using energy::StaticExclusionMask;
using namelist::FilesControls;
using namelist::RestraintControls;
using restraints::RestraintApparatus;
using structure::MdlMol;
using testing::StopWatch;
using topology::AtomGraph;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFileRole;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::detectCoordinateFileKind;
using trajectory::PhaseSpace;
using trajectory::TrajectoryFusion;

/// \brief Simple but central struct to store the primitive form of the collection of systems that
///        most STORMM applications will work with.  There are topologies, coordinate sets, and a
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
  /// \param policy_in      Response to bad user input, i.e. files of the wrong type
  /// \{
  SystemCache(ExceptionResponse policy_in = ExceptionResponse::DIE,
              MapRotatableGroups map_chemfe_rotators = MapRotatableGroups::NO,
              PrintSituation expectation_in = PrintSituation::OPEN_NEW);

  SystemCache(const FilesControls &fcon, const std::vector<RestraintControls> &rstcon,
              std::vector<MdlMol> *sdf_recovery,
              ExceptionResponse policy_in = ExceptionResponse::DIE,
              MapRotatableGroups map_chemfe_rotators = MapRotatableGroups::NO,
              PrintSituation expectation_in = PrintSituation::OPEN_NEW,
              StopWatch *timer_in = nullptr);

  SystemCache(const FilesControls &fcon, std::vector<MdlMol> *sdf_recovery,
              ExceptionResponse policy = ExceptionResponse::DIE,
              MapRotatableGroups map_chemfe_rotators = MapRotatableGroups::NO,
              PrintSituation expectation_in = PrintSituation::OPEN_NEW,
              StopWatch *timer_in = nullptr);

  SystemCache(const FilesControls &fcon, const std::vector<RestraintControls> &rstcon,
              ExceptionResponse policy_in = ExceptionResponse::DIE,
              MapRotatableGroups map_chemfe_rotators = MapRotatableGroups::NO,
              PrintSituation expectation_in = PrintSituation::OPEN_NEW,
              StopWatch *timer_in = nullptr);

  SystemCache(const FilesControls &fcon, ExceptionResponse policy = ExceptionResponse::DIE,
              MapRotatableGroups map_chemfe_rotators = MapRotatableGroups::NO,
              PrintSituation expectation_in = PrintSituation::OPEN_NEW,
              StopWatch *timer_in = nullptr);
  /// \}

  /// \brief The default copy and move constructors, copy and move assignment operators will be
  ///        effective for this object based entirely on Standard Template Library components.  
  /// \{
  SystemCache(const SystemCache &original) = default;
  SystemCache(SystemCache &&original) = default;
  SystemCache& operator=(const SystemCache &other) = default;
  SystemCache& operator=(SystemCache &&other) = default;
  /// \}
  
  /// \brief Get the number of systems
  int getSystemCount() const;

  /// \brief Get the number of topologies
  int getTopologyCount() const;

  /// \brief Get the number of unique labels across all systems
  int getLabelCount() const;

  /// \brief Get the number of different restraint apparatuses.  Many restraint apparatuses may be
  ///        unique in the sense that they have different stiffness and displacement settings,
  ///        even though they apply their restraints to the same atoms of different systems guided
  ///        by the same topology.
  int getRestraintCount() const;

  /// \brief Get the topology index of one of the coordinate sets contained in this cache.  This
  ///        will apply a bounds check to the coordinate index query.  This function should be used
  ///        to access topologies in the output of a getTopologyReference() or getTopologyPointer()
  ///        call when the full array is returned without function input arguments.  Calling either
  ///        getTopologyReference() or getTopologyPointer() with an integer argument implicitly
  ///        performs the operation of getting the appropriate topology list index, as would be
  ///        found by this function, and thus returns a reference or pointer to the appropriate
  ///        topology.
  ///
  /// \param coord_index  Index of the PhaseSpace entry object of interest
  int getSystemTopologyIndex(int coord_index) const;

  /// \brief Get the index of a coordinate set which provides an example of the system that one of
  ///        the topologies in the cache describes.
  ///
  /// \param topology_index  Index of the topology of interest
  int getCoordinateExample(int topology_index) const;

  /// \brief Get pointers to one or all of the unique topologies from within the compact topology
  ///        cache.
  ///
  /// Overloaded:
  ///   - Get a pointer to one topology
  ///   - Get a vector of pointers to all unique topologies in the SystemCache.
  ///
  /// \param topology_index  Index of the topology of interest
  /// \{
  const AtomGraph* getTopologyPointer(int topology_index) const;
  std::vector<const AtomGraph*> getTopologyPointer() const;
  /// \}

  /// \brief Get a reference to one of the unique topologies from within the topology cache, based
  ///        on an index into the cache itself.
  ///
  /// \param topology_index  Index of the topology of interes
  const AtomGraph& getTopologyReference(int topology_index) const;  
  
  /// \brief Get a pointer to a topology in the cache associated with a particular coordinate set.
  ///
  /// Overloaded:
  ///   - Get a pointer to the topology for a specific coordinate set (this will apply a bounds
  ///     check)
  ///   - Get pointers to topologies for all coordinate sets (this will not apply a bounds check)
  ///   - Const and non-const versions of each
  ///
  /// \param index  Index of the PhaseSpace entry of interest
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
  /// \param index  Index of the PhaseSpace entry of interest
  /// \{
  const AtomGraph& getSystemTopologyReference(int index) const;
  AtomGraph& getSystemTopologyReference(int index);
  /// \}
  
  /// \brief Get a pointer to a set of coordinates, velocities, and forces in the cache.
  ///
  /// Overloaded:
  ///   - Get a pointer to a specific PhaseSpace object (this will apply a bounds check)
  ///   - Get a pointer to the array for all PhaseSpace objects (this will not check bounds)
  ///   - Const and non-const versions of each
  ///
  /// \param index  Index of the requested PhaseSpace object
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
  /// \param index  Index of the system, the PhaseSpace object of interest
  /// \{
  const PhaseSpace& getCoordinateReference(int index) const;
  PhaseSpace& getCoordinateReference(int index);
  const std::vector<PhaseSpace>& getCoordinateReference() const;
  std::vector<PhaseSpace>& getCoordinateReference();
  /// \}

  /// \brief Get a pointer to the chemical features for a particular system.
  ///
  /// Overloaded:
  ///   - Return a const pointer to an object in a const SystemCache
  ///   - Return a non-const pointer to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const ChemicalFeatures* getFeaturesPointer(int index) const;
  ChemicalFeatures* getFeaturesPointer(int index);
  /// \}  

  /// \brief Get a reference to the chemical features for a particular system.
  ///
  /// Overloaded:
  ///   - Return a const reference to an object in a const SystemCache
  ///   - Return a non-const reference to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const ChemicalFeatures& getFeaturesReference(int index) const;
  ChemicalFeatures& getFeaturesReference(int index);
  /// \}

  /// \brief Get a pointer to the restraint apparatus for a particular system.
  ///
  /// Overloaded:
  ///   - Return a const pointer to an object in a const SystemCache
  ///   - Return a non-const pointer to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const RestraintApparatus* getRestraintPointer(int index) const;
  RestraintApparatus* getRestraintPointer(int index);
  /// \}

  /// \brief Get a reference to the restraint apparatus for a particular system.
  ///
  /// Overloaded:
  ///   - Return a const reference to an object in a const SystemCache
  ///   - Return a non-const reference to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const RestraintApparatus& getRestraintReference(int index) const;
  RestraintApparatus& getRestraintReference(int index);
  /// \}

  /// \brief Get a pointer to the static exclusion mask for a particular system.  These masks
  ///        will only have been calculated for systems with isolated boundary conditions.
  ///
  /// Overloaded:
  ///   - Return a const pointer to an object in a const SystemCache
  ///   - Return a non-const pointer to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const StaticExclusionMask* getSystemStaticMaskPointer(int index) const;
  StaticExclusionMask* getSystemStaticMaskPointer(int index);
  /// \}  

  /// \brief Get a reference to the static exclusion mask for a particular system.  These masks
  ///        will only have been calculated for systems with isolated boundary conditions.
  ///
  /// Overloaded:
  ///   - Return a const reference to an object in a const SystemCache
  ///   - Return a non-const reference to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const StaticExclusionMask& getSystemStaticMaskReference(int index) const;
  StaticExclusionMask& getSystemStaticMaskReference(int index);
  /// \}

  /// \brief Get a pointer to the forward exclusion mask for a particular system.  These masks
  ///        will only have been calculated for systems with periodic boundary conditions.
  ///
  /// Overloaded:
  ///   - Return a const pointer to an object in a const SystemCache
  ///   - Return a non-const pointer to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const ForwardExclusionMask* getSystemForwardMaskPointer(int index) const;
  ForwardExclusionMask* getSystemForwardMaskPointer(int index);
  /// \}  

  /// \brief Get a reference to the forward exclusion mask for a particular system.  These masks
  ///        will only have been calculated for systems with periodic boundary conditions.
  ///
  /// Overloaded:
  ///   - Return a const reference to an object in a const SystemCache
  ///   - Return a non-const reference to an object in a non-const SystemCache
  ///
  /// \param index  Index of the system of interest
  /// \{
  const ForwardExclusionMask& getSystemForwardMaskReference(int index) const;
  ForwardExclusionMask& getSystemForwardMaskReference(int index);
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

  /// \brief Get the input coordinate file name for a particular system.  This is not subject to
  ///        manipulations to differentiate the name as the output formats are.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  const std::string& getInputCoordinatesName(int system_index) const;
  
  /// \brief Get the name of the trajectory file associated with one of the systems.  The name
  ///        will be manipulated (by adding _## before the final '.', or at the end of the file
  ///        name if there is no '.').
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  std::string getTrajectoryName(int system_index) const;
  
  /// \brief Get the name of the checkpoint (restart) file associated with one of the systems.
  ///        This name will also be manipulated to differentiate systems making use of the same
  ///        label.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  std::string getCheckpointName(int system_index) const;

  /// \brief Get the label associated with a particular system in the cache.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  std::string getSystemLabel(int system_index) const;

  /// \brief Get the coordinate file type associated with a particular system's input.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  CoordinateFileKind getInputCoordinatesKind(const int system_index) const;
  
  /// \brief Get the coordinate file type associated with a particular system's trajectory.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  CoordinateFileKind getTrajectoryKind(const int system_index) const;

  /// \brief Get the coordinate file type associated with a particular system's checkpoint file.
  ///
  /// \param system_index  Index of the system from within the coordinates cache
  CoordinateFileKind getCheckpointKind(const int system_index) const;

  /// \brief When multiple systems' coordinates are to be combined into a particular file, the
  ///        file may be opened as a new file once, with whatever overwriting protocol, to write
  ///        the first set of coordinates.  Thereafter, the printing protocol must change: it
  ///        should be appended with other systems from the same label group.  This function will
  ///        determine how to modulate the printing protcol.
  ///
  /// Overloaded:
  ///   - Call with no arguments to get the baseline printing protocol, and that used for the
  ///     first frames written under any label group
  ///   - Call with the role of a coordinates file and the system index to get the modified
  ///     protocol specific to writing frames of that system
  ///
  /// \param purpose      The role that the file will play
  /// \param systm_index  A particular system for which to obtain the printing protocol, which may
  ///                     be modified from the original if the system shares the same label as any
  ///                     others
  /// \{
  PrintSituation getPrintingProtocol() const;
  PrintSituation getPrintingProtocol(CoordinateFileRole purpose, int system_index) const;
  /// \}

  /// \brief Set the printing protocol.
  ///
  /// \param expectation_in  The protocol to use in printing any files
  void setPrintingProtocol(const PrintSituation expectation_in);

private:

  /// The object will store the exception response behavior under which it was created, so that
  /// output irregularities may be handled in the same manner as input errors.
  ExceptionResponse policy;

  /// The SystemCache is built to organize a (potentially large) collection of systems for
  /// molecular mechanics calculations and then manage the output stream.  In a scenario where a
  /// program built on STORMM executes and then shuts down, this will mean file printing.  In
  /// a situation where a function based on STORMM takes inputs from a python environment and then
  /// returns results, it may still need to decide whether to overwrite or add to a growing array.
  /// This parameter indicates the SystemCache's way to handle existing data when writing it own.
  PrintSituation expectation;
  
  // Counts are kept as separate integers, to have a single authoritative number on some central
  // quantities in this systems cache.
  int system_count;       ///< An official record of the total number of systems in the cache
  int topology_count;     ///< The number of unique topologies
  int label_count;        ///< The number of unique labels
  int restraint_count;    ///< The number of unique restraint setups.  Each RestraintApparatus
                          ///<   object must pertain to a specific topology and each system may
                          ///<   reference but one RestraintApparatus, but many RestraintApparatus
                          ///<   objects may involve the same topology in order to guide different
                          ///<   systems.

  /// The question of whether to fuse trajectories and checkpoint files is pertinent to a program
  /// that is designed to run hundreds to thousands of calculations simultaneously.  The default
  /// behavior will be to fuse checkpoint files for systems grouped under the same label into a
  /// single file, if the provided format seems fit to accommodate that (in fact, there is no
  /// existing molecular dynamics file format that can handle multiple frames of both positions
  /// and velocities simultaneously, but for docking calculations and energy minimizations when
  /// the only output is the positions fusing the final states of many systems into a single
  /// file can be a great convenience).
  TrajectoryFusion file_merger_protocol;
  
  /// An array of all topologies to be read by the system: all free topologies and all topologies
  /// read as part of a MoleculeSystem.
  std::vector<AtomGraph> topology_cache;

  /// An array of all coordinate sets to be read by the system: all free coordinate sets and all
  /// coordinates read as part of a MoleculeSystem.
  std::vector<PhaseSpace> coordinates_cache;

  /// A collection of the unique labels across all systems
  std::vector<std::string> label_cache;
  
  /// Chemical features objects outline the important aspects of each system.
  std::vector<ChemicalFeatures> features_cache;

  /// Restraint collections supplement the topological energy surface in each system.
  std::vector<RestraintApparatus> restraints_cache;

  /// Static exclusion masks serve each topology, if appropriate
  std::vector<StaticExclusionMask> static_masks_cache;

  /// Forward exclusion masks serve each topology, if appropriate
  std::vector<ForwardExclusionMask> forward_masks_cache;

  /// The vector of all topology indices guiding each simulation.  This may contain repeats, if
  /// the various MoleculeSystem objects contain the same topology, but the list will be reduced
  /// when composing the synthesis objects.
  std::vector<int> topology_indices;

  /// The vector of all systems' indices into the array of unique labels.
  std::vector<int> label_indices;

  /// The number of systems sharing each label in the label cache
  std::vector<int> label_degeneracy;
  
  /// When systems share the same label, the degeneracy implies that each system can only be
  /// distinguished by a unique sub-index within that label group.
  std::vector<int> label_subindices;

  /// The vector of all restraint indices guiding each simulation.  Like topology_indices, this
  /// list may contain multiple systems pointing to the same restraint apparatus, especially if
  /// some systems are subject to no restraints.
  std::vector<int> restraint_indices;
  
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

  /// Names of coordinate input files for each system (for records keeping purposes)
  std::vector<std::string> system_input_coordinate_names;
  
  /// Names of trajectory files to write for each system
  std::vector<std::string> system_trajectory_names;

  /// Names of checkpoint files to write for each system
  std::vector<std::string> system_checkpoint_names;

  /// Labels applied to each system (whether from user input or auto-generated)
  std::vector<std::string> system_labels;

  /// File types for the various input coordinate files
  std::vector<CoordinateFileKind> system_input_coordinate_kinds;

  /// File types for the various output trajectory files
  std::vector<CoordinateFileKind> system_trajectory_kinds;

  /// File types for the various output checkpoint files
  std::vector<CoordinateFileKind> system_checkpoint_kinds;

  /// \brief Check that the requested index is within the number of systems present in the cache.
  ///
  /// \param index   The requested system index
  /// \param caller  Name of the calling function (for backtracing purposes after an error)
  void checkSystemBounds(int index, const char* caller) const;

  /// \brief Determine whether output files from multiple systems should be concatenated or fused.
  ///        This can happen in particular cases:
  ///        - The role of the files is checkpointing and the format provided is able to
  ///          accept multiple frames
  ///        - The role of the files is trajectories and the user has specified "concatenate files
  ///          from the same label"
  ///
  /// \param purpose       The role that the files written will play
  /// \param system_index  The index of the system about to be written--only if the system is
  ///                      grouped with others under a single, degenerate label will file fusion
  ///                      come into play
  bool determineFileMerger(CoordinateFileRole purpose, int system_index) const;
  
  /// \brief When multiple replicas are guided by one topology fall under the same system label,
  ///        it is not clear what to do when each of them is to be written to a trajectory or a
  ///        checkpoint file.  This function will examine the case of a particular system, with a
  ///        given name, and extend the name with a unique replica number to eliminate any
  ///        degeneracy in the
  ///
  /// \param fname_in      Name of the file to print.  In the event of a name collision with
  ///                      other systems under the same label, this name will be extended based on
  ///                      a system sub-index within that label group.  The extension will take
  ///                      the form of "_<sub-index>" just before the final '.' in the name, or at
  ///                      the end of the name if there is no '.' after the final OS directory
  ///                      separator.
  /// \param purpose       The role that the file will play, used to distinguish checkpoint files
  ///                      from trajectories
  /// \param system_index  System index out of all systems in the cache (this will be used to find
  ///                      the label index and its degeneracy)
  std::string nondegenerateName(const std::string &fname_in, CoordinateFileRole purpose,
                                int system_index) const;
};
  
} // namespace synthesis
} // namespace stormm

#endif
