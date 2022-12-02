// -*-c++-*-
#ifndef STORMM_PHASE_SPACE_H
#define STORMM_PHASE_SPACE_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Math/matrix_ops.h"
#include "Topology/atomgraph.h"
#include "trajectory_enumerators.h"
#include "write_frame.h"

namespace stormm {
namespace trajectory {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using math::computeBoxTransform;
using topology::AtomGraph;
using topology::UnitCellType;

/// \brief Collect constants and pointers to the components of a modifiable PhaseSpace object.
struct PhaseSpaceWriter {

  /// Constructor takes a large list of arguments passed in from the original PhaseSpace object
  PhaseSpaceWriter(const int natom_in, const UnitCellType unit_cell_in, double* xcrd_in,
                   double* ycrd_in, double* zcrd_in, double* umat_in, double* invu_in,
                   double* boxdim_in, double* xvel_in, double* yvel_in, double* zvel_in,
                   double* xfrc_in, double* yfrc_in, double* zfrc_in, double* xprv_in,
                   double* yprv_in, double* zprv_in, double* xnxt_in, double* ynxt_in,
                   double* znxt_in, double* vxprv_in, double* vyprv_in, double* vzprv_in,
                   double* vxnxt_in, double* vynxt_in, double* vznxt_in, double* fxprv_in,
                   double* fyprv_in, double* fzprv_in, double* fxnxt_in, double* fynxt_in,
                   double* fznxt_in);

  /// \brief Copy and move constructors.  The assignment operators are implicitly deleted.
  /// \{
  PhaseSpaceWriter(const PhaseSpaceWriter &original) = default;
  PhaseSpaceWriter(PhaseSpaceWriter &&original) = default;
  /// \}

  const int natom;               ///< Atom count for this system (still a constant)
  const UnitCellType unit_cell;  ///< The type of unit cell
  double* xcrd;                  ///< Cartesian X positions of all particles
  double* ycrd;                  ///< Cartesian Y positions of all particles
  double* zcrd;                  ///< Cartesian Z positions of all particles
  double* umat;                  ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  double* invu;                  ///< Transformation matrix to take coordinates into real space
  double* boxdim;                ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  double* xvel;                  ///< Cartesian X velocities of all particles
  double* yvel;                  ///< Cartesian Y velocities of all particles
  double* zvel;                  ///< Cartesian Z velocities of all particles
  double* xfrc;                  ///< Cartesian X forces acting on all particles
  double* yfrc;                  ///< Cartesian Y forces acting on all particles
  double* zfrc;                  ///< Cartesian Z forces acting on all particles
  double* xprv;                  ///< Prior Cartesian X positions of all particles
  double* yprv;                  ///< Prior Cartesian Y positions of all particles
  double* zprv;                  ///< Prior Cartesian Z positions of all particles
  double* xnxt;                  ///< Next, upcoming Cartesian X positions of all particles
  double* ynxt;                  ///< Next, upcoming Cartesian Y positions of all particles
  double* znxt;                  ///< Next, upcoming Cartesian Z positions of all particles 
  double* vxprv;                 ///< Prior Cartesian X velocities for all particles
  double* vyprv;                 ///< Prior Cartesian Y velocities for all particles
  double* vzprv;                 ///< Prior Cartesian Z velocities for all particles
  double* vxnxt;                 ///< Next (planned) Cartesian X velocities of particles
  double* vynxt;                 ///< Next (planned) Cartesian Y velocities of particles
  double* vznxt;                 ///< Next (planned) Cartesian Z velocities of particles
  double* fxprv;                 ///< Prior Cartesian X forces acting on all particles
  double* fyprv;                 ///< Prior Cartesian Y forces acting on all particles
  double* fzprv;                 ///< Prior Cartesian Z forces acting on all particles
  double* fxnxt;                 ///< Next (currently undergoing updates) forces acting on all
                                 ///<   particles in the Cartesian X direction
  double* fynxt;                 ///< Next forces on all particles in the Cartesian Y direction
  double* fznxt;                 ///< Next forces on all particles in the Cartesian Z direction
};

/// \brief Collect constants and pointers to the components of a read-only PhaseSpace object.
struct PhaseSpaceReader {

  /// The constructor takes a large list of arguments passed in from the original PhaseSpace
  /// object, or the cognate writer to make all of the associated pointers const.
  /// \{
  PhaseSpaceReader(const int natom_in, const UnitCellType unit_cell_in, const double* xcrd_in,
                   const double* ycrd_in, const double* zcrd_in, const double* umat_in,
                   const double* invu_in, const double* boxdim_in, const double* xvel_in,
                   const double* yvel_in, const double* zvel_in, const double* xfrc_in,
                   const double* yfrc_in, const double* zfrc_in, const double* xprv_in,
                   const double* yprv_in, const double* zprv_in, const double* xnxt_in,
                   const double* ynxt_in, const double* znxt_in, const double* vxprv_in,
                   const double* vyprv_in, const double* vzprv_in, const double* vxnxt_in,
                   const double* vynxt_in, const double* vznxt_in, const double* fxprv_in,
                   const double* fyprv_in, const double* fzprv_in, const double* fxnxt_in,
                   const double* fynxt_in, const double* fznxt_in);

  PhaseSpaceReader(const PhaseSpaceWriter &psw);
  /// \}
  
  /// \brief Copy and move constructors.  The assignment operators are implicitly deleted.
  /// \{
  PhaseSpaceReader(const PhaseSpaceReader &original) = default;
  PhaseSpaceReader(PhaseSpaceReader &&original) = default;
  /// \}
  
  const int natom;               ///< Atom count for this system (still a constant)
  const UnitCellType unit_cell;  ///< The type of unit cell
  const double* xcrd;            ///< Cartesian X positions of all particles
  const double* ycrd;            ///< Cartesian Y positions of all particles
  const double* zcrd;            ///< Cartesian Z positions of all particles
  const double* umat;            ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  const double* invu;            ///< Transformation matrix to take coordinates into real space
  const double* boxdim;          ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  const double* xvel;            ///< Cartesian X velocities of all particles
  const double* yvel;            ///< Cartesian Y velocities of all particles
  const double* zvel;            ///< Cartesian Z velocities of all particles
  const double* xfrc;            ///< Cartesian X forces acting on all particles
  const double* yfrc;            ///< Cartesian Y forces acting on all particles
  const double* zfrc;            ///< Cartesian Z forces acting on all particles
  const double* xprv;            ///< Prior Cartesian X positions of all particles
  const double* yprv;            ///< Prior Cartesian Y positions of all particles
  const double* zprv;            ///< Prior Cartesian Z positions of all particles
  const double* xnxt;            ///< Next, upcoming Cartesian X positions of all particles
  const double* ynxt;            ///< Next, upcoming Cartesian Y positions of all particles
  const double* znxt;            ///< Next, upcoming Cartesian Z positions of all particles
  const double* vxprv;           ///< Prior Cartesian X velocities for all particles
  const double* vyprv;           ///< Prior Cartesian Y velocities for all particles
  const double* vzprv;           ///< Prior Cartesian Z velocities for all particles
  const double* vxnxt;           ///< Next (planned) Cartesian X velocities of particles
  const double* vynxt;           ///< Next (planned) Cartesian Y velocities of particles
  const double* vznxt;           ///< Next (planned) Cartesian Z velocities of particles
  const double* fxprv;           ///< Prior Cartesian X forces acting on all particles
  const double* fyprv;           ///< Prior Cartesian Y forces acting on all particles
  const double* fzprv;           ///< Prior Cartesian Z forces acting on all particles
  const double* fxnxt;           ///< Next (currently undergoing updates) forces acting on all
                                 ///<   particles in the Cartesian X direction
  const double* fynxt;           ///< Next forces on all particles in the Cartesian Y direction
  const double* fznxt;           ///< Next forces on all particles in the Cartesian Z direction
};

/// \brief An object to complement a topology and hold positions, velocities, and forces of all
///        particles in a system.  This is not designed to be the most performant representation of
///        the system's structure.  Rather, it serves to hold a high-precision representation of a
///        single system and transport it between CPUs and high-performance accelerators.
class PhaseSpace {
public:

  /// \brief Construction of a phase space object, like a topology, is typically done from a file.
  ///
  /// Overloaded:
  ///   - Constructor fot an object with a number of atoms but no other information (assists in
  ///     delegation of initialization)
  ///   - Constructors for coordinate sets read from trajectory or restart files, with the option
  ///     of a specific frame number (if unspecified, the first frame is read).  A topology may
  ///     also be specified to check the atom count and sanity of the coordinates presented.
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param ag_reference  Topology with which to check atom count, proximity of bonded atoms, and
  ///                      sanity of bond angles (optional)
  /// \param caller        Name of the calling function (to help backtrace errors)
  /// \{
  PhaseSpace(int atom_count_in = 0, UnitCellType unit_cell_in = UnitCellType::NONE);

  PhaseSpace(const std::string &file_name_in,
             CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN, int frame_number = 0);

  PhaseSpace(const std::string &file_name_in, const AtomGraph &ag,
             CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
             int frame_number = 0);
  /// \}

  /// \brief Copy constructor handles assignment of internal POINTER-kind Hybrid objects
  ///
  /// \param original  The PhaseSpace object from which to make a deep copy
  PhaseSpace(const PhaseSpace &original);

  /// \brief Copy assignment operator likewise handles assignment of internal POINTER-kind Hybrid
  ///        objects
  ///  
  /// \param other     Another way to say original, in a different semantic context
  PhaseSpace& operator=(const PhaseSpace &other);

  /// \brief The move constructor prepares the original PhaseSpace object for destruction
  ///
  /// \param original  The PhaseSpace object from which to preserve content
  PhaseSpace(PhaseSpace &&original);

  /// \brief The move assignment operator looks much like the copy assignment operator.
  ///  
  /// \param other     Another way to say original, in a different semantic context
  PhaseSpace& operator=(PhaseSpace &&other);

  /// \brief Fill the object from information in some coordinate, restart, or trajectory file.
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param frame_number  The frame number to read (if the file is a trajectory, not a single
  ///                      point from the system's phase space)
  void buildFromFile(const std::string &file_name_in,
                     const CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                     int frame_number = 0);

  /// \brief Fill the object from information in three arrays.
  ///
  /// Overloaded:
  ///   - Fill from three C-style arrays
  ///   - Fill from three Standard Template Library vector objects
  ///
  /// \param xcrd        Cartesian X coordinates of positions, velocities, or forces
  /// \param ycrd        Cartesian Y coordinates of positions, velocities, or forces
  /// \param zcrd        Cartesian Z coordinates of positions, velocities, or forces
  /// \param kind        Type of coordinates coming in: fill the positions, velocities, or forces
  /// \param cycle_in    The point in the coordinate cycle to fill
  /// \param scale_bits  The number of bits after the decimal, applicable to fixed-precision
  ///                    representations of xcrd, ycrd, and zcrd (the box dimensions are always
  ///                    given as a double-precision array, in units of Angstroms)
  /// \param box_dims    Box dimensions, from which the tranformation matrices will be derived
  /// \{
  template <typename T>
  void fill(const T* xcrd, const T* ycrd, const T* zcrd,
            TrajectoryKind kind = TrajectoryKind::POSITIONS,
            CoordinateCycle cycle_in = CoordinateCycle::PRESENT, int scale_bits = 0,
            const double* box_dims = nullptr);

  template <typename T>
  void fill(const std::vector<T> &xcrd, const std::vector<T> &ycrd, const std::vector<T> &zcrd,
            TrajectoryKind kind = TrajectoryKind::POSITIONS,
            CoordinateCycle cycle_in = CoordinateCycle::PRESENT, int scale_bits = 0,
            const std::vector<double> &box_dims = {});
  /// \}

  /// \brief Get the name of the file associated with this object.
  std::string getFileName() const;

  /// \brief Get the number of atoms (particles, including virtual sites) in the object.
  int getAtomCount() const;

  /// \brief Get the unit cell type of the coordinate system
  UnitCellType getUnitCellType() const;

  /// \brief Get the time cycle stage indicating the arrays holding current coordinates.
  CoordinateCycle getCyclePosition() const;  
  
  /// \brief Get a pointer to the particle X, Y, or Z coordinates, velocities, or forces, on either
  ///        the host or device.  Use this when the entire abstract is unnecessary or would be
  ///        inefficient to retrieve.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///
  /// \param dim   Cartesian dimension of interest
  /// \param kind  Specify coordinates, velocities, or forces--anything that could be thought
  ///              of as a trajectory
  /// \param tier  Level at which to extract the data
  /// \{
  const double* getCoordinatePointer(CartesianDimension dim,
                                     TrajectoryKind kind = TrajectoryKind::POSITIONS,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getCoordinatePointer(CartesianDimension dim,
                               TrajectoryKind kind = TrajectoryKind::POSITIONS,
                               HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Get the coordinates returned in an X/Y/Z interlaced manner
  ///
  /// Overloaded:
  ///   - Get all coordinates
  ///   - Get coordinates for a range of atoms
  ///
  /// \param low_index   The lower atom index of a range
  /// \param high_index  The upper atom index of a range
  /// \param kind        Specify coordinates, velocities, or forces--anything that could be thought
  ///                    of as a trajectory
  /// \param tier        Level at which to extract the data
  /// \{
  std::vector<double>
  getInterlacedCoordinates(TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double>
  getInterlacedCoordinates(int low_index, int high_index,
                           TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a pointer to the box space transformation matrix that can track its evolution in
  ///        the PhaseSpace object's Hybrid data arrays, on either the host or device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///
  /// \param tier  Level at which to extract the data
  /// \{
  const double*
  getBoxSpaceTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSpaceTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get a pointer to the inverse (back to real space) transformation matrix that can track
  ///        its evolution in the PhaseSpace object's Hybrid data arrays, on either the host or
  ///        device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///
  /// \param tier  Level at which to extract the data
  /// \{
  const double* getInverseTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getInverseTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get a pointer to the box dimensions that can track its evolution in the PhaseSpace
  ///        object's Hybrid data arrays, on either the host or device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///
  /// \param tier  Level at which to extract the data
  /// \{
  const double* getBoxSizePointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSizePointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get the transformation matrix to take coordinates into box (fractional) space.  The
  ///        result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if STORMM is compiled to run on a GPU)
  std::vector<double> getBoxSpaceTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the transformation matrix to take coordinates from fractional space back into
  ///        real space.  The result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if STORMM is compiled to run on a GPU)
  std::vector<double> getInverseTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a pointer to the ARRAY-kind Hybrid object that holds the actual data.  Needed by
  ///        CoordinateFrame objects which want to be pointers into a PhaseSpace object.
  ///
  /// Overloaded:
  ///   - Get a const pointer to a const PhaseSpace object's data storage
  ///   - Get a non-const pointer to a mutable PhaeSpace object's data storage
  /// \{
  const Hybrid<double>* getStoragePointer() const;
  Hybrid<double>* getStoragePointer();
  /// \}

  /// \brief Initialize the forces (set them to zero)
  ///
  /// Overloaded:
  ///   - Update forces for an arbitrary point in the time cycle
  ///   - Update forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void initializeForces(CoordinateCycle orientation);
  void initializeForces();
  /// \}
  
  /// \brief Update the cycle position.
  ///
  /// Overloaded:
  ///   - Advance the cycle position based on its current setting (no input argument):
  ///     present >> future >> past >> present >> ...
  ///   - Set the cycle position to an arbitrary point
  ///
  /// \param  time_point  The point in the time cycle that shall become the PhaseSpace object's
  ///                     "present" coordinates.
  /// \{
  void updateCyclePosition();
  void updateCyclePosition(CoordinateCycle time_point);
  /// \}
  
  /// \brief Put the phase space data into a trajectory or checkpoint file.
  ///
  /// \param file_name     Name of the file to write
  /// \param current_time  Time progress of the simulation
  /// \param traj_kind     The type of trajectory to print (coordinates, velocities, or forces)
  /// \param output_kind   The format of the file to write (checkpoint files print position and
  ///                      velocity data by obligation, but trajectory files can contain either of
  ///                      these as well as forces)
  /// \param expectation   The condition in which the output file is expected to be found
  void exportToFile(const std::string &file_name, double current_time = 0.0,
                    TrajectoryKind traj_kind = TrajectoryKind::POSITIONS,
                    CoordinateFileKind output_kind = CoordinateFileKind::AMBER_INPCRD,
                    PrintSituation expectation = PrintSituation::UNKNOWN) const;
  
#ifdef STORMM_USE_HPC
  /// \brief Upload all information
  void upload();

  /// \brief Upload positional information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload particle positions for a specific point in the time cycle
  ///   - Upload particle positions for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void uploadPositions(CoordinateCycle orientation);
  void uploadPositions();
  /// \}

  /// \brief Upload current transformation matrices.
  void uploadTransformations();

  /// \brief Upload velocity information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload velocities for a specific point in the time cycle
  ///   - Upload velocities for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void uploadVelocities(CoordinateCycle orientation);
  void uploadVelocities();
  /// \}

  /// \brief Upload force information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload forces for a specific point in the time cycle
  ///   - Upload forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void uploadForces(CoordinateCycle orientation);
  void uploadForces();
  /// \}

  /// \brief Download all information
  void download();

  /// \brief Download positional information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download particle positions for a specific point in the time cycle
  ///   - Download particle positions for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void downloadPositions(CoordinateCycle orientation);
  void downloadPositions();
  /// \}

  /// \brief Download current transformation matrices.
  void downloadTransformations();

  /// \brief Download velocity information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download velocities for a specific point in the time cycle
  ///   - Download velocities for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void downloadVelocities(CoordinateCycle orientation);
  void downloadVelocities();
  /// \}

  /// \brief Download force information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download forces for a specific point in the time cycle
  ///   - Download forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (PAST, PRESENT, or FUTURE)
  /// \{
  void downloadForces(CoordinateCycle orientation);
  void downloadForces();
  /// \}
#endif
  /// \brief Get the abstract for this object, containing C-style pointers for the most rapid
  ///        access to any of its member variables.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const PhaseSpace object
  ///   - Get a writeable abstract from a mutable PhaseSpace object
  ///   - Get either object oriented with the past, present, or future positional arrays set as
  ///     holding the current coordinates (three such abstracts can be rotated over successive
  ///     cycles of dynamics to let the coordinates evolve, protected against race conditions,
  ///     without swapping the actual locations in memory)
  ///
  /// \param tier         Specify pointers on the host or device
  /// \param orientation  Arbitrarily selected point on the time cycle to have the reader or writer
  ///                     take as the current coordinates
  /// \{
  const PhaseSpaceReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  const PhaseSpaceReader data(CoordinateCycle orientation,
                              HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PhaseSpaceWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  PhaseSpaceWriter data(CoordinateCycle orientation,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

private:
  std::string file_name;               ///< Name of the file from which these coordinates (and
                                       ///<   perhaps velocities) derived.  Empty string indicates
                                       ///<   no file.
  int atom_count;                      ///< The number of atoms in the system
  UnitCellType unit_cell;              ///< The type of unit cell
  CoordinateCycle cycle_position;      ///< Indicates the place in the past >> present >> future
                                       ///<   cycle where the object is currently storing its
                                       ///<   relevant coordinates.  After a dynamics step, present
                                       ///<   becomes past, future becomes present, and what were
                                       ///<   the arrays holding past coordinates stand ready to
                                       ///<   accept the future configuration as it is assembled.
                                       ///<   With the work unit system, constraint updates might
                                       ///<   not occur atomically.  Three arrays, not just two as
                                       ///<   are used in other codes, are needed to protect
                                       ///<   against race conditions.
  Hybrid<double> x_coordinates;        ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;        ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;        ///< Cartesian Z coordinates of all particles

  // The actual nature of the following arrays, as well as the [x,y,z]_coordinates themselves, can
  // change based on the CoordinateCycle position.
  Hybrid<double> x_prior_coordinates;  ///< Previous step Cartesian X coordinates of all particles
  Hybrid<double> y_prior_coordinates;  ///< Previous step Cartesian Y coordinates of all particles
  Hybrid<double> z_prior_coordinates;  ///< Previous step Cartesian Z coordinates of all particles
  Hybrid<double> x_future_coordinates; ///< Upcoming step Cartesian X coordinates of all particles
  Hybrid<double> y_future_coordinates; ///< Upcoming step Cartesian Y coordinates of all particles
  Hybrid<double> z_future_coordinates; ///< Upcoming step Cartesian Z coordinates of all particles

  // Transformation matrices are given directly after the coordinates in the overall order of data
  Hybrid<double> box_space_transform;  ///< Matrix to transform coordinates into box space (3 x 3)
  Hybrid<double> inverse_transform;    ///< Matrix to transform coordinates into real space (3 x 3)
  Hybrid<double> box_dimensions;       ///< Three lengths and three angles defining the box
                                       ///<   (lengths are given in Angstroms, angles in radians)

  // Like coordinates, velocities and forces appear in PRESENT -> PAST -> FUTURE blocks
  Hybrid<double> x_velocities;         ///< Cartesian X velocities of all particles
  Hybrid<double> y_velocities;         ///< Cartesian Y velocities of all particles
  Hybrid<double> z_velocities;         ///< Cartesian Z velocities of all particles
  Hybrid<double> x_prior_velocities;   ///< Prior Cartesian X velocities of all particles
  Hybrid<double> y_prior_velocities;   ///< Prior Cartesian Y velocities of all particles
  Hybrid<double> z_prior_velocities;   ///< Prior Cartesian Z velocities of all particles
  Hybrid<double> x_future_velocities;  ///< Planned Cartesian X velocities of all particles
  Hybrid<double> y_future_velocities;  ///< Planned Cartesian Y velocities of all particles
  Hybrid<double> z_future_velocities;  ///< Planned Cartesian Z velocities of all particles
  Hybrid<double> x_forces;             ///< Cartesian X forces acting on all particles
  Hybrid<double> y_forces;             ///< Cartesian Y forces acting on all particles
  Hybrid<double> z_forces;             ///< Cartesian Z forces acting on all particles 
  Hybrid<double> x_prior_forces;       ///< Prior Cartesian X forces acting on all particles
  Hybrid<double> y_prior_forces;       ///< Prior Cartesian Y forces acting on all particles
  Hybrid<double> z_prior_forces;       ///< Prior Cartesian Z forces acting on all particles
  Hybrid<double> x_future_forces;      ///< Upcoming Cartesian X forces acting on all particles
  Hybrid<double> y_future_forces;      ///< Upcoming Cartesian Y forces acting on all particles
  Hybrid<double> z_future_forces;      ///< Upcoming Cartesian Z forces acting on all particles

  /// All of the above Hybrid objects are pointers into this single large array, segmented to hold
  /// each type of information with zero-padding to accommodate the HPC warp size.
  Hybrid<double> storage;

  /// \brief Allocate space for the object, based on a known number of atoms
  void allocate();
};

/// \brief Interpret the inverse transformation matrix to determine whether a unit cell exists, and
///        if it is rectilinear or not.
///
/// Overloaded:
///   - Accept a C-style pointer to the array of matrix elements
///   - Accept a Standard Template Library Vector or a Hybrid object
///
/// \param inverse_transform  The inverse transformation matrix
/// \{
UnitCellType determineUnitCellTypeByShape(const double* inverse_transform);
UnitCellType determineUnitCellTypeByShape(const std::vector<double> &inverse_transform);
UnitCellType determineUnitCellTypeByShape(const Hybrid<double> &inverse_transform);
/// \}

/// \brief Interlace three arrays of X, Y, and Z coordinates (i.e. positions, velocities, or
///        forces) into a result ordered { X(0), Y(0), Z(0), X(1), Y(1), ..., Y(N), Z(N) }.
///
/// \param xptr        Pointer to Cartesian X data
/// \param yptr        Pointer to Cartesian Y data
/// \param zptr        Pointer to Cartesian Z data
/// \param low_index   Starting index in each array
/// \param high_index  Upper limit of data to take from each array
std::vector<double> interlaceXYZ(const double* xptr, const double* yptr, const double* zptr,
                                 int low_index, int high_index);
  
} // namespace trajectory
} // namespace stormm

#include "phasespace.tpp"

#endif
