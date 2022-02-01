// -*-c++-*-
#ifndef OMNI_PHASE_SPACE_H
#define OMNI_PHASE_SPACE_H

#include "Accelerator/hybrid.h"
#include "Topology/atomgraph.h"
#include "Trajectory/trajectory_enumerators.h"

namespace omni {
namespace trajectory {

using cuda::Hybrid;
using cuda::HybridTargetLevel;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateFileKind;

/// \brief Collect constants and pointers to the components of a read-only PhaseSpace object.
struct PhaseSpaceReader {

  /// Constructor takes a large list of arguments passed in from the original PhaseSpace object
  PhaseSpaceReader(const int natom_in, const UnitCellType unit_cell_in, const double* xcrd_in,
                   const double* ycrd_in, const double* zcrd_in, const double* umat_in,
                   const double* invu_in, const double* boxdim_in, const double* xvel_in,
                   const double* yvel_in, const double* zvel_in, const double* xfrc_in,
                   const double* yfrc_in, const double* zfrc_in, const double* xprv_in,
                   const double* yprv_in, const double* zprv_in);

  const int natom;                ///< Atom count for this system (still a constant)
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
};

/// \brief Collect constants and pointers to the components of a modifiable PhaseSpace object.
struct PhaseSpaceWriter {

  /// Constructor takes a large list of arguments passed in from the original PhaseSpace object
  PhaseSpaceWriter(const int natom_in, const UnitCellType unit_cell_in, double* xcrd_in,
                   double* ycrd_in, double* zcrd_in, double* umat_in, double* invu_in,
                   double* boxdim_in, double* xvel_in, double* yvel_in, double* zvel_in,
                   double* xfrc_in, double* yfrc_in, double* zfrc_in, double* xprv_in,
                   double* yprv_in, double* zprv_in);

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
};

/// \brief An object to complement a topology and hold positions, velocities, and forces of all
///        particles in a system.  This is not designed to be the most performant representation of
///        the system's structure.  Rather, it serves to hold a high-precision representation of a
///        single system and transport it between CPUs and high-performance accelerators.
struct PhaseSpace {

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

  PhaseSpace(const std::string &file_name_in, CoordinateFileKind file_kind,
             int frame_number = 0);

  PhaseSpace(const std::string &file_name_in, CoordinateFileKind file_kind, const AtomGraph &ag, 
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

  /// \brief Take the default move assignment operator
  ///  
  /// \param other     Another way to say original, in a different semantic context
  PhaseSpace& operator=(PhaseSpace &&other);

  /// \brief Fill the object from information some coordinate, restart, or trajectory file.
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param frame_number  The frame number to read (if the file is a trajectory, not a single
  ///                      point from the system's phase space)
  void buildFromFile(const std::string &file_name_in, const CoordinateFileKind file_kind,
                     int frame_number = 0);

  /// \brief Get the name of the file associated with this object.
  std::string getFileName() const;

  /// \brief Get the number of atoms (particles, including virtual sites) in the object.
  int getAtomCount() const;

  /// \brief Get the unit cell type of the coordinate system
  UnitCellType getUnitCellType() const;

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

  /// \brief Initialize the forces (set them to zero)
  void initializeForces();

  /// \brief Get the transformation matrix to take coordinates into box (fractional) space.  The
  ///        result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if OMNI is compiled to run on a GPU)
  std::vector<double> getBoxSpaceTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the transformation matrix to take coordinates from fractional space back into
  ///        real space.  The result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if OMNI is compiled to run on a GPU)
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

#ifdef OMNI_USE_HPC
  /// \brief Upload all information
  void upload();

  /// \brief Upload current positional information
  void uploadPositions();

  /// \brief Upload current transformation matrices
  void uploadTransformations();

  /// \brief Upload velocity information
  void uploadVelocities();

  /// \brief Upload force information
  void uploadForces();

  /// \brief Upload force information
  void uploadNonbondedForces();

  /// \brief Upload previous positional information
  void uploadPriorPositions();

  /// \brief Download all information
  void download();

  /// \brief Download current positional information
  void downloadPositions();

  /// \brief Download current transformation matrices
  void downloadTransformations();

  /// \brief Download velocity information
  void downloadVelocities();

  /// \brief Download force information
  void downloadForces();

  /// \brief Download force information
  void downloadNonbondedForces();

  /// \brief Download previous positional information
  void downloadPriorPositions();
#endif
  /// \brief Get the abstract for this object, containing C-style pointers for the most rapid
  ///        access to any of its member variables.
  ///
  /// Overloaded:
  ///   - Get a read-only object
  ///   - Get a writeable object
  /// \{
  const PhaseSpaceReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PhaseSpaceWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

private:
  std::string file_name;              ///< Name of the file from which these coordinates (and
                                      ///<   perhaps velocities) derived.  Empty string indicates
                                      ///<   no file.
  int atom_count;                     ///< The number of atoms in the system
  UnitCellType unit_cell;             ///< The type of unit cell
  Hybrid<double> x_coordinates;       ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;       ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;       ///< Cartesian Z coordinates of all particles
  Hybrid<double> box_space_transform; ///< Matrix to transform coordinates into box space (3 x 3)
  Hybrid<double> inverse_transform;   ///< Matrix to transform coordinates into real space (3 x 3)
  Hybrid<double> box_dimensions;      ///< Three lengths and three angles defining the box (lengths
                                      ///<   are given in Angstroms, angles in radians)
  Hybrid<double> x_velocities;        ///< Cartesian X velocities of all particles
  Hybrid<double> y_velocities;        ///< Cartesian Y velocities of all particles
  Hybrid<double> z_velocities;        ///< Cartesian Z velocities of all particles
  Hybrid<double> x_forces;            ///< Cartesian X forces acting on all particles
  Hybrid<double> y_forces;            ///< Cartesian Y forces acting on all particles
  Hybrid<double> z_forces;            ///< Cartesian Z forces acting on all particles 
  Hybrid<double> x_prior_coordinates; ///< Previous step Cartesian X coordinates of all particles
  Hybrid<double> y_prior_coordinates; ///< Previous step Cartesian Y coordinates of all particles
  Hybrid<double> z_prior_coordinates; ///< Previous step Cartesian Z coordinates of all particles

  /// All of the above Hybrid objects are pointers into this single large array, segmented to hold
  /// each type of information with zero-padding to accommodate the HPC warp size.
  Hybrid<double> storage;

  /// \brief Allocate space for the object, based on a known number of atoms
  void allocate();
};

/// \brief Interpret the inverse transformation matrix to determine whether a unit cell exists, and
///        if it is rectilinear or not.
UnitCellType determineUnitCellTypeByShape(const double* inv_ptr);

} // namespace trajectory
} // namespace omni

#endif
