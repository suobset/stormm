// -*-c++-*-
#ifndef OMNI_PHASESPACE_SYNTHESIS_H
#define OMNI_PHASESPACE_SYNTHESIS_H

#ifdef OMNI_USE_HPC
#include "Accelerator/gpu_details.h"
#endif
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_util.h"
#include "Topology/atomgraph.h"
#include "Trajectory/barostat.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trajectory_enumerators.h"
#include "systemcache.h"

namespace omni {
namespace synthesis {

#ifdef OMNI_USE_HPC
using card::GpuDetails;
#endif
using card::Hybrid;
using card::HybridTargetLevel;
using diskutil::PrintSituation;
using numerics::default_globalpos_scale_bits;
using numerics::default_localpos_scale_bits;
using numerics::default_velocity_scale_bits;
using numerics::default_force_scale_bits;
using topology::UnitCellType;
using trajectory::Barostat;
using trajectory::BarostatKind;
using trajectory::CoordinateFileKind;
using trajectory::PhaseSpace;
using trajectory::Thermostat;
using trajectory::ThermostatKind;
using trajectory::TrajectoryKind;
  
/// \brief The reader for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisReader {

  /// \brief Constructor takes a straight list of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.
  PsSynthesisReader(int system_count_in, UnitCellType unit_cell_in,
                    ThermostatKind heat_bath_kind_in, BarostatKind piston_kind_in,
                    double time_step_in, const int* atom_starts_in, const int* atom_counts_in,
                    double gpos_scale_in, double lpos_scale_in, double vel_scale_in,
                    double frc_scale_in, int gpos_bits_in, int lpos_bits_in, int vel_bits_in,
                    int frc_bits_in, const llint* boxvecs_in, const double* umat_in,
                    const double* invu_in, const double* boxdims_in, const float* sp_umat_in,
                    const float* sp_invu_in, const float* sp_boxdims_in,
                    const longlong4* xyz_qlj_in, const llint* xvel_in, const llint* yvel_in,
                    const llint* zvel_in, const llint* xfrc_in, const llint* yfrc_in,
                    const llint* zfrc_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisReader(const PsSynthesisReader &original) = default;
  PsSynthesisReader(PsSynthesisReader &&other) = default;
  PsSynthesisReader& operator=(const PsSynthesisReader &original) = default;
  /// \}
  
  // System sizing information
  const int system_count;               ///< The number of independent systems
  const UnitCellType unit_cell;         ///< Type of unit cells (or none) each system resides in
  const ThermostatKind heat_bath_kind;  ///< The type of thermostat used throughout all systems
  const BarostatKind piston_kind;       ///< The type of barostat used throughout all systems
  const double time_step;               ///< Time step used in all simulations
  const int* atom_starts;               ///< Points at which each system starts in the atom list
  const int* atom_counts;               ///< Atom counts for all systems

  // Scaling factors: the PhaseSpaceSynthesis permits a customizable discretization of fixed-point
  // arithmetic.
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double lpos_scale;       ///< Local position coordinate scaling factor
  const double vel_scale;        ///< Velocity coordinate scaling factor
  const double frc_scale;        ///< Scaling factor for fixed-precision force accumulation
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  const double inv_lpos_scale;   ///< Inverse local coordinate scaling factor
  const double inv_vel_scale;    ///< Inverse velocity scaling factor
  const double inv_frc_scale;    ///< Inverse force scaling factor
  const float gpos_scale_f;      ///< Global position coordinate scaling factor
  const float lpos_scale_f;      ///< Local position coordinate scaling factor
  const float vel_scale_f;       ///< Velocity coordinate scaling factor
  const float frc_scale_f;       ///< Scaling factor for fixed-precision force accumulation
  const float inv_gpos_scale_f;  ///< Inverse global coordinate scaling factor
  const float inv_lpos_scale_f;  ///< Inverse local coordinate scaling factor
  const float inv_vel_scale_f;   ///< Inverse velocity scaling factor
  const float inv_frc_scale_f;   ///< Inverse force scaling factor
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const int lpos_bits;           ///< Local position coordinate bits after the decimal
  const int vel_bits;            ///< Velocity coordinate bits after the decimal
  const int frc_bits;            ///< Force component bits after the decimal
  
  // Pointers to the transformations and box vectors are likewise const--once created, this
  // object is valid for a system held in constant volume.
  const llint* boxvecs;     ///< Discretized box vectors
  const double* umat;       ///< Box (fractional) space transformation matrices, one per warp
  const double* invu;       ///< Inverse transformation matrices, one per warp
  const double* boxdims;    ///< Box dimensions (a, b, c, alpha, beta, gamma)
  const float* sp_umat;     ///< Single precision fractional coordinate transformation matrices
  const float* sp_invu;     ///< Single precision inverse transformation matrices
  const float* sp_boxdims;  ///< Single precision box dimensions

  // Pointers to the coordinate, velocity, and force data--these are mutable for accumulating
  // forces and letting a trajectory evolve.
  const longlong4* xyz_qlj;  ///< Non-wrapped coordinates with non-bonded property indices,
                             ///<   including non-bonded exclusions
  const llint* xvel;         ///< Cartesian X velocities
  const llint* yvel;         ///< Cartesian Y velocities
  const llint* zvel;         ///< Cartesian Z velocities
  const llint* xfrc;         ///< Discretized Cartesian X forces
  const llint* yfrc;         ///< Discretized Cartesian Y forces
  const llint* zfrc;         ///< Discretized Cartesian Z forces
};

/// \brief The writer for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisWriter {

  /// \brief Constructor takes a straight list of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.
  PsSynthesisWriter(int system_count_in, UnitCellType unit_cell_in,
                    ThermostatKind heat_bath_kind_in, BarostatKind piston_kind_in,
                    double time_step_in, const int* atom_starts_in, const int* atom_counts_in,
                    double gpos_scale_in, double lpos_scale_in, double vel_scale_in,
                    double frc_scale_in, int gpos_bits_in, int lpos_bits_in, int vel_bits_in,
                    int frc_bits_in, llint* boxvecs_in, double* umat_in, double* invu_in,
                    double* boxdims_in, float* sp_umat_in, float* sp_invu_in, float* sp_boxdims_in,
                    longlong4* xyz_qlj_in, llint* xvel_in, llint* yvel_in, llint* zvel_in,
                    llint* xfrc_in, llint* yfrc_in, llint* zfrc_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisWriter(const PsSynthesisWriter &original) = default;
  PsSynthesisWriter(PsSynthesisWriter &&other) = default;
  PsSynthesisWriter& operator=(const PsSynthesisWriter &original) = default;
  /// \}
  
  // Even in the writer, some information may not be altered.  This includes the simulation
  // conditions, the number of systems, and the number of atoms in each system.  The layout
  // of the data must not be corrupted.
  const int system_count;               ///< The number of independent systems
  const UnitCellType unit_cell;         ///< Type of unit cells (or none) each system resides in
  const ThermostatKind heat_bath_kind;  ///< The type of thermostat used throughout all systems
  const BarostatKind piston_kind;       ///< The type of barostat used throughout all systems
  const double time_step;               ///< Time step used in all simulations
  const int* atom_starts;               ///< Points at which each system starts in the atom list
  const int* atom_counts;               ///< Atom counts for all systems

  // Scaling factors: the PhaseSpaceSynthesis permits a customizable discretization of fixed-point
  // arithmetic.
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double lpos_scale;       ///< Local position coordinate scaling factor
  const double vel_scale;        ///< Velocity coordinate scaling factor
  const double frc_scale;        ///< Scaling factor for fixed-precision force accumulation
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  const double inv_lpos_scale;   ///< Inverse local coordinate scaling factor
  const double inv_vel_scale;    ///< Inverse velocity scaling factor
  const double inv_frc_scale;    ///< Inverse force scaling factor
  const float gpos_scale_f;      ///< Global position coordinate scaling factor
  const float lpos_scale_f;      ///< Local position coordinate scaling factor
  const float vel_scale_f;       ///< Velocity coordinate scaling factor
  const float frc_scale_f;       ///< Scaling factor for fixed-precision force accumulation
  const float inv_gpos_scale_f;  ///< Inverse global coordinate scaling factor
  const float inv_lpos_scale_f;  ///< Inverse local coordinate scaling factor
  const float inv_vel_scale_f;   ///< Inverse velocity scaling factor
  const float inv_frc_scale_f;   ///< Inverse force scaling factor
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const int lpos_bits;           ///< Local position coordinate bits after the decimal
  const int vel_bits;            ///< Velocity coordinate bits after the decimal
  const int frc_bits;            ///< Force component bits after the decimal

  // Pointers to the transformations and box vectors are mutable if the systems change volume.
  llint* boxvecs;     ///< Discretized box vectors
  double* umat;       ///< Box (fractional) space transformation matrices, one per warp
  double* invu;       ///< Inverse transformation matrices, one per warp
  double* boxdims;    ///< Box dimensions (a, b, c, alpha, beta, gamma)
  float* sp_umat;     ///< Single precision fractional coordinate transformation matrices
  float* sp_invu;     ///< Single precision inverse transformation matrices
  float* sp_boxdims;  ///< Single precision box dimensions

  // Pointers to the coordinate, velocity, and force data--these are mutable for accumulating
  // forces and letting a trajectory evolve.
  longlong4* xyz_qlj;  ///< Non-wrapped coordinates
  llint* xvel;         ///< Cartesian X velocities
  llint* yvel;         ///< Cartesian Y velocities
  llint* zvel;         ///< Cartesian Z velocities
  llint* xfrc;         ///< Discretized Cartesian X forces
  llint* yfrc;         ///< Discretized Cartesian Y forces
  llint* zfrc;         ///< Discretized Cartesian Z forces
};

/// \brief A fixed-precision representation of coordinates, velocities, and forces to manage a set
///        of simulations.  The time steps are stored in units of femtoseconds.  Coordinates are
///        stored as long long integers to utilize the 32-bit int pipeline for difference
///        computations, bypassing fp64 computations wherever possible to get the best throughput
///        on many visualization, AI-oriented, and gaming cards.
class PhaseSpaceSynthesis {
public:

  /// \brief The constructor works from a series PhaseSpace object, importing  The LabFrame is not
  ///        a POINTER-kind object of any sort.
  ///
  /// Overloaded:
  ///   - Take arrays of PhaseSpace objects and AtomGraph pointers
  ///   - Take a SystemCache object and unpack it
  ///
  /// \param ps            Input coordinates, velocities, and forces object
  /// \param time_step_in  The initial time step (will be converted to integer representation as
  ///                      described above, default 1.0fs which converts to llint(8))
  /// \param heat_bath_in  A thermostat to govern integration of the equations of motion
  /// \param piston_in     A barostat to govern box rescaling (if there is a unit cell)
  /// \{
  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list,
                      const std::vector<Thermostat> &heat_baths_in = { Thermostat() },
                      const std::vector<Barostat> &pistons_in = { Barostat() },
                      double time_step_in = 1.0,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const SystemCache &sysc,
                      const std::vector<Thermostat> &heat_baths_in,
                      const std::vector<Barostat> &pistons_in, double time_step_in = 1.0,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);
  /// \}

  /// \brief Copy and move constructors work much like their counterparts in the smaller
  ///        PhaseSpace object.
  ///
  /// \param original  The original PhaseSpaceSynthesis object
  /// \{
  PhaseSpaceSynthesis(const PhaseSpaceSynthesis &original);
  PhaseSpaceSynthesis(PhaseSpaceSynthesis &&original);
  /// \}

  /// \brief Copy and move assignment operators work much like their counterparts in the smaller
  ///        PhaseSpace object.
  ///
  /// \param other  The other PhaseSpaceSynthesis object
  /// \{
  PhaseSpaceSynthesis& operator=(const PhaseSpaceSynthesis &other) = delete;
  PhaseSpaceSynthesis& operator=(PhaseSpaceSynthesis &&other) = delete;
  /// \}
  
  /// \brief Get the reader or writer, as appropriate based on the const-ness of this object.
  ///
  /// \param tier  The level (host or device) at which to get the set of pointers
  /// \{
  const PsSynthesisReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PsSynthesisWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

#ifdef OMNI_USE_HPC
  /// \brief Get a special writer which allows the device to read and write to host-mapped data.
  ///        This form of the writer can be used in kernel calls that streamline download and
  ///        upload of specific systems within the PhaseSpaceSynthesis with the launch latency
  ///        of only a single kernel call.  This offers a huge reduction in bandwidth requirements
  ///        and lower latency than multiple cudaMemcpy calls.
  PsSynthesisWriter deviceViewToHostData();
  
  /// \brief Upload data to the device.
  ///
  /// Overloaded:
  ///   - Upload all data
  ///   - Upload data associated with a specific type of coordinates
  ///   - Upload data assoicated with a specific type of coordinates and specific systems
  ///
  /// \param kind                Choose from POSITIONS, VELOCITIES, FORCES (positional data will
  ///                            include box dimensions and is the only way to get box dimensions)
  /// \param system_index        The specific system to upload
  /// \param system_lower_bound  Upload systems in the range [ lower_bound, upper_bound )
  /// \param system_upper_bound  Upload systems in the range [ lower_bound, upper_bound )
  /// \{
  void upload();
  void upload(TrajectoryKind kind);
  void upload(TrajectoryKind kind, int system_index, const GpuDetails &gpu);
  void upload(TrajectoryKind kind, int system_lower_bound, int system_upper_bound,
              const GpuDetails &gpu);
  /// \}
  
  /// \brief Download data from the device
  ///
  /// Overloaded:
  ///   - Download all data
  ///   - Download data associated with a specific type of coordinates
  ///   - Download data assoicated with a specific type of coordinates and specific systems
  ///
  /// \param kind                Choose from POSITIONS, VELOCITIES, FORCES (positional data will
  ///                            include box dimensions and is the only way to get box dimensions)
  /// \param system_index        The specific system to download
  /// \param system_lower_bound  Download systems in the range [ lower_bound, upper_bound )
  /// \param system_upper_bound  Download systems in the range [ lower_bound, upper_bound )
  /// \{
  void download();
  void download(TrajectoryKind kind);
  void download(TrajectoryKind kind, int system_index, const GpuDetails &gpu);
  void download(TrajectoryKind kind, int system_lower_bound, int system_upper_bound,
                const GpuDetails &gpu);
  /// \}
#endif
  
  /// \brief Move particles according to a given set of forces, based on known masses, to complete
  ///        the coordinate update phase of Velocity-Verlet integration.  Also updates velocities
  ///        by half a time step.
  ///
  /// \param ag  System topology (for atomic masses)
  void velocityVerletCoordinateUpdate();

  /// \brief Update velocities by half a time step based on the newly computed forces in the
  ///        velocity Verlet integration scheme.
  ///
  /// \param ag  System topology (for atomic masses)
  void velocityVerletVelocityUpdate();

  /// \brief Assign velocities by an Andersen thermostat move
  ///
  /// \param ag  System topology (for atomic masses)
  void assignMaxwellVelocities();

  /// \brief Assign velocities by Berendsen coupling
  ///
  /// \param ag  System topology (for atomic masses)
  void berendsenThermocoupling();

  /// \brief Extract the phase space (plus forces) of a specific system within the synthesis.
  ///
  /// \param ps     Pointer to an allocated PhaseSpace object (i.e. the original) ready to accept
  ///               data from the synthesis (which may have evolved since it was first loaded)
  /// \param tier  The level (host or device) at which to get the data
  /// \param index  Index of the system of interest within the synthesis
  void extractPhaseSpace(PhaseSpace *ps, int index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Extract a specific type of coordinate from the synthesis into a pre-allocated space.
  ///
  /// \param ps        Pointer to an allocated PhaseSpace object (i.e. the original) ready to
  ///                  accept data from the synthesis (which may have evolved since it was first
  ///                  loaded)
  /// \param trajkind  Type of trajectory to copy
  /// \param tier      The level (host or device) at which to get the data
  /// \param index     Index of the system of interest within the synthesis
  void extractCoordinates(PhaseSpace *ps, int index,
                          TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                          HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Print a list of structures to a trajectory file.  Download will be performed when
  ///        calling this function, over the subset of relevant frames and data.
  ///
  /// \param system_indices  List of system coordinates / velocities / forces to print
  /// \param file_name       Name of the file to write, or base name of a set of files to write
  /// \param output_kind     The type of trajectory file to write
  /// \param expectation     The state that the output trajectory file is expected to be found in
  /// \param gpu             Specs of the GPU in use for the calculation (HPC mode only)
#ifdef OMNI_USE_HPC
  void printTrajectory(const std::vector<int> &system_indices, const std::string &file_name,
                       CoordinateFileKind output_kind, PrintSituation expectation,
                       const GpuDetails &gpu);
#else
  void printTrajectory(const std::vector<int> &system_indices, const std::string &file_name,
                       CoordinateFileKind output_kind, PrintSituation expectation);
#endif
  
private:
  int system_count;               ///< The number of systems to tend at once
  UnitCellType unit_cell;         ///< The types of unit cells.  All unit cells must exist in
                                  ///<   isolated boundary conditions, or all cells must exist in
                                  ///<   periodic boundary conditions with a rectilinear or
                                  ///<   triclinic unit cell.
  ThermostatKind heat_bath_kind;  ///< The type of thermostat that all systems will run.  Other
                                  ///<   details are variable between systems and stored in their
                                  ///<   own arrays.
  BarostatKind piston_kind;       ///< The pressure scaling method for all systems if the unit
                                  ///<   cells are periodic.  Additional details may vary between
                                  ///<   systems and are held in their own arrays.
  double time_step;               ///< The time step for every system, in femtoseconds

  // Scaling constants for fixed-precision coordinates (position, velocity) and also forces
  // in this PhaseSpaceSynthesis are fixed upon creation
  const double globalpos_scale;         ///< Global position coordinate scaling factor
  const double localpos_scale;          ///< Local position coordinate scaling factor
  const double velocity_scale;          ///< Velocity coordinate scaling factor
  const double force_scale;             ///< Scaling factor for fixed-precision force accumulation
  const double inverse_globalpos_scale; ///< Inverse global coordinate scaling factor
  const double inverse_localpos_scale;  ///< Inverse local coordinate scaling factor
  const double inverse_velocity_scale;  ///< Inverse velocity scaling factor
  const double inverse_force_scale;     ///< Inverse force scaling factor
  const int globalpos_scale_bits;       ///< Global position coordinate bits after the decimal
  const int localpos_scale_bits;        ///< Local position coordinate bits after the decimal
  const int velocity_scale_bits;        ///< Velocity coordinate bits after the decimal
  const int force_scale_bits;           ///< Force component bits after the decimal

  /// Starting positions for each system's stretch of atoms in xyz_qlj, (x,y,z)_velocities, and
  /// (x,y,z)_forces.  Atoms in each of those arrays will remain in their original orders, as
  /// presented in their respective topologies.
  Hybrid<int> atom_starts;

  ///< The numbers of atoms in each system
  Hybrid<int> atom_counts;

  /// Non-wrapped Cartesian X, Y, and Z coordinates of the system, fused with non-bonded properties
  /// for streamlined access in the majority of situations where they will be needed
  /// simultaneously.  This is its own array.
  Hybrid<longlong4> xyz_qlj;

  // These variables are POINTER-kind Hybrid objects targeting the llint_data array
  Hybrid<llint> x_velocities;   ///< Cartesian velocities of all particles in the X direction
  Hybrid<llint> y_velocities;   ///< Cartesian velocities of all particles in the Y direction
  Hybrid<llint> z_velocities;   ///< Cartesian velocities of all particles in the Z direction
  Hybrid<llint> x_forces;       ///< Cartesian forces acting on all particles in the X direction
  Hybrid<llint> y_forces;       ///< Cartesian forces acting on all particles in the Y direction
  Hybrid<llint> z_forces;       ///< Cartesian forces acting on all particles in the Z direction
  Hybrid<llint> box_vectors;    ///< Scaled real space transformation matrix--moving particles
                                ///<   between images by adding or subtracting multiples of these
                                ///<   vectors can be expeditious and keeps coordinate
                                ///<   representations consistent between the lab frame and the
                                ///<   primary unit cell.

  // The following are POINTER-kind Hybrid objects targeting the floating point data arrays
  Hybrid<double> box_space_transforms;   ///< Transformation matrices to take coordinates into
                                         ///<   box (fractional) space
  Hybrid<double> inverse_transforms;     ///< Transformation matrice back to go back to real space
  Hybrid<double> box_dimensions;         ///< Three box lengths and the angles in the planes normal
                                         ///<   to each axis in each system (synchronized with the
                                         ///<   transformation matrices--each update of the
                                         ///<   box_dimensions triggers an update of the
                                         ///<   transformations)
  Hybrid<float> sp_box_space_transforms; ///< Single precision fractional coordinate transforms
  Hybrid<float> sp_inverse_transforms;   ///< Single precision real space transforms
  Hybrid<float> sp_box_dimensions;       ///< Single precision box dimensions

  // Data arrays
  Hybrid<int> int_data;        ///< Counts of atoms and starting points for each system 
  Hybrid<llint> llint_data;    ///< The discretized data for all of phase space and forces
  Hybrid<double> double_data;  ///< Double-precision floating point transformations--these are the
                               ///<   standard for moving coordinates into a re-imaged
                               ///<   configuration for later computations.  Once particles are
                               ///<   in a re-imaged configuration for a given box size, they can
                               ///<   be manipulated in parallel with the coordinates in the lab
                               ///<   frame.  The unit cell transformations are computed in double
                               ///<   precision, then used to construct the long long int
                               ///<   box_vectors as scaled representation, then synchronized to
                               ///<   the box_vectors to make consistent representations for
                               ///<   calculating excursions from the unit cell boundaries.
  Hybrid<float> float_data;    ///< Single-precision floating point transformations--cannot be
                               ///<   trusted for actual coordinates, but as a guide on whether to
                               ///<   image a particle or move it between honeycomb pencils these
                               ///<   can be useful.

  // Pointer to the topology that describes this system
  const std::vector<AtomGraph*> topologies;

  /// \brief Allocate private array data
  void allocate(int atom_stride);
};

} // namespace trajectory
} // namespace omni

#endif
