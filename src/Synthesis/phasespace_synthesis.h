// -*-c++-*-
#ifndef STORMM_PHASESPACE_SYNTHESIS_H
#define STORMM_PHASESPACE_SYNTHESIS_H

#include "copyright.h"
#ifdef STORMM_USE_HPC
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

namespace stormm {
namespace synthesis {

#ifdef STORMM_USE_HPC
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
using trajectory::CoordinateCycle;
using trajectory::CoordinateFileKind;
using trajectory::PhaseSpace;
using trajectory::Thermostat;
using trajectory::ThermostatKind;
using trajectory::TrajectoryKind;
  
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
                    llint* xcrd_in, llint* ycrd_in, llint* zcrd_in, int* xcrd_ovrf_in,
                    int* ycrd_ovrf_in, int* zcrd_ovrf_in, llint* xvel_in, llint* yvel_in,
                    llint* zvel_in, int* xvel_ovrf_in, int* yvel_ovrf_in, int* zvel_ovrf_in,
                    llint* xfrc_in, llint* yfrc_in, llint* zfrc_in, int* xfrc_ovrf_in,
                    int* yfrc_ovrf_in, int* zfrc_ovrf_in, llint* xprv_in, llint* yprv_in,
                    llint* zprv_in, int* xprv_ovrf_in, int* yprv_ovrf_in, int* zprv_ovrf_in,
                    llint* xnxt_in, llint* ynxt_in,  llint* znxt_in, int* xnxt_ovrf_in,
                    int* ynxt_ovrf_in, int* znxt_ovrf_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisWriter(const PsSynthesisWriter &original) = default;
  PsSynthesisWriter(PsSynthesisWriter &&other) = default;
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
  llint* xcrd;     ///< Non-wrapped Cartesian X coordinates of all particles
  llint* ycrd;     ///< Non-wrapped Cartesian Y coordinates of all particles
  llint* zcrd;     ///< Non-wrapped Cartesian Z coordinates of all particles
  int* xcrd_ovrf;  ///< Non-wrapped Cartesian X coordinate overflow buffers
  int* ycrd_ovrf;  ///< Non-wrapped Cartesian Y coordinate overflow buffers
  int* zcrd_ovrf;  ///< Non-wrapped Cartesian Z coordinate overflow buffers
  llint* xvel;     ///< Cartesian X velocities
  llint* yvel;     ///< Cartesian Y velocities
  llint* zvel;     ///< Cartesian Z velocities
  int* xvel_ovrf;  ///< Cartesian X velocity overflow buffers
  int* yvel_ovrf;  ///< Cartesian Y velocity overflow buffers
  int* zvel_ovrf;  ///< Cartesian Z velocity overflow buffers
  llint* xfrc;     ///< Discretized Cartesian X forces
  llint* yfrc;     ///< Discretized Cartesian Y forces
  llint* zfrc;     ///< Discretized Cartesian Z forces
  int* xfrc_ovrf;  ///< Discretized Cartesian X force overflow buffers
  int* yfrc_ovrf;  ///< Discretized Cartesian Y force overflow buffers
  int* zfrc_ovrf;  ///< Discretized Cartesian Z force overflow buffers
  llint* xprv;     ///< Previous Cartesian X positions of particles
  llint* yprv;     ///< Previous Cartesian Y positions of particles
  llint* zprv;     ///< Previous Cartesian Z positions of particles
  int* xprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian X locations
  int* yprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian Y locations
  int* zprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian Z locations
  llint* xnxt;     ///< Future Cartesian X positions of particles
  llint* ynxt;     ///< Future Cartesian Y positions of particles
  llint* znxt;     ///< Future Cartesian Z positions of particles
  int* xnxt_ovrf;  ///< Overflow buffers for particles' future Cartesian X locations
  int* ynxt_ovrf;  ///< Overflow buffers for particles' future Cartesian Y locations
  int* znxt_ovrf;  ///< Overflow buffers for particles' future Cartesian Z locations
};

/// \brief The reader for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisReader {

  /// \brief The constructor can take lists of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.  Like the
  ///        CoordinateFrameReader, it can also convert its cognate writer by turning all of
  ///        the relevant pointers const.
  /// \{
  PsSynthesisReader(int system_count_in, UnitCellType unit_cell_in,
                    ThermostatKind heat_bath_kind_in, BarostatKind piston_kind_in,
                    double time_step_in, const int* atom_starts_in, const int* atom_counts_in,
                    double gpos_scale_in, double lpos_scale_in, double vel_scale_in,
                    double frc_scale_in, int gpos_bits_in, int lpos_bits_in, int vel_bits_in,
                    int frc_bits_in, const llint* boxvecs_in, const double* umat_in,
                    const double* invu_in, const double* boxdims_in, const float* sp_umat_in,
                    const float* sp_invu_in, const float* sp_boxdims_in,
                    const llint* xcrd_in, const llint* ycrd_in, const llint* zcrd_in,
                    const int* xcrd_ovrf_in, const int* ycrd_ovrf_in, const int* zcrd_ovrf_in,
                    const llint* xvel_in, const llint* yvel_in, const llint* zvel_in,
                    const int* xvel_ovrf_in, const int* yvel_ovrf_in, const int* zvel_ovrf_in,
                    const llint* xfrc_in, const llint* yfrc_in, const llint* zfrc_in,
                    const int* xfrc_ovrf_in, const int* yfrc_ovrf_in, const int* zfrc_ovrf_in,
                    const llint* xprv_in, const llint* yprv_in, const llint* zprv_in,
                    const int* xprv_ovrf_in, const int* yprv_ovrf_in, const int* zprv_ovrf_in,
                    const llint* xnxt_in, const llint* ynxt_in, const llint* znxt_in,
                    const int* xnxt_ovrf_in, const int* ynxt_ovrf_in, const int* znxt_ovrf_in);

  PsSynthesisReader(const PsSynthesisWriter &psyw);
  /// \}

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisReader(const PsSynthesisReader &original) = default;
  PsSynthesisReader(PsSynthesisReader &&other) = default;
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
  const llint* xcrd;     ///< Non-wrapped Cartesian X coordinates of all particles
  const llint* ycrd;     ///< Non-wrapped Cartesian Y coordinates of all particles
  const llint* zcrd;     ///< Non-wrapped Cartesian Z coordinates of all particles
  const int* xcrd_ovrf;  ///< Non-wrapped Cartesian X coordinate overflow buffers
  const int* ycrd_ovrf;  ///< Non-wrapped Cartesian Y coordinate overflow buffers
  const int* zcrd_ovrf;  ///< Non-wrapped Cartesian Z coordinate overflow buffers
  const llint* xvel;     ///< Cartesian X velocities
  const llint* yvel;     ///< Cartesian Y velocities
  const llint* zvel;     ///< Cartesian Z velocities
  const int* xvel_ovrf;  ///< Cartesian X velocity overflow buffers
  const int* yvel_ovrf;  ///< Cartesian Y velocity overflow buffers
  const int* zvel_ovrf;  ///< Cartesian Z velocity overflow buffers
  const llint* xfrc;     ///< Discretized Cartesian X forces
  const llint* yfrc;     ///< Discretized Cartesian Y forces
  const llint* zfrc;     ///< Discretized Cartesian Z forces
  const int* xfrc_ovrf;  ///< Discretized Cartesian X force overflow buffers
  const int* yfrc_ovrf;  ///< Discretized Cartesian Y force overflow buffers
  const int* zfrc_ovrf;  ///< Discretized Cartesian Z force overflow buffers
  const llint* xprv;     ///< Previous Cartesian X positions of particles
  const llint* yprv;     ///< Previous Cartesian Y positions of particles
  const llint* zprv;     ///< Previous Cartesian Z positions of particles
  const int* xprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian X locations
  const int* yprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian Y locations
  const int* zprv_ovrf;  ///< Overflow buffers for particles' previous Cartesian Z locations
  const llint* xnxt;     ///< Future Cartesian X positions of particles
  const llint* ynxt;     ///< Future Cartesian Y positions of particles
  const llint* znxt;     ///< Future Cartesian Z positions of particles
  const int* xnxt_ovrf;  ///< Overflow buffers for particles' future Cartesian X locations
  const int* ynxt_ovrf;  ///< Overflow buffers for particles' future Cartesian Y locations
  const int* znxt_ovrf;  ///< Overflow buffers for particles' future Cartesian Z locations
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
  ///   - Skip the specification of thermostats, barostats, and the time step (just use defaults),
  ///     but explicitly specify at least the global position scaling bit count and possibly other
  ///     bit counts
  ///
  /// \param ps_list       Array of input coordinates, velocities, and forces objects
  /// \param ag_list       Array of pointers to input topologies
  /// \param index_key     Indices of the given topology and coordinate objects to assemble into a
  ///                      larger list of systems to be held within the resulting
  ///                      PhaseSpaceSynthesis object.
  /// \param ps_index_key  Indices of the given coordinate objects to assemble into a larger list
  ///                      of systems to be held within the resulting PhaseSpaceSynthesis object.
  /// \param ag_index_key  Indices of the given topology objects to assemble into a larger list of
  ///                      systems to be held within the resulting PhaseSpaceSynthesis object.
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

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list, const std::vector<int> &index_key,
                      const std::vector<Thermostat> &heat_baths_in = { Thermostat() },
                      const std::vector<Barostat> &pistons_in = { Barostat() },
                      double time_step_in = 1.0,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list, const std::vector<int> &ps_index_key,
                      const std::vector<AtomGraph*> &ag_list, const std::vector<int> &ag_index_key,
                      const std::vector<Thermostat> &heat_baths_in = { Thermostat() },
                      const std::vector<Barostat> &pistons_in = { Barostat() },
                      double time_step_in = 1.0,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const SystemCache &sysc,
                      const std::vector<Thermostat> &heat_baths_in = { Thermostat() },
                      const std::vector<Barostat> &pistons_in = { Barostat() },
                      double time_step_in = 1.0,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list, int globalpos_scale_bits_in,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list, const std::vector<int> &index_key,
                      int globalpos_scale_bits_in,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list, const std::vector<int> &ps_index_key,
                      const std::vector<AtomGraph*> &ag_list, const std::vector<int> &ag_index_key,
                      int globalpos_scale_bits_in,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const SystemCache &sysc,
                      int globalpos_scale_bits_in,
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

  /// \brief Copy and move assignment operators are deleted.
  ///
  /// \param other  The other PhaseSpaceSynthesis object
  /// \{
  PhaseSpaceSynthesis& operator=(const PhaseSpaceSynthesis &other) = delete;
  PhaseSpaceSynthesis& operator=(PhaseSpaceSynthesis &&other) = delete;
  /// \}

  /// \brief Get the number of systems in the object.
  int getSystemCount() const;

  /// \brief Get the unit cell type.
  UnitCellType getUnitCellType() const;

  /// \brief Get the current position in the coordinate cycle.
  CoordinateCycle getCyclePosition() const;
  
  /// \brief Get the global position scaling bit count
  int getGlobalPositionBits() const;

  /// \brief Get the local position scaling bit count
  int getLocalPositionBits() const;

  /// \brief Get the velocity scaling bit count
  int getVelocityBits() const;

  /// \brief Get the force accumulation bit count
  int getForceAccumulationBits() const;
  
  /// \brief Get the reader or writer, as appropriate based on the const-ness of this object.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const object
  ///   - Get a writeable abstract for a mutable object
  ///   - Choose the stage in the time cycle that shall be deemed the relevant coordinates
  ///
  /// \param tier         The level (host or device) at which to get the set of pointers
  /// \param orientation  Stage in the object's time cycle to take as the current coordinates
  /// \{
  const PsSynthesisReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  const PsSynthesisReader data(CoordinateCycle orientation,
                               HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PsSynthesisWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  PsSynthesisWriter data(CoordinateCycle orientation,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Get a special writer which allows the device to read and write to host-mapped data.
  ///        This form of the writer can be used in kernel calls that streamline download and
  ///        upload of specific systems within the PhaseSpaceSynthesis with the launch latency
  ///        of only a single kernel call.  This offers a huge reduction in bandwidth requirements
  ///        and lower latency than multiple cudaMemcpy calls.
  ///
  /// Overloaded:
  ///   - Let the object automatically assign pointers to past, present, and future coordinates
  ///   - Specify which stage of the time cycle is to be taken as current
  ///
  /// \param orientation  Optional argument stipulating from which stage of the cycle to take the
  ///                     current coordinates, past coordinates, and pointers to any future
  ///                     coordinates.  If not specified the object's own cycle position will be
  ///                     used.
  /// \{
  PsSynthesisWriter deviceViewToHostData();
  PsSynthesisWriter deviceViewToHostData(CoordinateCycle orientation);
  /// \}
  
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
  void extractSystem(PhaseSpace *ps, int index,
                     HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Extract a specific type of coordinate from the synthesis into a pre-allocated space.
  ///
  /// \param ps        Pointer to an allocated PhaseSpace object (i.e. the original) ready to
  ///                  accept data from the synthesis (which may have evolved since it was first
  ///                  loaded)
  /// \param trajkind  Type of coordinates to copy
  /// \param tier      The level (host or device) at which to get the data
  /// \param index     Index of the system of interest within the synthesis
  void extractSystem(PhaseSpace *ps, int index,
                     TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                     HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Export a system's coordinates, velocities, and forces to a PhaseSpace object.
  ///
  /// \param index     Index of the system of interest within the synthesis
  /// \param tier      The level (host or device) at which to get the data
  PhaseSpace exportSystem(int index, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Export a system's coordinates, velocities, or forces to a compact CoordinateFrame
  ///        object.
  ///
  /// \param index     Index of the system of interest within the synthesis
  /// \param trajkind  Type of coordinates to copy
  /// \param tier      The level (host or device) at which to get the data
  CoordinateFrame exportCoordinates(int index, TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Initialize the forces within a PhaseSpaceSynthesis object.  This is the analog of the
  ///        eponymous function in the PhaseSpace object.
  ///
  /// Overloaded:
  ///   - Initialize one or more systems' forces on the host if STORMM is compiled for CPU only
  ///   - Initialize one or more systems' forces on either the host or the device if STORMM is
  ///     compiled for HPC
  ///
  /// \param gpu    Details of the GPU in use
  /// \param tier   The level (host or device) at which to initialize forces
  /// \param index  Index of the system of interest within the synthesis--if negative, all systems
  ///               will have their forces initialized
  /// \{
#ifdef STORMM_USE_HPC
  void initializeForces(const GpuDetails &gpu, HybridTargetLevel tier = HybridTargetLevel::HOST,
                        int index = -1);
#else
  void initializeForces(int index = -1);
#endif
  /// \}

  /// \brief Set the prior positions to the current forces, and initialize velocities to zero, as
  ///        part of the first step of conjugate gradient optimization.  This primes the system so
  ///        that the prior coordinates and velocities arrays can hold the prior forces and
  ///        temporary conjugate gradient memory.
  ///
  /// \param gpu    Details of the GPU in use
  /// \param tier   The level (host or device) at which to initialize vectors
#ifdef STORMM_USE_HPC
  void primeConjugateGradientCalculation(const GpuDetails &gpu,
                                         HybridTargetLevel tier = HybridTargetLevel::HOST);
#else
  void primeConjugateGradientCalculation();
#endif
  
  /// \brief Print a list of structures to a trajectory file.  Download will be performed when
  ///        calling this function, over the subset of relevant frames and data.
  ///
  /// \param system_indices  List of system coordinates / velocities / forces to print
  /// \param file_name       Name of the file to write, or base name of a set of files to write
  /// \param current_time    Current time progress of the group of simulations
  /// \param output_kind     The type of trajectory file to write
  /// \param expectation     The state that the output trajectory file is expected to be found in
  /// \param gpu             Specs of the GPU in use for the calculation (HPC mode only)
#ifdef STORMM_USE_HPC
  void printTrajectory(const std::vector<int> &system_indices, const std::string &file_name,
                       double current_time, CoordinateFileKind output_kind,
                       PrintSituation expectation, const GpuDetails &gpu);
#else
  void printTrajectory(const std::vector<int> &system_indices, const std::string &file_name,
                       double current_time, CoordinateFileKind output_kind,
                       PrintSituation expectation);
#endif
  
private:
  int system_count;               ///< The number of systems to tend at once
  UnitCellType unit_cell;         ///< The types of unit cells.  All unit cells must exist in
                                  ///<   isolated boundary conditions, or all cells must exist in
                                  ///<   periodic boundary conditions with a rectilinear or
                                  ///<   triclinic unit cell.
  CoordinateCycle cycle_position; ///< The stage of the time cycle (past >> present >> future) that
                                  ///<   holds the relevant, current coordinates
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

  // These variables are POINTER-kind Hybrid objects targeting the llint_data array
  Hybrid<llint> x_coordinates;         ///< Cartesian X coordinates of all particles
  Hybrid<llint> y_coordinates;         ///< Cartesian Y coordinates of all particles
  Hybrid<llint> z_coordinates;         ///< Cartesian Z coordinates of all particles
  Hybrid<int> x_coordinate_overflow;   ///< X direction coordinate overflow buffers
  Hybrid<int> y_coordinate_overflow;   ///< Y direction coordinate overflow buffers
  Hybrid<int> z_coordinate_overflow;   ///< Z direction coordinate overflow buffers
  Hybrid<llint> x_velocities;          ///< Cartesian X velocities of all particles
  Hybrid<llint> y_velocities;          ///< Cartesian Y velocities of all particles
  Hybrid<llint> z_velocities;          ///< Cartesian Z velocities of all particles
  Hybrid<int> x_velocity_overflow;     ///< Overflow buffers for velocities in the X direction
  Hybrid<int> y_velocity_overflow;     ///< Overflow buffers for velocities in the Y direction
  Hybrid<int> z_velocity_overflow;     ///< Overflow buffers for velocities in the Z direction
  Hybrid<llint> x_forces;              ///< Cartesian forces acting on particles in the X direction
  Hybrid<llint> y_forces;              ///< Cartesian forces acting on particles in the Y direction
  Hybrid<llint> z_forces;              ///< Cartesian forces acting on particles in the Z direction
  Hybrid<int> x_force_overflow;        ///< Force overflows acting on particles in the X direction
  Hybrid<int> y_force_overflow;        ///< Force overflows acting on particles in the Y direction
  Hybrid<int> z_force_overflow;        ///< Force overflows acting on particles in the Z direction
  Hybrid<llint> x_prior_coordinates;   ///< Prior Cartesian X coordinates of all particles
  Hybrid<llint> y_prior_coordinates;   ///< Prior Cartesian Y coordinates of all particles
  Hybrid<llint> z_prior_coordinates;   ///< Prior Cartesian Z coordinates of all particles
  Hybrid<int> x_prior_coord_overflow;  ///< Prior Cartesian X coordinates of all particles
  Hybrid<int> y_prior_coord_overflow;  ///< Prior Cartesian Y coordinates of all particles
  Hybrid<int> z_prior_coord_overflow;  ///< Prior Cartesian Z coordinates of all particles
  Hybrid<llint> x_future_coordinates;  ///< Prior Cartesian X coordinates of all particles
  Hybrid<llint> y_future_coordinates;  ///< Prior Cartesian Y coordinates of all particles
  Hybrid<llint> z_future_coordinates;  ///< Prior Cartesian Z coordinates of all particles
  Hybrid<int> x_future_coord_overflow; ///< Prior Cartesian X coordinates of all particles
  Hybrid<int> y_future_coord_overflow; ///< Prior Cartesian Y coordinates of all particles
  Hybrid<int> z_future_coord_overflow; ///< Prior Cartesian Z coordinates of all particles
  Hybrid<llint> box_vectors;           ///< Scaled real space transformation matrix--moving
                                       ///<   particles between images by adding or subtracting
                                       ///<   multiples of these vectors can be expeditious and
                                       ///<   keeps coordinate representations consistent between
                                       ///<   the lab frame and the primary unit cell.

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
  Hybrid<int> int_data;        ///< Counts of atoms and starting points for each system, plus
                               ///<   overflow data
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
  ///
  /// \param atom_stride  The total number of padded atoms in the system (sum over all individual
  ///                     systems with warp size padding in each of them)
  void allocate(size_t atom_stride);
};

} // namespace trajectory
} // namespace stormm

#endif
