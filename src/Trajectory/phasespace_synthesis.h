// -*-c++-*-
#ifndef OMNI_PHASESPACE_SYNTHESIS_H
#define OMNI_PHASESPACE_SYNTHESIS_H

#include "DataTypes/common_types.h"
#include "Accelerator/hybrid.h"
#include "Topology/atomgraph.h"
#include "barostat.h"
#include "phasespace.h"
#include "thermostat.h"

namespace omni {
namespace trajectory {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief The reader for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisReader {


  /// \brief Constructor takes a straight list of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.
  PsSynthesisReader(int system_count_in, UnitCellType unit_cell_in,
                    ThermostatKind heat_bath_kind_in, BarostatKind piston_kind_in,
                    double time_step_in, const int* atom_starts_in, const int* atom_counts_in,
                    const llint* boxvecs_in, const double* umat_in, const double* invu_in,
                    const double* boxdims_in, const float* sp_umat_in, const float* sp_invu_in,
                    const float* sp_boxdims_in, const longlong4* xyz_qlj_in, const double* xvel_in,
                    const double* yvel_in, const double* zvel_in, const llint* xfrc_in,
                    const llint* yfrc_in, const llint* zfrc_in);

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
  const double* xvel;        ///< Cartesian X velocities
  const double* yvel;        ///< Cartesian Y velocities
  const double* zvel;        ///< Cartesian Z velocities
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
                    const llint* boxvecs_in, const double* umat_in, const double* invu_in,
                    const double* boxdims_in, const float* sp_umat_in, const float* sp_invu_in,
                    const float* sp_boxdims_in, longlong4* xyz_qlj_in, double* xvel_in,
                    double* yvel_in, double* zvel_in, llint* xfrc_in, llint* yfrc_in,
                    llint* zfrc_in);

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
  longlong4* xyz_qlj;  ///< Non-wrapped coordinates
  double* xvel;        ///< Cartesian X velocities
  double* yvel;        ///< Cartesian Y velocities
  double* zvel;        ///< Cartesian Z velocities
  llint* xfrc;         ///< Discretized Cartesian X forces
  llint* yfrc;         ///< Discretized Cartesian Y forces
  llint* zfrc;         ///< Discretized Cartesian Z forces
};

/// \brief A fixed-precision representation of coordinates, velocities, and forces to manage a set
///        of simulations.  The time steps are stored in units of femtoseconds.  Coordinates are
///        stored as long long integers to utilize the 32-bit int pipeline for difference
///        computations, bypassing fp64 computations wherever possible to get the best throughput
///        on many visualization, AI-oriented, and gaming cards.
struct PhaseSpaceSynthesis {

  /// \brief The constructor works from a series PhaseSpace object, importing  The LabFrame is not
  ///        a POINTER-kind object of any sort.
  ///
  /// \param ps            Input coordinates, velocities, and forces object
  /// \param time_step_in  The initial time step (will be converted to integer representation as
  ///                      described above, default 1.0fs which converts to llint(8))
  /// \param heat_bath_in  A thermostat to govern integration of the equations of motion
  /// \param piston_in     A barostat to govern box rescaling (if there is a unit cell)
  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list, double time_step_in,
                      const std::vector<AtomGraph*> &ag_list,
                      const std::vector<Thermostat> &heat_baths_in,
                      const std::vector<Barostat> &pistons_in);

#ifdef OMNI_USE_HPC
  /// \brief Upload data to the device
  void upload();
  
  /// \brief Download data from the device
  void download();
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

  /// \brief Get the reader or writer, as appropriate based on the const-ness of this object.
  ///
  /// \param tier  The level (host or device) at which to get the set of pointers
  /// \{
  const PsSynthesisReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PsSynthesisWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

private:
  int system_count;              ///< The number of systems to tend at once
  UnitCellType unit_cell;        ///< The types of unit cells.  All unit cells must exist in
                                 ///<   isolated boundary conditions, or all cells must exist in
                                 ///<   periodic boundary conditions with a rectilinear or
                                 ///<   triclinic unit cell.
  ThermostatKind heat_bath_kind; ///< The type of thermostat that all systems will run.  Other
                                 ///<   details are variable between systems and stored in their
                                 ///<   own arrays.
  BarostatKind piston_kind;      ///< The pressure scaling method for all systems if the unit cells
                                 ///<   are periodic.  Additional details may vary between systems
                                 ///<   and are held in their own arrays.
  double time_step;              ///< The time step for every system, in femtoseconds

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
  Hybrid<double> x_velocities;  ///< Cartesian velocities of all particles in the X direction
  Hybrid<double> y_velocities;  ///< Cartesian velocities of all particles in the Y direction
  Hybrid<double> z_velocities;  ///< Cartesian velocities of all particles in the Z direction
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
};

} // namespace trajectory
} // namespace omni

#endif
