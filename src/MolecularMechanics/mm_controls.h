// -*-c++-*-
#ifndef OMNI_MM_CONTROLS_H
#define OMNI_MM_CONTROLS_H

#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Namelists/nml_dynamics.h"
#include "Synthesis/atomgraph_synthesis.h"

namespace omni {
namespace mm {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using namelist::default_dynamics_time_step;
using namelist::default_rattle_tolerance;
using synthesis::AtomGraphSynthesis;
  
/// \brief The C-style, always writeable abstract for the MolecularMechanicsControls object.  To
///        not be able to modify this object's contents would be nonsensical, as it is intended to
///        to keep counters of the simulation time step as well as force evaluation work units.
template <typename T> struct MMControlKit {

  /// \brief The constructor takes a straight list of values and pointers.  The step number is
  ///        left modifiable so that the object can be re-used over successive time steps.
  MMControlKit(int step_in, T dt_in, T rattle_tol_in, int* vwu_progress_in, int* nbwu_progress_in,
               int* pmewu_progress_in);

  /// \brief The usual copy and move constructors for an abstract apply here.  
  /// \{
  MMControlKit(const MMControlKit &original) = default;
  MMControlKit(MMControlKit &&original) = default;
  /// \}

  int step;             ///< The current simulation step
  const T dt;           ///< Simulation time step (only valid for dynamics)
  const T rattle_tol;   ///< Convergence tolerance for bond constraints
  int* vwu_progress;    ///< Progress counters for valence work units
  int* nbwu_progress;   ///< Progress counters for 
  int* pmewu_progress;  ///<
};

/// \brief A collection of contol data for molecular mechanics simulations, conveying the current
///        step number, progress counters through various work units, the time step, RATTLE
///        tolerance, and other critical parameters to guide calculations.
class MolecularMechanicsControls {
public:

  /// \brief The constructor can create an empty object with default parameters for the time step
  ///        and rattle tolerance, or accept user-specified values.
  MolecularMechanicsControls(double time_step_in = default_dynamics_time_step,
                             double rattle_tol_in = default_rattle_tolerance);

  /// \brief The copy constructor handles assignment of internal POINTER-kind Hybrid objects.
  ///
  /// \param original  The object to copy
  MolecularMechanicsControls(const MolecularMechanicsControls &original);

  /// \brief The copy assignment operator likewise handles assignment of internal POINTER-kind
  ///        Hybrid objects.
  ///
  /// \param other  Another way to say original, in a different semantic context
  MolecularMechanicsControls& operator=(const MolecularMechanicsControls &other);

  /// \brief The move constructor prepared the original object for destruction.
  ///
  /// \param original  The object from which to take contents
  MolecularMechanicsControls(MolecularMechanicsControls &&original);

  /// \brief The move assignment operator looks much like the copy assignment operator.
  ///
  /// \param other  The object from which to take contents
  MolecularMechanicsControls& operator=(MolecularMechanicsControls &&other);
  
  /// \brief Obtain a double-precision abstract for this object.
  MMControlKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Obtain a single-precision abstract for this object.
  MMControlKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Prime the work unit counters based on a particular GPU configuration.
  ///
  /// \param gpu  Description of the GPU available to the runtime process
  void primeWorkUnitCounters(const GpuDetails &gpu, const AtomGraphSynthesis &poly_ag);

  /// \brief Increment the step counter, moving the controls to a different progress counter.
  void incrementStep();
  
#ifdef OMNI_USE_HPC
  /// \brief Upload the object's contents to the device (needed so that CPU-primed work unit
  ///        counters can go into effect)
  void upload();

  /// \brief Download the object's contents from the device (useful for debugging)
  void download();
#else
                             
#endif
  
private:
  int step_number;             ///< The step counter for the simulation
  double time_step;            ///< Time step of the simulation, units of femtoseconds.  This may
                               ///<   be truncated with the precision level of the abstract.
  double rattle_tol;           ///< Rattle tolerance, again truncated with the precision level of
                               ///<   the abstract
  Hybrid<int> vwu_progress;    ///< Progress counters through valence work units, an array of two
                               ///<   times the number of lanes per warp so that the valence work
                               ///<   units kernel can perform resets of elements in this array
                               ///<   with maximum efficiency.  The kernel will work based on the
                               ///<   progress counter with index equal to the step number modulo
                               ///<   twice the warp size.  When the step counter is equal to the
                               ///<   a multiple of the warp size, one of the warps in the first
                               ///<   block of this kernel will reset counters in array indices
                               ///<   [ 0, warp size ).  When the step counter is equal to a
                               ///<   multiple of twice the warp size, one of the warps in this
                               ///<   kernel will reset the counters in array indices [ warp size,
                               ///<   2 * warp size ).
  Hybrid<int> nbwu_progress;   ///< Progress counters for non-bonded work units, arrayed in a
                               ///<   manner similar to the valence work units.
  Hybrid<int> pmewu_progress;  ///< Progress counters for long-ranged, mesh-based PME work units,
                               ///<   arrayed in a manner similar to the valence work units.
  Hybrid<int> progress_data;   ///< ARRAY-kind Hybrid object targeted by the above POINTER-kind
                               ///<   Hybrid objects
};

} // namespace mm
} // namespace omni

#include "mm_controls.tpp"

#endif
