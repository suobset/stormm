// -*-c++-*-
#ifndef OMNI_SCORECARD_H
#define OMNI_SCORECARD_H

#include "Constants/fixed_precision.h"
#include "Accelerator/hybrid.h"
#include "energy_enumerators.h"

namespace omni {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
using numerics::default_energy_scale_bits;

/// \brief Read-only abstract for the ScoreCard object.  This is needed more for completeness than
///        anything else, but could be useful in cases where virials are to be translated into
///        rescaling coefficients on the GPU and then immediately used to rescale positions in a
///        Berendsen barostatting situation.
struct ScoreCardReader {

  /// \brief The constructor is, as usual, a collection of the relevant constants and pointers.
  ScoreCardReader(int system_count_in, int data_stride_in, int sampled_step_count_in,
                  float nrg_scale_f_in, double nrg_scale_lf_in, float inverse_nrg_scale_f_in,
                  double inverse_nrg_scale_lf_in, const llint* instantaneous_accumulators_in,
                  const double* running_accumulators_in, const double* squared_accumulators_in,
                  const llint* time_series_in);

  const int system_count;                  ///< Number of independent systems tracked
  const int data_stride;                   ///< Size of the StateVariable enumerator rounded up
                                           ///<   to the nearest multiple of the HPC warp size
  const int sampled_step_count;            ///< The number of steps in the sample (a nice result
                                           ///<   is that updating this member variable in the
                                           ///<   parent ScoreCard right before making the
                                           ///<   abstract can send this data to the GPU)
  const float nrg_scale_f;                 ///< Conversion factor for fixed-precision accumulation
  const double nrg_scale_lf;               ///< Conversion factor for fixed-precision accumulation
  const float inverse_nrg_scale_f;         ///< Conversion for fixed-precision interpretation
  const double inverse_nrg_scale_lf;       ///< Conversion for fixed-precision interpretation
  const llint* instantaneous_accumulators; ///< State variables for each system
  const double* running_accumulators;      ///< Running sums of state variables for each system
  const double* squared_accumulators;      ///< Running squared sums of state variables for each
                                           ///<   system
  const llint* time_series;                ///< Details of energy values in each component's
                                           ///<   accumulator at each time step
};

/// \brief Writeable abstract for the ScoreCard object, useful for accumulating energies in many
///        kernels.
struct ScoreCardWriter {

  /// \brief The constructor is, as usual, a collection of the relevant constants and pointers.
  ScoreCardWriter(int system_count_in, int data_stride_in, int sampled_step_count_in,
                  float nrg_scale_f_in, double nrg_scale_lf_in, float inverse_nrg_scale_f_in,
                  double inverse_nrg_scale_lf_in, llint* instantaneous_accumulators_in,
                  double* running_accumulators_in, double* squared_accumulators_in,
                  llint* time_series_in);

  const int system_count;             ///< Number of independent systems tracked
  const int data_stride;              ///< Size of the StateVariable enumerator rounded up to the
                                      ///<   nearest multiple of the HPC warp size
  const int sampled_step_count;       ///< The number of steps in the sample
  const float nrg_scale_f;            ///< Conversion factor for fixed-precision accumulation
  const double nrg_scale_lf;          ///< Conversion factor for fixed-precision accumulation
  const float inverse_nrg_scale_f;    ///< Conversion for fixed-precision interpretation
  const double inverse_nrg_scale_lf;  ///< Conversion for fixed-precision interpretation
  llint* instantaneous_accumulators;  ///< State variables for each system
  double* running_accumulators;       ///< Running sums of state variables for each system
  double* squared_accumulators;       ///< Running squared sums of state variables for each system
  llint* time_series;                 ///< Details of energy values in each component's
                                      ///<   accumulator at each time step
};

/// \brief Track the energy components of a collection of systems in an HPC-capable array.  This
///        object uses the familiar trick of defining an enumerator (StateVariables) with a final
///        entry to indicate its total length, so that if more energy components need to be tracked
///        in the future the storage and indexing can automatically adjust with new entries.
class ScoreCard {
public:

  /// \brief The constructor requires only the number of systems.
  ///
  /// \param system_count_in    The number of systems to track
  /// \param capacity_in        The capacity to initially allocate for
  /// \param nrg_scale_bits_in  Number of bits after the decimal with which to store energy values
  ScoreCard(int system_count_in, int capacity_in = 16,
            int nrg_scale_bits_in = default_energy_scale_bits);

  /// \brief Get the number of systems that this object is tracking
  int getSystemCount() const;

  /// \brief Get the number of steps that have been sampled
  int getSampleSize() const;

  /// \brief Get the number of bits of fixed precision to which results are stored
  int getEnergyScaleBits() const;

  /// \brief Get the energy scaling factors in single- or double-precision floating point format
  /// \{
  template <typename T> T getEnergyScalingFactor() const;
  template <typename T> T getInverseEnergyScalingFactor() const;
  /// \}

#ifdef OMNI_USE_HPC
  /// \brief Upload data to the device (this could be useful in situations where the CPU is
  ///        required to compute additional energy terms, or other energy quantities are being
  ///        brought in from some module outside the typical OMNI molecular mechanics routines).
  void upload();

  /// \brief Download all data from the device (each of the report(...) functions below will call
  ///        the download method for the appropriate array in an HPC setting, as the data on the
  ///        device is assumed to be the most important).
  void download();
#endif

  /// \brief Get the appropriate abstract based on the const-ness of the abstract
  ///
  /// \param tier  Get pointers to data on the host or on the HPC accelerator device
  /// \{
  const ScoreCardReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  ScoreCardWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Reserve space for storing sampled energy component values.
  ///
  /// \param new_capacity  The new capacity to allocate for
  void reserve(int new_capacity);
  
  /// \brief Contribute a result into one of the instantaneous state variable accumulators.  This
  ///        is for CPU activity; the contributions will occur as part of each energy kernel using
  ///        pointers on the GPU.
  ///
  /// \param var           The state variable to which this contribution belongs, i.e. bond energy
  /// \param amount        Amount to contribute to said state variable (in fixed precision format)
  /// \param system_index  Index of the system (among a list of those being tracked) that the
  ///                      contrbution describes
  void contribute(StateVariable var, llint amount, int system_index = 0);

  /// \brief Increment the number of sampled steps.  This will automatically allocate additional
  ///        capacity if the sampled step count reaches the object's capacity.
  void incrementSampleCount();
  
  /// \brief Report the total energy for all systems.  Each result will be summed in the internal
  ///        fixed-point accumulation before conversion to real values in units of kcal/mol.
  std::vector<double> reportTotalEnergies();
  
  /// \brief Report the total energy for one system in the list, following the summation protocol
  ///        from reportTotalEnergies() above.
  ///
  /// \param system_index  Index of the system of interest
  double reportTotalEnergy(int system_index = 0);
  
  /// \brief Report instantaneous results in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \{
  std::vector<double> reportInstantaneousStates();
  std::vector<double> reportInstantaneousStates(int system_index);
  std::vector<double> reportInstantaneousStates(StateVariable aspect);
  double reportInstantaneousStates(StateVariable aspect, int system_index);
  /// \}

  /// \brief Report averaged results in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \{
  std::vector<double> reportAverageStates();
  std::vector<double> reportAverageStates(int system_index);
  std::vector<double> reportAverageStates(StateVariable aspect);
  double reportAverageStates(StateVariable aspect, int system_index);
  /// \}

  /// \brief Report averaged results in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \{
  std::vector<double> reportVarianceOfStates();
  std::vector<double> reportVarianceOfStates(int system_index);
  std::vector<double> reportVarianceOfStates(StateVariable aspect);
  double reportVarianceOfStates(StateVariable aspect, int system_index);
  /// \}

private:
  int system_count;                          ///< The number of systems in the collection (each
                                             ///<   system will get a separate set of accumulators
                                             ///<   and averages, although combined statistics can
                                             ///<   be collected).
  int data_stride;                           ///< The number of accumulators (i.e. bond, angle, ...
                                             ///<   dU/dLambda) per system.  This indicates the
                                             ///<   length of each system's subset of the data in
                                             ///<   the following arrays.
  int sampled_step_count;                    ///< The number of samples going into running averages
                                             ///<   and standard deviations
  int sample_capacity;                       ///< The maximum number of samples available to this
                                             ///<   ScoreCard (the time_series_accumulators array
                                             ///<   will be resized if the sampled_step_count
                                             ///<   reaches this capacity).
  
  /// Scaling factors for fixed-precision accumulation and interpretation.  The defaults from
  /// Constants/fixed_precision.h are good suggestions, but the precision that the code works in
  /// may offer a tradeoff between speed and stability.  For energy accumulation, the question is
  /// not what happens to the future of the dynamics, merely what values are reported and how
  /// accurate they might be.  However, in the interest of knowing quantities to very high
  /// precision these values can be altered at the user or developer's discretion.
  /// \{
  int nrg_scale_bits;
  float nrg_scale_f;
  double nrg_scale_lf;
  float inverse_nrg_scale_f;
  double inverse_nrg_scale_lf;
  /// \}

  Hybrid<llint> instantaneous_accumulators;  ///< Instantaneous accumulators for reporting the
                                             ///<   energy at one time step, accumulated in fixed
                                             ///<   precision after scaling by
                                             ///<   default_energy_scale_lf (see fixed_precision.h
                                             ///<   in Constants/ ).
  Hybrid<double> running_accumulators;       ///< Running sums of the energy and other state
                                             ///<   variables, collected over each sampled time
                                             ///<   step.  These are stored in double precision to
                                             ///<   guard against overflow in a long simulation.
                                             ///<   In an HPC setting, these will be computed on
                                             ///<   the device after preparing the instantaneous
                                             ///<   accumulators.
  Hybrid<double> squared_accumulators;       ///< Running sums of the squared energy components and
                                             ///<   other stat variables collected over each
                                             ///<   sampled time step.  In an HPC setting, these
                                             ///<   will be computed on the device after preparing
                                             ///<   the instantaneous accumulators.
  Hybrid<llint> time_series_accumulators;    ///< A history of values for each energy accumulator
                                             ///<   at each sampled time step.  This detailed
                                             ///<   array is resized as needed, or can be reserved.
};

} // namespace energy
} // namespace omni

#include "scorecard.tpp"

#endif
