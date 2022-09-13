// -*-c++-*-
#ifndef STORMM_THERMOSTAT_H
#define STORMM_THERMOSTAT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
namespace trajectory {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Enumerate the various thermostats available for simulations
enum class ThermostatKind {
  NONE, ANDERSEN, LANGEVIN, BERENDSEN
};

/// \brief The maximum depth that the random number cache can take, to guard against excessive
///        memory use.
constexpr int maximum_random_cache_depth = 15;

/// \brief The default cache depth for a thermostat will give high efficiency in the number of
///        write transactions issued for the 256-bit read required: eight bytes written for every
///        byte read.
constexpr int default_thermostat_cache_depth = 8;
  
/// \brief The default simulation temperature, chemistry's standard temperature and pressure (in
///        units of Kelvin).
constexpr double default_simulation_temperature = 298.15;

/// \brief The default random seed for thermostats.  This need not be the same random seed used
///        for any other application--in fact, having it be distinct guards against generators in
///        other parts of the program re-using the same sequence.
constexpr int default_thermostat_random_seed = 1329440765;
  
/// \brief Store the parameters for a simulation thermostat.  Includes Berendsen, Andersen, and
///        Langevin methods.  This class can be assembled like many of the control objects, i.e.
///        MinimizationControls, based on namelists.
class Thermostat {
public:

  /// \brief Constructors for the Thermostat object.  Any information can be added via setter
  ///        functions after constructing the object.
  ///
  /// Overloaded:
  ///   - Construct a blank thermostat that applies no regulation
  ///   - Construct a specific type of thermostat with default settings
  ///   - Accept a type of thermostat and a temperature to maintain
  ///   - Accept a temperature evolution profile
  ///   - Accept a total number of atoms and a compartmentalization plan
  ///
  /// \param kind_in           The type of thermostat to implement
  /// \param temperature_in    A flat temperature to apply at all times
  /// \param temperature_in
  /// \{
  Thermostat();
  Thermostat(ThermostatKind kind_in);
  Thermostat(ThermostatKind kind_in, double temperature_in);
  Thermostat(ThermostatKind kind_in, double init_temperature_in, double final_temperature_in,
             int initial_evolution_step_in, int final_evolution_step_in);
  Thermostat(ThermostatKind kind_in, int atom_count_in,
             const std::vector<int> &compartment_limits_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in, int initial_evolution_step_in,
             int final_evolution_step_in);
  /// \}

  /// \brief Get the kind of thermostat
  ThermostatKind getKind() const;

  /// \brief Get the total number of atoms that this thermostat can serve.
  int getAtomCount() const;

  /// \brief Get the step number according to this thermostat.  The thermostat will serve as an
  ///        official reference for the official simulation step number.
  int getStepNumber() const;

  /// \brief Get the random seed used to initialize Xoshiro256++ state vectors.
  int getRandomSeed() const;

  /// \brief Get the random number cache depth.  The quantity of random numbers that will be
  ///        produced and cached for each atom is three times this setting, as it implies values
  ///        for random perturbations in the Cartesian X, Y, and Z directions.
  int getRandomCacheDepth() const;

  /// \brief Get the step at which to begin changing the temperature from its initial value T(init)
  ///        to its final value T(final).  This is the last step that the thermostat will pull the
  ///        system towards a value of T(init).  Afterwards, the target will begin to change.
  int getInitialEvolutionStep() const;

  /// \brief Get the step at which the temperature evolution is expected to be complete.  This is
  ///        the step when the thermostat will begin pulling the system towards T(final), although
  ///        the system itself may or may not be there by this time.
  int getFinalEvolutionStep() const;

  /// \brief Indicate whether this thermostat applies a common target temperature to all atoms, or
  ///        (if otherwise) there is a list of target initial and final temperatures with shared
  ///        values among various subdivisions (compartments) of the atom set.
  bool commonTemperature() const;
  
  /// \brief Get the initial target temperature for this thermostat, in units of Kelvin.
  ///
  /// \param atom_index  Index of the atom of interest (only provided if there are different
  ///                    regulated temperatures for unique subsets of atoms)
  double getInitialTemperature(int atom_index = 0) const;

  /// \brief Get the end-point target temperature for this thermostat, in units of Kelvin.
  ///
  /// \param atom_index  Index of the atom of interest
  double getFinalTemperature(int atom_index = 0) const;

  /// \brief Get the current target temperature for this thermostat, in units of Kelvin.
  ///
  /// \param atom_index  Index of the atom of interest
  double getCurrentTemperatureTarget(int atom_index = 0) const;

  /// \brief Set the type of thermostat.
  ///   
  /// \param kind_in  The type of thermostat to use
  void setKind(ThermostatKind kind_in);

  /// \brief Set the number of atoms for which the thermostat is responsible.
  ///
  /// \param atom_count_in  The total number of atoms to assign to this thermostat, whether all
  ///                       part of one system or a padded number collected within a synthesis
  void setAtomCount(int atom_count_in);

  /// \brief Compartmentalize the thermostated atoms to maintain different temperatures across
  ///        distinct, or formally separate, groups.
  ///
  /// \param compartment_limits_in    The compartment boundaries, required to be monotonically
  ///                                 increasing and positive.  Starting with any number other than
  ///                                 zero will lead to the interpretation that there is an
  ///                                 implicit compartment from zero up to the first numbered
  ///                                 index.  Ending with any number less than the atom count
  ///                                 likewise implies a final compartment between the last stated
  ///                                 index and the end of the atom list.  No index in this list
  ///                                 may exceed the atom count.
  /// \param initial_temperatures_in  The initial temperatures for each compartment, i.e. system.
  ///                                 The length of this list must equal the number of compartments
  ///                                 determined from compartment_limits_in and the atom count.
  /// \param initial_temperatures_in  The final temperatures for each compartment.  The length of
  ///                                 this list must equal the number of compartments declined from
  ///                                 compartment_limits_in and the atom count.
  void setCompartments(const std::vector<int> &compartment_limits_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in);
  
  /// \brief Set the initial target temperature.
  ///
  /// Overloaded:
  ///   - Set the temperature to a single value (if the final temperatures are also consistent,
  ///     this will set common_temperature to TRUE and de-allocate any existing Hybrid data
  ///     related to unique temperature compartments).
  ///   - Set the temperatures of unique subsets of the atoms to a series of values based on an
  ///     existing compartmentalization
  ///   - Set the temperatures of unique subsets of the atoms to a series of values based on a new
  ///     compartmentalization
  ///
  /// \param init_temp_in     Temperature or temperatures that the thermostat shall start with, and
  ///                         continue to apply until the initial step in its evolution
  /// \param compartments_in  A new compartmentalization scheme (the same interpretation described
  ///                         in setCompartments() above applies)
  /// \{
  void setTemperature(double temp_init_in, double temp_final_in = -1.0);
  void setTemperature(const std::vector<double> &temp_init_in,
                      const std::vector<double> &temp_final_in);
  void setTemperature(const std::vector<double> &temp_init_in,
                      const std::vector<double> &temp_final_in,
                      const std::vector<int> &compartment_limits_in);
  /// \}

  /// \brief Set the step at which temperature evolution, away from T(init), shall begin.
  ///
  /// \param step_in  The step setting
  void setInitialStep(int step_in);
  
  /// \brief Set the step at which temperature evolution, away from T(final), shall begin.
  ///
  /// \param step_in  The step setting
  void setFinalStep(int step_in);
  
  /// \brief Set the random number cache depth.  This value is kept on a tight leash, as it can
  ///        easily lead to allocating too much memory.  A hundred-thousand atom system with cache
  ///        depth 8 leads to over 13.7MB of memory allocation.
  ///
  /// \param depth_in  The depth to take
  void setRandomCacheDepth(int depth_in);

  /// \brief Validate the initial or final target temperature.
  void validateTemperature(double temperature_in);

  /// \brief Validate the initial and final steps for applying the temperature evolution.
  ///
  /// \param step_in  The step number to validate (it will be checked for a positive value)
  void validateEvolutionWindow();

  /// \brief Validate the random number cache depth.  A maxmimum cache depth of 15 (45 total
  ///        random numbers per atom) is enforced, and tightened in cases of exceptionally large
  ///        systems.
  void validateRandomCacheDepth();
  
  /// \brief Initialize the random state vectors of a thermostat, whether for Langevin or Andersen
  ///        thermostating.  This passes a call on to a general-purpose function for seeding the
  ///        generators and will recruit the GPU if available to expedite the process.
  ///
  /// \param new_seed      A new random number seed to apply, if different from the one used
  ///                      during construction of this Thermostat object
  /// \param scrub_cycles  Number of cycles of Xoshiro256++ generation to run in order to ensure
  ///                      that the newly seeded generators produce high-quality results
  void initializeRandomStates(int new_seed = default_thermostat_random_seed,
                              int scrub_cycles = 25);

  /// \brief Fill the random number cache for a subset of the atoms.
  ///
  /// \param index_start  Starting index in the list of all atoms
  /// \param index_end    Upper bound of atoms for which to cache random numbers
  void fillRandomCache(int index_start, int index_end);

#ifdef STORMM_USE_HPC
  /// \brief Upload the thermostat's data from the host to the HPC device.  Because the GPU is
  ///        used to initialize the vector of random number states (random_sv), this array is
  ///        uploaded and downloaded for synchronization upon initialization.  However, the initial
  ///        and final temperature arrays covering each atom require explicit upload.
  void upload();

  /// \brief Download the thermostat's data from the HPC device to the host.  This is needed to
  ///        synchronize progress made on the GPU if any processes are to be carried out by the
  ///        CPU host.
  void download();
#endif
  
private:
  ThermostatKind kind;             ///< The type of thermostat
  int atom_count;                  ///< The total number of atoms, and thus the number of random
                                   ///<   state vectors, for which this thermostat is responsible
                                   ///<   (this times random_cache_depth gives the length of
                                   ///<   random_sv_bank)
  int step_number;                 ///< Number of the dynamics step, which should be consistent
                                   ///<   with any other record of the current simulation step.
                                   ///<   The thermostat, which will be included in any simulation
                                   ///<   even if set to NONE kind, is a good place to keep the
                                   ///<   official step number.
  int random_seed;                 ///< Seed for the first random state vector (position 0 in the
                                   ///<   random_sv_bank array)
  int random_cache_depth;          ///< Depth of the random cache
  int initial_evolution_step;      ///< The first step at which to initiate temperature evolution
  int final_evolution_step;        ///< Final step at which to temperature evolution is complete
  bool common_temperature;         ///< Flag to indicate that all particles use identical initial
                                   ///<   and final temperatures.
  double initial_temperature;      ///< Initial temperature to apply to all particles
  double final_temperature;        ///< Final temperature to apply to all particles

  /// Temperatures to apply from step 0 to the initiation of any requested evolution, across
  ///   various compartments of the simulation.  Different compartments, i.e. systems, or
  ///   different components of each system, may have different starting temperatures, but the
  ///   evolution must proceed along the same schedule for all systems.
  /// \{
  Hybrid<double> initial_temperatures;
  Hybrid<float> sp_initial_temperatures;
  /// \}
  
  /// Temperatures to apply from the end of any requested evolution until the end of the simulation
  /// \{
  Hybrid<double> final_temperatures;
  Hybrid<float> sp_final_temperatures;
  /// \}

  /// Bounds array for various compartments of the particles affected by this thermostat.  The name
  /// is more general than "atom_starts" or "system_bounds" because this object is designed to
  /// serve one or many systems, and could, in theory, be used to thermostat different parts of one
  /// or all systems at different values.  This array is not carried through to the GPU: it defines
  /// the values of the initial and final temperatures array, which are communicated to the GPU.
  std::vector<int> compartment_limits;  

  /// Xoshiro256++ state vectors for creating random numbers, one per atom of the simulation
  Hybrid<ullint4> random_sv;      

  /// \brief Allocate space for the random state vector and random number cache
  Hybrid<double> random_cache;
  Hybrid<float> sp_random_cache;

  /// \brief Allocate storage space for random number generator state vectors and cached random
  ///        numbers.  This will occur when the thermostat is set to Langevin or Andersen, or when
  ///        the atom count changes and the thermostat type is already Langevin or Andersen.
  void allocateRandomStorage();

  /// \brief Set the flag to indicate that a pair of common temperatures is in use for initial
  ///        and final targets of the thermostat across all atoms.  Setting the flag to FALSE will
  ///        trigger allocation of the member variables initial_temperatures and final_temperatures
  ///        to the proper sizes, and setting the flag to TRUE will de-allocate these arrays.
  ///
  /// \param setting_in  The state of temperature commonality
  void setCommonTemperature(bool setting_in);
};

/// \brief Return the name of the thermostat choice (an enumerator string conversion function)
///
/// \param kind  The type of thermostat
std::string getThermostatName(ThermostatKind kind);

} // namespace trajectory
} // namespace stormm

#endif
