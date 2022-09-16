#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Random/random.h"
#ifdef STORMM_USE_HPC
#include "Random/hpc_random.h"
#endif
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "thermostat.h"

namespace stormm {
namespace trajectory {

using card::HybridKind;
using card::HybridTargetLevel;
using math::roundUp;
using parse::realToString;
using parse::NumberFormat;
using random::fillRandomCache;
using random::initXoshiro256ppArray;
using random::RandomAlgorithm;
using random::RandomNumberKind;
using random::Xoshiro256ppGenerator;

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat() :
  kind{ThermostatKind::NONE}, atom_count{0}, padded_atom_count{0}, step_number{0},
    random_seed{default_thermostat_random_seed},
    random_cache_depth{default_thermostat_cache_depth},
    initial_evolution_step{0}, final_evolution_step{0}, common_temperature{true},
    initial_temperature{default_simulation_temperature},
    final_temperature{default_simulation_temperature},
    initial_temperatures{HybridKind::ARRAY, "tstat_init_temp"},
    sp_initial_temperatures{HybridKind::ARRAY, "tstat_init_tempf"},
    final_temperatures{HybridKind::ARRAY, "tstat_final_temp"},
    sp_final_temperatures{HybridKind::ARRAY, "tstat_final_tempf"},
    compartment_limits{},
    random_state_vector_xy{HybridKind::ARRAY, "tstat_rng_xycache"},
    random_state_vector_zw{HybridKind::ARRAY, "tstat_rng_zwcache"}
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;  
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in, const double temperature_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;
  initial_temperature = temperature_in;
  final_temperature = temperature_in;
  validateTemperature(initial_temperature);
  validateTemperature(final_temperature);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in, const double initial_temperature_in,
                       const double final_temperature_in, const int initial_evolution_step_in,
                       const int final_evolution_step_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;
  initial_temperature = initial_temperature_in;
  final_temperature = final_temperature_in;
  initial_evolution_step = initial_evolution_step_in;
  final_evolution_step = final_evolution_step_in;
  validateTemperature(initial_temperature);
  validateTemperature(final_temperature);
  validateEvolutionWindow();
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in, const int atom_count_in,
                       const std::vector<int> &compartment_limits_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;
  setAtomCount(atom_count_in);
  setCompartments(compartment_limits_in, initial_temperatures_in, final_temperatures_in);
  initial_evolution_step = initial_evolution_step_in;
  final_evolution_step = final_evolution_step_in;
  validateEvolutionWindow();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::allocateRandomStorage() {

  // Do not allocate extra memory which will not be used by certain thermostats.  When seeding
  // particle velocities, the process looks very much like applying an Andersen thermostat once at
  // the outset of the simulation.  In that case, it is faster to have a single generator, driven
  // by a single CPU thread, follow a similar strategy of creating up to 1024 generators and then
  // seeding a series of atoms' velocities with each of them.  That velocity initialization process
  // will get its own function, deferring to the array of generators for each atom in the case of
  // Andersen and Langevin thermostats, or using a temporary array of generators, outside of the
  // Thermostat object, in the cases of Berendesen or no thermostat.
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    random_state_vector_xy.resize(0);
    random_state_vector_zw.resize(0);
    random_cache.resize(0);
    sp_random_cache.resize(0);
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    {
      const size_t atom_count_zu = roundUp(atom_count, warp_size_int);
      const size_t rc_depth_zu = random_cache_depth * 3;
      random_state_vector_xy.resize(atom_count_zu);
      random_state_vector_zw.resize(atom_count_zu);
      random_cache.resize(atom_count_zu * rc_depth_zu);
      sp_random_cache.resize(atom_count_zu * rc_depth_zu);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setCommonTemperature(const bool setting_in) {
  common_temperature = setting_in;
  if (common_temperature) {
    initial_temperatures.resize(0);
    sp_initial_temperatures.resize(0);
    final_temperatures.resize(0);
    sp_final_temperatures.resize(0);
  }
  else {
    initial_temperatures.resize(atom_count);
    sp_initial_temperatures.resize(atom_count);
    final_temperatures.resize(atom_count);
    sp_final_temperatures.resize(atom_count);
  }
}
  
//-------------------------------------------------------------------------------------------------
ThermostatKind Thermostat::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getStepNumber() const {
  return step_number;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getRandomSeed() const {
  return random_seed;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getRandomCacheDepth() const {
  return random_cache_depth;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getInitialEvolutionStep() const {
  return initial_evolution_step;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getFinalEvolutionStep() const {
  return final_evolution_step;
}

//-------------------------------------------------------------------------------------------------
ullint4 Thermostat::getGeneratorState(const int atom_index, const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const ullint2 xy_state = random_state_vector_xy.readHost(atom_index);
      const ullint2 zw_state = random_state_vector_zw.readHost(atom_index);
      return { xy_state.x, xy_state.y, zw_state.x, zw_state.y };
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const ullint2 xy_state = random_state_vector_xy.readDevice(atom_index);
      const ullint2 zw_state = random_state_vector_zw.readDevice(atom_index);
      return { xy_state.x, xy_state.y, zw_state.x, zw_state.y };
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getCachedRandomResult(const PrecisionModel prec, const int atom_index,
                                         const int cache_row, const HybridTargetLevel tier) const {
  const size_t pos = (static_cast<size_t>(cache_row) * static_cast<size_t>(padded_atom_count)) +
                     static_cast<size_t>(atom_index);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      return random_cache.readHost(pos);
    case PrecisionModel::SINGLE:
      return sp_random_cache.readHost(pos);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      return random_cache.readDevice(pos);
    case PrecisionModel::SINGLE:
      return sp_random_cache.readDevice(pos);
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool Thermostat::commonTemperature() const {
  return common_temperature;
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getInitialTemperature(const int atom_index) const {
  if (common_temperature) {
    return initial_temperature;
  }
  else {
    return initial_temperatures.readHost(atom_index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getFinalTemperature(const int atom_index) const {
  if (common_temperature) {
    return final_temperature;
  }
  else {
    return final_temperatures.readHost(atom_index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getCurrentTemperatureTarget(const int atom_index) const {
  if (step_number <= initial_evolution_step) {
    return getInitialTemperature(atom_index);
  }
  else if (step_number <= final_evolution_step) {
    const double progress = static_cast<double>(step_number - initial_evolution_step) /
                            static_cast<double>(final_evolution_step - initial_evolution_step);
    if (common_temperature) {
      return ((1.0 - progress) * initial_temperature) + (progress * final_temperature);
    }
    else {
      return ((1.0 - progress) * initial_temperatures.readHost(atom_index)) +
             (progress * final_temperatures.readHost(atom_index));
    }
  }
  else {
    return getFinalTemperature(atom_index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const ThermostatReader<double> Thermostat::dpData(const HybridTargetLevel tier) const {
  return ThermostatReader<double>(kind, atom_count, padded_atom_count, step_number,
                                  random_cache_depth, initial_evolution_step, final_evolution_step,
                                  common_temperature, initial_temperature, final_temperature,
                                  initial_temperatures.data(tier), final_temperatures.data(tier),
                                  random_state_vector_xy.data(tier),
                                  random_state_vector_zw.data(tier), random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
ThermostatWriter<double> Thermostat::dpData(const HybridTargetLevel tier) {
  return ThermostatWriter<double>(kind, atom_count, padded_atom_count, step_number,
                                  random_cache_depth, initial_evolution_step, final_evolution_step,
                                  common_temperature, initial_temperature, final_temperature,
                                  initial_temperatures.data(tier), final_temperatures.data(tier),
                                  random_state_vector_xy.data(tier),
                                  random_state_vector_zw.data(tier), random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
const ThermostatReader<float> Thermostat::spData(const HybridTargetLevel tier) const {
  return ThermostatReader<float>(kind, atom_count, padded_atom_count, step_number,
                                 random_cache_depth, initial_evolution_step, final_evolution_step,
                                 common_temperature, initial_temperature, final_temperature,
                                 sp_initial_temperatures.data(tier),
                                 sp_final_temperatures.data(tier),
                                 random_state_vector_xy.data(tier),
                                 random_state_vector_zw.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
ThermostatWriter<float> Thermostat::spData(const HybridTargetLevel tier) {
  return ThermostatWriter<float>(kind, atom_count, padded_atom_count, step_number,
                                 random_cache_depth, initial_evolution_step, final_evolution_step,
                                 common_temperature, initial_temperature, final_temperature,
                                 sp_initial_temperatures.data(tier),
                                 sp_final_temperatures.data(tier),
                                 random_state_vector_xy.data(tier),
                                 random_state_vector_zw.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setKind(const ThermostatKind kind_in) {
  kind = kind_in;
  allocateRandomStorage();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setAtomCount(int atom_count_in) {
  atom_count = atom_count_in;
  padded_atom_count = roundUp(atom_count, warp_size_int);
  allocateRandomStorage();
}
  
//-------------------------------------------------------------------------------------------------
void Thermostat::setCompartments(const std::vector<int> &compartment_limits_in,
                                 const std::vector<double> &initial_temperatures_in,
                                 const std::vector<double> &final_temperatures_in) {

  // Check that there is an initial and final temperature for each compartment
  compartment_limits = compartment_limits_in;
  if (compartment_limits[0] != 0) {
    const std::vector<int> zero_vec(1, 0);
    compartment_limits.insert(compartment_limits.begin(), zero_vec.begin(), zero_vec.end());
  }
  else if (compartment_limits.back() != atom_count) {
    compartment_limits.push_back(atom_count);
  }

  // Check the sanity of the compartment demarcations
  const int n_comp = compartment_limits.size();
  for (int i = 1; i < n_comp; i++) {
    if (compartment_limits[i] <= compartment_limits[i - 1]) {
      rtErr("Compartment " + std::to_string(i) + " has invalid limits " +
            std::to_string(compartment_limits[i - 1]) + " - " +
            std::to_string(compartment_limits[i]) + ".", "Thermostat", "setCompartments");
    }
  }

  // Allocate arrays to hold individual atom temperatures.  This is done implicitly by marking
  // the fact that there are no longer a single, common initial or final temperatures.
  setCommonTemperature(false);
  double* ditemp_ptr = initial_temperatures.data();
  float* fitemp_ptr = sp_initial_temperatures.data();
  double* dftemp_ptr = final_temperatures.data();
  float* fftemp_ptr = sp_final_temperatures.data();
  for (int i = 0; i < n_comp; i++) {
    validateTemperature(initial_temperatures_in[i]);
    validateTemperature(final_temperatures_in[i]);
    for (int j = compartment_limits[i]; j < compartment_limits[i + 1]; j++) {
      ditemp_ptr[j] = initial_temperatures_in[i];
      fitemp_ptr[j] = initial_temperatures_in[i];
      dftemp_ptr[j] = final_temperatures_in[i];
      fftemp_ptr[j] = final_temperatures_in[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setTemperature(const double temp_init_in, const double temp_final_in) {
  setCommonTemperature(true);
  initial_temperature = temp_init_in;
  if (temp_final_in < 0.0) {
    final_temperature = initial_temperature;
  }
  validateTemperature(initial_temperature);
  validateTemperature(final_temperature);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setTemperature(const std::vector<double> &temp_init_in,
                                const std::vector<double> &temp_final_in) {
  
  // Fail immediately if the initial and final temperature vectors do not align.
  if (temp_init_in.size() != temp_final_in.size()) {
    rtErr("A thermostat cannot accept " + std::to_string(temp_init_in.size()) + " initial and " +
          std::to_string(temp_final_in.size()) + " final temperatures.", "Thermostat",
          "setTemperature");
  }
  if (common_temperature) {

    // If a common temperature is currently in effect, the only option is for the input vector to
    // account for every atom.  Otherwise produce an error.
    if (static_cast<int>(temp_init_in.size()) == atom_count) {
      setCommonTemperature(false);
      double* ditemp_ptr = initial_temperatures.data();
      float* fitemp_ptr = sp_initial_temperatures.data();
      double* dftemp_ptr = final_temperatures.data();
      float* fftemp_ptr = sp_final_temperatures.data();
      for (int i = 0; i < atom_count; i++) {
        ditemp_ptr[i] = temp_init_in[i];
        fitemp_ptr[i] = temp_init_in[i];
        dftemp_ptr[i] = temp_final_in[i];
        fftemp_ptr[i] = temp_final_in[i];
      }
    }
    else {
      rtErr("A thermostat with " + std::to_string(atom_count) + " atoms an no prior "
            "compartmentalization cannot accept series of " +
            std::to_string(temp_init_in.size()) + " temperatures.  Provide one initial and one "
            "final temperature per atom or a corresponding compartmentalization plan.",
            "Thermostat", "setInitialTemperature");
    }
  }
  else {

    // If there is already a compartmentalization plan, the input temperatures must meet it.
    if (temp_init_in.size() == compartment_limits.size() - 1LLU) {
      const size_t n_comp = temp_init_in.size();
      double* ditemp_ptr = initial_temperatures.data();
      float* fitemp_ptr = sp_initial_temperatures.data();
      double* dftemp_ptr = final_temperatures.data();
      float* fftemp_ptr = sp_final_temperatures.data();
      for (size_t i = 0LLU; i < n_comp; i++) {
        validateTemperature(temp_init_in[i]);
        validateTemperature(temp_final_in[i]);
        for (int j = compartment_limits[i]; j < compartment_limits[i + 1]; j++) {
          ditemp_ptr[j] = temp_init_in[i];
          fitemp_ptr[j] = temp_init_in[i];
          dftemp_ptr[j] = temp_final_in[i];
          fftemp_ptr[j] = temp_final_in[i];
        }
      }
    }
    else {
      rtErr("A thermostat with " + std::to_string(compartment_limits.size() - 1LLU) +
            " pre-existing groups of temperature-regulated atoms cannot accept lists of " +
            std::to_string(temp_init_in.size()) + " temperatures.  Provide one initial and one "
            "final temperature per group, or a corresponding compartmentalization plan.",
            "Thermostat", "setInitialTemperature");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setTemperature(const std::vector<double> &temp_init_in,
                                const std::vector<double> &temp_final_in,
                                const std::vector<int> &compartment_limits_in) {
  setCompartments(compartment_limits_in, temp_init_in, temp_final_in);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setInitialEvolutionStep(const int step_in) {
  initial_evolution_step = step_in;  
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setFinalEvolutionStep(const int step_in) {
  final_evolution_step = step_in;  
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setRandomCacheDepth(const int depth_in) {
  random_cache_depth = depth_in;
  validateRandomCacheDepth();
  allocateRandomStorage();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::incrementStep() {
  step_number += 1;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::decrementStep() {
  step_number -= 1;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateTemperature(const double temperature_in) {
  if (temperature_in < 0.0 || temperature_in >= 1.0e5) {
    rtErr("A temperature of " + realToString(temperature_in, 11, 4, NumberFormat::STANDARD_REAL) +
          " is not a sensible choice.", "Thermostat", "validateTemperature");
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateEvolutionWindow() {
  if (initial_evolution_step < 0) {
    initial_evolution_step = 0;
  }
  if (final_evolution_step <= initial_evolution_step) {
    final_evolution_step = initial_evolution_step;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateRandomCacheDepth() {
  if (random_cache_depth > maximum_random_cache_depth) {

    // A silent change occurs here.  The only effect that a bad cache depth could have is that the
    // simulation runs at an imperceptibly different rate than the one that might be expected.
    random_cache_depth = maximum_random_cache_depth;
  }
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    random_cache_depth = 0;
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::initializeRandomStates(const PrecisionModel prec, const int new_seed,
                                        const int scrub_cycles, const GpuDetails &gpu) {
  random_seed = new_seed;
#ifdef STORMM_USE_HPC
  initXoshiro256ppArray(&random_state_vector_xy, &random_state_vector_zw, new_seed, scrub_cycles,
                        gpu);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &random_cache,
                    padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                    RandomNumberKind::GAUSSIAN, 0, atom_count, gpu);
    break;
  case PrecisionModel::SINGLE:
    fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &sp_random_cache,
                    padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                    RandomNumberKind::GAUSSIAN, 0, atom_count, gpu);
    break;
  }
#else
  initXoshiro256ppArray(&random_state_vector_xy, &random_state_vector_zw, new_seed, scrub_cycles);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &random_cache,
                    padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                    RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  case PrecisionModel::SINGLE:
    fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &sp_random_cache,
                    padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                    RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void Thermostat::refresh(const size_t index_start, const size_t index_end,
                         const PrecisionModel mode) {

  // Return immediately with a warning if a developer tries to issue random numbers for thermostats
  // that do not use them.
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    rtWarn("The " + getThermostatName(kind) + " thermostat does not utilize random numbers.  No "
           "production will occur.");
    return;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    break;
  }
  switch (mode) {
  case PrecisionModel::DOUBLE:
    fillRandomCache(random_state_vector_xy.data(), random_state_vector_zw.data(),
                    random_cache.data(), atom_count, random_cache_depth * 3,
                    RandomAlgorithm::XOSHIRO_256PP, RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  case PrecisionModel::SINGLE:
    fillRandomCache(random_state_vector_xy.data(), random_state_vector_zw.data(),
                    sp_random_cache.data(), atom_count, random_cache_depth * 3,
                    RandomAlgorithm::XOSHIRO_256PP, RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void Thermostat::upload() {
  initial_temperatures.upload();
  sp_initial_temperatures.upload();
  final_temperatures.upload();
  sp_final_temperatures.upload();
  random_state_vector_xy.upload();
  random_state_vector_zw.upload();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::download() {
  initial_temperatures.download();
  sp_initial_temperatures.download();
  final_temperatures.download();
  sp_final_temperatures.download();
  random_state_vector_xy.download();
  random_state_vector_zw.download();
}
#endif

//-------------------------------------------------------------------------------------------------
std::string getThermostatName(const ThermostatKind kind) {
  switch (kind) {
  case ThermostatKind::NONE:
    return std::string("NONE");
  case ThermostatKind::ANDERSEN:
    return std::string("ANDERSEN");
  case ThermostatKind::LANGEVIN:
    return std::string("LANGEVIN");
  case ThermostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace stormm
