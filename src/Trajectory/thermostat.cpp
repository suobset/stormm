#include "copyright.h"
#include "Constants/scaling.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "thermostat.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat() :
    kind{ThermostatKind::NONE}, atom_count{0}, step_number{0}, random_seed{0},
    random_cache_depth{0}, initial_evolution_step{0}, final_evolution_step{0},
    common_temperature{true},
    initial_temperature{default_simulation_temperature},
    final_temperature{default_simulation_temperature},
    initial_temperatures{HybridKind::ARRAY, "tstat_init_temp"},
    sp_initial_temperatures{HybridKind::ARRAY, "tstat_init_tempf"},
    final_temperatures{HybridKind::ARRAY, "tstat_final_temp"},
    sp_final_temperatures{HybridKind::ARRAY, "tstat_final_tempf"},
    compartment_limits{}
    random_sv{HybridKind::ARRAY, "tstat_rng_cache"}
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;  
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
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in, const int atom_count_in,
                       const std::vector<int> &compartment_limits_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in) :
    Thermostat()
{
  // Additional initializations
  kind = kind_in;
  setAtomCount(atom_count_in);
  setCompartments(compartment_limits_in, initial_temperatures_in, final_temperatures_in);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::allocateRandomStorage() {
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    random_sv.resize(0);
    random_cache.resize(0);
    sp_random_cache.resize(0);
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    {
      const size_t atom_count_zu = atom_count;
      const size_t rc_depth_zu = random_cache_depth;
      random_sv.resize(atom_count_zu);
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
void Thermostat::setKind(const ThermostatKind kind_in) {
  kind = kind_in;
  allocateRandomStorage();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setAtomCount(int atom_count_in) {
  atom_count = atom_count_in;
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
  else if (compartment_limits.back() != atom_count_in) {
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
  (...)
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setFinalTemperature(const double temp_final_in) {
  if (common_temperature) {
    final_temperature = temp_final_in;    
  }
  else {

    // Test whether there is a consistent final temperature.  If so, set the common_temperature
    // flag and resize the Hybrid objects containing per-particle temperatures to zero.
    const size_t n_comp = compartment_limits.size();
    bool all_same = true;
    const Approx trial_temp(final_temperatures[0], constants::small);
    for (int i = 1; i < n_comp; i++) {
      all_same = (all_same && trial_temp.test(final_temperatures[i]));      
    }
    if (all_same) {
      setCommonTemperature(true);
    }
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void Thermostat::upload() {
  initial_temperatures.upload();
  sp_initial_temperatures.upload();
  final_temperatures.upload();
  sp_final_temperatures.upload();
  random_sv.upload();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::download() {
  initial_temperatures.download();
  sp_initial_temperatures.download();
  final_temperatures.download();
  sp_final_temperatures.download();
  random_sv.download();
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
