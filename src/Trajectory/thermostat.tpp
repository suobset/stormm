// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
ThermostatReader<T>::ThermostatReader(const ThermostatKind kind_in, const int atom_count_in,
                                      const int padded_atom_count_in, const int step_in,
                                      const int depth_in, const int init_evolution_in,
                                      const int end_evolution_in, const bool common_temperature_in,
                                      const T init_temperature_in, const T final_temperature_in,
                                      const T* init_temperatures_in,
                                      const T* final_temperatures_in, const ullint2* state_xy_in,
                                      const ullint2* state_zw_in, const T* cache_in) :
    kind{kind_in}, atom_count{atom_count_in}, padded_atom_count{padded_atom_count_in},
    step{step_in}, depth{depth_in}, init_evolution{init_evolution_in},
    end_evolution{end_evolution_in}, common_temperature{common_temperature_in},
    init_temperature{init_temperature_in}, final_temperature{final_temperature_in},
    init_temperatures{init_temperatures_in}, final_temperatures{final_temperatures_in},
    state_xy{state_xy_in}, state_zw{state_zw_in}, cache{cache_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ThermostatWriter<T>::ThermostatWriter(const ThermostatKind kind_in, const int atom_count_in,
                                      const int padded_atom_count_in, const int step_in,
                                      const int depth_in, const int init_evolution_in,
                                      const int end_evolution_in, const bool common_temperature_in,
                                      const T init_temperature_in,
                                      const T final_temperature_in, const T* init_temperatures_in,
                                      const T* final_temperatures_in, ullint2* state_xy_in,
                                      ullint2* state_zw_in, T* cache_in) :
    kind{kind_in}, atom_count{atom_count_in}, padded_atom_count{padded_atom_count_in},
    step{step_in}, depth{depth_in}, init_evolution{init_evolution_in},
    end_evolution{end_evolution_in}, common_temperature{common_temperature_in},
    init_temperature{init_temperature_in}, final_temperature{final_temperature_in},
    init_temperatures{init_temperatures_in}, final_temperatures{final_temperatures_in},
    state_xy{state_xy_in}, state_zw{state_zw_in}, cache{cache_in}
{}

} // namespace trajectory
} // namespace stormm
