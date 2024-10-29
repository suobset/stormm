// -*-c++-*-
#include <cmath>
#include "copyright.h"
#include "exchange_nexus.h"
#include "temp_distributions.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Namelists/nml_remd.h"
#include "Namelists/nml_dynamics.h"
#include "Potential/scorecard.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinate_series.h"

namespace stormm {
namespace sampling {

using symbols::boltzmann_constant;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using synthesis::StructureSource;
using synthesis::SynthesisCacheMap;
using synthesis::SystemGrouping;
using namelist::RemdControls;
using namelist::DynamicsControls;
using topology::AtomGraph;
using trajectory::CoordinateSeries;

//-------------------------------------------------------------------------------------------------
// TODO: Pass relevant &dynamics data in these constructors
ExchangeNexus::ExchangeNexus( const int system_count_in, const int total_swap_count,
                              const std::string &remd_type_in, const int frequency_swaps_in, 
                              const std::string &swap_store_in,
                              const std::string &temperature_dist_in,
                              const double exchange_probability_in, const double tolerance_in,
                              const int max_replicas_in, 
                              const double initial_temperature_in,
                              const double equilibrium_temperature_in,
                              const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                              const ScoreCard *sc_in) :
    system_count{system_count_in}, total_swaps{total_swap_count},
    remd_type{remd_type_in}, frequency_swaps{frequency_swaps_in},
    swap_store{swap_store_in}, temperature_dist{temperature_dist_in},
    exchange_probability{exchange_probability_in}, tolerance{tolerance_in}, 
    max_replicas{max_replicas_in},
    initial_temperature{initial_temperature_in}, 
    equilibrium_temperature{equilibrium_temperature_in},
    ps{const_cast<PhaseSpaceSynthesis*>(ps_in)},
    ag{const_cast<AtomGraphSynthesis*>(ag_in)},
    sc{const_cast<ScoreCard*>(sc_in)}
{}

//-------------------------------------------------------------------------------------------------
ExchangeNexus::ExchangeNexus( const int system_count_in, const int total_swap_count,
                              const std::string &remd_type_in, const int frequency_swaps_in, 
                              const std::string &swap_store_in,
                              const std::string &temperature_dist_in,
                              const double exchange_probability_in, const double tolerance_in,
                              const int max_replicas_in, 
                              const double initial_temperature_in,
                              const double equilibrium_temperature_in ) :
    system_count{system_count_in}, total_swaps{total_swap_count},
    remd_type{remd_type_in}, frequency_swaps{frequency_swaps_in},
    swap_store{swap_store_in}, temperature_dist{temperature_dist_in},
    exchange_probability{exchange_probability_in}, tolerance{tolerance_in}, 
    max_replicas{max_replicas_in},
    initial_temperature{initial_temperature_in}, 
    equilibrium_temperature{equilibrium_temperature_in}
{}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getPotentialEnergies(){
  return sc->reportInstantaneousStates(StateVariable::POTENTIAL_ENERGY);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getKineticEnergies(){
  return sc->reportInstantaneousStates(StateVariable::KINETIC);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getHamiltonian(){
  
    std::vector<double> potential_energy = getPotentialEnergies();
    std::vector<double> kinetic_energy = getKineticEnergies();
    std::vector<double> hamiltonian;
    if (kinetic_energy.size() >= potential_energy.size()) {
        for (size_t i = 0; i < kinetic_energy.size(); ++i) {
            hamiltonian.push_back(kinetic_energy[i] + potential_energy[i]);
        }
    }
    // TODO: Else condition, as well as making sure vector indices refer to same particle
    return hamiltonian;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getTempDistribution() {

  // TODO: Implement rigidity getters and setters and get this data out of AtomGraph and &dynamics
  const AtomGraph* cur_ag = ag->getSystemTopologyPointer(0);
  std::vector<double> temperatures = vanDerSpoel(exchange_probability, initial_temperature, 
                                                 equilibrium_temperature, 
                                                 cur_ag->getRigidWaterCount(),
                                                 cur_ag->getOrganicMoleculeCount(), tolerance, 
                                                 0, 0, 0, 0, 0, ExceptionResponse::WARN);
  return temperatures;
}

//-------------------------------------------------------------------------------------------------
double ExchangeNexus::getProbabilityAtReplica(const int replica, bool flag){
  // Check for energy, return back a probability
  // If a replica == phasespace; return current energies
  // Else, recalculate
  std::vector<double> hamiltonians = getHamiltonian();
  if(!flag) {
    double cur_energy = sc->reportInstantaneousStates(StateVariable::KINETIC, replica);
    return cur_energy;  
  }
  if(replica == max_replicas){
    return sc->reportInstantaneousStates(StateVariable::KINETIC, replica);
  }
  return sc->reportInstantaneousStates(StateVariable::KINETIC, replica+1);
}

//-------------------------------------------------------------------------------------------------
void ExchangeNexus::swapReplicas(int replica_a, double a_ener, int replica_b, double b_ener) {
  sc->contribute(StateVariable::KINETIC, b_ener, replica_a);
  sc->contribute(StateVariable::KINETIC, a_ener, replica_b);
}

//-------------------------------------------------------------------------------------------------
void ExchangeNexus::attemptSwaps(std::vector<double> cur_energy, std::vector<double> final_energy,
																 bool flag) {
    for (int i = 0; i < cur_energy.size(); ++i) {
        bool performSwap = false;
        int swapIndex = -1;

        // If flag is true
        if (flag) {
            if (cur_energy[i] < final_energy[i]) {
                if (i % 2 == 0 && i + 1 < cur_energy.size()) {
                    // Even i, attempt swap with i + 1
                    swapIndex = i + 1;
                    performSwap = true;
                } else if (i % 2 != 0 && i - 1 >= 0) {
                    // Odd i, attempt swap with i - 1
                    swapIndex = i - 1;
                    performSwap = true;
                }
            }
        }
        // If flag is false (inverse logic)
        else {
            if (cur_energy[i] >= final_energy[i]) {
                if (i % 2 != 0 && i + 1 < cur_energy.size()) {
                    // Odd i, attempt swap with i + 1
                    swapIndex = i + 1;
                    performSwap = true;
                } else if (i % 2 == 0 && i - 1 >= 0) {
                    // Even i, attempt swap with i - 1
                    swapIndex = i - 1;
                    performSwap = true;
                }
            }
        }

        // If swap is possible, invoke swap and record it
        if (performSwap && swapIndex != -1) {
            swapReplicas(i, final_energy[i], swapIndex, final_energy[swapIndex]);

            // Record the swap information
            SwapRecord record;
            record.index1 = i;
            record.index2 = swapIndex;
            record.successful = true;
            swapHistory.push_back(record);
        }
    }
}

//-------------------------------------------------------------------------------------------------
void ExchangeNexus::setPhaseSpaceSynthesis(const PhaseSpaceSynthesis* ps_in){
  ps = const_cast<PhaseSpaceSynthesis*>(ps_in);
}

//-------------------------------------------------------------------------------------------------
void ExchangeNexus::setAtomGraphSynthesis(const AtomGraphSynthesis* ag_in){
  ag = const_cast<AtomGraphSynthesis*>(ag_in);
}

//-------------------------------------------------------------------------------------------------
void ExchangeNexus::setScoreCard(const ScoreCard* sc_in){
  sc = const_cast<ScoreCard*>(sc_in);
}

} // namespace sampling
} // namespace stormm
