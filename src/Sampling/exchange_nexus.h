// -*-c++-*-
#ifndef STORMM_EXCHANGE_NEXUS_H
#define STORMM_EXCHANGE_NEXUS_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
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

using chemistry::PhaseSpaceSynthesis;
using energy::ScoreCard;
using energy::StateVariable;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
using trajectory::TrajectoryKind;

struct SwapRecord {
  int index1;
  int index2;
  bool successful;
};

/// \brief
class ExchangeNexus {
public:

  /// \brief Main constructor that takes in all details, plus an existing PhaseSpaceSynthesis,
  ///        AtomGraphSynthesis, and a ScoreCard object.
  ExchangeNexus(  int system_count_in, int total_swap_count,
                  const std::string &remd_type_in, int frequency_swaps_in, 
                  const std::string &swap_store_in, const std::string &temperature_dist_in,
                  double exchange_probability_in, double tolerance_in, int max_replicas_in, 
                  double initial_temperature_in, double equilibrium_temperature_in,
                  const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                  const ScoreCard *sc_in);

  /// \brief Constructor that takes in all data, but no *Synthesis objects.
  ExchangeNexus(  int system_count_in, int total_swap_count,
                  const std::string &remd_type_in, int frequency_swaps_in, 
                  const std::string &swap_store_in, const std::string &temperature_dist_in,
                  double exchange_probability_in, double tolerance_in, int max_replicas_in, 
                  double initial_temperature_in, double equilibrium_temperature_in);
  

  /// \brief Copy and Move Constructors
  /// \{
  ExchangeNexus(const ExchangeNexus &original) = default;
  ExchangeNexus(ExchangeNexus &&other) = default;
  /// /}
  
  /// \brief Function to return the potential energies of all particles in the system
  ///         in a vector.
  std::vector<double> getPotentialEnergies();
  
  /// \brief Function to return the kinetic energies of all particles in the system
  ///         in a vector.
  std::vector<double> getKineticEnergies();
  
  /// \brief Function to calculate the Hamiltonian of each particle in the current system
  ///
  std::vector<double> getHamiltonian();
  
  
  /// \brief Function to return the temperature distributions by invoking one of
  ///         temp_distributions.h and passing the correct algorithm in.
  ///
  ///  \param  initial_temperature      The temperature at the bottom of the range in which a REMD
  ///                                   simulation has to be conducted in.
  ///  \param  equilibrium_temperature  The temperature at the top of the range in which a REMD
  ///                                   simulation has to be conducted in.
  ///  \param  temperature_dist         The algorithm to follow when calculating a temperature
  ///                                   distribution.
  ///  \param  exchange_probability     The user input probability at which to perform a swap.
  ///  \param  cur_ag                   The AtomGraph system for the current REMD.
  std::vector<double> getTempDistribution();

  /// \brief Function to return the probability of a system at a given replica.
  ///
  /// \param ps       The PhaseSpace for which we want to check the probability
  /// \param replica  The replica at which we want to check the probability
  double getProbabilityAtReplica(const int replica, bool flag);

  /// \brief Function to attempt swapping of replicas from one index to another
  ///
  /// \param replica_a      The first of the two replicas to swap
  /// \param a_ener         The energy of the CURRENT pre-swap replica_a
  /// \param replica_b      The second of the two replicas to swap
  /// \param b_ener         The energy of the CURRENT pre-swap replica b
  void swapReplicas(int replica_a, double a_ener, int replica_b, double b_ener);

  /// \brief Function to attempt replica swaps based on the current and final states
  ///        of each replica
  ///
  /// \param cur_energy  		Vector with the energies in the current state
  /// \param final_energy   Vector with the energies in the swapped state
  /// \param flag           Boolean value to indicate swapping of odd or even replicas
  void attemptSwaps(std::vector<double> cur_energy, std::vector<double> final_energy, bool flag);
  
  /// \brief Set a PhaseSpaceSynthesis in this object if it was created without one
  void setPhaseSpaceSynthesis(const PhaseSpaceSynthesis* ps_in);

  /// \brief Set an AtomGraphSynthesis in this object if it was created without one
  void setAtomGraphSynthesis(const AtomGraphSynthesis* ag_in);

  /// \brief Set a ScoreCard in this object if it was created without one
  void setScoreCard(const ScoreCard* sc_in);
  
private:  
  int system_count;
  int total_swaps;
  std::string remd_type;
  int frequency_swaps;
  int max_replicas;
  std::string swap_store;
  std::string temperature_dist;
  double initial_temperature;
  double equilibrium_temperature;
  double exchange_probability;
  double tolerance;
  PhaseSpaceSynthesis *ps;
  AtomGraphSynthesis *ag;
  ScoreCard *sc;
  std::vector<SwapRecord> swapHistory;
};

/// \brief A class to conduct Replica Exchange Molecular Dynamics; taking in constraints
///        across multiple simulations. Takes in data from the current running synthesis
///        and user inputs regarding temperatures, probabilities, and temperature distribution
///        algorithms.

} // namespace sampling
} // namespace stormm

#endif
