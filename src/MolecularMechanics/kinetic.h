// -*-c++-*-
#ifndef STORMM_KINETIC_H
#define STORMM_KINETIC_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace mm {

using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using energy::ScoreCard;
using energy::StateVariable;
using numerics::velocity_scale_nonoverflow_bits;
using structure::ApplyConstraints;
using symbols::boltzmann_constant;
using symbols::boltzmann_constant_f;
using symbols::gafs_to_kcal;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateCycle;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
using trajectory::Thermostat;
using trajectory::ThermostatReader;

/// \brief Compute the kinetic energy of a simulation.
///
/// Overloaded:
///   - Supply coordinates in various representations
///   - Provide const pointers to original objects, or pass by const references
///
/// \param xvel              Velocities of all particles along the Cartesian X axis
/// \param yvel              Velocities of all particles along the Cartesian Y axis
/// \param zvel              Velocities of all particles along the Cartesian Z axis
/// \param xvel_ovrf         Overflow bits for extreme, 95-bit fixed-precision representations of
///                          particle velocities along the Cartesian X axis
/// \param yvel_ovrf         Overflow bits for particle velocities along the Cartesian Y axis
/// \param zvel_ovrf         Overflow bits for particle velocities along the Cartesian Z axis
/// \param natom             The number of atoms in the system
/// \param psw               Contains the velocities of all particles.  It is expected that this
///                          abstract will have been taken such that the relevant velocities are
///                          found in the member array variables vxalt, vyalt, and vzalt, based on
///                          the first half velocity-Verlet integration time step being contributed
///                          on top of the velocities xvel, yvel, and zvel which were in effect as
///                          of the moment forces were computed.  This is the same principle seen
///                          in the Amber restart file, where the velocities are a half step
///                          "ahead" of the positions.
/// \param ps                Positions, velocities, and forces acting on all particles
/// \param orientation       Point in the time cycle at which the relevant velocities can be found
/// \param sc                The energy-tracking object, complete with the energy precision model
/// \param masses            Array of masses for all particles
/// \param cdk               Contains the mass of each particle in the simulation
/// \param poly_auk          Holds the masses of atoms in a synthesis of systems
/// \param ag                Topological details of the simulation, including masses of particles
/// \param system_index      The index of the system within sc for which to compute kinetic energy
/// \param nrg_scale_factor  Scaling factor determining the precision of energy accumulation
/// \param inv_vel_scale     Inverse velocity scaling factor taking fixed-precision velocities into
///                          Angstroms per femtosecond
/// \{
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                        const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                        const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                        Tcalc inv_vel_scale = 1.0);

template <typename Tcoord, typename Tcalc>
void evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel, ScoreCard *sc,
                       const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                       const ChemicalDetailsKit &cdk, int system_index = 0,
                       Tcalc inv_vel_scale = 1.0);
  
void evalKineticEnergy(const PhaseSpaceWriter *psw, ScoreCard *sc, const ChemicalDetailsKit &cdk,
                       int system_index);

void evalKineticEnergy(const PhaseSpaceWriter &psw, ScoreCard *sc, const ChemicalDetailsKit &cdk,
                       int system_index);

void evalKineticEnergy(const PhaseSpace *ps, CoordinateCycle orientation, ScoreCard *sc,
                       const AtomGraph *ag, int system_index = 0);

void evalKineticEnergy(const PhaseSpace &ps, CoordinateCycle orientation, ScoreCard *sc,
                       const AtomGraph &ag, int system_index = 0);

void evalKineticEnergy(const PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                       int system_index = 0);

void evalKineticEnergy(const PhaseSpace &ps, ScoreCard *sc, const AtomGraph &ag,
                       int system_index = 0);

template <typename Tmass, typename Tmass2, typename Tmass4, typename Tcalc>
void evalKineticEnergy(const PsSynthesisReader &poly_psr, ScoreCard *sc,
                       const SyAtomUpdateKit<Tmass, Tmass2, Tmass4> &poly_auk);
/// \}

/// \brief Compute the rescaled kinetic energy for a REMD Simulation
///        based on PhaseSpaceSynthesis
///
/// Overloaded:
///   - Supply coordinates in various representations
///   - Provide const pointers to original objects, or pass by const references
///
/// \param xvel              Velocities of all particles along the Cartesian X axis
/// \param yvel              Velocities of all particles along the Cartesian Y axis
/// \param zvel              Velocities of all particles along the Cartesian Z axis
/// \param xvel_ovrf         Overflow bits for extreme, 95-bit fixed-precision representations of
///                          particle velocities along the Cartesian X axis
/// \param yvel_ovrf         Overflow bits for particle velocities along the Cartesian Y axis
/// \param zvel_ovrf         Overflow bits for particle velocities along the Cartesian Z axis
/// \param masses            Array of masses for all particles
/// \param natom             The number of atoms in the system
/// \param nrg_scale_factor  Scaling factor determining the precision of energy accumulation
/// \param inv_vel_scale     Inverse velocity scaling factor taking fixed-precision velocities into
///                          Angstroms per femtosecond
/// \{
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalRescaledKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                                const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                                const Tcalc inv_vel_scale);
/// \}

/// \brief Compute the Particle Momentum & Velocity Rescalings, based on the Rescaled Kinetic
///        Energy Above
///
/// Overloaded:
///   - Supply coordinates in various representations
///   - Provide const pointers to original objects, or pass by const references
///
/// \param xvel              Velocities of all particles along the Cartesian X axis
/// \param yvel              Velocities of all particles along the Cartesian Y axis
/// \param zvel              Velocities of all particles along the Cartesian Z axis
/// \param xvel_ovrf         Overflow bits for extreme, 95-bit fixed-precision representations of
///                          particle velocities along the Cartesian X axis
/// \param yvel_ovrf         Overflow bits for particle velocities along the Cartesian Y axis
/// \param zvel_ovrf         Overflow bits for particle velocities along the Cartesian Z axis
/// \param masses            Array of masses for all particles
/// \param natom             The number of atoms in the system
/// \param inv_vel_scale     Inverse velocity scaling factor taking fixed-precision velocities into
///                          Angstroms per femtosecond
/// \{
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalMomenta(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                  const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                  const Tmass* masses, const int natom, const Tcalc inv_vel_scale);
/// \}

/// \brief Compute the temperature of the system.  It is expected that the kinetic energy will
///        have already been computed.
///
/// Overloaded:
///   - Provide an abstract or the original topology itself
///   - Operate on one system or many
///
/// Descriptions of input parameters follow from the foundational evalKineticEnergy() above, in
/// addition to:
///
/// \param ndof         The number of degrees of freedom in the system
/// \param cdk          In addition to masses, this object holds the number of degrees of freedom,
///                     with or without constraints, for the system as a whole
/// \param poly_auk     Also contains the number of degrees of freedom for each system in the
///                     synthesis of systems
/// \param use_cnst     Indicate whether geometric constraints are in effect (TRUE) or not (FALSE),
///                     critical for determining the number of degrees of freedom
/// \param ke_computed  Indicate whether the kinetic energy has already been computed in the energy
///                     tracking object's instantaneous buffers
/// \{
template <typename Tcoord, typename Tmass, typename Tcalc>
Tcalc computeTemperature(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                         const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                         const Tmass* masses, int natom, int ndof, Tcalc nrg_scale_factor,
                         Tcalc inv_vel_scale);

template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace *ps, const AtomGraph *ag, ApplyConstraints use_cnst);

template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace &ps, const AtomGraph &ag, ApplyConstraints use_cnst);

template <typename Tmass, typename Tmass2, typename Tmass4, typename Tcalc>
void computeTemperature(const PsSynthesisReader &poly_psr, ScoreCard *sc,
                        const SyAtomUpdateKit<Tmass, Tmass2, Tmass4> &poly_auk,
                        const ThermostatReader<Tcalc> &tstr, bool ke_computed = false);
  
double computeTemperature(const PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                          ApplyConstraints use_cnst, int system_index,
                          PrecisionModel prec = PrecisionModel::SINGLE);

double computeTemperature(const PhaseSpaceSynthesis &ps, const AtomGraphSynthesis &ag,
                          ApplyConstraints use_cnst, int system_index,
                          PrecisionModel prec = PrecisionModel::SINGLE);

void computeTemperature(ScoreCard *sc, const ChemicalDetailsKit &cdk, bool cnst, int index = 0);

void computeTemperature(ScoreCard *sc, const AtomGraph &ag, bool cnst, int index = 0);

void computeTemperature(ScoreCard *sc, const AtomGraph *ag, bool cnst, int index = 0);
/// \}

/// \brief Initiate the Momenta Rescale for a given REMD Simulation
///
/// \param PsSynthesis      The given PhaseSpace Synthesis for the current instance
/// \param AgSynthesis      The given AtomGraph Synthesis for the current instance
/// \param ScoreCard        Current ScoreCard tracking energies for all particles
/// \param swap_indices     Vector with indices that have been swapped in a REMD simulation
///                         for which the momenta have to be rescaled
/// \param temp             Vector with the temperature distributions during REMD simulation
template <typename Tcalc>
std::vector<llint> initiateMomentaRescale(const PhaseSpaceSynthesis *PsSynthesis, 
                                          const AtomGraphSynthesis *AgSynthesis, 
                                          const ScoreCard *sc, std::vector<int> swap_indices, 
                                          std::vector<double> temp);

/// \brief Initiate the Kinetic Energy Rescale for a given REMD Simulation
///
/// \param PsSynthesis      The given PhaseSpace Synthesis for the current instance
/// \param AgSynthesis      The given AtomGraph Synthesis for the current instance
/// \param ScoreCard        Current ScoreCard tracking energies for all particles
template <typename Tcalc>
std::vector<llint> initiateKineticEnergyRescale(const PhaseSpaceSynthesis *PsSynthesis, 
                                                const AtomGraphSynthesis *AgSynthesis, 
                                                const ScoreCard *sc);

} // namespace mm
} // namespace stormm

#include "kinetic.tpp"

#endif
