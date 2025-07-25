// -*-c++-*-
#ifndef STORMM_INTEGRATION_H
#define STORMM_INTEGRATION_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Potential/cellgrid.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using energy::CellGrid;
using energy::CellGridReader;
using numerics::force_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using symbols::boltzmann_constant;
using symbols::boltzmann_constant_f;
using symbols::kcal_to_gafs;
using symbols::kcal_to_gafs_f;
using synthesis::AtomGraphSynthesis;
using synthesis::SyAtomUpdateKit;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;

/// \brief Advance particle velocities according to the velocity-Verlet first half update.
///
/// Overloaded:
///   - Provide a single system
///   - Provide a synthesis of systems
///   - Provide the original objects or their abstracts
///
/// \param xvel:  Velocities of all particles along the Cartesian X axis
/// \param yvel:  Velocities of all particles along the Cartesian Y axis
/// \param zvel:  Velocities of all particles along the Cartesian Z axis
/// \param xfrc:  Forces acting on all particles along the Cartesian X axis
/// \param yfrc:  Forces acting on all particles along the Cartesian Y axis
/// \param zfrc:  Forces acting on all particles along the Cartesian Z axis
/// \param psw:   Abstract for a single system's coordinates, containing positions, velocities, and
///               forces in double-precision real representations
/// \param ps:    A single system's coordinates, containing positions, velocities, and forces in
///               double-precision real representations
/// \param cdk:   Contains masses of all particles
/// \param ag:    Topology details for a single system
/// \param tstr:  Abstract for the thermostat
/// \param tst:   Thermostat configured for the system or systems of interest
/// \{
template <typename Tcoord, typename Tcalc>
void velocityVerletVelocityUpdate(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                  const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                  int natom, const Tcalc* masses, Tcoord* vxalt, Tcoord* vyalt,
                                  Tcoord* vzalt, const ThermostatReader<Tcalc> &tstr,
                                  const int* xvel_ovrf = nullptr, const int* yvel_ovrf = nullptr,
                                  const int* zvel_ovrf = nullptr, const int* xfrc_ovrf = nullptr,
                                  const int* yfrc_ovrf = nullptr, const int* zfrc_ovrf = nullptr,
                                  int* vxalt_ovrf = nullptr, int* vyalt_ovrf = nullptr,
                                  int* vzalt_ovrf = nullptr, int atom_offset = 0,
                                  Tcalc vel_scale_factor = 1.0, Tcalc frc_scale_factor = 1.0);
  
void velocityVerletVelocityUpdate(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                                  const ThermostatReader<double> &tstr);

void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst);

void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst);

template <typename T, typename T2, typename T4>
void velocityVerletVelocityUpdate(PsSynthesisWriter *poly_psw,
                                  const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                                  const ThermostatReader<T> &tstr);

void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                                  const Thermostat *tst,
                                  PrecisionModel prec = PrecisionModel::SINGLE);

void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                                  const Thermostat &tst,
                                  PrecisionModel prec = PrecisionModel::SINGLE);

template <typename Tcoord, typename Tacc, typename Tcoord4> 
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg,
                                  const AtomGraphSynthesis &poly_ag, const Thermostat &tst,
                                  PrecisionModel prec = PrecisionModel::SINGLE);

template <typename Tcoord, typename Tacc, typename Tcoord4> 
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg_qq,
                                  const CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> &cg_lj,
                                  const AtomGraphSynthesis &poly_ag, const Thermostat &tst,
                                  PrecisionModel prec = PrecisionModel::SINGLE);
/// \}

/// \brief Advance particle velocities according to the velocity-Verlet second half update, then
///        advance coordinates according to the new velocities.  Overloading and descriptions of
///        input parameters follow from velocityVerletVelocityUpdate(), above, in addition to:
///
/// \param xcrd               Original positions of all particles along the Cartesian X axis
/// \param ycrd               Original positions of all particles along the Cartesian Y axis
/// \param zcrd               Original positions of all particles along the Cartesian Z axis
/// \param xalt               Developing positions of all particles along the X axis
/// \param yalt               Developing positions of all particles along the Y axis
/// \param zalt               Developing positions of all particles along the Z axis
/// \param gpos_scale_factor  Scaling factor for converting global positions to the fixed-precision
///                           representation
/// \param vel_scale_factor   Scaling factor for covnerting particle velocities to the
///                           fixed-precision representation
/// \param frc_scale_factor   Scaling factor for accumulating (and applying) forces in
///                           fixed-precision
/// \{
template <typename Tcoord, typename Tcalc>
void velocityVerletCoordinateUpdate(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                    const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                    int natom, const Tcalc* masses, Tcoord* xalt, Tcoord* yalt,
                                    Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt, Tcoord* vzalt,
                                    const ThermostatReader<Tcalc> &tstr,
                                    const int* xcrd_ovrf = nullptr, const int* ycrd_ovrf = nullptr,
                                    const int* zcrd_ovrf = nullptr, const int* xfrc_ovrf = nullptr,
                                    const int* yfrc_ovrf = nullptr, const int* zfrc_ovrf = nullptr,
                                    int* xalt_ovrf = nullptr, int* yalt_ovrf = nullptr,
                                    int* zalt_ovrf = nullptr, int* vxalt_ovrf = nullptr,
                                    int* vyalt_ovrf = nullptr, int* vzalt_ovrf = nullptr,
                                    int atom_offset = 0, Tcalc gpos_scale_factor = 1.0,
                                    Tcalc vel_scale_factor = 1.0, Tcalc frc_scale_factor = 1.0);
  
void velocityVerletCoordinateUpdate(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                                    const ThermostatReader<double> &tstr);

void velocityVerletCoordinateUpdate(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst);

void velocityVerletCoordinateUpdate(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst);

template <typename T, typename T2, typename T4>
void velocityVerletCoordinateUpdate(PsSynthesisWriter *poly_psw,
                                    const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                                    const ThermostatReader<T> &tstr);

void velocityVerletCoordinateUpdate(PhaseSpaceSynthesis *poly_ps,
                                    const AtomGraphSynthesis *poly_ag, const Thermostat *tst);

void velocityVerletCoordinateUpdate(PhaseSpaceSynthesis *poly_ps,
                                    const AtomGraphSynthesis &poly_ag, const Thermostat &tst);
/// \}

} // namespace structure
} // namespace stormm

#include "integration.tpp"

#endif
