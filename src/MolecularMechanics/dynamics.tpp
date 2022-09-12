// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Scorecard dynamics(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, Tforce* xfrc, Tforce* yfrc,
                   Tforce* zfrc, Tcoord* xnxt, Tcoord* ynxt, Tcoord* znxt, Thermostat *heat_bath,
                   const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                   const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar, const ConstraintKit<Tcalc> &bcr,
                   const VirtualSiteKit<Tcalc> &vsk, const StaticExclusionMaskReader &ser,
                   const DynamicsControls &dyncon, const int nrg_scale_bits,
                   const Tcalc gpos_factor, const Tcalc force_factor) {

  // Loop for the requested number of cycles.  One calculation
  Scorecard sc(1, dyncon.getStepCount() + 1, nrg_scale_bits);
  for (int step = 0; step < dyncon.getStepCount(); step++) {

    // Evaluate the force and energy for a system in vacuum with isolated boundary conditions
    evalNonbValeRestMM<Tcoord, Tforce,
                       Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                              UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk, nbk,
                                              ser, rar, EvaluateForce::YES, 0, step);
    transmitVirtualSiteForces<Tcalc, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                            UnitCellType::NONE, vsk);

    // Commit the energy, all components (energy computations are obligatory in CPU functions)
    sc.commit(StateVariable::ALL_STATES);
    sc.incrementSampleCount();

    // Move particles, placing them in the {x,y,z}nxt arrays.
    

    // Apply constraints
  }
}
