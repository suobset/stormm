// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_CONSTANTS_H
#define STORMM_ATOMGRAPH_CONSTANTS_H

#include "copyright.h"

namespace stormm {
namespace topology {

/// \brief Use a rounding tolerance to decide whether to smear a uniform charge over all particles
///        in the system in order to achieve an integer charge on the system as a whole.
constexpr double default_charge_rounding_tol = 0.001;

/// \brief Use the lowest bit of a 32-bit float representing the range [1.0, 2.0) as the default
///        input for the smoothCharges function, ensuring that all charges in the range (-2e, +2e)
///        will be expressed in increments that a floating point number can represent exactly.
///        The hexfloat would read 0x1p-23.
constexpr double default_charge_precision_inc = 1.1920928955078125E-7;

/// \brief Default 1:4 non-bonded screening factors are taken from various major codes.
/// \{
constexpr double amber_default_elec14_screen = 1.2;
constexpr double amber_default_vdw14_screen = 2.0;
constexpr double charmm_default_elec14_screen = 1.0;
constexpr double charmm_default_vdw14_screen = 1.0;
constexpr double glycam_default_elec14_screen = 1.0;
constexpr double glycam_default_vdw14_screen = 1.0;
/// \}
  
} // namespace topology
} // namespace stormm

#endif
