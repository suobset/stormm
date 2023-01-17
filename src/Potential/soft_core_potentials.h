// -*-c++-*-
#ifndef STORMM_SOFT_CORE_POTENTIALS_H
#define STORMM_SOFT_CORE_POTENTIALS_H

#include <cmath>
#include "copyright.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace energy {

/// \brief Compute the electrostatic contribution for an interaction damped to become a quadratic
///        function with a relative maximum at -1.0 Angstroms (the location of the maximum ensures
///        that the function derivative has a finite, nonzero value even as the inter-particle
///        distance approaches zero).  If the inter-particle distance is too close to zero, no
///        force will be logged.
///
///
/// \param r               Distance between the two particles
/// \param clash_distance  Absolute distance between any two particles, below which they would be
///                        declared to clash
/// \param qiqj            Product of the two particles' charges and Coulomb's constant, attenuated
///                        by any other prefactors such as 1:4 scaling terms.
/// \param ele_contrib     The energy contribution from this specific interaction.  This will be
///                        accumulated and returned.
/// \param fmag            The magnitude of the force between two particles.  This will be
///                        accumulated and returned.  Providing a nullptr will cause the force
///                        calculation to be skipped.
template <typename Tcalc>
void quadraticCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                                 Tcalc *ele_contrib, Tcalc *fmag);

/// \brief Compute the Lennard-Jones contribution for an interactions damped to become a quartic
///        function with a relative maximum at -1.0 Angstroms (the location of the maximum ensures
///        that the function derivative has a finite, nonzero value even as the inter-particle
///        distance approaches zero).
///
/// \param r            Distance between the two particles
/// \param clash_ratio  The minimum ratio of the distance between any two particles and the pair's
///                     Lennard-Jones sigma parameter.  If the inter-particle distance is too low,
///                     the soft-core potential will engage.
/// \param lja          The Lennard-Jones B parameter, pre-scaled by any attenuation constants such
///                     as 1:4 prefactors
/// \param ljb          The Lennard-Jones B parameter, pre-scaled by any attenuation constants
/// \param vdw_contrib  The energy contribution from this specific interaction.  This will be
///                     accumulated and returned.
/// \param fmag         The magnitude of the force between two particles.  This will be accumulated
///                     and returned.  Providing a nullptr will skip the force calculation.
template <typename Tcalc>
void quarticCoreLennardJones(const Tcalc r, const Tcalc clash_ratio, const Tcalc lja,
                             const Tcalc ljb, Tcalc *vdw_contrib, Tcalc *fmag);
                             
} // namespace energy
} // namespace stormm

#include "soft_core_potentials.tpp"

#endif
