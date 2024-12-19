// -*-c++-*-
#ifndef STORMM_NAMELIST_COMMON_H
#define STORMM_NAMELIST_COMMON_H

#include "copyright.h"
#include "Potential/energy_enumerators.h"
#include "namelist_emulator.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

/// \brief The default cutoffs for electrostatic and van-der Waals interactions.  These apply only
///        in the case of periodic dynamics and will be superceded by information in &pppm
///        namelists if present in the user input.
/// \{
constexpr double default_electrostatic_cutoff = 8.0;
constexpr double default_van_der_waals_cutoff = 10.0;
/// \}

/// \brief The minimum cutoffs for electrostatic and van-der Waals interactions.  These apply only
///        in the case of periodic dynamics and will be superceded by information in &pppm
///        namelists if present in the user input.
/// \{
constexpr double minimum_elec_cutoff = 0.0;
constexpr double minimum_vdw_cutoff = 4.5;
/// \}

/// \brief The default method for making van-der Waals interactions vanish at the cutoff.
constexpr char default_vdw_cutoff_style[] = "cutoff";

using energy::VdwSumMethod;
  
/// \brief A free function to encapsulate the addition of common cutoff keywords to a namelist
///        emulator object.  This function will add the following keywords to the control block:
///        cut, elec_cut, vdw_cut, vdw_style.  The specific cutoffs for eletrostatic or van-der
///        Waals interactions will override a general setting for either cutoff, if including the
///        function addRangedInterpretation(), below.  The vdw_style keyword works off of the
///        VdwSumMethod enumerator.
///
/// \param t_nml  The namelist object to modify
void addRangedInteractionControls(NamelistEmulator *t_nml);
 
/// \brief A free function to encapsulate the interprertation of common cutoff keywords to a
///        namelist emulator object: cut, elec_cut, vdw_cut, vdw_style.  See the above function
///        addRangedInteractionControls() for further details.
///
/// \param electrostatic_cutoff  The electrostatic cutoff for the simulation
/// \param van_der_waals_cutoff  The van-der Waals cutoff for the simulation
/// \param vdw_style             Manner in which to make van-der Waals interactions vanish at the
///                              chosen cutoff
/// \param t_nml                 The namelist object to interpret
/// \param policy                Course of action to take if bad input is encountered
void addRangedInteractionInterpretation(double *electrostatic_cutoff, double* van_der_waals_cutoff,
                                        VdwSumMethod *vdw_style, const NamelistEmulator &t_nml,
                                        ExceptionResponse policy);

} // namelist namelist
} // namepace stormm

#endif
