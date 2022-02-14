// -*-c++-*-
#ifndef OMNI_NML_DYNAMICS_H
#define OMNI_NML_DYNAMICS_H

#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

/// \brief Default values for molecular dynamics
/// \{
constexpr int default_dynamics_nstlim = 100;
/// \}
  
/// \brief Produce a namelist for specifying molecular dynamics directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the energy
///        minimization input.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line);

} // namespace namelist
} // namespace omni

#endif
