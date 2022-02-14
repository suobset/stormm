// -*-c++-*-
#ifndef OMNI_NML_MINIMIZE_H
#define OMNI_NML_MINIMIZE_H

#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

/// \brief Default values for energy minimization
/// \{
constexpr int default_minimize_maxcyc  = 200;
constexpr int default_minimize_ncyc    = 50;
constexpr double default_minimize_dx0  = 0.01;
constexpr double default_minimize_drms = 0.0001;
/// \}
  
/// \brief Produce a namelist for specifying energy minimization directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the
///        molecular dynamics input, obviating the need for the imin setting found in the general
///        &cntrl namelist of sander and pmemd.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line);

} // namespace namelist
} // namespace omni

#endif
