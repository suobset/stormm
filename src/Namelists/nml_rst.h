// -*-c++-*-
#ifndef OMNI_NML_RST_H
#define OMNI_NML_RST_H

#include "input.h"
#include "namelist_emulator.h"
#include "Parsing/textfile.h"

namespace omni {
namespace namelists {

using parse::NamelistEmulator;

/// \brief Produce a namelist for specifying an NMR restraint, equivalent to the eponymous namelist
///        in sander of pmemd.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
NamelistEmulator rstInput(const TextFile &tf, int *start_line);

} // namespace namelists
} // namespace omni

#endif
