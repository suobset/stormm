// -*-c++-*-
#ifndef OMNI_NML_RST_H
#define OMNI_NML_RST_H

#include "input.h"
#include "namelist_emulator.h"
#include "Parsing/textfile.h"

namespace omni {
namespace namelist {

/// \brief Produce a namelist for specifying basic input and output files, which can take the place
///        of a great deal of command line input in the Amber pmemd and sander programs.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this will
///                    wrap back to the beginning of the file in search of a unique &files
///                    namelist)
NamelistEmulator filesInput(const TextFile &tf, int *start_line);

} // namespace namelist
} // namespace omni

#endif
