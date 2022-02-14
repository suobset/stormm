// -*-c++-*-
#ifndef OMNI_NML_RANDOM_H
#define OMNI_NML_RANDOM_H

#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

/// \brief Default values for molecular dynamics
/// \{
constexpr int default_random_seed    = 30965871;
constexpr int default_random_streams = 1048576;
constexpr int default_random_stride  = 64;
constexpr int default_random_warmup  = 96;
/// \}
  
/// \brief Produce a namelist for specifying random number generation protocols.  This subsumes the
///        igseed keyword found in the &cntrl namelist of the Amber sander program, and much more.
///        OMNI random number streams created by any GPU can be reproduced exactly on the CPU.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
NamelistEmulator randomInput(const TextFile &tf, int *start_line);

/// \brief Check the random number seed to ensure that it is not too small.
///
/// \param igseed         The provided random seed
/// \param warmup_cycles  Number of warmup cycles slated to process when creating the first
///                       random number stream generator (modified if necessary and returned)
void validateRandomSeed(int igseed, int *warmup_cycles);

} // namespace namelist
} // namespace omni

#endif
