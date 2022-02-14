// -*-c++-*-
#ifndef OMNI_NML_SOLVENT_H
#define OMNI_NML_SOLVENT_H

#include "input.h"
#include "namelist_emulator.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph_enumerators.h"

namespace omni {
namespace namelist {

using topology::ImplicitSolventModel;

/// \brief Default values for the implicit solvent model
/// \{
constexpr int default_solvent_igb = 0;
constexpr double default_solvent_rgbmax = 0.0;
constexpr double default_solvent_intdiel = 1.0;
constexpr double default_solvent_extdiel = 78.5;
constexpr double default_solvent_saltcon = 0.0;
/// \}

/// \brief Translate the numerical input for the implicit solvent model into one of the recognized
///        models available in OMNI.
///
/// \param igb_val  Numerical value of the implicit solvent model, as read from the parent namelist
ImplicitSolventModel extractImplicitSolventModel(int igb_val);

/// \brief Produce a namelist for defining the implicit solvent model, replicating various inputs
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the
///        molecular dynamics input as well as the Particle-Mesh Ewald namelist.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
NamelistEmulator solventInput(const TextFile &tf, int *start_line);

} // namespace namelist
} // namespace omni

#endif
