#include "Topology/atomgraph_enumerators.h"
#include "namelist_element.h"
#include "nml_solvent.h"

namespace omni {
namespace namelist {

using topology::ImplicitSolventModel;
  
//-------------------------------------------------------------------------------------------------
ImplicitSolventModel extractImplicitSolventModel(const int igb_val) {
  switch (igb_val) {
  case 0:
  case 6:
    return ImplicitSolventModel::NONE;
  case 1:
    return ImplicitSolventModel::HCT_GB;
  case 2:
    return ImplicitSolventModel::OBC_GB;
  case 5:
    return ImplicitSolventModel::OBC_GB_II;
  case 7:
    return ImplicitSolventModel::NECK_GB;
  case 8:
    return ImplicitSolventModel::NECK_GB_II;
  default:
    rtErr("Unrecognized implicit solvent model, igb = " + std::to_string(igb_val) + ".",
          "extractImplicitSolventModel");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator solventInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("solvent", CaseSensitivity::AUTOMATIC, ExceptionResponse::DIE,
                         "Wraps directives needed to define the solvent model of a molecular "
                         "mechanical system.  This is not the same as the explicit water model in "
                         "use, and generally refers to implicit solvent methods.");
  t_nml.addKeyword(NamelistElement("igb", NamelistType::INTEGER,
                                   std::to_string(default_solvent_igb)));
  t_nml.addKeyword(NamelistElement("intdiel", NamelistType::REAL,
                                   std::to_string(default_solvent_intdiel)));
  t_nml.addKeyword(NamelistElement("extdiel", NamelistType::REAL,
                                   std::to_string(default_solvent_extdiel)));
  t_nml.addKeyword(NamelistElement("saltcon", NamelistType::REAL,
                                   std::to_string(default_solvent_saltcon)));
  t_nml.addKeyword(NamelistElement("rgbmax", NamelistType::REAL,
                                   std::to_string(default_solvent_rgbmax)));
  t_nml.addHelp("igb", "Definition of the Generalized Born model to use.  Numerical settings "
                "follow those found in the eponymous keyword of Amber sander's &cntrl namelist: "
                "0 or 6 (no Generalized Born model is used), 1 (Hawkins, Cramer, Truhlar), 2 or 5 "
                "(distinct models contributed by Onufriev, Bashford, and Case), 7 or 8 (distinct "
                "models contributed by Mongan, the \"neck\" GB method).");
  t_nml.addHelp("intdiel", "Internal dielectric constant.  No value other than the default of "
                "1.0 has been well tested.");
  t_nml.addHelp("extdiel", "External dielectric constant for the surrounding continuum solvent.");
  t_nml.addHelp("saltcon", "Concentration of 1-1 monovalent ions in the continuum solvent.");
  t_nml.addHelp("rgbmax", "Maximum distance at which two particles can contribute to each other's "
                "effective Born radii.  Not to be confused with the non-bonded cutoff, which is "
                "the distance past which actual interactions between particles are discarded.");
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.  All calls to this function should
  // proceed in consecutive calls, to make use of the updates to start_line and avoid reading any
  // instance of this namelist twice or skipping instances of it in the search for some other
  // namelist.  An alternative is to keep an independent counter to track progress through the
  // input file in search for &rst namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::NO, tf.getLineCount());

  // Checks on the input
  extractImplicitSolventModel(t_nml.getIntValue("igb"));
  
  return t_nml;
}

} // namespace namelist
} // namespace omni
