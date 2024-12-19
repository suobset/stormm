#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "namelist_common.h"

namespace stormm {
namespace namelist {

using energy::translateVdwSumMethod;
using parse::realToString;
using parse::NumberFormat;
  
//-------------------------------------------------------------------------------------------------
void addRangedInteractionControls(NamelistEmulator *t_nml) {
  t_nml->addKeyword("cut", NamelistType::REAL, std::to_string(default_van_der_waals_cutoff));
  t_nml->addKeyword("elec_cut", NamelistType::REAL, std::to_string(default_electrostatic_cutoff));
  t_nml->addKeyword("vdw_cut", NamelistType::REAL, std::to_string(default_van_der_waals_cutoff));
  t_nml->addKeyword("vdw_tail", NamelistType::STRING, std::string(default_vdw_cutoff_style));
  t_nml->addHelp("cut", "The inter-particle distance at which to begin neglecting pairwise, "
                 "particle-particle interactions, in units of Angstroms.  If given, this unifying "
                 "parameter will take precedence over separate specifications for the "
                 "electrostatic or Lennard-Jones (van-der Waals) cutoffs.");
  t_nml->addHelp("elec_cut", "The inter-particle distance at which to begin discounting "
                 "electrostatic interactions, in units of Angstroms.  This applies to all methods "
                 "for evaluating the electrostatic potential.");
  t_nml->addHelp("vdw_cut", "The inter-particle distance at which to begin discounting "
                 "van-der Waals interactions, in units of Angstroms.  This applies to all methods "
                 "for evaluating the van-der Waals potential).");
  t_nml->addHelp("vdw_tail", "The manner in which van-der Waals interactions vanish at the "
                 "chosen cutoff.  Valid choices include \"cutoff\", \"truncation\", \"smooth\", "
                 "\"cubic\", \"infinite\", or \"pme\".");
}

//-------------------------------------------------------------------------------------------------
void addRangedInteractionInterpretation(double *electrostatic_cutoff, double *van_der_waals_cutoff,
                                        VdwSumMethod *vdw_style, const NamelistEmulator &t_nml,
                                        const ExceptionResponse policy) {
  t_nml.assignVariable(electrostatic_cutoff, "elec_cut");
  t_nml.assignVariable(van_der_waals_cutoff, "vdw_cut");
  t_nml.assignVariable(electrostatic_cutoff, "cut");
  t_nml.assignVariable(van_der_waals_cutoff, "cut");
  try {
    *vdw_style = translateVdwSumMethod(t_nml.getStringValue("vdw_tail"));
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The van-der Waals tail style \"" + t_nml.getStringValue("vdw_style") + "\" does not "
            "translate into any known method.", "addRangedInteractionInterpretation");
    case ExceptionResponse::WARN:
      rtWarn("The van-der Waals tail style \"" + t_nml.getStringValue("vdw_style") + "\" does not "
             "translate into any known method.  The van-der Waals interactions will be truncated "
             "at the cutoff.", "addRangedInteractionInterpretation");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *vdw_style = VdwSumMethod::CUTOFF;
  }

  // Check additional inputs
  if (*electrostatic_cutoff < minimum_elec_cutoff) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The electrostatic cutoff (" +
            realToString(*electrostatic_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ") cannot be "
            "set below a minimum value of " +
            realToString(minimum_elec_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ".",
            "addRangedInteractionInterpretation");
    case ExceptionResponse::WARN:
      rtWarn("The electrostatic cutoff (" +
             realToString(*electrostatic_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ") cannot "
             "be set below a minimum value of " +
             realToString(minimum_elec_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ".  This "
             "minimum value will be applied.", "addRangedInteractionInterpretation");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *electrostatic_cutoff = minimum_elec_cutoff;
  }
  if (*van_der_waals_cutoff < minimum_vdw_cutoff) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A van-der Waals cutoff of " +
            realToString(*van_der_waals_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " is too "
            "short to permit accurate force and energy computations.", "DynamicsControls",
            "addRangedInteractionInterpretation");
    case ExceptionResponse::WARN:
      rtWarn("A van-der Waals cutoff of " +
             realToString(*van_der_waals_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " is too "
             "short to permit accurate force and energy computations.  The minimum value of " +
             realToString(minimum_vdw_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " will be "
             "applied.", "addRangedInteractionInterpretation");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *van_der_waals_cutoff = minimum_vdw_cutoff;
  }
}

} // namespace namelist
} // namepace stormm
