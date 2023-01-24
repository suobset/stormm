#include "copyright.h"
#include "Namelists/input.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "nml_conformer.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using namelist::InputStatus;
using namelist::NamelistElement;
using namelist::NamelistType;
using namelist::readNamelist;
using errors::rtErr;
using errors::rtWarn;
using parse::CaseSensitivity;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
using parse::WrapTextSearch;

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const ExceptionResponse policy_in,
                                     const WrapTextSearch wrap) :
    policy{policy_in}, core_atom_mask{std::string("")}, core_data_item{std::string("")},
    core_rk2{0.0}, core_rk3{0.0}, core_r2{0.0}, core_r3{0.0}, anchor_conformation{std::string("")},
    sample_chirality{false}, sample_cis_trans{false}, prevent_hbonds{false},
    running_states{default_conf_running_states},
    final_states{default_conf_final_states},
    rotation_samples{default_conf_rotation_samples},
    rotatable_bond_limit{default_conf_max_rotatable_bonds},
    max_seeding_attempts{default_conf_max_seeding_attempts},
    system_trials{default_conf_max_system_trials},
    rmsd_tolerance{default_conf_rmsd_tolerance}
{}

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const TextFile &tf, int *start_line, bool *found_nml,
                                     const ExceptionResponse policy_in,
                                     const WrapTextSearch wrap) :
    ConformerControls(policy_in)
{
  const NamelistEmulator t_nml = conformerInput(tf, start_line, found_nml, policy, wrap);
  if (t_nml.getKeywordStatus("core_mask") != InputStatus::MISSING) {
    core_atom_mask = t_nml.getStringValue("core_mask", "atoms");
    core_data_item = t_nml.getStringValue("core_mask", "data_item");

    // Get the rk2 restraint value
    if (t_nml.getKeywordStatus("core_mask", "rk2") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "rk2");
    }
    else if (t_nml.getKeywordStatus("core_mask", "repulsion") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "repulsion");
    }
    else if (t_nml.getKeywordStatus("core_mask", "stiffness") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "stiffness");
    }

    // Get the rk3 restraint value
    if (t_nml.getKeywordStatus("core_mask", "rk3") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "rk3");
    }
    else if (t_nml.getKeywordStatus("core_mask", "attraction") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "attraction");
    }
    else if (t_nml.getKeywordStatus("core_mask", "stiffness") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "stiffness");
    }

    // Get the r2 restraint value
    if (t_nml.getKeywordStatus("core_mask", "r2") != InputStatus::MISSING) {
      core_r2 = t_nml.getRealValue("core_mask", "r2");
    }
    else if (t_nml.getKeywordStatus("core_mask", "demand") != InputStatus::MISSING) {
      core_r2 = t_nml.getRealValue("core_mask", "demand");
    }

    // Get the r3 restraint value
    if (t_nml.getKeywordStatus("core_mask", "r3") != InputStatus::MISSING) {
      core_r3 = t_nml.getRealValue("core_mask", "r3");
    }
    else if (t_nml.getKeywordStatus("core_mask", "grace") != InputStatus::MISSING) {
      core_r3 = t_nml.getRealValue("core_mask", "grace");
    }
  }
  if (t_nml.getKeywordStatus("anchor_conf") != InputStatus::MISSING) {
    anchor_conformation = t_nml.getStringValue("anchor_conf");
  }
  sample_chirality = strcmpCased(t_nml.getStringValue("sample_chirality"), "true",
                                 CaseSensitivity::NO);
  sample_cis_trans = strcmpCased(t_nml.getStringValue("sample_cis_trans"), "true",
                                 CaseSensitivity::NO);
  prevent_hbonds = strcmpCased(t_nml.getStringValue("prevent_hbonds"), "true",
                               CaseSensitivity::NO);
  running_states = t_nml.getIntValue("running_states");
  final_states = t_nml.getIntValue("final_states");
  rotation_samples = t_nml.getIntValue("rotation_samples");
  rotatable_bond_limit = t_nml.getIntValue("max_rotatable_bonds");
  max_seeding_attempts = t_nml.getIntValue("max_seeding_attempts");
  system_trials = t_nml.getIntValue("max_system_trials");
  rmsd_tolerance = t_nml.getRealValue("rmsd_tol");

  // Validate input
  validateSampleChirality(t_nml.getStringValue("sample_chirality"));
  validateSampleCisTrans(t_nml.getStringValue("sample_cis_trans"));
  validatePreventHBonds(t_nml.getStringValue("prevent_hbonds"));
  validateStateCounts();
}

//-------------------------------------------------------------------------------------------------
const std::string& ConformerControls::getCoreAtomMask() const {
  return core_atom_mask;
}

//-------------------------------------------------------------------------------------------------
const std::string& ConformerControls::getCoreDataItemName() const {
  return core_data_item;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreRK2Value() const {
  return core_rk2;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreRK3Value() const {
  return core_rk3;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreR2Value() const {
  return core_r2;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreR3Value() const {
  return core_r3;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::sampleChirality() const {
  return sample_chirality;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::sampleCisTrans() const {
  return sample_cis_trans;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::preventHydrogenBonding() const {
  return prevent_hbonds;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRunningStateCount() const {
  return running_states;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getFinalStateCount() const {
  return final_states;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRotationSampleCount() const {
  return rotation_samples;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRotatableBondLimit() const {
  return rotatable_bond_limit;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getMaxSeedingAttempts() const {
  return max_seeding_attempts;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getSystemTrialCount() const {
  return system_trials;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getRMSDTolerance() const {
  return rmsd_tolerance;
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateCoreRestraint() const {
  if (core_r2 < 0.0 || core_r2 > core_r3) {
    rtErr("Displacement limits of " + realToString(core_r2, 7, 4, NumberFormat::STANDARD_REAL) +
          " to " + realToString(core_r3, 7, 4, NumberFormat::STANDARD_REAL) + " are invalid.",
          "ConformerControls", "validateCoreRestraint");
  }
  if (core_rk2 < 0.0) {
    rtErr("A repulsive restraint penalty of " +
          realToString(core_rk2, 7, 4, NumberFormat::STANDARD_REAL) + " kcal/mol-A^2 is "
          "unacceptable.", "ConformerControls", "validateCoreRestraint");
  }
  if (core_rk3 < 0.0) {
    rtErr("A restraint penalty of " + realToString(core_rk3, 7, 4, NumberFormat::STANDARD_REAL) +
          " kcal/mol-A^2 is unacceptable.", "ConformerControls", "validateCoreRestraint");
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSampleChirality(const std::string &directive) const {
  if (strcmpCased(directive, "true") == false && strcmpCased(directive, "false") == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The sample_chirality keyword accepts input of 'true' or 'false'.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleChirality");
    case ExceptionResponse::WARN:
      rtWarn("The sample_chirality keyword accepts input of 'true' or 'false'.  Input " +
             directive + " is invalid.  Chirality will not be sampled.", "ConformerControls",
             "validateSampleChirality");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSampleCisTrans(const std::string &directive) const {
  if (strcmpCased(directive, "true") == false && strcmpCased(directive, "false") == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The sample_cis_trans keyword accepts input of 'true' or 'false'.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleCisTrans");
    case ExceptionResponse::WARN:
      rtWarn("The sample_cis_trans keyword accepts input of 'true' or 'false'.  Input " +
             directive + " is invalid.  Cis- and trans- states of molecules will not be sampled.",
             "ConformerControls", "validateSampleCisTrans");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validatePreventHBonds(const std::string &directive) const {
  if (strcmpCased(directive, "true") == false && strcmpCased(directive, "false") == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The prevent_hbonds keyword accepts input of 'true' or 'false' to apply restraints "
            "that will inhibit hydrogen bond formation within the molecules.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleCisTrans");
    case ExceptionResponse::WARN:
      rtWarn("The prevent_hbonds keyword accepts input of 'true' or 'false' to apply restraints "
            "that will inhibit hydrogen bond formation within the molecules.  Input " +
             directive + " is invalid.  Intramolecular hydrogen bonds will be free to form.",
             "ConformerControls", "validatePreventHBonds");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateStateCounts() {
  if (running_states > active_states_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of running states (" + std::to_string(running_states) + ") cannot exceed "
            "the maximum allowed number of " + std::to_string(active_states_limit) + ".",
            "ConformerControls", "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of running states (" + std::to_string(running_states) + ") exceeds "
             "the maximum allowed number of " + std::to_string(active_states_limit) + " and will "
             "be reduced to limit memory usage.", "ConformerControls", "validateStateCounts");
      running_states = active_states_limit;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (system_trials > active_states_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of energy minimizations for a single system (" +
            std::to_string(system_trials) + ") cannot exceed the maximum allowed number of " +
            std::to_string(active_states_limit) + ".", "ConformerControls", "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of energy minimizations for a single system (" +
             std::to_string(system_trials) + ") cannot exceed the maximum allowed number of " +
             std::to_string(active_states_limit) + " and will be reduced.", "ConformerControls",
             "validateStateCounts");
      system_trials = active_states_limit;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (final_states > system_trials) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of final states (" + std::to_string(final_states) + ") cannot exceed the "
            "number of energy minimizations evaluated per system (" +
            std::to_string(system_trials) + ").", "ConformerControls", "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of final states (" + std::to_string(final_states) + ") cannot exceed the "
             "number of energy minimizations evaluated per system (" +
             std::to_string(system_trials) + ") and will be reduced.", "ConformerControls",
             "validateStateCounts");
      final_states = system_trials;
      break;      
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (rmsd_tolerance < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("A negative positional root mean squared deviation criterion to distinguish unique "
             "conformers is non-sensical.  A value of zero will be applied.", "ConformerControls",
             "validateStateCounts");
      rmsd_tolerance = 0.0;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator conformerInput(const TextFile &tf, int *start_line, bool *found,
                                const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("conformer", CaseSensitivity::AUTOMATIC, policy,
                         "Collects instructions for conformer sampling in STORMM.");
  const std::string core_help("Define a subset of atoms that will be prevented from moving, or "
                              "constrained to move in proximity to their initial positions.");
  const std::vector<std::string> core_keys_help = {
    "A data item in a system's SD file to seek out, which will contain one of the following: a "
    "list of atom indices, a list of atom names, or an atom mask specific to that file.  If "
    "the coordinates are not present in an SD file or the SD file does not contain such a data "
    "item, no core atoms will be defined.",
    "An atom mask defining the core atoms.  If given, this will supercede any data items in SD "
    "files and apply to all systems in the calculation.",
    "The restraint penalty, in kcal/mol-Angstrom^2, defining a harmonic repulsive potential that "
    "pushes core atoms away from their initial coordinates.",
    "Alias for 'rk2' in this context.  If both are supplied, the value of 'rk2' will take "
    "precedence.",
    "The harmonic restraint stiffness, in units of kcal/mol-Angstrom^2, preventing atoms from "
    "wandering away from their initial positions.",
    "Alias for 'rk3'.  If both are defined, the value of 'rk3' will take precedence.",
    "Alias for 'rk2' and 'rk3', in units of kcal/mol-Angstrom^2 with a default of 16.0.  This "
    "value will apply to both stiffness constants, but specifying either 'rk2', 'repulsion', "
    "'rk3', or 'attraction' will override the effect.",
    "The distance, in Angstroms, at which to stop applying a repulsive potential pushing atoms "
    "away from their initial positions.",
    "Alias for 'r2', with a default value of 0.0 Angstroms.  If 'r2' is supplied, that value will "
    "take precedence.",
    "The distance, in Angstroms, at which to begin applying an attractive potential that keeps "
    "atoms from wandering away from their initial positions.",
    "Alias for 'r3', with a default value of 0.0 Angstroms.  If 'r3' is supplied, that value will "
    "take precedence." };
  t_nml.addKeyword(NamelistElement("core_mask",
                                   { "data_item", "atoms", "rk2", "repulsion", "rk3", "attraction",
                                     "stiffness", "r2", "demand", "r3", "grace" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), "16.0", std::string(""),
                                     "0.0", std::string(""), "0.0", "0.0" },
                                   DefaultIsObligatory::NO,
                                   InputRepeats::NO, core_help, core_keys_help,
                                   { SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL }));
  t_nml.addKeyword(NamelistElement("anchor_conf", NamelistType::STRING, std::string("")));
  t_nml.addKeyword(NamelistElement("sample_chirality", NamelistType::STRING,
                                   std::string(default_conf_chirality)));
  t_nml.addKeyword(NamelistElement("sample_cis_trans", NamelistType::STRING,
                                   std::string(default_conf_cis_trans)));
  t_nml.addKeyword(NamelistElement("prevent_hbonds", NamelistType::STRING,
                                   std::string(default_conf_stop_hbonds)));
  t_nml.addKeyword(NamelistElement("running_states", NamelistType::INTEGER,
                                   std::to_string(default_conf_running_states)));
  t_nml.addKeyword(NamelistElement("final_states", NamelistType::INTEGER,
                                   std::to_string(default_conf_final_states)));
  t_nml.addKeyword(NamelistElement("rotation_samples", NamelistType::INTEGER,
                                   std::to_string(default_conf_rotation_samples)));
  t_nml.addKeyword(NamelistElement("max_rotatable_bonds", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_rotatable_bonds)));
  t_nml.addKeyword(NamelistElement("max_seeding_attempts", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_seeding_attempts)));
  t_nml.addKeyword(NamelistElement("max_system_trials", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_system_trials)));
  t_nml.addKeyword(NamelistElement("rmsd_tol", NamelistType::REAL,
                                   std::to_string(default_conf_rmsd_tolerance)));
  t_nml.addHelp("core_mask", "Atom mask for common core atoms.  These atoms will be held in a "
                "rigid configuration during energy minimization and other sampling operations.");
  t_nml.addHelp("sample_chirality", "Sample chiral states of identifiable chiral centers.  "
                "Specify 'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("sample_cis_trans", "Sample cis and trans states of double bonds.  Specify "
                "'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("prevent_hbonds", "A quick way to have STORMM prevent hydrogen bonding between "
                "donors and acceptor atoms that it can identify in the molecule(s).  This will "
                "establish a restraint ensemble for each case with default parameters to prevent "
                "donor and acceptor pairs from coming too close.");
  t_nml.addHelp("running_states", "Number of energy-minimizations to carry out at one time in "
                "order to generate a smaller set of final states.");
  t_nml.addHelp("final_states", "Number of final, energy-minimized states to accept as unique "
                "conformers.");
  t_nml.addHelp("rotation_samples", "Number of samples to apply to each rotatable bond.  "
                "Locations for the samples will be chosen based on the detected minima along each "
                "bond's rotation profile.");
  t_nml.addHelp("max_rotatable_bonds", "The maximum number of rotatable bonds to explicitly "
                "sample.  This quickly runs into a combinatorial problem, but there is a "
                "guardrail in the max_system_trials keyword.");
  t_nml.addHelp("max_seeding_attempts", "If a conformer's initial configuration contains a clash "
                "between atoms (their van-der Waals radii are violated), randomization of the "
                "configuration will occur for this number of attempts.  If, after exhausting this "
                "provision, a stable configuration still cannot be found, the input configuration "
                "will be accepted as the initial coordinates for subsequent energy "
                "minimizations.");
  t_nml.addHelp("max_system_trials", "The maximum number of trials that will be made for each "
                "system.  Explicit sampling of chirality, cis-trans isomers, and then rotatable "
                "bonds will proceed in that priority, but the maximum number of sampled states "
                "will be cut off at this value.");
  t_nml.addHelp("rmsd_tol", "Positional, mass-weighted root-mean squared deviation between "
                "conformers required to declare uniqueness.");
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
