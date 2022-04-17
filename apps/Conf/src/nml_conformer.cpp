#include "nml_conformer.h"

#include "../../../src/Namelists/input.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"

namespace conf_app {
namespace user_input {

using omni::constants::CaseSensitivity;
using omni::namelist::NamelistElement;
using omni::namelist::NamelistType;
using omni::namelist::readNamelist;
using omni::errors::rtErr;
using omni::errors::rtWarn;
using omni::parse::CaseSensitivity;
using omni::parse::strcmpCased;
using omni::parse::WrapTextSearch;

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const ExceptionResponse policy_in) :
    policy{policy_in}, common_atom_mask{std::string("")}, sample_chirality{false},
    sample_cis_trans{false}, prevent_hbonds{false},
    running_states{default_conf_running_states},
    final_states{default_conf_final_states},
    rotation_samples{default_conf_rotation_samples},
    rotatable_bond_limit{default_conf_max_rotatable_bonds},
    system_trials{default_conf_max_system_trials},
    rmsd_tolerance{default_conf_rmsd_tolerance}
{}

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const TextFile &tf, int *start_line,
                                     const ExceptionResponse policy_in) :
    ConformerControls(policy_in)
{
  const NamelistEmulator t_nml = conformerInput(tf, start_line, policy);
  common_atom_mask = t_nml.getStringValue("commonmask");
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
  system_trials = t_nml.getIntValue("max_system_trials");
  rmsd_tolerance = t_nml.getRealValue("rmsd_tol");

  // Validate input
  validateSampleChirality(t_nml.getStringValue("sample_chirality"));
  validateSampleCisTrans(t_nml.getStringValue("sample_cis_trans"));
  validatePreventHBonds(t_nml.getStringValue("prevent_hbonds"));
  validateStateCounts();
}

//-------------------------------------------------------------------------------------------------
std::string ConformerControls::getCommonAtomMask() const {
  return common_atom_mask;
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
int ConformerControls::getSystemTrialCount() const {
  return system_trials;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getRmsdTolerance() const {
  return rmsd_tolerance;
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
NamelistEmulator conformerInput(const TextFile &tf, int *start_line,
                                const ExceptionResponse policy) {
  NamelistEmulator t_nml("conformer", CaseSensitivity::AUTOMATIC, policy,
                         "Collects instructions for conformer sampling in OMNI.");
  t_nml.addKeyword(NamelistElement("commonmask", NamelistType::STRING, "MISSING"));
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
  t_nml.addKeyword(NamelistElement("max_system_trials", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_system_trials)));
  t_nml.addKeyword(NamelistElement("rmsd_tol", NamelistType::REAL,
                                   std::to_string(default_conf_rmsd_tolerance)));
  t_nml.addHelp("commonmask", "Atom mask for common core atoms.  These atoms will be held in a "
                "rigid configuration during energy minimization and other sampling operations.");
  t_nml.addHelp("sample_chirality", "Sample chiral states of identifiable chiral centers.  "
                "Specify 'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("sample_cis_trans", "Sample cis and trans states of double bonds.  Specify "
                "'yes' / 'true' to sample or 'no' / 'false' to decline.");
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
  t_nml.addHelp("max_system_trials", "The maximum number of trials that will be made for each "
                "system.  Explicit sampling of chirality, cis-trans isomers, and then rotatable "
                "bonds will proceed in that priority, but the maximum number of sampled states "
                "will be cut off at this value.");
  t_nml.addHelp("rmsd_tol", "Positional, mass-weighted root-mean squared deviation between "
                "conformers required to declare uniqueness.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount());
  return t_nml;
}

} // namespace user_input
} // namespace conf_app
