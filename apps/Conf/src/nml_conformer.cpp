#include "nml_conformer.h"

#include "../../../src/Namelists/input.h"

namespace conf_app {
namespace user_input {

using omni::constants::CaseSensitivity;
using omni::namelist::NamelistElement;
using omni::namelist::NamelistType;
using omni::namelist::readNamelist;
using omni::parse::WrapTextSearch;

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const ExceptionResponse policy_in) :
    policy{policy_in}, common_atom_mask{std::string("")}, sample_chirality{false},
    sample_cis_trans{false}, prevent_hbonds{false},
    running_states{default_conf_running_states},
    final_states{default_conf_final_states},
    rotation_samples{default_conf_rotation_samples},
    rotatable_bond_limit{default_conf_max_rotatable_bonds},
    system_trial_limit{default_conf_max_system_trials},
    rmsd_tolerance{default_conf_rmsd_tolerance}
{}

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const TextFile &tf, int *start_line,
                                     const ExceptionResponse policy_in) :
  ConformerControls(policy_in)
{
  const NamelistEmulator t_nml = conformerInput(tf, start_line, policy);
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
int ConformerControls::getSystemTrialLimit() const {
  return system_trial_limit;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getRmsdTolerance() const {
  return rmsd_tolerance;
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
