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
NamelistEmulator conformerInput(const TextFile &tf, int *start_line,
                                const ExceptionResponse policy) {
  NamelistEmulator t_nml("conformer", CaseSensitivity::AUTOMATIC, policy,
                         "Collects instructions for conformer sampling in OMNI.");
  t_nml.addKeyword(NamelistElement("commonmask", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("rotation_samples", NamelistType::INTEGER,
                                   std::to_string(default_conf_rotation_samples)));
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
  t_nml.addKeyword(NamelistElement("rmsd_tol", NamelistType::REAL,
                                   std::to_string(default_conf_rmsd_tolerance)));
  t_nml.addHelp("commonmask", "Atom mask for common core atoms.  These atoms will be held in a "
                "rigid configuration during energy minimization and other sampling operations.");
  t_nml.addHelp("rotation_samples", "Number of samples to make about each rotatable bond.");
  t_nml.addHelp("sample_chirality", "Sample chiral states of identifiable chiral centers.  "
                "Specify 'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("sample_cis_trans", "Sample cis and trans states of double bonds.  Specify "
                "'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("final_states", "Number of final, energy-minimized states to accept as unique "
                "conformers.");
  t_nml.addHelp("rmsd_tol", "Positional, mass-weighted root-mean squared deviation between "
                "conformers required to declare uniqueness.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount());
  return t_nml;
}

} // namespace user_input
} // namespace conf_app
