#include "namelist_element.h"
#include "nml_restraint.h"

namespace omni {
namespace namelist {

//-------------------------------------------------------------------------------------------------
NamelistEmulator restraintInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("restraint", CaseSensitivity::AUTOMATIC, ExceptionResponse::DIE,
                         "Replicates the Amber NMR restraint namelist within OMNI.");
  t_nml.addKeyword(NamelistElement("iat1", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("iat2", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("iat3", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("iat4", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("nstep1", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("nstep2", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("ninc", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("r1", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r2", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r3", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r4", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r1a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r2a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r3a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r4a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("rk2", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("rk3", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("rk2a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("rk3a", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("mask_i", NamelistType::STRING, ""));
  t_nml.addKeyword(NamelistElement("mask_j", NamelistType::STRING, ""));
  t_nml.addKeyword(NamelistElement("mask_k", NamelistType::STRING, ""));
  t_nml.addKeyword(NamelistElement("mask_l", NamelistType::STRING, ""));
  t_nml.addHelp("iat1", "The first atom in the restraint (this or mask_i is required)");
  t_nml.addHelp("iat2", "The second atom in the restraint (this or mask_j is required)");
  t_nml.addHelp("iat3", "The third atom in the restraint (optional, will convert a distance "
                "restraint into an angle restraint if provided)");
  t_nml.addHelp("iat4", "The fourth atom in the restraint (optional, will convert an angle "
                "bending restraint into a dihedral angle restraint if provided)");
  t_nml.addHelp("nstep1", "Step number at which to begin applying the restraint, with "
                "displacement parameters r1, r2, r3, and r4, plus stiffness parameter rk2 and "
                "rk3.");
  t_nml.addHelp("nstep2", "Step number by which to finish applying the restraint, with "
                "displacement parameters r1a, r2a, r3a, and r4a, plus stiffness parameter rk2a "
                "and rk3a (if r[1-4]a and rk[2,3]a are unspecified, they are taken to be the same "
                "as their initial counterparts and the restraint does not scale between nstep1 "
                "and nstep2).");
  t_nml.addHelp("ninc", "Method by which to modulate the restraint parameters between r[1-4] and "
                "r[1-4]a, as well as between rk[2,3] and rk[2,3]a, between steps nstep1 and "
                "nstep2.");
  t_nml.addHelp("r1", "Leftmost distance at which the left-hand harmonic restraint linearizes");
  t_nml.addHelp("r2", "Middle distance parameter, to the left of which the left-hand harmonic "
                "restraint applies.  The left-hand harmonic restraint evaluates as "
                "k2 * (R - r2)^2 for r1 <= R <= r2.");
  t_nml.addHelp("r3", "Middle distance parameter, to the left of which the potential is flat and "
                "to the right of which the right-hand harmonic restraint is applied.  The "
                "potential is flat, force 0, for a distance R such that r2 <= R <= r3.");
  t_nml.addHelp("r4", "Rightmost distance parameter, defining the right-hand harmonic restraint "
                "on the interval r3 <= R <= r4.  To the right of r4 the right-hand harmonic "
                "restraint linearizes.");
  t_nml.addHelp("rk2", "Stiffness constant for the left-hand harmonic restraint.");
  t_nml.addHelp("rk2", "Stiffness constant for the right-hand harmonic restraint.");
  t_nml.addHelp("r1a", "If specified along with nstep1, nstep2, and ifvari > 0, this specifies "
                "the final value of the left-most restraint distance parameter.");
  t_nml.addHelp("r2a", "Final value of the right-most limit of the left-hand harmonic restraint, "
                "given ifvari > 0 and appropriate nstep values.");
  t_nml.addHelp("r3a", "Final value of the left-most limit of the right-hand harmonic restraint, "
                "given ifvari > 0 and appropriate nstep values.");
  t_nml.addHelp("r4a", "Final value of the right-most limit of the right-hand harmonic restraint "
                "(beyond which it linearizes), given ifvari > 0 and appropriate nstep values.");
  t_nml.addHelp("rk2a", "Final value of the left-hand harmonic restraint stiffness, given ifvari "
                "> 0 and appropriate nstep values.");
  t_nml.addHelp("rk3a", "Final value of the right-hand harmonic restraint stiffness, given ifvari "
                "> 0 and appropriate nstep values.");
  t_nml.addHelp("mask_i", "Ambmask string evaluating to the first atom in the restraint (may be "
                "given in place of iatm1, but one of these is required)");
  t_nml.addHelp("mask_j", "Ambmask string evaluating to the second atom in the restraint (may be "
                "given in place of iatm2, but one of these is required)");
  t_nml.addHelp("mask_k", "Ambmask string evaluating to the third atom in the restraint (may be "
                "given in place of iatm3, and either will convert a distance restraint into an "
                "angle bending restraint)");
  t_nml.addHelp("mask_l", "Ambmask string evaluating to the fourth atom in the restraint (may be "
                "given in place of iatm4, and either will convert a distance restraint into an "
                "angle bending restraint)");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.  All calls to this function should
  // proceed in consecutive calls, to make use of the updates to start_line and avoid reading any
  // instance of this namelist twice or skipping instances of it in the search for some other
  // namelist.  An alternative is to keep an independent counter to track progress through the
  // input file in search for &rst namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::NO, tf.getLineCount());

  return t_nml;
}

} // namespace namelist
} // namespace omni
