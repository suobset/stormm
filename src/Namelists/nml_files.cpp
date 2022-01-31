#include "nml_rst.h"

namespace omni {
namespace namelists {

using parse::NamelistElement;

//-------------------------------------------------------------------------------------------------
NamelistEmulator filesInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("files", CaseSensitivity::AUTO, ExceptionResponse::DIE, "Collects file "
                         "names for OMNI programs, offloading work that would otherwise require "
                         "command-line arguments.");
  t_nml.addKeyword(NamelistElement("-p", NamelistType::STRING, "prmtop", DefaultIsObligatory::NO,
                                   InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("-c", NamelistType::STRING, "inpcrd", DefaultIsObligatory::NO,
                                   InputRepeats::YES));
  const std::string sys_help;
  sys_help = "Expression for a complete system, linking a topology file explicitly to a starting "
             "coordinates file, with the option of that coordinates file being a trajectory with "
             "more than one frame.  This keyword provides a means to read more than one frame "    
             "from a trajectory starting coordinates file, if the frame_end subkey is given and "
             "greater than frame_start.  All starting coordinates will be paired to the same "
             "topology object.  Like several other specifiers in this namelist, this keyword is "
             "repeatable.";
  const std::vector<std::string> sys_keys_help = { "Topology file", "Starting coordinates file"
  t_nml.addKeyword(NamelistElement("-sys", { "-p", "-c", "frame_start", "frame_end", "kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::STRING },
                                   { "prmtop", "inpcrd", "0", "AMBER_ASCII_RST" },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sys_help,
                                   { "Topology file", "Starting coordinates file",
                                     "Starting frame (if the coordinates are a trajectory)",
                                     "Ending frame (if the coordinates are a trajectory).  If "
                                     "unspecified, only the starting frame will be read.  "
                                     "Otherwise, distinct systems will be made for the given "
                                     "topology and every frame between frame_start and frame_end.",
                                     "Type of coordinates file (if unspecified, the type will be "
                                     "detected automatically." }));
  t_nml.addKeyword(NamelistElement("-o", NamelistType::STRING, "mdout"));
  t_nml.addKeyword(NamelistElement("-x", NamelistType::STRING, "mdcrd"));
  t_nml.addKeyword(NamelistElement("-warn", NamelistType::STRING, "warnings"));
  t_nml.addKeyword(NamelistElement("-error", NamelistType::STRING, "errors"));
  t_nml.addHelp("-p", "System topology file.  Repeatable for multiple systems.  Also accepts "
                "regular expressions.");
  t_nml.addHelp("-c", "Input coordinates file.  Repeatable for multiple systems.  Also accepts "
                "regular expressions.  Input coordinate files will be matched to topology files "
                "using atom counts and sanity of valence parameters, if free -c and -p parameters "
                "are provided.  Otherwise, use the system keyword and its subkeys to tie specific "
                "sets of starting coordinates to each topology.");
  t_nml.addHelp("-o", "Output diagnostics file, equivalent to mdout from Amber's sander program.  "
                "Reports for all systems will be included in this file.  Output can become "
                "voluminous for multiple systems if all dump their details into one file.  Use "
                "the \"outfmt\" keyword with setting \"INDIVIDUAL\" to obtain separate files for "
                "each system along with a \".master\" output file providing details of the entire "
                "run.");
  t_nml.addHelp("-x", "Trajectory output file (base name) for each system.  The actual name of "
                "each output file will be \"yyy_(number).zzz\", where \"yyy\" is any part of the "
                "-x string value preceding the final dot (.), \"_(number)\" is the number of the "
                "system in some internal list (key printed in the master output file).");
  t_nml.addHelp("-warn", "Warnings reported for the run, collecting results from all systems.");
  t_nml.addHelp("-error", "Errors reported for the run, collecting results from all systems.");

  // There is expected to be one unique &files namelist in a given input file.  Seek it out by
  // wrapping back to the beginning of the input file if necessary.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTestSearch::YES, tf.getLineCount());
  
  return t_nml;
}
  
