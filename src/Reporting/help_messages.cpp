#include <cstring>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/nml_conformer.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_report.h"
#include "Namelists/nml_restraint.h"
#include "Namelists/nml_solvent.h"
#include "Parsing/parse.h"
#include "Parsing/textfile.h"
#include "help_messages.h"

namespace stormm {
namespace display {

using constants::CaseSensitivity;
using namelist::conformerInput;
using namelist::dynamicsInput;
using namelist::filesInput;
using namelist::minimizeInput;
using namelist::reportInput;
using namelist::restraintInput;
using namelist::solventInput;
using namelist::SubkeyRequirement;
using parse::strcmpCased;
using parse::TextFile;
using parse::TextOrigin;

//-------------------------------------------------------------------------------------------------
bool detectHelpSignal(const int argc, const char* argv[], const HelpSignalKind kind,
                      const std::vector<std::string> &help_words) {
  const int nh_words = help_words.size();
  switch (kind) {
  case HelpSignalKind::NO_ARGS:
  case HelpSignalKind::NO_ARGS_ONLY:
  case HelpSignalKind::KEYWORD:
    if (argc == 1) {
      return true;
    }
    else {
      return false;
    }
    break;
  case HelpSignalKind::KEYWORD_ONLY:
    break;
  }
  switch (kind) {
  case HelpSignalKind::NO_ARGS_ONLY:
    break;
  case HelpSignalKind::NO_ARGS:
  case HelpSignalKind::KEYWORD:
  case HelpSignalKind::KEYWORD_ONLY:  
    for (int i = 0; i < argc; i++) {
      if (nh_words == 0) {
        if (strcmpCased(argv[i], "help", CaseSensitivity::NO) ||
            strcmpCased(argv[i], "-help", CaseSensitivity::NO) ||
            strcmpCased(argv[i], "--help", CaseSensitivity::NO)) {
          return true;
        }
      }
      else {
        for (int j = 0; j < nh_words; j++) {
          if (strcmpCased(argv[i], help_words[j], CaseSensitivity::YES)) {
            return true;
          }
        }
      }
    }
    break;
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const int argc, const char* argv[],
                         const std::vector<std::string> &module_name) {
  const int n_names = module_name.size();
  bool found = false;
  for (int i = 0; i < n_names; i++) {
    for (int j = 0; j < argc; j++) {
      if (strcmp(argv[j], module_name[i].c_str()) == 0) {
        found = (found || displayNamelistHelp(module_name[i]));
      }
    }
  }
  return found;
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const int argc, const char* argv[], const std::string &module_name) {
  return displayNamelistHelp(argc, argv, std::vector<std::string>(1, module_name));
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const std::string &module_name) {
  const TextFile tf(std::string(""), TextOrigin::RAM);
  int start_line = 0;
  if (strcmpCased(module_name, "&files", CaseSensitivity::YES)) {
    const std::vector<SubkeyRequirement> sys_keyword_reqs(1, SubkeyRequirement::REQUIRED);
    const NamelistEmulator t_nml = filesInput(tf, &start_line, sys_keyword_reqs);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&restraint", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = restraintInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&solvent", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = solventInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&minimize", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = minimizeInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&conformer", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = conformerInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&dynamics", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = restraintInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&report", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = reportInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else {
    rtErr("No namelist " + module_name + " is known in the STORMM libraries.",
          "displayNamelistHelp");
  }
  __builtin_unreachable();
}

} // namespace display
} // namespace stormm
