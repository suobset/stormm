#include <cstring>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/namelist_inventory.h"
#include "Namelists/nml_conformer.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_precision.h"
#include "Namelists/nml_random.h"
#include "Namelists/nml_remd.h"
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
using namelist::ffmorphInput;
using namelist::filesInput;
using namelist::minimizeInput;
using namelist::namelist_inventory;
using namelist::precisionInput;
using namelist::randomInput;
using namelist::remdInput;
using namelist::reportInput;
using namelist::restraintInput;
using namelist::solventInput;
using namelist::KeyRequirement;
using parse::strcmpCased;
using parse::TextFile;
using parse::TextOrigin;

//-------------------------------------------------------------------------------------------------
bool detectHelpSignal(const int argc, const char* argv[], const HelpSignalKind kind,
                      const std::vector<std::string> &help_words) {
  const int nh_words = help_words.size();

  // If the only way to get help is by supplying no arguments, return TRUE or FALSE based on
  // whether there are no command line arguments.  If the help is triggered by a keyword, but not
  // keyword only, assume that running the program with no arguments is also intended to trigger
  // the help messages.
  switch (kind) {
  case HelpSignalKind::NO_ARGS_ONLY:
    return (argc == 1);
  case HelpSignalKind::NO_ARGS:    
  case HelpSignalKind::KEYWORD:
    if (argc == 1) {
      return true;
    }
    break;
  case HelpSignalKind::KEYWORD_ONLY:
    break;
  }

  // If there are keywords to parse, cases where the program is supposed to print help messages
  // only in the case that no command line arguments are supplied would already have receive a
  // signal of FALSE.  However, if the program is set to provide help messages with a particular
  // command line argument, parse each argument here.
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
                         const std::vector<std::string> &module_name,
                         const std::vector<NamelistToken> &custom_namelists) {
  const int n_names = module_name.size();
  bool found = false;
  for (int i = 0; i < n_names; i++) {
    const std::string module_less_amp = (module_name[i][0] == '&') ? module_name[i].substr(1) :
                                                                     module_name[i];
    for (int j = 0; j < argc; j++) {
      if (strcmp(argv[j], module_name[i].c_str()) == 0 ||
          strcmp(argv[j], module_less_amp.c_str()) == 0) {
        found = (found || displayNamelistHelp(module_name[i], custom_namelists));
      }
    }
  }
  return found;
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const int argc, const char* argv[], const std::string &module_name,
                         const std::vector<NamelistToken> &custom_namelists) {
  return displayNamelistHelp(argc, argv, std::vector<std::string>(1, module_name),
                             custom_namelists);
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const std::string &module_name,
                         const std::vector<NamelistToken> &custom_namelists) {
  const TextFile tf(std::string(""), TextOrigin::RAM);
  const std::string search_name = (module_name[0] == '&') ? module_name : "&" + module_name;
  const size_t nmlcount = namelist_inventory.size();
  for (size_t i = 0; i < nmlcount; i++) {
    int start_line = 0;
    bool found = false;
    if (strcmpCased(search_name, namelist_inventory[i].getTitle(), CaseSensitivity::YES)) {
      const NamelistEmulator t_nml = namelist_inventory[i].invoke(tf, &start_line, &found);
      t_nml.printHelp();
      return true;
    }
  }
  if (custom_namelists.size() > 0) {
    const size_t custom_nmlcount = custom_namelists.size();
    for (size_t i = 0; i < custom_nmlcount; i++) {
      int start_line = 0;
      bool found = false;
      if (strcmpCased(search_name, custom_namelists[i].getTitle(), CaseSensitivity::YES)) {
        const NamelistEmulator t_nml = custom_namelists[i].invoke(tf, &start_line, &found);
        t_nml.printHelp();
        return true;
      }
    }
  }
  return false;
}

} // namespace display
} // namespace stormm
