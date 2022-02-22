#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "command.h"

using omni::synthesis::SystemCache;
using conf_app::user_input::UserSettings;

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv);

  // Read topologies and coordinate files
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior());

}
