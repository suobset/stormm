#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/polynumeric.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Topology/amber_prmtop_util.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "command.h"

using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::errors::rtWarn;
using omni::diskutil::osSeparator;
using conf_app::user_input::UserSettings;

//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv);

  // Read topologies and coordinate files
  for (int i = 0; i < ui.getSystemCount(); i++) {
    
  }
}
