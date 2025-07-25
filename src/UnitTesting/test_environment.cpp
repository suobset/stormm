#include <algorithm>
#include <cstring>
#include "copyright.h"
#include "FileManagement/directory_util.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "test_environment.h"

namespace stormm {
namespace testing {

using diskutil::DrivePathType;
using diskutil::findStormmPath;
using diskutil::getDrivePathType;
using diskutil::makePathAbsolute;
using diskutil::stormmMkdir;
using diskutil::stormmBatchRmdir;
using diskutil::stormmRmdir;
using diskutil::openOutputFile;
using diskutil::osSeparator;
using diskutil::removeFile;
using diskutil::separatePath;
using namelist::InputStatus;
using namelist::NamelistEmulator;
using namelist::NamelistType;
using parse::uppercase;
using parse::NumberFormat;

//-------------------------------------------------------------------------------------------------
TestEnvironment::TestEnvironment(const int argc, const char* argv[], CommandLineParser *clip,
                                 const TmpdirStatus tmpdir_required,
                                 const ExceptionResponse policy) :
    user_mods{"user_mods", "Relevant information taken in from the command line"},
    verbose_level{TestVerbosity::COMPACT},
    general_tolerance{1.0e-4},
    random_seed{827493},
    stormm_home_path{""},
    stormm_source_path{""},
    tmpdir_path{""},
    tmpdir_created{false},
    tmpdir_writeable{true},
    files_created{},
    directories_created{},
    remove_files{true},
    do_intensive_tests{true},
    snapshot_behavior{SnapshotOperation::COMPARE},
    display_time{false}
{
  // Fill out the command-line keywords relevant to testing
  user_mods.suppressHelpOnNoArgs();
  user_mods.addStandardAmberInputs("-ig_seed");
  NamelistEmulator *t_nml = user_mods.getNamelistPointer();
  t_nml->addKeyword("-tolscale", NamelistType::REAL, std::to_string(1.0));
  t_nml->addHelp("-tolscale", "Scaling factor to be applied on the tolerance for unit tests");
  t_nml->addKeyword("-verbose", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-verbose", "The verbosity to apply in reporting test results");
  t_nml->addKeyword("-stormmhome", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-stormmhome", "Custom path to the STORMM home directory");
  t_nml->addKeyword("-stormmsrc", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-stormmsrc", "Custom path to the STORMM source directory");
  t_nml->addKeyword("-tmpdir", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-tmpdir", "Custom path to the STORMM temporary directory");
  t_nml->addKeyword("-keep_files", NamelistType::BOOLEAN);
  t_nml->addHelp("-keep_files", "Specify to retain files generated by the test cases in the "
                 "temporary directory");
  t_nml->addKeyword("-snapshot", NamelistType::BOOLEAN);
  t_nml->addHelp("-snapshot", "Specify to rebuild snapshot files when remaking test cases");
  t_nml->addKeyword("-timings", NamelistType::BOOLEAN);
  t_nml->addHelp("-timings", "Request timings data from the tests");
  t_nml->addKeyword("-intensive", NamelistType::BOOLEAN);
  t_nml->addHelp("-intensive", "Trigger intensive unit tests for certain programs");
  if (clip != nullptr) {
    user_mods.coordinateWithPartner(clip);
  }
  user_mods.parseUserInput(argc, argv);  

  // Translate various command-line arguments
  display_time = t_nml->getBoolValue("-timings");
  if (t_nml->getBoolValue("-snapshot")) {
    snapshot_behavior = SnapshotOperation::SNAPSHOT;
  }
  
  // Loop over command line input and try to parse instructions
  bool cli_vbs = (t_nml->getKeywordStatus("-verbose") == InputStatus::USER_SPECIFIED);
  bool cli_tol = (t_nml->getKeywordStatus("-tolscale") == InputStatus::USER_SPECIFIED);
  bool cli_rng = (t_nml->getKeywordStatus("-ig_seed") == InputStatus::USER_SPECIFIED);
  bool cli_dir = (t_nml->getKeywordStatus("-tmpdir") == InputStatus::USER_SPECIFIED);
  bool cli_ohm = (t_nml->getKeywordStatus("-stormmhome") == InputStatus::USER_SPECIFIED);
  bool cli_osc = (t_nml->getKeywordStatus("-stormmsrc") == InputStatus::USER_SPECIFIED);
  do_intensive_tests = t_nml->getBoolValue("-intensive");
  
  // Seek environment variables, if there are no command-line arguments to supersede them
  const char* vbs_char = std::getenv("STORMM_VERBOSE");
  const char* tol_char = std::getenv("STORMM_TEST_TOL");
  const char* rng_char = std::getenv("STORMM_PRNG_SEED");
  const char* dir_char = std::getenv("STORMM_TMPDIR");
  const char* ohm_char = std::getenv("STORMM_HOME");
  const char* osc_char = std::getenv("STORMM_SOURCE");
  std::string vbs_str = (vbs_char == nullptr) ? "" : std::string(vbs_char);
  std::string tol_str = (tol_char == nullptr) ? "" : std::string(tol_char);
  std::string rng_str = (rng_char == nullptr) ? "" : std::string(rng_char);
  std::string dir_str = (dir_char == nullptr) ? "" : std::string(dir_char);
  std::string ohm_str = (ohm_char == nullptr) ? "" : std::string(ohm_char);
  std::string osc_str = (osc_char == nullptr) ? "" : std::string(osc_char);
  if (cli_vbs) {
    vbs_str = t_nml->getStringValue("-verbose");
  }
  if (cli_tol) {
    tol_str = std::to_string(t_nml->getRealValue("-tolscale"));
  }
  if (cli_rng) {
    rng_str = std::to_string(t_nml->getIntValue("-ig_seed"));
  }
  if (cli_dir) {
    dir_str = t_nml->getStringValue("-tmpdir");
  }
  if (cli_ohm) {
    ohm_str = t_nml->getStringValue("-stormmhome");
  }
  if (cli_osc) {
    osc_str = t_nml->getStringValue("-stormmsrc");
  }
  
  // Translate the critical arguments, whatever their source
  if (vbs_str.size() > 0) {
    if (vbs_str == "FULL") {
      verbose_level = TestVerbosity::FULL;
    }
    else if (vbs_str == "COMPACT") {
      verbose_level = TestVerbosity::COMPACT;
    }
    else if (vbs_str == "FAIL" || vbs_str == "FAILURE_ONLY") {
      verbose_level = TestVerbosity::FAILURE_ONLY;
    }
    else {
      rtWarn("Invalid setting for environment variable STORMM_VERBOSE (" + vbs_str + ").  The "
             "default value of COMPACT will be taken instead.");
    }
  }
  if (tol_str.size() > 0) {
    if (verifyNumberFormat(tol_str.c_str(), NumberFormat::SCIENTIFIC) ||
        verifyNumberFormat(tol_str.c_str(), NumberFormat::STANDARD_REAL) ||
        verifyNumberFormat(tol_str.c_str(), NumberFormat::INTEGER)) {
      general_tolerance = strtod(tol_str.c_str(), nullptr);
    }
    else {
      rtWarn("Invalid setting for environment variable STORMM_TEST_TOL (" + tol_str + ").  The "
             "default value of 1.0e-4 will be used instead.");
    }
  }
  if (rng_str.size() > 0) {
    if (verifyNumberFormat(tol_str.c_str(), NumberFormat::INTEGER)) {
      random_seed = strtol(tol_str.c_str(), nullptr, 10);
    }
    else {
      rtWarn("Invalid integer format for environment variable STORMM_PRNG_SEED (" + rng_str +
             ").  The default value of 827493 will be taken instead.", "TestEnvironment");
    }
  }
  if (dir_str.size() > 0) {
    tmpdir_path = dir_str;
  }
  if (ohm_str.size() > 0) {
    stormm_home_path = ohm_str;
  }
  if (osc_str.size() > 0) {
    stormm_source_path = osc_str;
  }

  // ${STORMM_HOME} needs to be understood from the command line or the environment.  The default
  // is a last-ditch effort to infer it from the current working directory or the executable path.
  char tmp_c[512];
  const std::string exe_path(argv[0]);
  if (stormm_home_path.size() == 0 || stormm_source_path.size() == 0) {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
    if (GetCurrentDirectory(512, tmp_c) == nullptr) {
      rtErr("Encountered getcwd() error during initialization.", "TestEnvironment.");
    }
#else
    if (getcwd(tmp_c, 512) == nullptr) {
      rtErr("Encountered getcwd() error during initialization.", "TestEnvironment.");
    }
#endif
  }
  if (stormm_home_path.size() == 0) {
    switch (verbose_level) {
    case TestVerbosity::FULL:
      rtWarn("${STORMM_HOME} is not set in the shell environment.  Define this variable for best "
             "results.", "TestEnvironment");
      break;
    case TestVerbosity::COMPACT:
    case TestVerbosity::FAILURE_ONLY:
      break;
    }
    const std::string known_program = std::string("test") + osSeparator() + "bin" + osSeparator() +
                                      "test_hybrid";
    stormm_home_path = findStormmPath(exe_path, known_program);
    if (stormm_home_path.size() == 0) {
      stormm_home_path = findStormmPath(std::string(tmp_c) + osSeparator() + exe_path,
                                        known_program);
    }
    if (stormm_home_path.size() > 0) {
      switch (verbose_level) {
      case TestVerbosity::FULL:
      case TestVerbosity::COMPACT:
        rtAlert("${STORMM_HOME} is taken to be " + stormm_home_path, "TestEnvironment");
        break;
      case TestVerbosity::FAILURE_ONLY:
        break;
      }
    }
  }

  // Also try to find ${STORMM_SOURCE}
  if (stormm_source_path.size() == 0) {
    switch (verbose_level) {
    case TestVerbosity::FULL:
      rtWarn("${STORMM_SOURCE} is not set in the shell environment.  Define this variable for "
             "best results.", "TestEnvironment");
      break;
    case TestVerbosity::COMPACT:
    case TestVerbosity::FAILURE_ONLY:
      break;
    }
    const std::string known_code = std::string("test") + osSeparator() + "Parsing" +
      osSeparator() + "test_parse.cpp";
    stormm_source_path = findStormmPath(exe_path, known_code);
    if (stormm_source_path.size() == 0) {
      stormm_source_path = findStormmPath(std::string(tmp_c) + osSeparator() + exe_path,
                                          known_code);
    }
    if (stormm_source_path.size() > 0) {
      switch (verbose_level) {
      case TestVerbosity::FULL:
      case TestVerbosity::COMPACT:
        rtAlert("${STORMM_SOURCE} is taken to be " + stormm_source_path, "TestEnvironment");
        break;
      case TestVerbosity::FAILURE_ONLY:
        break;
      }
    }
  }

  // Finally, find ${STORMM_TMPDIR}. Determine whether to remove any temporary, created files.
  bool tmpdir_defined = (tmpdir_path.size() > 0);
  if (tmpdir_path.size() == 0) {
    tmpdir_path = stormm_source_path + osSeparator() + "tmp";
  }
  remove_files = (t_nml->getBoolValue("-keep_files") == false);

  // Create the temporary directory for these tests
  DrivePathType tmpdir_exist = getDrivePathType(tmpdir_path);
  switch (tmpdir_exist) {
  case DrivePathType::FILE:
    if (tmpdir_required == TmpdirStatus::REQUIRED) {
      const std::string advice = (tmpdir_defined) ?
        std::string("  Define STORMM_TMPDIR in the environment or specify "
                    "-tmpdir in the arguments list to this executable.") :
        std::string("");
      rtWarn("STORMM temporary directory path " + tmpdir_path + " already exists and is a file." +
             advice, "TestEnvironment");
    }
    tmpdir_writeable = false;
    break;
  case DrivePathType::DIRECTORY:
    {
      std::string test_filename(tmpdir_path + osSeparator() + "test.test");
      while (getDrivePathType(test_filename) == DrivePathType::FILE) {
        test_filename += ".x";
      }
      try {
        std::ofstream foutp = openOutputFile(test_filename);
      }
      catch (std::runtime_error) {
        if (tmpdir_required == TmpdirStatus::REQUIRED) {
          rtWarn("STORMM temporary directory path " + tmpdir_path + " is not writeable.  Bad "
                 "permissions are likely to blame.", "TestEnvironment.");
        }
        tmpdir_writeable = false;
      }
      if (getDrivePathType(test_filename) == DrivePathType::FILE) {
        removeFile(test_filename);
      }
    }
    break;
  case DrivePathType::REGEXP:
    try {
      tmpdir_created = true;
      directories_created = stormmMkdir(tmpdir_path);
      std::reverse(directories_created.begin(), directories_created.end());
    }
    catch (std::runtime_error) {
      tmpdir_created = false;
      tmpdir_writeable = false;
    }
    break;
  }

  // There is no need to free the memory allocated by std::getenv() above.  Those char*'s are
  // merely pointers to memory pertaining to environment variables that already existed prior
  // to calling std::getenv(), and will be free'd by whatever machinery created the memory.
}

//-------------------------------------------------------------------------------------------------
TestEnvironment::TestEnvironment(const int argc, const char* argv[],
                                 const TmpdirStatus tmpdir_required,
                                 const ExceptionResponse policy) :
    TestEnvironment(argc, argv, nullptr, tmpdir_required, policy)
{}

//-------------------------------------------------------------------------------------------------
TestEnvironment::TestEnvironment(const int argc, const char* argv[],
                                 const ExceptionResponse policy) :
    TestEnvironment(argc, argv, nullptr, TmpdirStatus::NOT_REQUIRED, policy)
{}

//-------------------------------------------------------------------------------------------------
TestEnvironment::~TestEnvironment() {
  
  // Remove files
  if (remove_files) {
    const int n_files = files_created.size();
    for (int i = 0; i < n_files; i++) {
      DrivePathType pth_type = getDrivePathType(files_created[i]);
      if (pth_type == DrivePathType::FILE) {
        remove(files_created[i].c_str());
      }
    }
  }

  // Remove the temporary directory and any others created along the way
  const int n_dirs = directories_created.size();
  stormmBatchRmdir(directories_created);
}

//-------------------------------------------------------------------------------------------------
TestVerbosity TestEnvironment::getVerbosity() const {
  return verbose_level;
}

//-------------------------------------------------------------------------------------------------
double TestEnvironment::getTolerance() const {
  return general_tolerance;
}

//-------------------------------------------------------------------------------------------------
int TestEnvironment::getRandomSeed() const {
  return random_seed;
}

//-------------------------------------------------------------------------------------------------
std::string TestEnvironment::getStormmHomePath() const {
  return stormm_home_path;
}

//-------------------------------------------------------------------------------------------------
std::string TestEnvironment::getStormmSourcePath() const {
  return stormm_source_path;
}

//-------------------------------------------------------------------------------------------------
std::string TestEnvironment::getTemporaryDirectoryPath() const {
  return tmpdir_path;
}

//-------------------------------------------------------------------------------------------------
bool TestEnvironment::getTemporaryDirectoryAccess() const {
  return tmpdir_writeable;
}

//-------------------------------------------------------------------------------------------------
bool TestEnvironment::doIntensiveTests() const {
  return do_intensive_tests;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> TestEnvironment::getListOfFilesCreated() const {
  return files_created;
}

//-------------------------------------------------------------------------------------------------
SnapshotOperation TestEnvironment::takeSnapshot() const {
  return snapshot_behavior;
}

//-------------------------------------------------------------------------------------------------
bool TestEnvironment::getDisplayTimingsOrder() const {
  return display_time;
}

//-------------------------------------------------------------------------------------------------
void TestEnvironment::setFileRemoval(const bool remove_files_in) {
  remove_files = remove_files_in;
}

//-------------------------------------------------------------------------------------------------
void TestEnvironment::logFileCreated(const std::string &path) {
  std::string abs_path = makePathAbsolute(path);
  switch (getDrivePathType(path)) {
  case DrivePathType::FILE:
    break;
  case DrivePathType::DIRECTORY:
    rtErr("The file " + abs_path + " is in fact a directory.", "TestEnvironment",
          "logFileCreated");
  case DrivePathType::REGEXP:
    rtWarn("The file " + abs_path + " does not exist at the time it was catalogged.",
           "TestEnvironment", "logFileCreated");
    break;
  }
  files_created.push_back(abs_path);
}

//-------------------------------------------------------------------------------------------------
void TestEnvironment::logDirectoryCreated(const std::string &path) {
  std::string abs_path = makePathAbsolute(path);
  switch (getDrivePathType(abs_path)) {
  case DrivePathType::FILE:
    rtErr("The directory " + abs_path + " is in fact a file.", "TestEnvironment",
          "logDirectoryCreated");
  case DrivePathType::DIRECTORY:
    break;
  case DrivePathType::REGEXP:
    rtWarn("The directory " + abs_path + " does not exist at the time it was catalogged.",
           "TestEnvironment", "logDirectoryCreated");
    break;
  }
  directories_created.push_back(abs_path);
}

//-------------------------------------------------------------------------------------------------
void TestEnvironment::logDirectoryCreated(const std::vector<std::string> &paths) {
  const int n_dirs = paths.size();
  for (int i = 0; i < n_dirs; i++) {
    std::string abs_path = makePathAbsolute(paths[i]);
    switch (getDrivePathType(abs_path)) {
    case DrivePathType::FILE:
      rtErr("The directory " + abs_path + " is in fact a file.", "TestEnvironment",
            "logDirectoryCreated");
    case DrivePathType::DIRECTORY:
      break;
    case DrivePathType::REGEXP:
      rtWarn("The directory " + abs_path + " does not exist at the time it was catalogged.",
             "TestEnvironment", "logDirectoryCreated");
      break;
    }
    directories_created.push_back(abs_path);
  }
}

} // namespace testing
} // namespace stormm
