#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/FileManagement/directory_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::errors::rtWarn;
using stormm::review::stormmSplash;
using namespace stormm::diskutil;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  const char osc = osSeparator();

  // Section 1
  section("Directory search and listing");

  // Section 2
  section("Directory creation");

  // Section 3
  section("File creation");

  // Section 4
  section("File destruction");

  // Section 5
  section("Directory destruction");

  // Section 6
  section("File name manipulation");
  
  // Test the directory and regular expression search features
  section(1);
  const std::string base_dir = oe.getStormmSourcePath() + osc + "test" + osc + "FileManagement";
  const std::string xm_dir = base_dir + osc + "example";
  TestPriority xm_priority = (getDrivePathType(xm_dir) == DrivePathType::DIRECTORY) ?
                             TestPriority::CRITICAL : TestPriority::ABORT;
  const std::vector<std::string> xm_fi = listFilesInPath(xm_dir, SearchStyle::RECURSIVE);
  check(xm_fi.size(), RelationalOperator::EQUAL, 14, "The number of files found in " + xm_dir +
	" does not meet expectations.", xm_priority);
  const std::string xm_regexp = xm_dir + osc + ".*" + osc + "some_text.txt";
  const std::vector<std::string> xm_st_nr_fi = listFilesInPath(xm_regexp,
                                                               SearchStyle::NON_RECURSIVE);
  check(xm_st_nr_fi.size(), RelationalOperator::EQUAL, 7, "The number of files named "
        "\"some_text.txt\" found in the first sub-level of " + xm_dir + " does not meet "
        "expectations.", xm_priority);
  const std::string xm_regexp2 = xm_dir + osc + ".*" + osc + ".*" + osc + "some_text.txt";
  const std::vector<std::string> xm_st_r_fi = listFilesInPath(xm_regexp2, SearchStyle::RECURSIVE);
  check(xm_st_r_fi.size(), RelationalOperator::EQUAL, 1, "The number of files named "
        "\"some_text.txt\" found in all sub-levels of " + xm_dir + " does not meet expectations.",
        xm_priority);

  // Try creating some directories (inside the tmpdir), then make some files in it using various
  // opening protocols.
  section(2);
  TestPriority tmpdir_priority = (oe.getTemporaryDirectoryAccess()) ? TestPriority::CRITICAL :
                                                                      TestPriority::ABORT;
  if (tmpdir_priority == TestPriority::ABORT) {
    rtWarn("An unwriteable temporary directory implies that some of the subsequent tests must be "
           "skipped.  Check the $STORMM_TMPDIR environment variable and set it to a directory "
           "where you have write permissions.", "test_file_management");
  }
  std::string new_folder = oe.getTemporaryDirectoryPath() + osc + "trial";
  while (getDrivePathType(new_folder) == DrivePathType::DIRECTORY) {
    new_folder += "_x";
  }
  const std::string nested_folder = new_folder + osc + "nest" + osc + "sblev";
  std::vector<std::string> my_folders;
  if (oe.getTemporaryDirectoryAccess()) {

    // This effort, appending the lists of created directories to a list for removal later,
    // could be accomplished by dumping them on the test environment stack for automatic removal
    // at the end of the program, but I want to test the directory removal explicitly.
    std::vector<std::string> additions = stormmMkdir(new_folder);
    my_folders.insert(my_folders.end(), additions.begin(), additions.end());
    additions = stormmMkdir(nested_folder);
    my_folders.insert(my_folders.end(), additions.begin(), additions.end());
  }
  check(getDrivePathType(new_folder) == DrivePathType::DIRECTORY, "Failed to create a folder "
	"inside the temporary directory.", tmpdir_priority);
  CheckResult nested_dir_works = check(getDrivePathType(nested_folder) == DrivePathType::DIRECTORY,
                                       "Failed to create a nested folder inside the temporary "
                                       "directory.", tmpdir_priority);
  section(3);
  const std::string simple_file_name = nested_folder + osc + "scribble_file.txt";
  TestPriority fi_priority;
  if (nested_dir_works == CheckResult::SUCCESS) {
    fi_priority = TestPriority::CRITICAL;
    std::ofstream foutp = openOutputFile(simple_file_name, PrintSituation::OPEN_NEW, "Open a "
                                         "simple file de novo to test basic output capabilities.");
    const std::string buffer("Sequence of characters.\n");
    foutp.write(buffer.c_str(), buffer.size());
    foutp.close();
  }
  else {
    fi_priority = TestPriority::ABORT;
  }
  check(getDrivePathType(simple_file_name) == DrivePathType::FILE, "The file " + simple_file_name +
	" was not created as expected.", fi_priority);
  std::string before, after;
  const std::string path1("./trj/abc/pos/my.crd");
  splitPath(path1, &before, &after);
  check(before, RelationalOperator::EQUAL, "./trj/abc/pos/my", "Path splitting with " + path1 +
        " fails.");
  check(after, RelationalOperator::EQUAL, "crd", "Path splitting with " + path1 + " fails.");
  const std::string path2("/trj/abc/pos/mycrd");
  splitPath(path2, &before, &after);
  check(before, RelationalOperator::EQUAL, "/trj/abc/pos/mycrd", "Path splitting with " + path2 +
        " fails.");
  check(after, RelationalOperator::EQUAL, "", "Path splitting with " + path2 + " fails.");
  const std::string path3(".crd");
  splitPath(path3, &before, &after);
  check(before, RelationalOperator::EQUAL, "", "Path splitting with " + path3 + " fails.");
  check(after, RelationalOperator::EQUAL, "crd", "Path splitting with " + path3 + " fails.");

  // There is a possibility that the temporary directory could not be created, which would
  // negate a number of tests but also prevent the successful creation or reading of a TextFile
  // object.  The TextFile object requires the name of a valid file for construction, so it
  // must get protection of its own to ensure that trying to create it does not cause the test
  // program to exit with an uncaught exception.  In order to run tests on the TextFile object
  // without failing to track them in the overall tally, put a placeholder test on the other
  // side of the protective condition that will log as having been skipped if the real but
  // protected test cannot be run.
  if (nested_dir_works == CheckResult::SUCCESS) {
    std::ofstream foutp = openOutputFile(simple_file_name, PrintSituation::APPEND, "Open a simple "
                                         "file for appending more text.");
    const std::string buffer("Another set of characters.\n");
    foutp.write(buffer.c_str(), buffer.size());
    foutp.close();
    TextFile tf(simple_file_name);
    check(tf.getLineCount(), RelationalOperator::EQUAL, 2, "The number of lines in " +
          simple_file_name + " does not match expectations.", fi_priority);
  }
  else {
    check(true, "Placeholder test", fi_priority);
  }
  CHECK_THROWS_SOFT(std::ofstream foutp =
                    openOutputFile(simple_file_name, PrintSituation::OPEN_NEW, "Reopen a simple "
                                   "file using the same name as one that already exists."),
                    "A new file was opened under OPEN_NEW despite one already existing with the "
                    "same name.", fi_priority);

  // Test file removal
  section(4);
  if (nested_dir_works == CheckResult::SUCCESS) {
    removeFile(simple_file_name);
  }
  check(getDrivePathType(simple_file_name) == DrivePathType::REGEXP, "The file " +
        simple_file_name + " was not removed.");
  CHECK_THROWS_SOFT(removeFile(simple_file_name, ExceptionResponse::DIE), "The file " +
                    simple_file_name + " no longer exists, but attempting to remove it did not "
                    "throw an error.", fi_priority);

  // Test directory cleanup
  section(5);
  if (nested_dir_works == CheckResult::SUCCESS) {
    stormmBatchRmdir(my_folders);
  }
  check(getDrivePathType(nested_folder) == DrivePathType::REGEXP, "The nested directory " +
        nested_folder + " was not removed.", tmpdir_priority);
  check(getDrivePathType(new_folder) == DrivePathType::REGEXP, "The directory " +
        new_folder + " was not removed.", tmpdir_priority);
  CHECK_THROWS_SOFT(stormmRmdir(nested_folder, ExceptionResponse::DIE), "Attempted to remove a "
                    "directory that should no longer exist.", tmpdir_priority);

  // Test one more thing: create a very deep directory, then scramble the directories to remove
  // along its path to ensure that stormmBatchRmdir can still take care of them.
  section(2);
  const std::string super_nested_folder = nested_folder + osc + "a" + osc + "z" + osc + "tbn";
  if (oe.getTemporaryDirectoryAccess()) {
    stormmMkdir(new_folder);
    my_folders = stormmMkdir(super_nested_folder);
  }
  check(my_folders.size(), RelationalOperator::EQUAL, 5, "The number of new directories needed "
	"to create " + super_nested_folder + " within " + new_folder + " is incorrect.",
	tmpdir_priority);
  section(5);
  if (oe.getTemporaryDirectoryAccess() && my_folders.size() >= 5) {
    std::swap(my_folders[0], my_folders[4]);
    std::swap(my_folders[1], my_folders[3]);
    std::swap(my_folders[2], my_folders[0]);
    stormmBatchRmdir(my_folders);
    stormmRmdir(new_folder);
  }
  check(getDrivePathType(new_folder) == DrivePathType::REGEXP, "The directory " +
        new_folder + " was not removed after creating a deep nested folder within it and "
        "scrambling the tree of subdirectories.", tmpdir_priority);

  // Test file name mechanics
  section(6);
  const std::string jungle = oe.getStormmSourcePath() + osc + "its" + osc + "a" + osc + "jungle" +
                             osc + "out" + osc + "there.book";
  check(getBaseName(jungle), RelationalOperator::EQUAL, std::string("there.book"), "The base "
        "name of a file path is not returned correctly.");
  const std::string sea_creature = substituteNameExtension("sea.turtle", "snake");
  check(sea_creature, RelationalOperator::EQUAL, std::string("sea.snake"), "Substitution of a "
        "file name extension does not work as expected.");
  const std::string island = substituteNameExtension("bikini", "atoll");
  check(island, RelationalOperator::EQUAL, std::string("bikini.atoll"), "Substitution of a "
        "file name extension does not work as expected when no extension is initially present.");
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return countGlobalTestFailures();
}
