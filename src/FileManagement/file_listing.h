// -*-c++-*-
#ifndef STORMM_FILE_LISTING_H
#define STORMM_FILE_LISTING_H

#include <vector>
#include <string>

namespace stormm {
namespace diskutil {

/// \brief Enumerate the possible types of directory search (this is equivalent to whether there is
///        a -r modifier after ls).
enum class SearchStyle {
  RECURSIVE, NON_RECURSIVE
};

/// \brief Enumerate whether a path can be determined to be a file, directory, or, barring any
///        other explanation, a regular expression.
enum class DrivePathType {
  FILE, DIRECTORY, REGEXP
};

/// \brief An enumerator to make note of the operating systems
enum class OperatingSystem {
  LINUX, UNIX, WINDOWS, MAC_OS
};

/// Pre-processor directives determine the separators to use in path and file concatenation
/// \{
#if defined(__linux__)
constexpr OperatingSystem detected_os = OperatingSystem::LINUX;
#elif defined(unix) || defined(__unix) || defined(__unix__)
constexpr OperatingSystem detected_os = OperatingSystem::UNIX;
#elif defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
constexpr OperatingSystem detected_os = OperatingSystem::WINDOWS;
#elif defined(__APPLE__) || defined(__MACH__)
constexpr OperatingSystem detected_os = OperatingSystem::MAC_OS;
#endif
/// \}

/// \brief Produce the correct file path separator.
char osSeparator();

/// \brief Test whether a path is a file, directory, or must instead be a regular expression.
///
/// \param path  The path of interest
DrivePathType getDrivePathType(const std::string &path);

/// \brief Given a path that has been established to be a directory, list all files within it.
///        Recursively descend into subdirectories, if requested.
///
/// \param dir_path  The path to search
/// \param r_option  Option to use recursion or not
std::vector<std::string> listDirectory(const std::string &path,
                                       SearchStyle r_option = SearchStyle::NON_RECURSIVE,
                                       DrivePathType entity_kind = DrivePathType::FILE);

/// \brief Detect whether a path appears to be absolute, based on the existence of a leading slash
///        on Unix, Linux, and Mac OS systems or one of the drive letter codes on Windows.
///
/// \param path  The path of the file or directory
bool pathIsAbsolute(const std::string &path);

/// \brief Make a path into an absolute path, modified as appropriate for the system architecture.
///
/// \param path  The name of the file or directory
std::string makePathAbsolute(const std::string &path);

/// \brief Get the normative path by removing any trailing slashes.
///
/// \param path  The path to normalize
std::string getNormPath(const std::string &path);

/// \brief Get the base name of a path by taking everything after the last OS-specific separator.
///
/// \param path  The path from which to extract a base name
std::string getBaseName(const std::string &path);

/// \brief Extract the various directory levels in a given path
///
/// \param path  The path to analyze
std::vector<std::string> separatePath(const std::string &path);

/// \brief Split a path into parts before and after the final dot (".")
///
/// \param path    The original path
/// \param before  The part of the path before the final dot, might be an empty string (returned)
/// \param after   The part of the path after the final dot, might be an empty string (returned)
void splitPath(const std::string &path, std::string *before, std::string *after);
  
/// \brief Given a path that has been established to be a regular expression, list all files it
///        could describe, recursively descend into subdirectories if requested.  The recursion
///        only kicks in once the regular expression has been interpreted into a path that is,
///        itself, a directory.  Otherwise, "A/*/B/[a-z]*" will find all subdirectories of A/ such
///        that they contain their own subdirectory B/ and within each such B/ all files beginning
///        with a lowercase letter or, if recursion is activated, descend recursively into any and
///        all subdirectories of B/ whose name begins with a lowercase letter listing all files.
///        The search string "A/*/B/[a-z]*" would not find a file A/foo.txt, but it would find
///        A/bar/B/foo.txt, and if recursion were activated it would find A/bar/B/C/D/E/foo.txt.
///
/// \param regexp_path  The regular expression to evaluate
/// \param r_option     Option to use recursion or not
std::vector<std::string> listFilesInPath(const std::string &regexp_path, SearchStyle r_option);

/// \brief Try to find the STORMM home or source directories based on a hypothetical path.  A known
///        file within the desired STORMM directory is provided to help identify the correct root.
///
/// \param path       The path that may contain clues as to the STORMM home directory
/// \param extension  Name of the file that must be found in some partial path when that partial
///                   path is the correct root
std::string findStormmPath(const std::string &path, const std::string &extension);

} // namespace parse
} // namespace stormm

#endif
