#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <regex>
#include <string.h>
#include "file_listing.h"
#include "Reporting/error_format.h"

namespace omni {
namespace diskutil {

using errors::rtErr;

//-------------------------------------------------------------------------------------------------
char osSeparator() {
  switch (detected_os) {
  case (OperatingSystem::LINUX):
  case (OperatingSystem::UNIX):
  case (OperatingSystem::MAC_OS):
    return '/';
  case (OperatingSystem::WINDOWS):
    return '\\';
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
DrivePathType getDrivePathType(const std::string &path) {

  // Try getting the status.  If this fails, the path is probably a regular expression.
  struct stat path_stat;
  if (stat(path.c_str(), &path_stat) == 0) {
    if (S_ISREG(path_stat.st_mode) == 1) {
      return DrivePathType::FILE;
    }
    else if (S_ISDIR(path_stat.st_mode) == 1) {
      return DrivePathType::DIRECTORY;
    }
  }

  // Interpret anything else as a regular expression
  return DrivePathType::REGEXP;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> listDirectory(const std::string &dir_path, const SearchStyle r_option,
                                       const DrivePathType entity_kind) {
  
  // Attempt to open the directory
  std::vector<std::string> ls_result;
  DIR *dir;
  if ((dir = opendir(dir_path.c_str())) == NULL) {
    return ls_result;
  }
  struct dirent *ent;
  while ((ent = readdir(dir)) != NULL) {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
      continue;
    }
    
    // Test whether the result at hand is, itself, a directory.
    std::string nested_dir_path = dir_path + osSeparator() + std::string(ent->d_name);
    DrivePathType item_type = getDrivePathType(nested_dir_path);
    if (item_type == DrivePathType::DIRECTORY) {
      switch (r_option) {
      case SearchStyle::RECURSIVE:
        {
          std::vector<std::string> nested_result = listDirectory(nested_dir_path, r_option,
                                                                 entity_kind);
          ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
          break;
        }
      case SearchStyle::NON_RECURSIVE:
        if (entity_kind == DrivePathType::DIRECTORY) {
          ls_result.push_back(nested_dir_path);
        }
        break;
      }
    }
    else if (item_type == DrivePathType::FILE) {
      if (entity_kind == DrivePathType::FILE) {
        std::string new_name(ent->d_name);
        ls_result.push_back(dir_path + osSeparator() + new_name);
      }
    }
  }
  closedir(dir);

  return ls_result;
}

//-------------------------------------------------------------------------------------------------
bool pathIsAbsolute(const std::string &path) {

  // Assume that a zero-length path is a relative path to the current directory
  if (path.size() == 0) {
    return false;
  }
  const char sep_char = osSeparator();
  switch (detected_os) {
  case (OperatingSystem::LINUX):
  case (OperatingSystem::UNIX):
  case (OperatingSystem::MAC_OS):
    return (path[0] == sep_char);
  case (OperatingSystem::WINDOWS):
    return (path.size() > 2 && path[0] >= 'A' && path[0] <= 'Z' && path[1] == ':' &&
            path[2] == sep_char);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string makePathAbsolute(const std::string &path) {

  // Return immediately if there is nothing to be done
  if (pathIsAbsolute(path)) {
    return path;
  }

  // Detect the current directory
  char tmp_c[512];
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  if (GetCurrentDirectory(512, tmp_c) == nullptr) {
    rtErr("Encountered getcwd() error.", "makePathAbsolute");
  }
#else
  if (getcwd(tmp_c, 512) == nullptr) {
    rtErr("Encountered getcwd() error.", "makePathAbsolute");
  }
#endif
  const std::string program_cwd(tmp_c);
  return program_cwd + osSeparator() + path;
}

//-------------------------------------------------------------------------------------------------
std::string getNormPath(const std::string &path) {
  const char sep_char = osSeparator();
  int i = path.size() - 1;
  while (i >= 0 && path[i] == sep_char) {
    i--;
  }
  return path.substr(0, i + 1);
}

//-------------------------------------------------------------------------------------------------
std::string getBaseName(const std::string &path) {
  const char sep_char = osSeparator();
  const int plength = path.size();
  int last_separator = 0;
  for (int i = 0; i < plength; i++) {
    if (path[i] == sep_char) {
      last_separator = i + 1;
    }
  }
  return path.substr(last_separator, plength - last_separator);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> separatePath(const std::string &path) {

  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  const bool absolute_path = (path[0] == sep_char);

  // Break the regular expression into a series of directory levels
  std::vector<std::string> levels;
  int i = 0;
  const int p_length = path.size();
  while (i < p_length) {

    // Take in characters until the next separator
    std::string new_level;
    while (i < p_length) {
      char t_char = path[i];
      i++;
      if (t_char != sep_char) {
        new_level += t_char;
      }
      else {
        break;
      }
    }
    if (new_level.size() > 0) {
      levels.push_back(new_level);
    }
  }

  return levels;
}

//-------------------------------------------------------------------------------------------------
void splitPath(const std::string &path, std::string *before, std::string *after) {
  const int pathsize = path.size();
  int dotpos = pathsize;
  for (int i = pathsize - 1; i >= 0; i--) {
    if (path[i] == '.') {
      dotpos = i;
      break;
    }
  }
  *after = (dotpos < pathsize) ? path.substr(dotpos + 1, pathsize - dotpos - 1) : "";
  *before = (dotpos > 0) ? path.substr(0, dotpos) : "";
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> listFilesInPath(const std::string &regexp_path, SearchStyle r_option) {
  std::vector<std::string> ls_result;

  // Detect a regular file or directory
  DrivePathType td = getDrivePathType(regexp_path);
  switch (td) {
  case DrivePathType::FILE:
    ls_result.push_back(regexp_path);
    return ls_result;
  case DrivePathType::DIRECTORY:
    {
      std::vector<std::string> nested_result = listDirectory(regexp_path, r_option);
      ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
      return ls_result;
    }
  case DrivePathType::REGEXP:
    break;
  }
    
  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  const bool absolute_path = (regexp_path[0] == sep_char);
  
  // Break the regular expression into a series of directory levels
  const std::vector<std::string> levels = separatePath(regexp_path);

  // For each level, look at the appropriate directory and search for matching subdirectories
  const int n_levels = levels.size();  
  std::string partial_path = (absolute_path) ? "" : ".";
  std::regex complete_expr(regexp_path);
  for (int i = 0; i < n_levels; i++) {
    std::string test_path = partial_path + sep_char + levels[i];
    DrivePathType pd = getDrivePathType(test_path);
    switch (pd) {
    case DrivePathType::FILE:
      if (i == n_levels - 1) {
        if (regex_match(test_path, complete_expr)) {
          ls_result.push_back(test_path);
        }
      }
      break;
    case DrivePathType::DIRECTORY:
      partial_path = test_path;
      break;
    case DrivePathType::REGEXP:
      {
        std::regex level_expr(levels[i]);
        std::vector<std::string> test_paths = listDirectory(partial_path,
                                                            SearchStyle::NON_RECURSIVE,
                                                            DrivePathType::FILE);
        std::vector<std::string> additional_paths = listDirectory(partial_path,
                                                                  SearchStyle::NON_RECURSIVE,
                                                                  DrivePathType::DIRECTORY);
        test_paths.insert(test_paths.end(), additional_paths.begin(), additional_paths.end());
        const int n_paths = test_paths.size();
        for (int j = 0 ; j < n_paths; j++) {
          
          // Test that the deepest level of this path conforms to the regular expression at hand.
          if (regex_match(getBaseName(getNormPath(test_paths[j])), level_expr)) {

            // This path is valid.  If there are more levels to explore, continue adding to the
            // partial path.  Otherwise, if the path is a file, add it to the list, or if it is
            // a directory, list its files, using recursion as requested.
            std::string refined_path = test_paths[j];
            for (int k = i + 1; k < n_levels; k++) {
              refined_path += sep_char + levels[k];
            }
            std::vector<std::string> nested_result = listFilesInPath(refined_path, r_option);
            ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
          }
        }

        // Continue incrementing the partial path
        partial_path = test_path;
      }
      break;
    }
  }
  
  return ls_result;
}

//-------------------------------------------------------------------------------------------------
std::string findOmniPath(const std::string &path, const std::string &extension) {
  const std::string abs_path = pathIsAbsolute(path) ? path : makePathAbsolute(path);
  const std::vector<std::string> pieces = separatePath(abs_path);
  const int n_piece = pieces.size();
  const char sep_char = osSeparator();
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  std::string partial_path("");
#else
  std::string partial_path(1, sep_char);
#endif
  for (int i = 0; i < n_piece; i++) {
    const std::string guess_filename = partial_path + sep_char + extension;
    if (getDrivePathType(guess_filename) == DrivePathType::FILE) {
      return partial_path;
    }
    if (i > 0) {
      partial_path += sep_char;
    }
    partial_path += pieces[i];
  }

  // If this point is reached, a good guess as to ${OMNI_HOME} could not be made 
  return std::string("");
}

} // namespace diskutil
} // namespace omni
