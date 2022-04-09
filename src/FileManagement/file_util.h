// -*-c++-*-
#ifndef OMNI_FILE_UTIL_H
#define OMNI_FILE_UTIL_H

#include <fstream>
#include <iostream>
#include <string>
#include "Constants/behavior.h"

namespace omni {
namespace diskutil {

using constants::ExceptionResponse;

/// \brief Enumerate the situations that might be encountered when writing a trajectory file
enum class PrintSituation {
  OPEN_NEW,   ///< Expect to open a new file, and only open a new file if no file by the given name
              ///<   already exists
  APPEND,     ///< Append an existing file, or open a new file if no file by the given name exists
  OVERWRITE,  ///< Always open a new file, overwriting any existing file by the same name
  UNKNOWN     ///< Option to let the program decide what printing behavior to use, in lieu of user
              ///<   or developer input
};

/// \brief Format of the file (just a fancy way of coding "binary" or "ascii")
enum class DataFormat {
  ASCII,   ///< Ascii, formatted and human-readable output
  BINARY   ///< Binary, direct-to-disk copy of numbers
};

/// \brief Open a file for output writing.  This encapsulates error messages in the event that
///        the file cannot be opened as expected.
///
/// \param filename     Name of the file to write
/// \param expectation  Conditions to look for in order to successfully open the file for writing
/// \param style        The kind of file to write, ascii or binary
/// \param description  Description of the file, to help identify problems if it cannot be written
std::ofstream openOutputFile(const std::string &filename,
                             PrintSituation expectation = PrintSituation::OPEN_NEW,
                             const std::string &description = std::string(""),
                             DataFormat style = DataFormat::ASCII);

/// \brief Remove a file.  This encapsulates error messages in the event that the file cannot be
///        removed as requested.
///
/// \param filename  Name of the file to remove
/// \param policy    Indicates what to do if the file cannot be removed for some reason
int removeFile(const std::string &filename, ExceptionResponse policy = ExceptionResponse::WARN);

} // namespace diskutil
} // namespace omni

#endif
