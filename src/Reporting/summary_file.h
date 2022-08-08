// -*-c++-*-
#ifndef STORMM_SUMMARY_FILE_H
#define STORMM_SUMMARY_FILE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"

namespace stormm {
namespace review {

/// \brief Default width for STORMM output files intended for human parsing.
constexpr int default_output_file_width = 80;
  
/// \brief Find the ideal format width for an output stream.  This decides based on whether the
///        output is going to the terminal (std::cout) or a file.  In the former case, STORMM
///        will try to print to the terminal's full width.  In the latter case, STORMM will use a
///        preset width given by the default_output_file_width constant.
int findFormatWidth(const std::ostream *foutp);
  
/// \brief Print the "splash screen" to stdout or to a file.  The splash screen provides the
///        authorship, origin of the software, and customary disclaimer.
///
/// \param foutp  Writeable file stream to accept the output
void stormmSplash(std::ostream *foutp = &std::cout);

/// \brief Summarize the work done in a particular STORMM run.
///
/// \param description  Content of the header to write, a description of the program or activity
/// \param foutp        Writeable file stream for the output file
/// \param file_width   Width of the file to write (default -1 will choose the default file width)
void summaryHeader(const std::string &description, std::ofstream *foutp, int file_width = -1);
 
} // namespace review
} // namespace stormm


#endif
