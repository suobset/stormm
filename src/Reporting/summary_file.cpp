#include <sys/ioctl.h>
#include <unistd.h>
#include "display.h"
#include "error_format.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using display::terminalHorizontalRule;
using errors::RTMessageKind;
using errors::terminalFormat;

//-------------------------------------------------------------------------------------------------
int findFormatWidth(const std::ostream *foutp) {
  int recommended_width = default_output_file_width;
  if (foutp->rdbuf() == std::cout.rdbuf()) {
    struct winsize console_dims;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);
    recommended_width = (console_dims.ws_col > 2) ? console_dims.ws_col - 1 :
                                                    default_output_file_width;
  }
  return recommended_width;  
}

//-------------------------------------------------------------------------------------------------
void stormmSplash(std::ostream *foutp) {
  int recommended_width = findFormatWidth(foutp);
  terminalHorizontalRule("+", "+", recommended_width, foutp);
  std::string buffer("Copyright 2022, Roivant Sciences\n\n");
  buffer += "Permission is hereby granted, free of charge, to any person obtaining a copy of "
            "this software and associated documentation files (the \"Software\"), to deal in the "
            "Software without restriction, including without limitation the rights to use, copy, "
            "modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, "
            "and to permit persons to whom the Software is furnished to do so, subject to the "
            "following conditions:\n\n"
            "The above copyright notice and this permission notice shall be included in all "
            "copies or substantial portions of the Software.\n\n"
            "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
            "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR "
            "A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT "
            "HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF "
            "CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE "
            "OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.";
  buffer = terminalFormat(buffer, nullptr, nullptr, 0, 0, 0, recommended_width,
                          RTMessageKind::ERROR);
  buffer += "\n";
  foutp->write(buffer.c_str(), buffer.size());
  terminalHorizontalRule("+", "+", recommended_width, foutp);
}

//-------------------------------------------------------------------------------------------------
void summaryHeader(const std::string &description, std::ofstream *foutp) {

  // Lead with the STORMM authorship and disclaimer
  stormmSplash(foutp);
  foutp->write("\n", 1);

  // Write the description, whether of a program itself or the program feature in use
  int recommended_width = findFormatWidth(foutp);
  terminalHorizontalRule("+", "+", recommended_width, foutp);
  foutp->write(description.c_str(), description.size());
  terminalHorizontalRule("+", "+", recommended_width, foutp);
}
  
} // namespace review
} // namespace stormm
