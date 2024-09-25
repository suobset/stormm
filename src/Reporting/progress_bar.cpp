#include <unistd.h>
#include <sys/ioctl.h>
#include "copyright.h"
#include "error_format.h"
#include "summary_file.h"
#include "progress_bar.h"

namespace stormm {
namespace reporting {

using review::default_output_file_width;
  
//-------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(int n, bool showbar, std::ostream& out) :
    progress(0), nCycles(n), lastPercent(0), doShowBar(showbar), 
    updateIsCalled(false), doneChar("#"), todoChar(" "), 
    openingBracketChar("["), closingBracketChar("]"), output(&out)
{
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::initialize(int n, bool showbar, std::ostream& out) {
  progress = 0;
  nCycles = n;
  lastPercent = 0;
  doShowBar = showbar;
  updateIsCalled = false;
  doneChar = "#";
  todoChar = " ";
  openingBracketChar = "[";
  closingBracketChar = "]";
  output = &out;
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::reset() {
  progress = 0;
  updateIsCalled = false;
  lastPercent = 0;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setIterations(int iter) {
  if (iter <= 0) {
    rtErr("The number of iterations must be positive.", "ProgressBar", "setIterations");
  }
  nCycles = iter;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setDoneChar(const std::string& sym) {
  doneChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTodoChar(const std::string& sym) {
  todoChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOpeningBracketChar(const std::string& sym) {
  openingBracketChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setClosingBracketChar(const std::string& sym) {
  closingBracketChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::showBar(bool flag) {
  doShowBar = flag;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOutputStream(std::ostream& stream) {
  output = &stream;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTerminalWidth(int width) {
  if (width > 0) {
    termWidth = width;
  } else {
    rtErr("Terminal width must be positive.", "ProgressBar", "setTerminalWidth");
  }
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::update() {
    if (nCycles == 0) {
        rtErr("The number of cycles has not been set.", "ProgressBar", "update");
    }
    if (!updateIsCalled) {
        updateIsCalled = true;
    }

    // Calculate the percentage of progress
    const int percent = progress * 100 / (nCycles - 1);
    
    if (doShowBar) {
        int barWidth = termWidth - (openingBracketChar.length() + closingBracketChar.length() + 7);
        barWidth = std::max(barWidth, 10);  // Ensure minimum bar width

        // Preallocate the string to avoid dynamic allocations
        std::string progressBarString;
        progressBarString.resize(termWidth);

        // Place the carriage return character to move to the start of the line
        progressBarString[0] = '\r';

        // Fill in the opening bracket character
        size_t position = 1;
        for (char c : openingBracketChar) {
            progressBarString[position++] = c;
        }

        // Calculate the number of done and todo characters
        int numDoneChars = percent * barWidth / 100;

        // Fill the done characters in the string
        for (int i = 0; i < numDoneChars; ++i) {
            progressBarString[position++] = doneChar[0];
        }

        // Fill the remaining part of the bar with todo characters
        for (int i = numDoneChars; i < barWidth; ++i) {
            progressBarString[position++] = todoChar[0];
        }

        // Fill in the closing bracket character
        for (char c : closingBracketChar) {
            progressBarString[position++] = c;
        }

        // Add a space and the percentage display
        progressBarString[position++] = ' ';
        std::string percentString = std::to_string(percent) + "%";
        for (char c : percentString) {
            progressBarString[position++] = c;
        }

        // Resize the string to the actual length
        progressBarString.resize(position);

        // Output the entire string at once, overwriting the current line
        *output << progressBarString << std::flush;
    } else {
        // If the bar is hidden, only show the percentage.
        std::string percentString = std::to_string(percent) + "%";
        *output << "\r" << percentString << std::flush;
    }

    // Update internal counters
    lastPercent = percent;
    ++progress;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::updateTerminalWidth() {
  if (isatty(fileno(stdout))) {
    struct winsize ws;
    if (ioctl(fileno(stdout), TIOCGWINSZ, &ws) == 0) {
      termWidth = ws.ws_col;
    }
    else {
      termWidth = default_output_file_width;
    }
  }
  else {
    termWidth = default_output_file_width;
  }
}

} // namespace reporting
} // namespace stormm
