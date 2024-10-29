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
ProgressBar::ProgressBar(int cycle_count_in, bool show_bar_in, std::ostream& output_in) :
    progress{0}, cycle_count{cycle_count_in}, last_percent{0}, show_bar{show_bar_in}, 
    update_called{false},
    finished_mark{"#"},
    todo_mark{" "}, 
    prgb_open_bracket{"["},
    prgb_close_bracket{"]"},
    output{&output_in},
    terminal_width{default_output_file_width},
    bar_contents{}
{
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
const std::string& ProgressBar::getState() const {
  return bar_contents;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::initialize(int cycle_count_in, bool show_bar_in, std::ostream& output_in) {
  progress = 0;
  cycle_count = cycle_count_in;
  last_percent = 0;
  show_bar = show_bar_in;
  update_called = false;
  finished_mark = "#";
  todo_mark = " ";
  prgb_open_bracket = "[";
  prgb_close_bracket = "]";
  output = &output_in;
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::reset() {
  progress = 0;
  update_called = false;
  last_percent = 0;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setIterations(int iter) {
  if (iter <= 0) {
    rtErr("The number of iterations must be positive.", "ProgressBar", "setIterations");
  }
  cycle_count = iter;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setDoneChar(const std::string& sym) {
  finished_mark = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTodoChar(const std::string& sym) {
  todo_mark = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOpeningBracketChar(const std::string& sym) {
  prgb_open_bracket = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setClosingBracketChar(const std::string& sym) {
  prgb_close_bracket = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::showBar(bool flag) {
  show_bar = flag;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOutputStream(std::ostream& stream) {
  output = &stream;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTerminalWidth(int width) {
  if (width > 0) {
    terminal_width = width;
    allocateBarContents();
  } else {
    rtErr("Terminal width must be positive.", "ProgressBar", "setTerminalWidth");
  }
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::allocateBarContents() {
    // Calculate the width of the progress bar area
    int barWidth = terminal_width - (prgb_open_bracket.length() + prgb_close_bracket.length() + 7); 
    barWidth = std::max(barWidth, 10);  // Ensure minimum width

    // Resize the contents of bar_contents to fit terminal width
    bar_contents.clear();
    bar_contents.resize(terminal_width);

    // Initialize with the structure: [###     ] 50%
    size_t pos = 0;

    // Opening bracket
    for (char c : prgb_open_bracket) {
        bar_contents[pos++] = c;
    }

    // Fill with initial todo marks (space by default)
    for (int i = 0; i < barWidth; ++i) {
        bar_contents[pos++] = todo_mark[0];
    }

    // Closing bracket
    for (char c : prgb_close_bracket) {
        bar_contents[pos++] = c;
    }

    // Add percentage placeholder (spaces initially)
    bar_contents[pos++] = ' ';
    for (int i = 0; i < 4; ++i) {  // Reserve 3 characters for "XXX%" and a space
        bar_contents[pos++] = ' ';
    }

    // This preps the entire progress bar with brackets and placeholders
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::update() {
    if (cycle_count == 0) {
        rtErr("The number of cycles has not been set.", "ProgressBar", "update");
    }

    if (!update_called) {
        update_called = true;
    }

    // Calculate the percentage of progress
    const int percent = progress * 100 / (cycle_count - 1);

    // Only update if the percentage has changed
    if (percent != last_percent) {
        int barWidth = terminal_width - (prgb_open_bracket.length() + prgb_close_bracket.length() + 7);
        int numDoneChars = percent * barWidth / 100;
        int prevDoneChars = last_percent * barWidth / 100;

        // Update only the difference between the old and new percentage
        if (show_bar) {
            for (int i = prevDoneChars; i < numDoneChars; ++i) {
                bar_contents[prgb_open_bracket.length() + i] = finished_mark[0];
            }
            for (int i = numDoneChars; i < barWidth; ++i) {
                bar_contents[prgb_open_bracket.length() + i] = todo_mark[0];
            }

            // Update the percentage display
            std::string percentString = std::to_string(percent) + "%";
            size_t percentPos = bar_contents.size() - 4;
            for (char c : percentString) {
                bar_contents[percentPos++] = c;
            }

            // Dump the fully updated bar_contents to the output
            *output << "\r" << bar_contents << std::flush;
        } else {
            // If bar is hidden, only show the percentage
            *output << "\r" << std::to_string(percent) + "%" << std::flush;
        }

        last_percent = percent;  // Update the last_percent to the current one
    }

    // Increment progress
    ++progress;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::updateTerminalWidth() {
    if (isatty(fileno(stdout))) {
        struct winsize ws;
        if (ioctl(fileno(stdout), TIOCGWINSZ, &ws) == 0) {
            terminal_width = ws.ws_col;
        }
        else {
            terminal_width = default_output_file_width;
        }
    }
    else {
        terminal_width = default_output_file_width;
    }

    // Re-allocate the bar_contents to match new terminal width
    allocateBarContents();

    // Refill completion symbols based on last_percent
    int barWidth = terminal_width - (prgb_open_bracket.length() + prgb_close_bracket.length() + 7);
    int numDoneChars = last_percent * barWidth / 100;

    for (int i = 0; i < numDoneChars; ++i) {
        bar_contents[prgb_open_bracket.length() + i] = finished_mark[0];
    }

    // Update the percentage text
    std::string percentString = std::to_string(last_percent) + "%";
    size_t percentPos = bar_contents.size() - 4; // Assuming the last 4 chars are reserved for "XXX%"
    for (char c : percentString) {
        bar_contents[percentPos++] = c;
    }
}

} // namespace reporting
} // namespace stormm