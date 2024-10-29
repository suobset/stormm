// -*-c++-*-
#ifndef STORMM_PROGRESSBAR_H
#define STORMM_PROGRESSBAR_H

#include <iostream>
#include <ostream>
#include <string>
#include <stdexcept>

namespace stormm {
namespace reporting {

class ProgressBar {
public:
  /// \brief The constructor for a ProgressBar.
  ProgressBar(int cycle_count_in = 0, bool show_bar_in = true,
              std::ostream &output_in = std::cout);

  /// \brief Get the current state of the progress bar, for inspection.
  const std::string& getState() const;
  
  /// \brief  Initializes a new Progress Bar from scratch, taking in user-specified information.
  ///         To be used when we want to do a hard-refresh on the existing ProgressBar object.
  ///
  /// \param  out     Determines which C standard output stream to use. Default is cout, can
  ///                 use cerr or C++ methods such as printf.
  void initialize(int n = 0, bool show_bar_in = true, std::ostream &out = std::cout);

  /// \brief  Resets the current instance of the ProgressBar object (sets percentage to 0).
  ///         To be used before a new loop, with every other setting remaining the same.
  void reset();
  
  /// \brief  Set a new number of iterations for an existing ProgressBar object.
  ///
  /// \param iter The new number of iterations for the ProgressBar object
  void setIterations(int iter);

  /// \brief  Set a new "done" char for a ProgressBar object.  This is the char that appears when
  ///         a percentage is done in the ProgressBar.
  ///
  /// \param sym  The string to use for a Done char (default is '#')
  void setDoneChar(const std::string &sym);

  /// \brief  Set a new "todo" char for a ProgressBar object.  This is the char that populates the
  ///         remaining of the ProgressBar object.
  ///
  /// \param sym  The string to use for a Todo char (default is ' ')
  void setTodoChar(const std::string &sym);

  /// \brief  Set a new opening bracket for a ProgressBar object.  This is the char during the
  ///         start of a ProgressBar object.
  ///
  /// \param sym  The string to use for an opening char of the bar (default is '[')
  void setOpeningBracketChar(const std::string &sym);
  
  /// \brief  Set a new closing bracket for a ProgressBar object.  This is the char to use for a
  ///         closing char of the bar (default is ']')
  void setClosingBracketChar(const std::string &sym);

  /// \brief  Toggles the visibility of the ProgressBar. Default is true.
  ///
  /// \param flag  Boolean value to show/hide the ProgressBar. If hidden, the percentage value is
  ///              still shown.
  void showBar(bool flag = true);
  
  /// \brief  Function to set the output stream of the current ProgressBar object.  Default is
  ///         cout, can be changed to any std::ostream &object.
  ///
  /// \param stream   The standard output stream to use for rendering the progress bar.
  void setOutputStream(std::ostream &stream);

  /// \brief  Allows the user to set a custom terminal width for the ProgressBar.
  ///         If this function is called, it will override the dynamically inferred width.
  ///
  /// \param width  The new width of the terminal in characters
  void setTerminalWidth(int width);

  /// \brief  Function to allocate the bar_contents string for the ProgressBar object.
  ///         This is called when the bar is first initialized or when 
  ///         the terminal width is changed.
  void allocateBarContents();

  /// \brief  Function to update the ProgressBar, incrementing the number of iterations by 1,
  ///         calculating the appropriate percentage, and rendering the ProgressBar in the
  ///         terminal.
  void update();
  
private:

  int progress;                   ///< The amount of iterations that have been completed by 
                                  ///< the ProgressBar
  int cycle_count;                ///< The number of iterations the ProgressBar has to go 
                                  ///< through in total
  int last_percent;               ///< The percentage value that the ProgressBar calculated so far
  bool show_bar;                  ///< Boolean value determining if we want to show the 
                                  ///< visual ProgressBar
  bool update_called;             ///< Boolean flag to determine if the bar has been rendered 
                                  ///< at least once
  std::string finished_mark;      ///< String that appears when a percentage is done
  std::string todo_mark;          ///< String that appears for the rest of the ProgressBar
  std::string prgb_open_bracket;  ///< String at the start of a ProgressBar
  std::string prgb_close_bracket; ///< String at the end of a ProgressBar
  std::ostream* output;           ///< C std output stream throuth which to render the bar
  int terminal_width;             ///< The detected width of the terminal screen
  std::string bar_contents;       ///< The current state of the progress bar.  Updates that edit
                                  ///<   the state of the progress bar will alter a minimal number
                                  ///<   of characters in this string, and any display of the
                                  ///<   progress bar will dump the entire contents of this string
                                  ///<   to the chosen output mode.
  
  /// \brief  Function to dynamically infer terminal width to draw an appropriate ProgressBar.
  ///         This function is called implicitly and does not require user input.
  void updateTerminalWidth();
};

} // namespace reporting
} // namespace stormm

#endif // STORMM_PROGRESSBAR_H

