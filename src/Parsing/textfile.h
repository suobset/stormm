// -*-c++-*-
#ifndef OMNI_TEXTFILE_H
#define OMNI_TEXTFILE_H

#include <string>
#include <vector>

namespace omni {
namespace parse {

/// \brief Many searches in TextFile objects will begin at a particular line.  If the query is not
///        found, it may be necessary to wrap the search back to the beginning and continue until
///        the original starting line.  This will indicate whether to do that.
enum class WrapTextSearch {
  NO, YES
};

/// \brief Differentiate between text data originating on a disk and in RAM
enum class TextOrigin {
  DISK, RAM
};

/// \brief Structure for translating a text file into a compact, rapidly parsable vector of
///        characters in CPU RAM.  The struct contains two nested struct definitions, for a Reader
///        and a Writer.  Private objects of these structs can then be accessed with eponymous
///        data() getter functions, returning the appropriate kind of access depending on the
///        const-ness of the TextFile object itself.
struct TextFile {

  /// Abstract for read-only access
  struct Reader {
    int line_count;
    const int* line_limits;
    const char* text;
    const std::string file_name;
  };

  /// \brief Constructor for taking an ascii file or a very long, formatted string and transforming
  ///        it into a std::vector of characters with line limits recorded.
  ///
  /// Overloaded:
  ///   - Construct an empty object after taking no arguments
  ///   - Construct a complete object from a named file (throws an exception if the file does not
  ///     exist)
  ///
  /// \param file_name   Name of the input file
  /// \param source      Origin of the text--disk or RAM
  /// \param content     Content for the TextFile, and perhaps later an ASCII text file, to hold
  /// \param caller      (Optional) name of the calling function
  TextFile();
  TextFile(const std::string &file_name, TextOrigin source = TextOrigin::DISK,
           const std::string &content = std::string(""),
           const std::string &caller = std::string(""));

  /// \brief Default destructor
  ~TextFile() = default;

  /// \brief Get the name of the original file.
  std::string getFileName() const;

  /// \brief Get the line count of a text file after converting it to a character vector in memory.
  int getLineCount() const;

  /// \brief Get one of the line limits of a text file converted to a character vector in memory.
  ///
  /// \param index   The line index.  The array contains one more indices than the text file had
  ///                lines, to allow the end of the last line to be determined.
  int getLineLimits(int index) const;

  /// \brief Get one character of a text file after converting it to a character vector in memory.
  ///
  /// \param index   The character index, as ascertained by line limits and some offset
  char getText(int index) const;

  /// \brief Get a char* to a specific index in the object
  ///
  /// \param index   The character index, as ascertained by line limits and some offset
  const char* getTextPointer(int index) const;

  /// \brief Get an abstract of a text file's CPU-RAM representation, for ease of use.
  ///
  /// \param tf   The text file, as read into memory
  const Reader& data() const;

private:

  /// Name of the file that was read
  std::string orig_file;

  /// The number of lines detected in the file
  int line_count;

  /// Limits for each line's text in the concatenated character array
  std::vector<int> line_limits;

  /// The text, sans carriage returns (see line limits to determine their locations)
  std::vector<char> text;

  /// Read-only access; write access to a struct like this, via a convenient pointer struct,
  /// could be very problematic.  Will the line limits stay the same?  What if the file needs
  /// to grow?  Editing and growth of a TextFile will need to be done by public member functions
  /// of the actual object, not an abstract.
  Reader readable;

  /// \brief Set the name of the corresponding file based on the characteristics of several
  ///        constructor input variables.
  ///
  /// \param file_name   Name of the input file (if the source is DISK--otherwise if the source
  ///                    is RAM and the content is blank, then this is assumed to be content and
  ///                    a blank corresponding file name is returned)
  /// \param source      Origin of the text--disk or RAM
  /// \param content     Content for the TextFile, and perhaps later an ASCII text file, to hold
  std::string setFileName(const std::string &file_name, TextOrigin source,
                          const std::string &content);

  /// \brief Break a large, formatted string into separate lines based on carriage returns.
  ///
  /// \param text_in  The text to parse.  See setFileName and above in the TextFile documentation
  ///                 to understand how this could be either of two constructor inputs.
  void linesFromString(const std::string &text_in);
};

} // namespace parse
} // namespace omni

#endif
