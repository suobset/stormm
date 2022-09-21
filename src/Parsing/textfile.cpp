#include <iostream>
#include <fstream>
#include "copyright.h"
#include "textfile.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
TextFileReader::TextFileReader(const int line_count_in, const int* line_limits_in,
                               const char* text_in, const std::string file_name_in) :
  line_count{line_count_in}, line_limits{line_limits_in}, text{text_in}, file_name{file_name_in}
{}

//-------------------------------------------------------------------------------------------------
TextFile::TextFile() :
    orig_file{std::string("")},
    line_count{0},
    line_limits{std::vector<int>(1, 0)},
    text{std::vector<char>()}
{}

//-------------------------------------------------------------------------------------------------
TextFile::TextFile(const std::string &file_name, const TextOrigin source,
                   const std::string &content, const std::string &caller) :
  orig_file{setFileName(file_name, source, content)},
  line_count{0},
  line_limits{},
  text{}
{
  switch (source) {
  case TextOrigin::DISK:
    {
      // Start reading from a file.  The content variable will be ignored.
      std::ifstream finp;
      finp.open(file_name.c_str());
      if (finp.is_open() == false) {
        if (caller.size() == 0) {
          rtErr(file_name + " was not found.", "TextFile");
        }
        else {
          rtErr(file_name + " was not found when called from " + caller + ".", "TextFile");
        }
      }
      std::string line;
      int total_chars = 0;
      int line_counter = 0;
      text.resize(0);
      line_limits.resize(0);
      line_limits.push_back(total_chars);
      while (std::getline(finp, line)) {
        const int line_length = line.size();
        for (int i = 0; i < line_length; i++) {
          text.push_back(line[i]);
        }
        total_chars += line_length;
        line_limits.push_back(total_chars);
        line_counter++;
      }
      line_count = line_counter;

      // Close input file
      finp.close();
    }
    break;
  case TextOrigin::RAM:
    {
      // If ther content variable is blank, this implies that the file_name variable is the
      // actual content.  Either way, scan the content for line breaks (carriage returns) and
      // assign a file name as appropriate.
      if (content.size() > 0) {
        linesFromString(content);
      }
      else {
        linesFromString(file_name);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::getFileName() const {
  return orig_file;
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLineCount() const {
  return line_count;
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLineLimits(const int index) const {
  return line_limits[index];
}

//-------------------------------------------------------------------------------------------------
char TextFile::getChar(const int index) const {
  return text[index];
}

//-------------------------------------------------------------------------------------------------
const char* TextFile::getTextPointer(const int index) const {
  return &text[index];
}

//-------------------------------------------------------------------------------------------------
const TextFileReader TextFile::data() const {
  return TextFileReader(line_count, line_limits.data(), text.data(), orig_file);
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::extractString(const int line_number, const int start_pos,
                                    const int string_length) const {
  const int actual_length = checkAvailableLength(line_number, start_pos, string_length);
  std::string buffer;
  buffer.resize(actual_length);
  const char* text_ptr = &text[line_limits[line_number] + start_pos];
  for (int i = 0; i < actual_length; i++) {
    buffer[i] = text_ptr[i];
  }
  return buffer;
}

//-------------------------------------------------------------------------------------------------
char4 TextFile::extractChar4(const int line_number, const int start_pos,
                             const int string_length) const {
  const int actual_length = checkAvailableLength(line_number, start_pos, string_length);
  const char* text_ptr = &text[line_limits[line_number] + start_pos];
  return { (actual_length > 0) ? text_ptr[0] : ' ', (actual_length > 1) ? text_ptr[1] : ' ',
           (actual_length > 2) ? text_ptr[2] : ' ', (actual_length > 3) ? text_ptr[3] : ' ' };
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::setFileName(const std::string &file_name, const TextOrigin source,
                                  const std::string &content) {
  switch (source) {
  case TextOrigin::DISK:
    return file_name;
  case TextOrigin::RAM:
    return (content.size() > 0) ? file_name : std::string("");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void TextFile::linesFromString(const std::string &text_in) {
  int n_br = 0;
  const int n_char = text_in.size();
  for (int i = 0; i < n_char; i++) {
    n_br += (text_in[i] == '\n');
  }
  line_count = n_br + (text_in[text_in.size() - 1] != '\n');
  line_limits.resize(line_count + 1);
  text.resize(n_char - n_br + 1);
  int n_tx = 0;
  n_br = 0;
  for (int i = 0; i < n_char; i++) {
    text[n_tx] = text_in[i];
    if (text_in[i] == '\n') {
      n_br++;
      line_limits[n_br] = n_tx;
    }
    else {
      n_tx++;
    }
  }

  // Ensure that the final line limit is catalogged
  line_limits[line_count] = n_tx;
}

//-------------------------------------------------------------------------------------------------
int TextFile::checkAvailableLength(const int line_number, const int start_pos,
                                   const int string_length) const {
  if (line_number >= line_count || line_number < 0) {
    rtErr("The text file originating in " + orig_file + " has " + std::to_string(line_count) +
          " lines and cannot return text from line " + std::to_string(line_number) + ".",
          "TextFile", "checkAvailableLength");
  }
  const int available_chars = line_limits[line_number + 1] - line_limits[line_number];
  if (start_pos >= available_chars || start_pos + string_length >= available_chars) {
    rtErr("Line " + std::to_string(line_number) + " of text file " + orig_file + " has " +
          std::to_string(available_chars) + " characters (requested starting position " +
          std::to_string(start_pos) + ", length " + std::to_string(string_length) + ").",
          "TextFile", "checkAvailableLength");
  }
  return ((string_length < 0) ? available_chars - start_pos : string_length);
}

} // namespace parse
} // namespace stormm
