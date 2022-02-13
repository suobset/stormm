#include <iostream>
#include <fstream>
#include "textfile.h"
#include "Reporting/error_format.h"

namespace omni {
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
char TextFile::getText(const int index) const {
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
  line_count = n_br;
  line_limits.resize(n_br + 1);
  text.resize(n_char - n_br);
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
  line_limits[line_count] = n_char - n_br;
}

} // namespace parse
} // namespace omni
