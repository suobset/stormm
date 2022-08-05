#include <sys/ioctl.h>
#include <unistd.h>
#include "display.h"

namespace stormm {
namespace display {

//-------------------------------------------------------------------------------------------------
void terminalHorizontalRule(const std::string &left_corner, const std::string &right_corner,
                            const int width) {

  // Obtain the console size
  struct winsize console_dims;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);
  int n_dash = (width == 0) ? console_dims.ws_col - 1 : width;
  n_dash -= left_corner.size() + right_corner.size();
  std::string hrule(n_dash, '-');
  hrule = left_corner + hrule + right_corner;
  printf("%s\n", hrule.c_str());
}

} // namespace display
} // namespace stormm
