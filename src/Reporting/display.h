// -*-c++-*-
#ifndef STORMM_DISPLAY_H
#define STORMM_DISPLAY_H

#include <fstream>
#include <iostream>
#include <string>

namespace stormm {
namespace display {
  
/// \brief Print a horizontal bar of --- across the screen, using a corner of the developer's
///        choice (default '+')
void terminalHorizontalRule(const std::string &left_corner = std::string("+"),
                            const std::string &right_corner = std::string("+"),
                            int width = 0, std::ostream *foutp = &std::cout);

} // namespace display
} // namespace stormm

#endif
