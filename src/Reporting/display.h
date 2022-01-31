// -*-c++-*-
#ifndef OMNI_DISPLAY_H
#define OMNI_DISPLAY_H

#include <string>

namespace omni {
namespace display {

/// \brief Print a horizontal bar of --- across the screen, using a corner of the developer's
///        choice (default '+')
void terminalHorizontalRule(const std::string &left_corner = std::string("+"),
                            const std::string &right_corner = std::string("+"),
                            int width = 0);

} // namespace display
} // namespace omni

#endif
