// -*-c++-*-
#ifndef STORMM_WRITE_SDF_H
#define STORMM_WRITE_SDF_H

#include <vector>
#include "copyright.h"
#include "coordinateframe.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {
    
/// \brief Write one or more frames to a trajectory file in a format-agnostic manner.  This
///        function accepts text strings, with all of the information for the various lines of the
///        file concatenated together with the requisite carriage returns.
///
/// Overloaded:
///   - Accept a single TextFile object, assumed to be one frame
///   - Accept a vector of TextFile objects, each assumed to be one frame
///   - Accept one or more TextFile objects with a corresponding number of CoordinateFrame or
///     PhaseSpace objects, or an appropriately sized PhaseSpaceSynthesis object and vectors of
///     line and column indices at which the coordinates exist.  In this case, the coordinates of
///     the TextFile objects will be replaced with those in the coordinate objects before writing
///     and the TextFile objects themselves will be unaltered.
/// \{
void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const PhaseSpace &ps, int total_width, int fraction,
                const std::vector<int> &crd_lines, int x_start, int y_start, int z_start);

void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf,
                const CoordinateFrame &cf, int total_width, int fraction,
                const std::vector<std::vector<int>> &crd_lines, int x_start, int y_start,
		int z_start);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
