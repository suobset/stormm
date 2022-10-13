// -*-c++-*-
#ifndef STORMM_WRITE_SDF_H
#define STORMM_WRITE_SDF_H

#include <vector>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinateframe.h"
#include "coordinate_swap_plan.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using diskutil::PrintSituation;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
  
/// \brief Write one or more frames to a trajectory file in a format-agnostic manner.  This
///        function accepts text strings, with all of the information for the various lines of the
///        file concatenated together with the requisite carriage returns.
///
/// Overloaded:
///   - Accept a single TextFile object, assumed to be one frame
///   - Accept a vector of TextFile objects, each assumed to be one frame
///   - Accept one or more TextFile objects with a corresponding number of CoordinateFrame or
///     PhaseSpace objects, an appropriately sized PhaseSpaceSynthesis object, or abstracts
///     thereof, and vectors of line and column indices at which the coordinates exist.  In this
///     case, the coordinates of the TextFile objects will be replaced with those in the coordinate
///     objects before writing and the TextFile objects themselves will be unaltered.
/// \{
void writeFrame(std::ofstream *foutp, const TextFile &tf);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf);

void writeFrame(std::ofstream *foutp, const TextFile &tf, const double* xcrd, const double* ycrd,
                const double* zcrd, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const double* xcrd, const double* ycrd, const double* zcrd,
                const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const PhaseSpaceReader &psr, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateFrameReader &cfr, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateFrame &cf, const CoordinateSwapPlan &excision);

void writeFrame(std::ofstream *foutp, const std::vector<TextFile> &tf_list,
                const PsSynthesisReader &poly_psr,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices);

void writeFrame(const std::string &filename, PrintSituation expectation,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
