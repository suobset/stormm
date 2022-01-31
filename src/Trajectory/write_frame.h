// -*-c++-*-
#ifndef OMNI_WRITE_FRAME_H
#define OMNI_WRITE_FRAME_H

#include <string>
#include <vector>
#include "FileManagement/file_util.h"
#include "trajectory_enumerators.h"

namespace omni {
namespace trajectory {

using diskutil::PrintSituation;

/// \brief Write a coordinate trajectory file or input coordinates file.  Each of the overloaded
///        versions of the function feeds into the base case, which uses pointers to the data
///        at hand, with a verified dimension.
///
/// Overloaded:
///   - Takes double pointers and the number of atoms
///   - Takes std::vector<double> objects
///   - Takes a Coordinates object based on Hybrid objects
///
/// \param filename        Name of the file to write
/// \param file_kind       The type of coordinate (or restart) file to write
/// \param expectation     Dictates writing behavior based on the presence or absence of any
///                        existing file with the same name
/// \param title           Title of the coordinate file, to place on the first line
/// \param natom           The number of atoms in the system
/// \param x_crd           Vector / array of Cartesian x coordinates for all atoms
/// \param y_crd           Vector / array of Cartesian y coordinates for all atoms
/// \param z_crd           Vector / array of Cartesian z coordinates for all atoms
/// \param box_dimensions  Six-element vector of box dimensions (can be obtained from a
///                        transformation matrix)
/// \param time_point      The time point to write at the top of a restart file
/// \{
void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                int natom, const double* x_crd, const double* y_crd, const double* z_crd,
                const double* x_vel, const double* y_vel, const double* z_vel,
                const double* box_dimensions);

void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                const std::vector<double> &x_crd, const std::vector<double> &y_crd,
                const std::vector<double> &z_crd, const std::vector<double> &x_vel,
                const std::vector<double> &y_vel, const std::vector<double> &z_vel,
                const std::vector<double> &box_dimensions);
/// \}

} // namespace trajectory
} // namespace omni

#endif
