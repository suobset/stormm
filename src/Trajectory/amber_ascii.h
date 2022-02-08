// -*-c++-*-
#ifndef OMNI_AMBER_ASCII_H
#define OMNI_AMBER_ASCII_H

#include "Constants/behavior.h"
#include "Accelerator/hybrid.h"
#include "Parsing/textfile.h"

namespace omni {
namespace trajectory {

using constants::ExceptionResponse;
using card::Hybrid;
using parse::TextFile;

/// \brief Obtain just the atom count from an Amber ASCII-format restart (or inpcrd) file.
///
/// \param tf  The text file to read from (must already be parsed and in RAM, as other routines
///            will go back and read other parts of it--reading the file twice from disk would
///            be wasteful)
int getAmberRestartAtomCount(const TextFile &tf);

/// \brief Obtain the atom count and time from an Amber ASCII-format restart file.
///
/// \param tf          The text file to read from (must already be parsed and in RAM, as other
///                    routines will go back and read other parts of it--reading the file twice
///                    from disk would be wasteful)
/// \param atom_count  The atom count (returned)
/// \param start_time  The initial time (returned)
void getAmberRestartAtomCountAndTime(const TextFile &tf, int *atom_count, double *start_time,
                                     ExceptionResponse policy = ExceptionResponse::WARN);

/// \brief Read coordinates from an Amber inpcrd or ASCII-format restart file.
///
/// \param tf                   The text file to read from
/// \param x_coordinates        Cartesian X coordinates of all particles (returned)
/// \param y_coordinates        Cartesian Y coordinates of all particles (returned)
/// \param z_coordinates        Cartesian Z coordinates of all particles (returned)
/// \param box_space_transform  Transformation to take coordinates into fractional (box) space
/// \param inverse_transform    Transformation to take fractional coordinates into real space
/// \param box_dimensions       Box dimensions (taken directly from the file)
void getAmberInputCoordinates(const TextFile &tf, Hybrid<double> *x_coordinates,
                              Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                              Hybrid<double> *box_space_transform,
                              Hybrid<double> *inverse_transform, Hybrid<double> *box_dimensions);

/// \brief Read velocities from an Amber ASCII-format restart file.
///
/// \param tf            The text file to read from
/// \param x_velocities  Cartesian X velocities of all particles (returned)
/// \param y_velocities  Cartesian Y velocities of all particles (returned)
/// \param z_velocities  Cartesian Z velocities of all particles (returned)
void getAmberRestartVelocities(const TextFile &tf, Hybrid<double> *x_velocities,
                               Hybrid<double> *y_velocities,  Hybrid<double> *z_velocities);

/// \brief Read coordinates from an Amber ASCII-format .crd trajectory.
///
/// \param tf             The text file to read from
/// \param x_coordinates  Cartesian X coordinates of all particles (returned)
/// \param y_coordinates  Cartesian Y coordinates of all particles (returned)
/// \param z_coordinates  Cartesian Z coordinates of all particles (returned)
/// \param box_lengths    Box lengths (excluding angles).  Left unchanged if there is no box
///                       information.
/// \param frame_number   The frame number to seek in the file.  Raises an exception if the end of
///                       the file is reached before finding this frame.
void readAmberCrdFormat(const TextFile &tf, Hybrid<double> *x_coordinates,
                        Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                        Hybrid<double> *box_space_transform, Hybrid<double> *inverse_transform,
                        int frame_number = 0);

} // namespace trajectory
} // namespace omni

#endif
