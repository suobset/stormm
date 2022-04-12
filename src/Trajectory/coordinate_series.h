// -*-c++-*-
#ifndef OMNI_COORDINATE_SERIES_H
#define OMNI_COORDINATE_SERIES_H

#include "Constants/scaling.h"
#include "trajectory_enumerators.h"
#include "coordinateframe.h"
#include "phasespace.h"

namespace omni {
namespace trajectory {

/// \brief Store the coordinates and box information for a series of frames, in one of several
///        levels of precision.  Individual frames can be extracted into CoordinateFrame, and
///        PhaseSpace, and PhaseSpaceSynthesis objects, or new frames can be added from
///        CoordinateFrame objects.  This object is not the CoordinateFrame equivalent of a
///        PhaseSpaceSynthesis object, however: its purpose is to collect only the coordinates of
///        many frames of a single system at the appropriate level of precision.
template <typename T>  class CoordinateSeries {
public:

  /// \brief There are several options for constructing this collection of coordinate frames.
  ///
  /// Overloaded:
  ///   - Allocate to hold a given number of atoms and frames
  ///   - Create from any of the coordinate file formats (restart and input coordinate formats
  ///     will create only one frame, but this limit can be increased later)
  ///   - From an existing PhaseSpace or CoordinateFrame object (as a pointer or a copy, with the
  ///     option to make many copies immediately)
  ///
  /// \param natom_in          The number of atoms expected
  /// \param nframe_in         Initial number of frames to allocate for
  /// \param unit_cell_in      The type of unit cell to prepare for (this can be modified after
  ///                          creating the object)
  /// \param file_name_in      File to read from
  /// \param frame_numbers_in  Frame numbers of the file to read (default all frames)
  /// \param replica_count_in  The number of times to replicate a series of one of more frames
  ///                          read from a file (this is useful for immediately making many copies
  ///                          of a restrart or input coordinates file, which has only one frame)
  /// \param atom_count_in     The number of atoms to expect (critical if the input file is an
  ///                          Amber .crd format trajectory, otherwise can be left at zero to be
  ///                          filled in when the file is read)
  /// \param ps                Pre-existing object with a complete coordinate set to use as a
  ///                          template.  When constructing from a pre-existing PhaseSpace object,
  ///                          nframe_in indicates a number of copies to allocate for and create.
  /// \param cf                Pre-existing object with a complete coordinate set to use as a
  ///                          template.  When constructing from a pre-existing CoordinateFrame,
  ///                          nframe_in indicates a number of copies to allocate for and create.
  /// \{
  explicit CoordinateSeries(int natom_in = 0, int nframe_in = 0,
                            UnitCellType unit_cell_in = UnitCellType::NONE);
  explicit CoordinateSeries(const std::string &file_name_in,
                            CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                            const std::vector<int> &frame_numbers_in = {},
                            int replica_count_in = 1, int atom_count_in = 0);
  explicit CoordinateSeries(const TextFile &tf, CoordinateFileKind file_kind,
                            const std::vector<int> &frame_numbers_in = {},
                            int replica_count_in = 1, int atom_count_in = 0);
  explicit CoordinateSeries(PhaseSpace *ps, int nframe_in);
  explicit CoordinateSeries(const PhaseSpace &ps, int nframe_in);
  explicit CoordinateSeries(CoordinateFrame *cf, int nframe_in);
  explicit CoordinateSeries(const CoordinateFrame &cf, int nframe_in);
  /// \}

private:
  std::string file_name;                ///< File from which the series was derived, if applicable
  int atom_count;                       ///< Number of atoms in each frame
  int frame_count;                      ///< Total number of frames currently allocated
  Hybrid<int> frame_numbers;            ///< Frame numbers of the different frames derived from
                                        ///<   some trajectory, if applicable
  Hybrid<T> x_coordinates;              ///< Cartesian X coordinates of all particles, with each
                                        ///<   frame's coordinates padded by the warp size
  Hybrid<T> y_coordinates;              ///< Cartesian Y coordinates of all particles 
  Hybrid<T> z_coordinates;              ///< Cartesian Z coordinates of all particles
  Hybrid<double> box_space_transforms;  ///< Matrices to transform each frame into box space
  Hybrid<double> inverse_transforms;    ///< Matrix to transform each frame into real space
  Hybrid<double> box_dimensions;        ///< Lengths and angles defining each frame's unit cell
                                        ///<   (lengths are given in Angstroms, angles in radians)
  Hybrid<T> particle_data;              ///< Array targeted by the x_, y_, and z_coordinates
                                        ///<   POINTER-kind Hybrid objects
  Hybrid<double> box_data;              ///< Array targeted by box transform and dimension
                                        ///<   POINTER-kind Hybrid objects
};
 
} // namespace trajectory 
} // namespace omni

#endif
