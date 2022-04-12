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
  /// \param file_name         File to read from
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
  explicit CoordinateSeries(const std::string &file_name,
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

  /// \brief Import coordinates from a CoordinateFrame or PhaseSpace object.  The original object
  ///        must have the same number of atoms as the CoordinateSeries itself, or else a range of
  ///        atoms within the original coordinate object must be specified that fits the
  ///        CoordinateSeries.  The default behavior is to push the new coordinates to the back of
  ///        the list, but any frame index within the bounds of the current list may also be
  ///        specified.
  ///
  /// Overloaded:
  ///   - Accept all types of single-frame coordinate objects
  ///   - Accept all atoms or a subset of the atoms that fits the current atom count of the
  ///     CoordinateSeries
  /// \param cfr          Coordinates to import.  The CoordinateFrameReader is the most basic
  ///                     object available for importation.  Both CoordinateFrame and PhaseSpace
  ///                     objects can create CoordinateFrameReaders or writers, and the writers
  ///                     can be const-ified into readers.
  /// \param cfw          Coordinates to import
  /// \param cf           Coordinates to import
  /// \param ps           Coordinates to import
  /// \param frame_index  Index of the frame into which the coordinates should be imported.  The
  ///                     default value of -1 adds the new coordinates to the end of the list.
  void importCoordinateSet(const CoordinateFrameReader &cfr, int atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const CoordinateFrameReader &cfr, int frame_index = -1);
  void importCoordinateSet(const CoordinateFrameWriter &cfw, int atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const CoordinateFrameWriter &cfw, int frame_index = -1);
  void importCoordinateSet(const CoordinateFrame &cf, int atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const CoordinateFrame &cf, int frame_index = -1);
  void importCoordinateSet(const CoordinateFrame *cf, int atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const CoordinateFrame *cf, int frame_index = -1);
  void importCoordinateSet(const PhaseSpace &ps, int atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const PhaseSpace &ps, int frame_index = -1);
  void importCoordinateSet(const PhaseSpace *ps, nt atom_start, int atom_end,
                           int frame_index = -1);
  void importCoordinateSet(const PhaseSpace *ps, int frame_index = -1);
  
  /// \brief 
  
private:
  int atom_count;                       ///< Number of atoms in each frame.  This number is also
                                        ///<   modifiable with a special call to the resize()
                                        ///<   member function, and between frames the space for
                                        ///<   atoms is padded by the warp size, but there is no
                                        ///<   concept of an atom capacity in the same way that
                                        ///<   there is a frame capacity.
  int frame_count;                      ///< Total number of frames currently in the object
  int frame_capacity;                   ///< Total frame capacity of the object
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
