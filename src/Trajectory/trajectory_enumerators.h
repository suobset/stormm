// -*-c++-*-
#ifndef OMNI_TRAJECTORY_ENUMERATORS_H
#define OMNI_TRAJECTORY_ENUMERATORS_H

#include <string>

namespace omni {
namespace trajectory {

/// \brief Options for the type of coordinate file to write
enum class CoordinateFileKind {
  AMBER_CRD,         ///< The awful .crd format trajectory file, 10 x %8.3f per line, with or
                     ///<   without box coordinates
  AMBER_INPCRD,      ///< The still-used input coordinates format, as output by tleap, 6 x %12.7f
                     ///<   per line, with or without box coordinates
  AMBER_ASCII_RST,   ///< The deprecated ascii formatted Amber restart file format, 6 x %12.7f
                     ///<   per line, with velocities and with or without box coordinates (a
                     ///<   restart file with no velocities is an INPUT_COORDINATES file)
  AMBER_NETCDF,      ///< The binary Amber NetCDF trajectory format
  AMBER_NETCDF_RST,  ///< The binary Amber NetCDF restart format
  UNKNOWN            ///< The coordinate file kind is not (yet) understood
};

/// \brief An enumerator to track different Amber NetCDF variable identifiers.  The associated
///        translation function produces the strings of interest.
enum class AncdfVariable {
  NCFRAME,         ///< 
  NCSPATIAL,       ///<
  NCATOM,          ///<
  NCCELL_SPATIAL,  ///<
  NCCELL_LENGTHS,  ///<
  NCCELL_ANGULAR,  ///<
  NCCELL_ANGLES,   ///<
  NCCOORDS,        ///<
  NCVELO,          ///<
  NCTEMPERATURE,   ///<
  NCTIME,          ///<
  NCLABEL          ///<
};

/// \brief Enumerate the kind of information that a trajectory can contain
enum class TrajectoryKind {
  POSITIONS,   // Positions of particles
  VELOCITIES,  // Velocities of particles
  FORCES       // Total forces on all particles (the two force arrays are summed in cases where
               //   the forces are partitioned)
};
  
/// \brief Translate the CoordinateFileKind enumeration.
///
/// \param cfkind  The enumerator instance of interest
std::string getCoordinateFileKindName(CoordinateFileKind cfkind);

/// \brief Detect various coordinate file types.
///
/// \param file_name  Name of the file to test
/// \param caller     Name of the calling function (optional)
CoordinateFileKind detectCoordinateFileKind(const std::string &file_name,
                                            const std::string &caller = std::string(""));
  
/// \brief Translate the AncdfVariable enumeration.  This provides keywords to serve as landmarks
///        in an Amber binary trajectory or restart file.
///
/// \param key  The enumerator instance of interest
std::string getCoordinateFileKindName(AncdfVariable key);

} // namespace trajectory
} // namespace omni

#endif
