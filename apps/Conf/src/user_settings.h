// -*-c++-*-
#ifndef OMNI_USER_SETTINGS_H
#define OMNI_USER_SETTINGS_H

#include "../../../src/Trajectory/trajectory_enumerators.h"

namespace conf_app {
namespace user_input {

using omni::trajectory::CoordinateFileKind;

/// \brief Object to encapsulate a system, a coupled set of coordinates and a single topology.
///        This can be used to enforce conf.omni to read certain coordinates in the context of a
///        particular topology, even when other viable topologies might be available.
struct MoleculeSystem {

  /// \brief The constructor can make a blank system or automatically fill values.
  ///
  /// Overloaded:
  ///   - Create a blank object
  ///   - Create an object with a complete list of files and read settings
  /// \{
  MoleculeSystem();
  MoleculeSystem(const std::string &topology_file_in, const std::string &coordinate_file_in,
                 int frame_start_in, int frame_end_in, CoordinateFileKind coordinate_kind_in);
  /// \}

private:
  std::string topology_file_name;
  std::string coordinate_file_name;
  int frame_start;
  int frame_end;
  CoordinateFileKind coordinate_kind;
}

/// \brief Object to hold general user input data, including file names or regular expressions for
///        topology and coordinate files, energy minimization settings, analysis protocols, and
///        output controls.
struct UserSettings {

private:
  std::vector<std::string> topology_file_names;
  std::vector<std::string> coordinate_file_names;
  std::vector<std::string>
};

} // namespace user_input
} // namespace conf_app

#endif
