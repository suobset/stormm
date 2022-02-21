#include <cstdlib>
#include <string>
#include <vector>
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Namelists/namelist_emulator.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/polynumeric.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Potential/valence_potential.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "user_settings.h"

namespace conf_app {
namespace user_input {

using omni::diskutil::DrivePathType;
using omni::diskutil::getBaseName;
using omni::diskutil::getDrivePathType;
using omni::diskutil::splitPath;
using omni::energy::evaluateBondTerms;
using omni::energy::evaluateAngleTerms;
using omni::errors::rtErr;
using omni::errors::rtWarn;
using omni::namelist::MoleculeSystem;
using omni::parse::NumberFormat;
using omni::parse::TextOrigin;
using omni::parse::verifyNumberFormat;
using omni::energy::ScoreCard;
using omni::topology::ValenceKit;
using omni::trajectory::CoordinateFrameReader;
using omni::trajectory::detectCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const int argc, const char* argv[]) :
    policy{ExceptionResponse::DIE},
    input_file{std::string(default_conformer_input_file)},
    common_core_mask{std::string("")},
    file_io_input{}, line_min_input{}, solvent_input{}, prng_input{},
    topology_cache{},
    initial_coordinates_cache{},
    topology_indices{}
{
  // Local variables to store command line arguments
  int cval_igseed = 0;
  std::string cval_report_file;
  std::string cval_conf_file_name;
  std::vector<std::string> cval_topology_file_names;
  std::vector<std::string> cval_coordinate_file_names;
  
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile  = false;
  bool cli_igseed   = false;
  bool cli_report   = false;
  bool cli_confname = false;
  for (int i = 1; i < argc; i++) {
    if (i < argc - 1 && strcmp(argv[i], "-i") == 0) {
      input_file = std::string(argv[i + 1]);
      cli_inpfile = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-igseed") == 0) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER) == false) {
        rtErr("The random seed " + std::string(argv[i + 1]) + " is not valid.", "UserSettings");
      }
      cval_igseed = atoi(argv[i + 1]);
      cli_igseed = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-p") == 0) {
      std::string tmp_top(argv[i + 1]);
      cval_topology_file_names.push_back(tmp_top);
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-c") == 0) {
      std::string tmp_crd(argv[i + 1]);
      cval_coordinate_file_names.push_back(tmp_crd);
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-o") == 0) {
      cval_report_file = std::string(argv[i + 1]);
      cli_report = true;
      i++;
    }
    else if (i < argc - 1 && strcmp(argv[i], "-xname") == 0) {
      cval_conf_file_name = std::string(argv[i + 1]);
      cli_confname = true;
      i++;
    }
    else if (strcmp(argv[i], "-warn") == 0) {
      policy = ExceptionResponse::WARN;
    }
    else if (strcmp(argv[i], "-silent") == 0) {
      policy = ExceptionResponse::SILENT;
    }
    else {
      rtErr("Command line argument " + std::string(argv[i]) + " was not recognized.",
            "UserSettings");
    }
  }

  // Process the input file.  Take only the first instance of each namelist, as found by searching
  // from the beginning.
  if (getDrivePathType(input_file) != DrivePathType::FILE) {
    const std::string descriptor = (cli_inpfile) ? std::string("user specified") :
                                                   std::string("default");
    rtErr("The " + descriptor + " input file " + input_file + " was not found or could not be "
          "read.", "UserSettings");
  }
  TextFile inp_tf(input_file, TextOrigin::DISK, "Input deck for conformer.omni", "UserSettings");
  int start_line = 0;
  file_io_input = FilesControls(inp_tf, &start_line, policy);
  start_line = 0;
  line_min_input = MinimizeControls(inp_tf, &start_line, policy);
  start_line = 0;
  solvent_input = SolventControls(inp_tf, &start_line, policy);
  start_line = 0;
  prng_input = RandomControls(inp_tf, &start_line, policy);

  // Superimpose, or contribute, command line directives
  if (cli_igseed) {
    prng_input.setRandomSeed(cval_igseed);    
  }
  if (cli_report) {
    file_io_input.setReportFileName(cval_report_file);
  }
  if (cli_confname) {
    file_io_input.setGeneralTrajectoryFileName(cval_conf_file_name);
  }
  if (cval_topology_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_topology_file_names.size(); i++) {
      file_io_input.addFreeTopologyName(cval_topology_file_names[i]);
    }
  }
  if (cval_coordinate_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_coordinate_file_names.size(); i++) {
      file_io_input.addFreeCoordinateName(cval_coordinate_file_names[i]);
    }
  }
  
  // Read all free topologies and free coordinate sets, then determine which free topology
  // matches which free coordinate set.
  int n_free_top = file_io_input.getFreeTopologyCount();
  topology_cache.reserve(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    try {
      topology_cache.push_back(AtomGraph(file_io_input.getFreeTopologyName(i)));
    }
    catch (std::runtime_error) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The format of topology " + file_io_input.getFreeTopologyName(i) + " could not be "
              "understood.", "UserSettings");
      case ExceptionResponse::WARN:
        rtWarn("The format of topology " + file_io_input.getFreeTopologyName(i) + " could not be "
               "understood.  The file will be skipped.", "UserSettings");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  n_free_top = topology_cache.size();
  int n_free_crd = file_io_input.getFreeCoordinatesCount();
  initial_coordinates_cache.reserve(n_free_crd);
  std::vector<int> initial_coordinates_frame_count;
  initial_coordinates_frame_count.reserve(n_free_crd);
  std::vector<CoordinateFileKind> initial_coordinates_kind;
  initial_coordinates_kind.reserve(n_free_crd);
  for (int i = 0; i < n_free_crd; i++) {
    const std::string crd_name = file_io_input.getFreeCoordinateName(i);
    const CoordinateFileKind kind = detectCoordinateFileKind(crd_name);
    switch (kind) {
    case CoordinateFileKind::AMBER_CRD:
      rtErr("Amber .crd format trajectories cannot be read without a means of providing the "
            "atom count.  Use the -sys keyword to couple such trajectories to a specific "
            "topology file so that multiple frames can be read from such a format.",
            "UserSettings");
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::AMBER_ASCII_RST:
      initial_coordinates_cache.push_back(CoordinateFrame(crd_name, kind));
      initial_coordinates_frame_count.push_back(1);
      initial_coordinates_kind.push_back(kind);
      break;
    case CoordinateFileKind::AMBER_NETCDF:
    case CoordinateFileKind::AMBER_NETCDF_RST:
      break;
    case CoordinateFileKind::UNKNOWN:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The format of " + crd_name + " could not be understood.", "UserSettings");
      case ExceptionResponse::WARN:
        rtWarn("The format of " + crd_name + " could not be understood.  The file will be skipped",
               "UserSettings");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  n_free_crd = initial_coordinates_cache.size();

  // CHECK
  printf("There were %d free topologies and %d free coordinate sets.\n", n_free_top, n_free_crd);
  // END CHECK
  
  // Filter the unique topologies and match them to systems.  List the atom counts of each.
  int max_match = 0;
  std::vector<int> topology_atom_counts(n_free_top);
  std::vector<int> coordinate_atom_counts(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    topology_atom_counts[i] = topology_cache[i].getAtomCount();
  }
  for (int i = 0; i < n_free_crd; i++) {
    coordinate_atom_counts[i] = initial_coordinates_cache[i].getAtomCount();
  }

  // Make a table of the unique topology atom counts.  Track which coordinate sets might
  // be tied to each group of topologies.
  std::vector<bool> topology_covered(n_free_top, false);
  std::vector<bool> coordinates_covered(n_free_crd, false);
  std::vector<int> topology_series(n_free_top);
  std::vector<int> coordinate_series(n_free_top);
  std::vector<int> unique_topology_sizes;
  std::vector<int> unique_topology_size_bounds(1, 0);
  std::vector<int> unique_coordinate_size_bounds(1, 0);
  int top_series_counter = 0;
  int crd_series_counter = 0;
  for (int i = 0; i < n_free_top; i++) {
    if (topology_covered[i]) {
      continue;
    }

    // First, scan to see if any other topologies share the same number of atoms.
    topology_covered[i] = true;
    int n_samesize_topology = 1;
    topology_series[top_series_counter] = i;
    top_series_counter++;
    const int iatom_count = topology_atom_counts[i];
    for (int j = i + 1; j < n_free_top; j++) {
      if (topology_atom_counts[j] == iatom_count) {
        topology_covered[j] = true;
        topology_series[top_series_counter] = j;
        top_series_counter++;
      }
    }
    for (int j = 0; j < n_free_crd; j++) {
      if (coordinates_covered[j] == false && coordinate_atom_counts[j] == iatom_count) {
        coordinates_covered[j] = true;
        coordinate_series[crd_series_counter] = j;
        crd_series_counter++;
      }
    }
    unique_topology_sizes.push_back(iatom_count);
    unique_topology_size_bounds.push_back(top_series_counter);
    unique_coordinate_size_bounds.push_back(crd_series_counter);
  }

  // Check that all coordinates were covered
  std::vector<std::string> orphan_coordinates;
  for (int i = 0; i < n_free_crd; i++) {
    if (coordinates_covered[i] == false) {
      orphan_coordinates.push_back(initial_coordinates_cache[i].getFileName());
    }
  }
  const int n_orphan = orphan_coordinates.size();
  if (n_orphan > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      {
        std::string orphan_errmsg;
        if (n_orphan > 6) {
          for (int i = 0; i < 3; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
          orphan_errmsg += "  (... more files, list truncated ...)\n";
          for (int i = n_orphan - 3; i < n_orphan; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
        }
        else {
          for (int i = 0; i < 6; i++) {
            orphan_errmsg += "  " + orphan_coordinates[i] + '\n';
          }
        }
        rtWarn("A total of " + std::to_string(orphan_coordinates.size()) + " free coordinate sets "
               "do not correspond to any of the provided topologies, due to atom count "
               "mismatches.\n\nOrphan topologies:\n" + orphan_errmsg, "UserSettings");
      }
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  
  // Test each coordinate set with respect to topologies of the appropriate size.  Evaluate bond
  // and angle energy, and take the best scoring result as the indicator of which topology
  // describes which coordinate set.  Make new systems based on each coordinate set.
  const int n_unique_sizes = unique_topology_sizes.size();
  ScoreCard sc(1);
  std::string trajectory_base, trajectory_ext, restart_base, restart_ext;
  splitPath(file_io_input.getTrajectoryFileName(), &trajectory_base, &trajectory_ext);
  splitPath(file_io_input.getCheckpointFileName(), &restart_base, &restart_ext);
  for (int i = 0; i < n_unique_sizes; i++) {

    // Loop over all coordinates in this size group.  Try interpreting them with each topology.
    // Store each result as a unique MoleculeSystem and expand the list.
    const int j_llim = unique_coordinate_size_bounds[i];
    const int j_hlim = unique_coordinate_size_bounds[i + 1];
    const int k_llim = unique_topology_size_bounds[i];
    const int k_hlim = unique_topology_size_bounds[i + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      const int icrdj = coordinate_series[j];
      const CoordinateFrameReader cfr(initial_coordinates_cache[icrdj].data());
      int best_topology;
      double min_bondang_e, min_bond_e, min_angl_e;
      for (int k = k_llim; k < k_hlim; k++) {
        const int tpk = topology_series[k];
        const ValenceKit<double> vk = topology_cache[tpk].getDoublePrecisionValenceKit();
        const double bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
        const double angl_e = evaluateAngleTerms(vk, cfr, &sc, 0);
        const double bondang_e = bond_e + angl_e;
        if (k == k_llim || bondang_e < min_bondang_e) {
          best_topology = tpk;
          min_bondang_e = bondang_e;
          min_bond_e = bond_e;
          min_angl_e = angl_e;
        }
      }

      // Construct the trajectory and restart file names for this system based on generic paths.
      // Hack a solution in the odd event that the user has stuffed their file differentiation
      // behind the final dot.
      const std::string orig_crd_file = initial_coordinates_cache[icrdj].getFileName();
      std::string orig_base, orig_ext;
      splitPath(getBaseName(orig_crd_file), &orig_base, &orig_ext);
      const std::string trajectory_middle = (orig_base.size() > 0) ? orig_base : orig_ext;

      // Add this pair to the list of systems
      file_io_input.addSystem(MoleculeSystem(topology_cache[best_topology].getFileName(),
                                             initial_coordinates_cache[icrdj].getFileName(),
                                             trajectory_base + trajectory_middle + trajectory_ext,
                                             restart_base + trajectory_middle + restart_ext, 0,
                                             initial_coordinates_frame_count[icrdj], 1,
                                             initial_coordinates_kind[icrdj],
                                             file_io_input.getOutputCoordinateFormat(),
                                             file_io_input.getCheckpointFormat()));
    }
  }

  // 
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator conformerInput(const TextFile &tf, int *start_line,
                                const ExceptionResponse policy) {
}

} // namespace user_input
} // namespace conf_app
