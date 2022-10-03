#include "copyright.h"
#include "systemcache.h"
#include "FileManagement/file_listing.h"
#include "Math/summation.h"
#include "Parsing/parse.h"
#include "Potential/scorecard.h"
#include "Potential/valence_potential.h"
#include "MoleculeFormat/mdl_mol_format.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"

namespace stormm {
namespace synthesis {

using diskutil::getBaseName;
using diskutil::splitPath;
using energy::evaluateBondTerms;
using energy::evaluateAngleTerms;
using energy::ScoreCard;
using math::prefixSumInPlace;
using math::PrefixSumType;
using namelist::MoleculeSystem;
using parse::findStringInVector;
using structure::MdlMolObj;
using structure::readStructureDataFile;
using topology::UnitCellType;
using topology::ValenceKit;
using trajectory::detectCoordinateFileKind;
  
//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache() :
    system_count{0}, topology_cache{}, coordinates_cache{}, features_cache{}, restraints_cache{},
    static_masks_cache{}, forward_masks_cache{}, topology_indices{}, example_indices{},
    topology_cases{}, topology_case_bounds{}, system_trajectory_names{}, system_checkpoint_names{},
    system_labels{}, system_trajectory_kinds{}, system_checkpoint_kinds{}
{}

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const FilesControls &fcon, const std::vector<RestraintControls> &rstcon,
                         const ExceptionResponse policy,
                         const MapRotatableGroups map_chemfe_rotators, StopWatch *timer_in) :
    SystemCache()
{
  // Read all free topologies, using a try-catch block to filter out things that may not work.
  int n_free_top = fcon.getFreeTopologyCount();
  topology_cache.reserve(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    try {
      topology_cache.push_back(AtomGraph(fcon.getFreeTopologyName(i)));
    }
    catch (std::runtime_error) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The format of topology " + fcon.getFreeTopologyName(i) +
              " could not be understood.", "UserSettings");
      case ExceptionResponse::WARN:
        rtWarn("The format of topology " + fcon.getFreeTopologyName(i) +
               " could not be understood.  The file will be skipped.", "UserSettings");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  n_free_top = topology_cache.size();
  
  // Read all free coordinate sets, using file type detection to filter the list
  int n_free_crd = fcon.getFreeCoordinatesCount();
  std::vector<CoordinateFrame> tmp_coordinates_cache;
  tmp_coordinates_cache.reserve(n_free_crd);
  std::vector<int> tmp_coordinates_frame_count;
  tmp_coordinates_frame_count.reserve(n_free_crd);
  std::vector<CoordinateFileKind> tmp_coordinates_kind;
  tmp_coordinates_kind.reserve(n_free_crd);
  for (int i = 0; i < n_free_crd; i++) {
    const std::string crd_name = fcon.getFreeCoordinateName(i);
    const CoordinateFileKind kind = detectCoordinateFileKind(crd_name);
    switch (kind) {
    case CoordinateFileKind::AMBER_CRD:
      rtErr("Amber .crd format trajectories cannot be read without a means of providing the "
            "atom count.  Use the -sys keyword to couple such trajectories to a specific "
            "topology file so that multiple frames can be read from such a format.",
            "UserSettings");
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::AMBER_ASCII_RST:
      tmp_coordinates_cache.push_back(CoordinateFrame(crd_name, kind));
      tmp_coordinates_frame_count.push_back(1);
      tmp_coordinates_kind.push_back(kind);
      break;
    case CoordinateFileKind::AMBER_NETCDF:
    case CoordinateFileKind::AMBER_NETCDF_RST:
      break;
    case CoordinateFileKind::SDF:
      {
        const std::vector<MdlMolObj> all_frames = readStructureDataFile(crd_name);
        const size_t n_entries = all_frames.size();
        for (size_t i = 0; i < n_entries; i++) {
          tmp_coordinates_cache.push_back(all_frames[i].exportCoordinateFrame());
          tmp_coordinates_frame_count.push_back(1);
          tmp_coordinates_kind.push_back(kind);
        }
      }
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
  n_free_crd = tmp_coordinates_cache.size();

  // Filter the unique topologies and match them to systems.  List the atom counts of each.
  int max_match = 0;
  std::vector<int> topology_atom_counts(n_free_top);
  std::vector<int> coordinate_atom_counts(n_free_crd);
  for (int i = 0; i < n_free_top; i++) {
    topology_atom_counts[i] = topology_cache[i].getAtomCount();
  }
  for (int i = 0; i < n_free_crd; i++) {
    coordinate_atom_counts[i] = tmp_coordinates_cache[i].getAtomCount();
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
      orphan_coordinates.push_back(tmp_coordinates_cache[i].getFileName());
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

  // Pull the original list of systems out of the &files namelist
  std::vector<MoleculeSystem> sysvec;
  int nsys = fcon.getSystemDefinitionCount();
  sysvec.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    sysvec.push_back(fcon.getSystem(i));
    const std::string traj_str = sysvec[i].getTrajectoryFileName();
    if (traj_str.size() == 0LLU) {
      sysvec[i].setTrajectoryFileName(fcon.getTrajectoryFileName());
    }
    const std::string chkp_str = sysvec[i].getCheckpointFileName();
    if (chkp_str.size() == 0LLU) {
      sysvec[i].setCheckpointFileName(fcon.getCheckpointFileName());
    }
  }
  
  // Test each coordinate set with respect to topologies of the appropriate size.  Evaluate
  // valence term energy, and take the best scoring result as the indicator of which topology
  // describes which coordinate set.  Make new systems based on each coordinate set.
  const int n_unique_sizes = unique_topology_sizes.size();
  ScoreCard sc(1);
  std::string trajectory_base, trajectory_ext, restart_base, restart_ext;
  splitPath(fcon.getTrajectoryFileName(), &trajectory_base, &trajectory_ext);
  splitPath(fcon.getCheckpointFileName(), &restart_base, &restart_ext);
  std::vector<bool> topology_in_use(n_free_top, false);
  int n_paired_systems = 0;
  for (int i = 0; i < n_unique_sizes; i++) {

    // Loop over all coordinates in this size group.  Try interpreting them with each topology.
    // Store each result as a unique MoleculeSystem and expand the list.
    const int j_llim = unique_coordinate_size_bounds[i];
    const int j_hlim = unique_coordinate_size_bounds[i + 1];
    const int k_llim = unique_topology_size_bounds[i];
    const int k_hlim = unique_topology_size_bounds[i + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      const int icrdj = coordinate_series[j];
      const CoordinateFrameReader cfr(tmp_coordinates_cache[icrdj].data());
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
      const std::string orig_crd_file = tmp_coordinates_cache[icrdj].getFileName();
      std::string orig_base, orig_ext;
      splitPath(getBaseName(orig_crd_file), &orig_base, &orig_ext);
      std::string trajectory_middle("_");
      trajectory_middle += (orig_base.size() > 0) ? orig_base + "." : orig_ext + ".";
      
      // Add this pair to the list of systems
      const std::string syslabel = std::string("PairedSystem") + std::to_string(n_paired_systems);
      sysvec.push_back(MoleculeSystem(topology_cache[best_topology].getFileName(),
                                      tmp_coordinates_cache[icrdj].getFileName(),
                                      trajectory_base + trajectory_middle + trajectory_ext,
                                      restart_base + trajectory_middle + restart_ext, syslabel, 0,
                                      tmp_coordinates_frame_count[icrdj], 1,
                                      tmp_coordinates_kind[icrdj],
                                      fcon.getOutputCoordinateFormat(),
                                      fcon.getCheckpointFormat()));
      n_paired_systems++;

      // Note that the topology is used
      topology_in_use[best_topology] = true;
    }
  }

  // Loop back over systems and make sure that each has a unique restart file (even if these files
  // are not to be written).  Systems with similar trajectory file names must all share the same
  // original topology.
  nsys = sysvec.size();
  bool name_collision;
  do {
    name_collision = false;
    std::vector<std::string> unique_restart_names;
    std::vector<int> name_copies;
    unique_restart_names.reserve(nsys);
    name_copies.reserve(nsys);
    int n_unique = 0;
    for (int i = 0; i < nsys; i++) {
      const int pos = findStringInVector(unique_restart_names, sysvec[i].getCheckpointFileName());
      if (pos != n_unique) {
        name_copies[pos] += 1;
        name_collision = true;
      }
      else {
        unique_restart_names.push_back(sysvec[i].getCheckpointFileName());
        name_copies.push_back(1);
        n_unique++;
      }
    }
    if (name_collision) {
      for (int i = 0; i < n_unique; i++) {
        if (name_copies[i] == 1) {
          continue;
        }
        const std::string overused_name = unique_restart_names[i];
        std::string overused_base, overused_ext;
        splitPath(overused_name, &overused_base, &overused_ext);
        const bool base_substantial = (overused_base.size() > 0LLU);
        const bool ext_substantial = (overused_ext.size() > 0LLU);
        const bool both_substantial = (base_substantial && ext_substantial);
        int overused_count = 0;
        for (int j = 0; j < nsys; j++) {
          if (sysvec[j].getCheckpointFileName() == overused_name) {
            if (both_substantial) {
              sysvec[j].setCheckpointFileName(overused_base + "_" +
                                              std::to_string(overused_count) + "." + overused_ext);
            }
            else if (base_substantial) {
              sysvec[j].setCheckpointFileName(overused_base + "_" +
                                              std::to_string(overused_count));
            }
            else {
              sysvec[j].setCheckpointFileName(overused_ext + "_" +
                                              std::to_string(overused_count));
            }
            overused_count++;
          }
        }
      }
    }
  } while (name_collision);
  do {
    name_collision = false;
    std::vector<std::string> unique_trajectory_names;
    std::vector<std::string> unique_trajectory_basis;
    std::vector<int> divergent_topologies;
    unique_trajectory_names.reserve(nsys);
    unique_trajectory_basis.reserve(nsys);
    divergent_topologies.reserve(nsys);
    int n_unique = 0;
    for (int i = 0; i < nsys; i++) {
      const int pos = findStringInVector(unique_trajectory_names,
                                         sysvec[i].getTrajectoryFileName());
      if (pos != n_unique) {

        // Two systems name the same trajectory, but is that a problem?  Not if the topologies
        // describing each system are the same.
        if (unique_trajectory_basis[pos] != sysvec[i].getTopologyFileName()) {
          divergent_topologies[pos] += 1;
          name_collision = true;
        }
      }
      else {
        unique_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
        unique_trajectory_basis.push_back(sysvec[i].getTopologyFileName());
        divergent_topologies.push_back(1);
        n_unique++;
      }
    }
    if (name_collision) {
      for (int i = 0; i < n_unique; i++) {
        if (divergent_topologies[i] == 1) {
          continue;
        }
        const std::string overused_name = unique_trajectory_names[i];
        std::string overused_base, overused_ext;
        splitPath(overused_name, &overused_base, &overused_ext);
        const bool base_substantial = (overused_base.size() > 0LLU);
        const bool ext_substantial = (overused_ext.size() > 0LLU);
        const bool both_substantial = (base_substantial && ext_substantial);
        int overused_count = 0;
        std::vector<std::string> new_trajectory_basis;
        std::vector<std::string> new_trajectory_names;
        for (int j = 0; j < nsys; j++) {
          if (sysvec[j].getTrajectoryFileName() != overused_name) {
            continue;
          }
          const std::string jtop_name = sysvec[j].getTopologyFileName();
          if (jtop_name == unique_trajectory_basis[i]) {
            continue;
          }
          bool new_basis_found = false;
          for (int k = 0; k < overused_count; k++) {
            if (jtop_name == new_trajectory_basis[k]) {
              sysvec[j].setTrajectoryFileName(new_trajectory_names[k]);
              new_basis_found = true;
            }
          }
          if (new_basis_found == false) {
            std::string ntraj_name;
            if (both_substantial) {
              ntraj_name = overused_base + "_" + std::to_string(overused_count) + "." +
                           overused_ext;
            }
            else if (base_substantial) {
              ntraj_name = overused_base + "_" + std::to_string(overused_count);
            }
            else {
              ntraj_name = overused_ext + "_" + std::to_string(overused_count);
            }
            sysvec[j].setTrajectoryFileName(ntraj_name);
            new_trajectory_names.push_back(ntraj_name);
            new_trajectory_basis.push_back(jtop_name);
            overused_count++;
          }
        }
      }      
    }
  } while (name_collision);

  // Loop back over the systems (now representing all entries, including the paired free topologies
  // and coordinate sets).  If the topology has already been read, don't read it again.  Read
  // coordinates (and perhaps velocities, if available) into phase space objects.
  std::vector<std::string> current_topology_holdings;
  current_topology_holdings.reserve(n_free_top);
  int cache_track = 0;
  for (int i = 0; i < n_free_top; i++) {
    if (topology_in_use[i]) {
      current_topology_holdings.push_back(topology_cache[cache_track].getFileName());
      cache_track++;
    }
    else {
      topology_cache.erase(topology_cache.begin() + cache_track);
    }
  }
  for (int i = 0; i < nsys; i++) {
    int top_idx = findStringInVector(current_topology_holdings, sysvec[i].getTopologyFileName());
    bool topology_ok = false;
    if (top_idx >= current_topology_holdings.size()) {
      try {
        topology_cache.push_back(AtomGraph(sysvec[i].getTopologyFileName()));
        current_topology_holdings.push_back(sysvec[i].getTopologyFileName());
        top_idx = topology_cache.size() - 1LLU;
        topology_ok = true;
      }
      catch (std::runtime_error) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The format of topology " + sysvec[i].getTopologyFileName() +
                " for system " + std::to_string(i) + " could not be understood.", "UserSettings");
        case ExceptionResponse::WARN:
          rtWarn("The format of topology " + sysvec[i].getTopologyFileName() +
                 " could not be understood.  System " + std::to_string(i + 1) +
                 " will be skipped.", "UserSettings");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
    }
    else {
      topology_ok = true;
    }
    if (topology_ok) {
      switch (detectCoordinateFileKind(sysvec[i].getInputCoordinateFileName(), "SystemCache")) {
      case CoordinateFileKind::AMBER_CRD:
        break;
      case CoordinateFileKind::AMBER_INPCRD:
      case CoordinateFileKind::AMBER_ASCII_RST:
        coordinates_cache.push_back(PhaseSpace(sysvec[i].getInputCoordinateFileName(),
                                               sysvec[i].getInputCoordinateFileKind()));
        topology_indices.push_back(top_idx);
        system_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
        system_checkpoint_names.push_back(sysvec[i].getCheckpointFileName());
        system_labels.push_back(sysvec[i].getLabel());
        system_count += 1;
        break;
      case CoordinateFileKind::AMBER_NETCDF:
      case CoordinateFileKind::AMBER_NETCDF_RST:
        break;
      case CoordinateFileKind::SDF:
        {
          // If the user has specified the frame range [ 0, -1 ), take that as a special case
          // indicating that "all frames" shall be read and produce no warning or error.
          // Otherwise, check the specified range against the file's actual contents.
          const int fr_init  = sysvec[i].getStartingFrame();
          const int fr_final = sysvec[i].getFinalFrame();
          const std::vector<MdlMolObj> frame_selection = (fr_init == 0 && fr_final == -1) ?
            readStructureDataFile(sysvec[i].getInputCoordinateFileName()) :
            readStructureDataFile(sysvec[i].getInputCoordinateFileName(), fr_init, fr_final);
          const int jlim = frame_selection.size();
          for (int j = 0; j < jlim; j++) {

            // Apply a simple check to help ensure that the data in the SDF applies to the stated
            // topology.  SD files can contain different molecules, but if it passes this test
            // then there is good reason to trust the user.
            const int top_atom_count = topology_cache[top_idx].getAtomCount();
            if (frame_selection[j].getAtomCount() != top_atom_count) {
              const std::string base_err("The atom count in frame " + std::to_string(j) +
                                         " of SD file " +
                                         getBaseName(sysvec[i].getInputCoordinateFileName()) +
                                         " (" + std::to_string(frame_selection[j].getAtomCount()) +
                                         ") does not agree with the atom count in topology " +
                                         getBaseName(sysvec[i].getTopologyFileName()) + " (" +
                                         std::to_string(top_atom_count) + ").");
              switch (policy) {
              case ExceptionResponse::DIE:
                rtErr(base_err, "SystemCache");
              case ExceptionResponse::WARN:
                rtErr(base_err + "  This frame will be skipped.", "SystemCache");
                break;
              case ExceptionResponse::SILENT:
                break;
              }
            }
            coordinates_cache.push_back(frame_selection[j].exportPhaseSpace());
            topology_indices.push_back(top_idx);
            system_trajectory_names.push_back(sysvec[i].getTrajectoryFileName());
            system_checkpoint_names.push_back(sysvec[i].getCheckpointFileName());
            system_labels.push_back(sysvec[i].getLabel());
            system_trajectory_kinds.push_back(sysvec[i].getTrajectoryFileKind());
            system_checkpoint_kinds.push_back(sysvec[i].getCheckpointFileKind());
            system_count += 1;
          }
        }
        break;
      case CoordinateFileKind::UNKNOWN:
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The format of coordinate file " + sysvec[i].getTopologyFileName() +
                " for system " + std::to_string(i) + " could not be understood.",
                "UserSettings");
        case ExceptionResponse::WARN:
          rtWarn("The format of coordinate file " + sysvec[i].getTopologyFileName() +
                 " could not be understood.  System " + std::to_string(i + 1) +
                 " will be skipped.", "UserSettings");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        break;
      }
    }
  }
  
  // Collect examples of all systems, taking the first system in the list to use each topology
  // as the example of coordinates for that topology.
  const int top_count = topology_cache.size();
  const int system_count = coordinates_cache.size();
  example_indices.resize(top_count, -1);
  for (int i = 0; i < system_count; i++) {
    const int top_idx = topology_indices[i];
    if (example_indices[top_idx] == -1) {
      example_indices[top_idx] = i;
    }
  }

  // Make a list of all the systems making use of each topology, plus a bounds list to navigate it
  topology_cases.resize(system_count);
  topology_case_bounds.resize(top_count + 1, 0);
  for (int i = 0; i < system_count; i++) {
    topology_case_bounds[topology_indices[i]] += 1;
  }
  prefixSumInPlace<int>(&topology_case_bounds, PrefixSumType::EXCLUSIVE, "SystemCache");
  for (int i = 0; i < system_count; i++) {
    const int case_idx = topology_case_bounds[topology_indices[i]];
    topology_cases[case_idx] = i;
    topology_case_bounds[topology_indices[i]] = case_idx + 1;
  }
  for (int i = top_count; i > 0; i--) {
    topology_case_bounds[i] = topology_case_bounds[i - 1];
  }
  topology_case_bounds[0] = 0;

  // Create ChemicalFeatures objects to pair with each system--the majority of the work in
  // computing the ChemicalFeatures lies in ring detection and drawing the Lewis structure to
  // determine formal charges and bond orders, aspects which are invariant to coordinates.
  // Compute the ChemicalFeatures for each topology and then reassess the chiralities for each
  // system.
  features_cache.reserve(system_count);
  for (int i = 0; i < system_count; i++) {

    // Determine whether chemical features have been computed for this topology already.
    // If not, create them anew.  If so, copy the features and recompute chiral orientations.
    const int top_idx = topology_indices[i];
    if (i > example_indices[top_idx]) {
      features_cache.emplace_back(features_cache[example_indices[top_idx]]);
      features_cache[i].findChiralOrientations(CoordinateFrameReader(coordinates_cache[i]));
    }
    else {
      features_cache.emplace_back(&topology_cache[top_idx],
                                  CoordinateFrameReader(coordinates_cache[i]), map_chemfe_rotators,
                                  300.0, timer_in);
    }
  }

  // Create Exclusion masks of the appropriate kind
  std::vector<bool> need_static_mask(top_count, false);
  std::vector<bool> need_forward_mask(top_count, false);
  for (int i = 0; i < system_count; i++) {
    const int top_idx = topology_indices[i];
    switch (coordinates_cache[i].getUnitCellType()) {
    case UnitCellType::NONE:
      need_static_mask[top_idx] = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      need_forward_mask[top_idx] = true;
      break;
    }
  }
  static_masks_cache.reserve(top_count);
  forward_masks_cache.reserve(top_count);
  for (int i = 0; i < top_count; i++) {
    if (need_static_mask[i]) {
      static_masks_cache.emplace_back(&topology_cache[i]);
    }
    else {
      static_masks_cache.emplace_back(nullptr);
    }
    if (need_forward_mask[i]) {
      forward_masks_cache.emplace_back(&topology_cache[i]);
    }
    else {
      forward_masks_cache.emplace_back(nullptr);
    }
  }

  // Use the chemical features objects to make restraint apparatuses associated with the various
  // labels found in &restraint namelists.
  const int nrst_nml = rstcon.size();

  // CHECK
  printf("Detected %d &restraint namelists in the input.\n", nrst_nml);
  // END CHECK
  
  int nrst_labels = 0;
  std::vector<std::string> tmp_rst_labels;
  for (int i = 0; i < nrst_nml; i++) {

    // CHECK
    printf("Apply to %s\n", rstcon[i].getSystemLabel().c_str());
    // END CHECK
  }
}

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const FilesControls &fcon, const ExceptionResponse policy,
                         const MapRotatableGroups map_chemfe_rotators, StopWatch *timer_in) :
    SystemCache(fcon, {}, policy, map_chemfe_rotators, timer_in)
{}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCount() const {
  return topology_cache.size();
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemTopologyIndex(const int coord_index) const {
  if (coord_index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(coord_index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getSystemTopologyIndex");
  }
  return topology_indices[coord_index];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getCoordinateExample(const int topology_index) const {
  return example_indices[topology_index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* SystemCache::getTopologyPointer(const int topology_index) const {
  if (topology_index >= static_cast<int>(topology_cache.size())) {
    rtErr("Index " + std::to_string(topology_index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyPointer");
  }
  return &topology_cache[topology_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<const AtomGraph*> SystemCache::getTopologyPointer() const {
  const size_t ntop = topology_cache.size();
  std::vector<const AtomGraph*> result(ntop);
  for (size_t i = 0; i < ntop; i++) {
    result[i] = &topology_cache[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& SystemCache::getTopologyReference(const int topology_index) const {
  if (topology_index >= static_cast<int>(topology_cache.size())) {
    rtErr("Index " + std::to_string(topology_index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyReference");
  }
  return topology_cache[topology_index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* SystemCache::getSystemTopologyPointer(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getSystemTopologyPointer");
  }
  const AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraph* SystemCache::getSystemTopologyPointer(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getSystemTopologyPointer");
  }
  AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
std::vector<const AtomGraph*> SystemCache::getSystemTopologyPointer() const {
  const size_t nsys = coordinates_cache.size();
  std::vector<const AtomGraph*> result(nsys);
  const AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &tp_data[topology_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> SystemCache::getSystemTopologyPointer() {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &tp_data[topology_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> SystemCache::getSystemTopologyPointerCC() const {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  const AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = const_cast<AtomGraph*>(&tp_data[topology_indices[i]]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& SystemCache::getSystemTopologyReference(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getSystemTopologyReference");
  }
  return topology_cache[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraph& SystemCache::getSystemTopologyReference(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getSystemTopologyReference");
  }
  return topology_cache[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace* SystemCache::getCoordinatePointer(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getCoordinatePointer");
  }
  const PhaseSpace* crd_cache_ptr = coordinates_cache.data();
  return &crd_cache_ptr[index];
}

//-------------------------------------------------------------------------------------------------
PhaseSpace* SystemCache::getCoordinatePointer(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getCoordinatePointer");
  }
  PhaseSpace* crd_cache_ptr = coordinates_cache.data();
  return &crd_cache_ptr[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<const PhaseSpace*> SystemCache::getCoordinatePointer() const {
  const PhaseSpace* ps_ptr = coordinates_cache.data();
  const int nsys = coordinates_cache.size();
  std::vector<const PhaseSpace*> result(nsys);
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &ps_ptr[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<PhaseSpace*> SystemCache::getCoordinatePointer() {
  PhaseSpace* ps_ptr = coordinates_cache.data();
  const int nsys = coordinates_cache.size();
  std::vector<PhaseSpace*> result(nsys);
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &ps_ptr[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace& SystemCache::getCoordinateReference(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getCoordinateReference");
  }
  return coordinates_cache[index];
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& SystemCache::getCoordinateReference(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getCoordinateReference");
  }
  return coordinates_cache[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<PhaseSpace>& SystemCache::getCoordinateReference() const {
  return coordinates_cache;
}

//-------------------------------------------------------------------------------------------------
std::vector<PhaseSpace>& SystemCache::getCoordinateReference() {
  return coordinates_cache;
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures* SystemCache::getFeaturesPointer(const int index) const {
  if (index >= static_cast<int>(features_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(features_cache.size()) + ".", "SystemCache", "getFeaturesPointer");
  }
  return &features_cache[index];
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures* SystemCache::getFeaturesPointer(const int index) {
  if (index >= static_cast<int>(features_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(features_cache.size()) + ".", "SystemCache", "getFeaturesPointer");
  }
  return &features_cache[index];
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures& SystemCache::getFeaturesReference(const int index) const {
  if (index >= static_cast<int>(features_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(features_cache.size()) + ".", "SystemCache", "getFeaturesReference");
  }
  return features_cache[index];
}

//-------------------------------------------------------------------------------------------------
ChemicalFeatures& SystemCache::getFeaturesReference(const int index) {
  if (index >= static_cast<int>(features_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(features_cache.size()) + ".", "SystemCache", "getFeaturesReference");
  }
  return features_cache[index];
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus* SystemCache::getRestraintPointer(const int index) const {
  if (index >= static_cast<int>(restraints_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(restraints_cache.size()) + ".", "SystemCache", "getRestraintPointer");
  }
  return &restraints_cache[index];
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus* SystemCache::getRestraintPointer(const int index) {
  if (index >= static_cast<int>(restraints_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(restraints_cache.size()) + ".", "SystemCache", "getRestraintPointer");
  }
  return &restraints_cache[index];
}

//-------------------------------------------------------------------------------------------------
const RestraintApparatus& SystemCache::getRestraintReference(const int index) const {
  if (index >= static_cast<int>(restraints_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(restraints_cache.size()) + ".", "SystemCache", "getRestraintReference");
  }
  return restraints_cache[index];
}

//-------------------------------------------------------------------------------------------------
RestraintApparatus& SystemCache::getRestraintReference(const int index) {
  if (index >= static_cast<int>(restraints_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(restraints_cache.size()) + ".", "SystemCache", "getRestraintReference");
  }
  return restraints_cache[index];
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMask* SystemCache::getSystemStaticMaskPointer(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getStaticMaskPointer");
  }
  const int top_idx = topology_indices[index];
  return &static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMask* SystemCache::getSystemStaticMaskPointer(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getStaticMaskPointer");
  }
  const int top_idx = topology_indices[index];
  return &static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const StaticExclusionMask& SystemCache::getSystemStaticMaskReference(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getStaticMaskReference");
  }
  const int top_idx = topology_indices[index];
  return static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
StaticExclusionMask& SystemCache::getSystemStaticMaskReference(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getStaticMaskReference");
  }
  const int top_idx = topology_indices[index];
  return static_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const ForwardExclusionMask* SystemCache::getSystemForwardMaskPointer(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getForwardMaskPointer");
  }
  const int top_idx = topology_indices[index];
  return &forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
ForwardExclusionMask* SystemCache::getSystemForwardMaskPointer(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getForwardMaskPointer");
  }
  const int top_idx = topology_indices[index];
  return &forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
const ForwardExclusionMask& SystemCache::getSystemForwardMaskReference(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getForwardMaskReference");
  }
  const int top_idx = topology_indices[index];
  return forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
ForwardExclusionMask& SystemCache::getSystemForwardMaskReference(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache",
          "getForwardMaskReference");
  }
  const int top_idx = topology_indices[index];
  return forward_masks_cache[top_idx];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCaseCount(const int topology_index) const {
  return topology_case_bounds[topology_index + 1] - topology_case_bounds[topology_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SystemCache::getTopologicalCases(const int topology_index) const {
  const int llim = topology_case_bounds[topology_index];
  const int hlim = topology_case_bounds[topology_index + 1];
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = topology_cases[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::getSystemTrajectoryName(const int system_index) const {
  if (system_index < 0 || system_index >= static_cast<int>(system_trajectory_names.size())) {
    rtErr("Index " + std::to_string(system_index) + " is invalid for an array of length " +
          std::to_string(system_trajectory_names.size()) + ".", "SystemCache",
          "getTrajectoryFileName");
  }
  return system_trajectory_names[system_index];
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::getSystemCheckpointName(const int system_index) const {
  if (system_index < 0 || system_index >= static_cast<int>(system_checkpoint_names.size())) {
    rtErr("Index " + std::to_string(system_index) + " is invalid for an array of length " +
          std::to_string(system_checkpoint_names.size()) + ".", "SystemCache",
          "getCheckpointFileName");
  }
  return system_checkpoint_names[system_index];
}

//-------------------------------------------------------------------------------------------------
std::string SystemCache::getSystemLabel(const int system_index) const {
  if (system_index < 0 || system_index >= static_cast<int>(system_labels.size())) {
    rtErr("Index " + std::to_string(system_index) + " is invalid for an array of length " +
          std::to_string(system_labels.size()) + ".", "SystemCache", "getSystemLabel");
  }
  return system_labels[system_index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getSystemTrajectoryKind(const int system_index) const {
  if (system_index < 0 || system_index >= static_cast<int>(system_trajectory_kinds.size())) {
    rtErr("Index " + std::to_string(system_index) + " is invalid for an array of length " +
          std::to_string(system_trajectory_kinds.size()) + ".", "SystemCache",
          "getSystemTrajectoryKind");
  }
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFileKind SystemCache::getSystemCheckpointKind(const int system_index) const {
  if (system_index < 0 || system_index >= static_cast<int>(system_checkpoint_kinds.size())) {
    rtErr("Index " + std::to_string(system_index) + " is invalid for an array of length " +
          std::to_string(system_checkpoint_kinds.size()) + ".", "SystemCache",
          "getSystemCheckpointKind");
  }
}
  
} // namespace synthesis
} // namespace stormm
