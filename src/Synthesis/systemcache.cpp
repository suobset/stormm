#include "systemcache.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "Potential/scorecard.h"
#include "Potential/valence_potential.h"
#include "Topology/atomgraph_abstracts.h"

namespace omni {
namespace synthesis {

using omni::diskutil::getBaseName;
using omni::diskutil::splitPath;
using omni::energy::evaluateBondTerms;
using omni::energy::evaluateAngleTerms;
using omni::energy::ScoreCard;
using omni::namelist::MoleculeSystem;
using omni::parse::findStringInVector;
using omni::topology::ValenceKit;
using omni::trajectory::detectCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache() :
    topology_cache{},
    coordinates_cache{},
    topology_indices{},
    example_indices{}
{}

//-------------------------------------------------------------------------------------------------
SystemCache::SystemCache(const FilesControls &fcon, const ExceptionResponse policy) :
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
  }
  
  // Test each coordinate set with respect to topologies of the appropriate size.  Evaluate bond
  // and angle energy, and take the best scoring result as the indicator of which topology
  // describes which coordinate set.  Make new systems based on each coordinate set.
  const int n_unique_sizes = unique_topology_sizes.size();
  ScoreCard sc(1);
  std::string trajectory_base, trajectory_ext, restart_base, restart_ext;
  splitPath(fcon.getTrajectoryFileName(), &trajectory_base, &trajectory_ext);
  splitPath(fcon.getCheckpointFileName(), &restart_base, &restart_ext);
  std::vector<bool> topology_in_use(n_free_top, false);
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
      const std::string trajectory_middle = (orig_base.size() > 0) ? orig_base : orig_ext;

      // Add this pair to the list of systems
      sysvec.push_back(MoleculeSystem(topology_cache[best_topology].getFileName(),
                                      tmp_coordinates_cache[icrdj].getFileName(),
                                      trajectory_base + trajectory_middle + trajectory_ext,
                                      restart_base + trajectory_middle + restart_ext, 0,
                                      tmp_coordinates_frame_count[icrdj], 1,
                                      tmp_coordinates_kind[icrdj],
                                      fcon.getOutputCoordinateFormat(), 
                                      fcon.getCheckpointFormat()));

      // Note that the topology is used
      topology_in_use[best_topology] = true;
    }
  }

  // Loop back over the systems (now representing all entries, including the paired free topologies
  // and coordinate sets).  If the topology has already been read, don't read it again.  Read
  // coordinates (and perhaps velocities, if available) into phase space objects.
  std::vector<std::string> current_topology_holdings;
  current_topology_holdings.reserve(n_free_top);
  for (int i = 0; i < n_free_top; i++) {
    if (topology_in_use[i]) {
      current_topology_holdings.push_back(topology_cache[i].getFileName());
    }
    else {
      topology_cache.erase(topology_cache.begin() + i);
    }
  }
  nsys = sysvec.size();
  for (int i = 0; i < nsys; i++) {
    int top_idx = findStringInVector(current_topology_holdings, sysvec[i].getTopologyFileName());
    bool topology_ok = false;
    if (top_idx >= current_topology_holdings.size()) {
      try {
        topology_cache.push_back(AtomGraph(sysvec[i].getTopologyFileName()));
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
      bool coordinates_ok = false;
      try {
        coordinates_cache.push_back(PhaseSpace(sysvec[i].getInputCoordinateFileName(),
                                               sysvec[i].getInputCoordinateFileKind()));
        coordinates_ok = true;
      }
      catch (std::runtime_error) {
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
      }
      if (coordinates_ok) {
        topology_indices.push_back(top_idx);
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
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getSystemCount() const {
  return coordinates_cache.size();
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyCount() const {
  return topology_cache.size();
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getTopologyIndex(const int coord_index) const {
  if (coord_index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(coord_index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyIndex");
  }
  return topology_indices[coord_index];
}

//-------------------------------------------------------------------------------------------------
int SystemCache::getCoordinateExample(const int topology_index) const {
  return example_indices[topology_index];
}
                                      
//-------------------------------------------------------------------------------------------------
const AtomGraph* SystemCache::getTopologyPointer(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyPointer");
  }
  const AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraph* SystemCache::getTopologyPointer(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyPointer");
  }
  AtomGraph* tp_cache_ptr = topology_cache.data();
  return &tp_cache_ptr[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
std::vector<const AtomGraph*> SystemCache::getTopologyPointer() const {
  const size_t nsys = coordinates_cache.size();
  std::vector<const AtomGraph*> result(nsys);
  const AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &tp_data[topology_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> SystemCache::getTopologyPointer() {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = &tp_data[topology_indices[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> SystemCache::getTopologyPointerCC() const {
  const size_t nsys = coordinates_cache.size();
  std::vector<AtomGraph*> result(nsys);
  const AtomGraph* tp_data = topology_cache.data();
  for (size_t i = 0; i < nsys; i++) {
    result[i] = const_cast<AtomGraph*>(&tp_data[topology_indices[i]]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& SystemCache::getTopologyReference(const int index) const {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyReference");
  }
  return topology_cache[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
AtomGraph& SystemCache::getTopologyReference(const int index) {
  if (index >= static_cast<int>(coordinates_cache.size())) {
    rtErr("Index " + std::to_string(index) + " is invalid for an array of length " +
          std::to_string(coordinates_cache.size()) + ".", "SystemCache", "getTopologyReference");
  }
  return topology_cache[topology_indices[index]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph>& SystemCache::getTopologyReference() const {
  return topology_cache;
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph>& SystemCache::getTopologyReference() {
  return topology_cache;
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

} // namespace synthesis
} // namespace omni
