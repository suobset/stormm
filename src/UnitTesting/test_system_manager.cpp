#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "test_system_manager.h"

namespace stormm {
namespace testing {

using constants::getEnumerationName;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using diskutil::getBaseName;
using diskutil::osSeparator;
  
//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager() :
    topology_base{std::string("")}, topology_extn{std::string("")},
    coordinate_base{std::string("")}, coordinate_extn{std::string("")}, 
    fault_response{ExceptionResponse::WARN}, fault_found{false},
    topology_names{}, coordinate_names{}, topology_exist{}, topology_success{}, coordinate_exist{},
    coordinate_success{}, all_topologies{}, all_coordinates{}
{}
  
//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::string &topology_base_in,
                                     const std::string &topology_extn_in,
                                     const std::vector<std::string> &topology_names_in,
                                     const std::string &coordinate_base_in,
                                     const std::string &coordinate_extn_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager()
{
  topology_base = topology_base_in;
  topology_extn = topology_extn_in;
  coordinate_base = coordinate_base_in;
  coordinate_extn = coordinate_extn_in;
  const size_t n_tops = topology_names_in.size();
  const size_t n_crds = coordinate_names_in.size();
  if (n_tops != n_crds) {
    rtErr("The number of topology (" + std::to_string(n_tops) + ") and coordinate (" +
          std::to_string(n_crds) + ") file names must be equal.", "TestSystemManager");
  }
  system_count = n_tops;
  topology_names.reserve(system_count);
  coordinate_names.reserve(system_count);
  const char osc = osSeparator();
  for (int i = 0; i < system_count; i++) {
    topology_names.emplace_back(topology_base_in + osc + topology_names_in[i] + '.' +
                                topology_extn_in);
    coordinate_names.emplace_back(coordinate_base_in + osc + coordinate_names_in[i] + '.' +
                                  coordinate_extn_in);
  }
  fault_response = fault_response_in;
  all_go_response = all_go_response_in;
  topology_exist.resize(system_count);
  topology_success.resize(system_count);
  coordinate_exist.resize(system_count);
  coordinate_success.resize(system_count);
  all_topologies.reserve(system_count);
  for (int i = 0; i < system_count; i++) {
    topology_exist[i] = (getDrivePathType(topology_names[i]) == DrivePathType::FILE);
    coordinate_exist[i] = (getDrivePathType(coordinate_names[i]) == DrivePathType::FILE);
    if (topology_exist[i]) {
      try {
        all_topologies.emplace_back(topology_names[i], ExceptionResponse::DIE);
        topology_success[i] = true;
      }
      catch (std::runtime_error) {
        all_topologies.emplace_back();
        topology_success[i] = false;
        fault_found = true;
      }
    }
    else {
      all_topologies.emplace_back();
      topology_success[i] = false;
      fault_found = true;
    }
    if (coordinate_exist[i]) {
      try {
        all_coordinates.emplace_back(coordinate_names[i]);
        coordinate_success[i] = true;
      }
      catch (std::runtime_error) {
        all_coordinates.emplace_back();
        coordinate_success[i] = false;
        fault_found = true;
      }      
    }
    else {
      all_coordinates.emplace_back();
      coordinate_success[i] = false;
      fault_found = true;
    }
    if (topology_success[i] && coordinate_success[i]) {
      if (all_topologies[i].getAtomCount() != all_coordinates[i].getAtomCount()) {
        rtErr("Atom counts in topology " + getBaseName(topology_names[i]) + " (" +
              std::to_string(all_topologies[i].getAtomCount()) + ") and coordinate file " +
              getBaseName(coordinate_names[i]) + " (" +
              std::to_string(all_coordinates[i].getAtomCount()) + ") do not agree.",
              "TestSystemManager");
      }
    }
  }

  // Report errors
  if (fault_found) {
    std::string file_errors;
    for (int i = 0; i < system_count; i++) {
      if (topology_exist[i] == false) {
        file_errors += "  - " + getBaseName(topology_names[i]) + " (does not exist)\n";
      }
      else if (topology_success[i] == false) {
        file_errors += "  - " + getBaseName(topology_names[i]) + " (unreadable)\n";
      }
      if (coordinate_exist[i] == false) {
        file_errors += "  - " + getBaseName(coordinate_names[i]) + " (does not exist)\n";
      }
      else if (coordinate_success[i] == false) {
        file_errors += "  - " + getBaseName(coordinate_names[i]) + " (unreadable)\n";
      }
    }
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Some critical files were not found or were unreadable.  These include:\n" +
            file_errors + "This is unrecoverable.", "TestSystemManager");
    case ExceptionResponse::WARN:
      rtWarn("Some critical files were not found or were unreadable.  These include:\n" +
             file_errors + "Subsequent tests will be handled with priority " +
             getEnumerationName(getTestingStatus()) + ".", "TestSystemManager");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::string &topology_base_in,
                                     const std::vector<std::string> &topology_names_in,
                                     const std::string &coordinate_base_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager(topology_base_in, std::string(""), topology_names_in, coordinate_base_in,
                      std::string(""), coordinate_names_in, policy, fault_response_in,
                      all_go_response_in)
{}

//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::vector<std::string> &topology_names_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager(std::string(""), std::string(""), topology_names_in, std::string(""),
                      std::string(""), coordinate_names_in, policy, fault_response_in,
                      all_go_response_in)
{}

//-------------------------------------------------------------------------------------------------
void TestSystemManager::checkIndexing(const int index, const char* caller) const {
  if (index < 0 || index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "TestSystemManager", caller);
  }
}

//-------------------------------------------------------------------------------------------------
int TestSystemManager::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyBasePath() const {
  return topology_base;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateBasePath() const {
  return coordinate_base;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyExtension() const {
  return topology_extn;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateExtension() const {
  return coordinate_extn;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyFile(const int index) const {
  checkIndexing(index, "getTopologyFile");
  return topology_names[index];
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateFile(const int index) const {
  checkIndexing(index, "getCoordinateFile");
  return coordinate_names[index];
}

//-------------------------------------------------------------------------------------------------
TestPriority TestSystemManager::getTestingStatus() const {
  return (fault_found) ? fault_response : all_go_response;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace TestSystemManager::exportPhaseSpace(const int index) const {
  checkIndexing(index, "exportPhaseSpace");
  return all_coordinates[index];
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame TestSystemManager::exportCoordinateFrame(const int index) const {
  checkIndexing(index, "exportCoordinateFrame");
  return CoordinateFrame(all_coordinates[index]);
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace& TestSystemManager::viewCoordinates(const int index) const {
  checkIndexing(index, "viewCoordinates");
  return all_coordinates[index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& TestSystemManager::getTopologyReference(const int index) const {
  checkIndexing(index, "viewTopology");
  return all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph>& TestSystemManager::getTopologyReference() const {
  return all_topologies;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* TestSystemManager::getTopologyPointer(const int index) const {
  checkIndexing(index, "viewTopology");
  return &all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> TestSystemManager::getTopologyPointer() {
  std::vector<AtomGraph*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = &all_topologies[i];
  }
  return result;
}

} // namespace testing
} // namespace stormm
