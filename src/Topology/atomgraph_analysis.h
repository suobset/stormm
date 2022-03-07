// -*-c++-*-
#ifndef OMNI_ATOMGRAPH_ANALYSIS_H
#define OMNI_ATOMGRAPH_ANALYSIS_H

#include <vector>
#include "atomgraph.h"

namespace omni {
namespace topology {

/// \brief Apply a series of checks to identify the water model in use within a topology.  If there
///        is no water, the model will be set to NONE.
///
/// \param ag  The topology to analyze
WaterModel identifyWaterModel(const AtomGraph &ag);

/// \brief Identify the virtual site types present in a topology.
///
/// Overloaded:
///   - Take a const pointer to the topology
///   - Take a const reference to the topology
///   - Operate on the vector of virtual site types with a trusted indication of its length
/// 
/// \param ag        The topology of interest
/// \param vs_types  Array containing the frame types of all virtual sites in the system
/// \param nsite     The number of virtual sites in the system
/// \{
std::string listVirtualSiteFrameTypes(const AtomGraph &ag);

std::string listVirtualSiteFrameTypes(const AtomGraph *ag);

std::string listVirtualSiteFrameTypes(const int* vs_types, int nsite);
/// \}

/// \brief The essential work of the selectRotatingAtoms() function below, abstracted to allow
///        more flexibility in the way the function is called.
///
/// \param nbk     Non-bonded details (needed for exclusion lists) of the original topology
/// \param cdk     Chemical details (here, atom names) of the original topology
/// \param atom_i  The root atom of the rotatable bond
/// \param atom_j  The second atom of the rotatable bond.  This can be taken as the center of
///                coordinates for rotating the local system, if convenient.
std::vector<int> mapRotatingGroup(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                  const int atom_i, const int atom_j);
  
/// \brief Select atoms for rotation, given a rotatable bond found in a molecule's chemical
///        features.
///
/// Overloaded:
///   - Take a const reference to an AtomGraph object
///   - Take a const pointer to an AtomGraph object
///
/// \param ag      System topology
/// \param atom_i  The root atom of the rotatable bond
/// \param atom_j  The second atom of the rotatable bond.  This can be taken as the center of
///                coordinates for rotating the local system, if convenient.
/// \{
std::vector<int> selectRotatingAtoms(const AtomGraph &ag, int atom_i, int atom_j);
std::vector<int> selectRotatingAtoms(const AtomGraph *ag, int atom_i, int atom_j);
/// \}

  
} // namespace topology
} // namespace omni

#endif
