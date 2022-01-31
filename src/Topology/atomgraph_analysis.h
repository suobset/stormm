// -*-c++-*-
#ifndef OMNI_ATOMGRAPH_ANALYSIS_H
#define OMNI_ATOMGRAPH_ANALYSIS_H

#include "atomgraph.h"

namespace omni {
namespace topology {

/// \brief Apply a series of checks to identify the water model in use within a topology.  If there
///        is no water, the model will be set to NONE.
///
/// \param ag  The topology to analyze
WaterModel identifyWaterModel(const AtomGraph &ag);

} // namespace topology
} // namespace omni

#endif
