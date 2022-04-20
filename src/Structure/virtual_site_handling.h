// -*-c++-*-
#ifndef OMNI_VIRTUAL_SITE_PLACEMENT_H
#define OMNI_VIRTUAL_SITE_PLACEMENT_H

#include <cmath>
include "Math/vector_ops.h"
#include "local_arrangement.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace omni {
namespace structure {

using math::dot;
using math::project;
using math::crossProduct;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::VirtualSiteKind;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Reference function for placing virtual sites, using double-precision math throughout.
///
/// Overloaded:
///   - Take the raw coordinate and box transformation pointers plus a virtual sites abstract
///   - Take a modifiable PhaseSpace object and the system topology
///   - Take a modifiable CoordinateFrame object and the system topology
///   - Accept the topology by const reference or by const pointer
///
/// \param xcrd       Cartesian X coordinates of all atoms
/// \param ycrd       Cartesian Y coordinates of all atoms
/// \param zcrd       Cartesian Z coordinates of all atoms
/// \param umat       Transformation matrix taking coordinates into unit cell fractional space
/// \param invu       Transformation matrix taking coordinates back to real space
/// \param unit_cell  The system's unit cell type
/// \param ps         Coordinates of the system as a mutable PhaseSpace object
/// \param cf         Coordinates of the system as a mutable CoordinateFrame object
/// \param vsk        Virtual sites details abstracted from the original topology
/// \param ag         System topology containing virtual site specifications        
/// \{
void placeVirtualSites(double* xcrd, double* ycrd, double* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<double> &vsk);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag);
/// \}

/// \brief Transmit forces on virtual sites to their frame atoms, using double-precision math
///        throughout.
///
/// Overloaded:
///   - Pass the topology by const pointer or by const reference
///   - Accept the virtual site abstract directly
///
/// \param ps   Coordinates and forces of the system
/// \param ag   System topology containing virtual site specifications        
/// \param vsk  Virtual site abstract (obtained from the original topology)
/// \{
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag);

void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag);

void transmitVirtualSiteForces(PhaseSpace *ps, const VirtualSiteKit<double> &vsk);
/// \}
  
} // namespace structure
} // namespace omni

#endif
