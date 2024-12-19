// -*-c++-*-
#ifndef STORMM_VIRTUAL_SITE_PLACEMENT_H
#define STORMM_VIRTUAL_SITE_PLACEMENT_H

#include <cmath>
#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Structure/local_arrangement.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"
#include "local_arrangement.h"

namespace stormm {
namespace structure {

using numerics::force_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::operator+;
using stmath::dot;
using stmath::project;
using stmath::crossProduct;
using synthesis::AtomGraphSynthesis;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::SyAtomUpdateKit;
using synthesis::SyValenceKit;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::VirtualSiteKind;
using topology::VirtualSiteKit;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Given the coordinate arrays, indices of the parent atom and virtual site, and the
///        Cartesian X, Y, and Z displacements of the virtual site from its parent atom, place the
///        virtual site relative to its parent atom.  This is a common series of operations among
///        all virtual site frames, and this function encapsulates it for int95_t coordinates.
///
/// \param vsite_atom   Index of the virtual site in the atom list, indicative of where to find its
///                     coordinates in the arrays xcrd, ycrd, and zcrd
/// \param parent_atom  The first frame atom, to which the virtual site will be considered to be in
///                     a 1:1 bonded relationship (the virtual site is, for purposes of creating an
///                     exclusions list, a part of the parent atom)
/// \param xdisp        Displacement of the virtual site from its parent atom along the Cartesian
///                     X axis
/// \param ydisp        Displacement of the virtual site along the Cartesian Y axis
/// \param zdisp        Displacement of the virtual site along the Cartesian Z axis
/// \param xcrd         Cartesian X coordinates of all particles
/// \param ycrd         Cartesian Y coordinates of all particles
/// \param zcrd         Cartesian Z coordinates of all particles
/// \param xcrd_ovrf    Array of overflow bits for the Cartesian X coordinates
/// \param ycrd_ovrf    Array of overflow bits for the Cartesian Y coordinates
/// \param zcrd_ovrf    Array of overflow bits for the Cartesian Z coordinates
template <typename Tcoord>
void drawVirtualSite(const int vsite_atom, const int parent_atom, const int95_t xdisp,
                     const int95_t ydisp, const int95_t zdisp, Tcoord *xcrd, Tcoord *ycrd,
                     Tcoord *zcrd, int* xcrd_ovrf, int* ycrd_ovrf, int* zcrd_ovrf);

/// \brief Basic function for placing a virtual site given coordinates and indices of the particle,
///        its frame, and parameters.
///
/// Overloaded:
///   - Take templated pointers to raw coordinate and box transformation matrices, plus a virtual
///     sites abstract from the topology matching the precision of the transformation matrices
///   - Take modifiable forms of any of the coordinate objects or their abstracts (most will imply
///     a specific coordinate representation)
///
/// \param xcrd               Cartesian X coordinates of all particles
/// \param ycrd               Cartesian Y coordinates of all particles
/// \param zcrd               Cartesian Z coordinates of all particles
/// \param umat               Transformation matrix taking coordinates into fractional space
/// \param invu               Transformation matrix taking coordinates back to real space
/// \param unit_cell          The system's unit cell type
/// \param vsite_atom         Index of the virtual site in the atom list, indicative of where to
///                           find its coordinates in the arrays xcrd, ycrd, and zcrd
/// \param parent_atom        The first frame atom, to which the virtual site will be considered to
///                           be in a 1:1 bonded relationship (the virtual site is, for purposes
///                           of creating an exclusions list, a part of the parent atom)
/// \param frame2_atom        The second frame atom
/// \param frame3_atom        The third frame atom
/// \param frame4_atom        The fourth frame atom
/// \param frame_type         The type of frame for the virtual site
/// \param gpos_scale_factor  Scaling factor to convert real coordinates into a fixed-precision
///                           representation (the inverse is computed internally)
/// \param xcrd_ovrf          Optional array of overflow bits for the Cartesian X coordinates
/// \param ycrd_ovrf          Optional array of overflow bits for the Cartesian Y coordinates
/// \param zcrd_ovrf          Optional array of overflow bits for the Cartesian Z coordinates
template <typename Tcoord, typename Tcalc>
void placeVirtualSite(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                      const double* invu, UnitCellType unit_cell, int vsite_atom, int parent_atom,
                      int frame2_atom, int frame3_atom, int frame4_atom,
                      VirtualSiteKind frame_type, Tcalc frame_d1, Tcalc frame_d2, Tcalc frame_d3,
                      Tcalc gpos_scale_factor = 1.0, int* xcrd_ovrf = nullptr,
                      int* ycrd_ovrf = nullptr, int* zcrd_ovrf = nullptr);

/// \brief Function for placing virtual site in a system or a collection of systems.
///
/// Overloaded:
///   - Take templated pointers to raw coordinate and box transformation matrices, plus a virtual
///     sites abstract from the topology matching the precision of the transformation matrices
///   - Take modifiable forms of any of the coordinate objects or their abstracts (most will imply
///     a specific coordinate representation)
///
/// Descriptions of input parameters follow from the individual placement function
/// placeVirtualSite(), above, in addition to:
///
/// \param ps        Coordinates of the system as a mutable PhaseSpace object
/// \param psw       Coordinates of the system as a mutable PhaseSpace abstract
/// \param cf        Coordinates of the system as a mutable CoordinateFrame object
/// \param cfw       Coordinates of the system as a mutable CoordinateFrame abstract
/// \param poly_ps   Coordinates and forces describing a a synthesis of systems
/// \param poly_psw  Abstract for a synthesis of system coordinates
/// \param vsk       Virtual sites details abstracted from the original topology
/// \param ag        System topology containing virtual site specifications        
/// \param poly_vk   Valence work unit details in an abstract from the topology of a synthesis of
///                  many systems
/// \param poly_auk  Virtual site instructions (along with other details of atom movement) for work
///                  units in a synthesis of systems
/// \{
template <typename Tcoord, typename Tcalc>
void placeVirtualSites(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<Tcalc> &vsk, Tcalc gpos_scale_factor = 1.0);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag, CoordinateCycle affix,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag, CoordinateCycle affix,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

template <typename Tcalc>
void placeVirtualSites(PhaseSpaceWriter *psw, const VirtualSiteKit<Tcalc> vsk,
                       bool edit_alternate_image = true);

template <typename Tcalc>
void placeVirtualSites(CoordinateFrameWriter *cfw, const VirtualSiteKit<Tcalc> vsk);

template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void placeVirtualSites(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                       const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk);

void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       CoordinateCycle affix, PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       CoordinateCycle affix, PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);

void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       PrecisionModel prec = PrecisionModel::DOUBLE);
/// \}

/// \brief Transmit forces on virtual sites to their frame atoms, using double-precision math
///        throughout.
///
/// Overloaded:
///   - Place a single virtual site based on generic coordinate and force arrays with generic
///     details of the frame
///   - Take templated pointers to raw coordinate and box transformation matrices, plus a virtual
///     sites abstract from the topology matching the precision of the transformation matrices in a
///     single system or a synthesis of many systems
///   - Take modifiable forms of appropriate coordinate objects or their abstracts (each will imply
///     a specific coordinate representation)
///
/// \param xcrd                Cartesian X coordinates of all particles
/// \param ycrd                Cartesian Y coordinates of all particles
/// \param zcrd                Cartesian Z coordinates of all particles
/// \param xfrc                Cartesian X forces acting on all particles
/// \param yfrc                Cartesian Y forces acting on all particles
/// \param zfrc                Cartesian Z forces acting on all particles
/// \param umat                Transformation matrix taking coordinates into fractional space
/// \param invu                Transformation matrix taking coordinates back to real space
/// \param unit_cell           The system's unit cell type
/// \param ps                  Coordinates of the system as a mutable PhaseSpace object
/// \param cf                  Coordinates of the system as a mutable CoordinateFrame object
/// \param vsk                 Virtual sites details abstracted from the original topology
/// \param ag                  System topology containing virtual site specifications
/// \param gpos_scale_factor   Scaling factor to convert real coordinates into a fixed-precision
///                            representation (the inverse is computed internally)
/// \param force_scale_factor  Scaling factor to convert real-valued force components into a
///                            fixed-precision representation (the inverse is computed internally)
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, UnitCellType unit_cell, int vsite_atom,
                               int parent_atom, int frame2_atom, int frame3_atom, int frame4_atom,
                               VirtualSiteKind frame_type, Tcalc frame_d1, Tcalc frame_d2,
                               Tcalc frame_d3, Tcalc gpos_scale_factor = 1.0,
                               Tcalc force_scale_factor = 1.0, const int* xcrd_ovrf = nullptr,
                               const int* ycrd_ovrf = nullptr, const int* zcrd_ovrf = nullptr,
                               int* xfrc_ovrf = nullptr, int* yfrc_ovrf = nullptr,
                               int* zfrc_ovrf = nullptr);
  
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, const UnitCellType unit_cell,
                               const VirtualSiteKit<Tcalc> &vsk, Tcalc gpos_scale_factor = 1.0,
                               Tcalc force_scale_factor = 1.0);

void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag,
                               PrecisionModel prec = PrecisionModel::DOUBLE);

void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag,
                               PrecisionModel prec = PrecisionModel::DOUBLE);

template <typename Tcalc>
void transmitVirtualSiteForces(PhaseSpaceWriter *psw, const VirtualSiteKit<Tcalc> vsk);

template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void transmitVirtualSiteForces(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                               const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk);

void transmitVirtualSiteForces(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                               PrecisionModel prec = PrecisionModel::DOUBLE);

void transmitVirtualSiteForces(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                               PrecisionModel prec = PrecisionModel::DOUBLE);

/// \}
  
} // namespace structure
} // namespace stormm

#include "virtual_site_handling.tpp"

#endif
