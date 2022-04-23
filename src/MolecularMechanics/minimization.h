// -*-c++-*-
#ifndef OMNI_MINIMIZATION_H
#define OMNI_MINIMIZATION_H

#include "Namelists/nml_minimize.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace omni {
namespace mm {

using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using namelist::MinimizeControls;
using restraints::RestraintApparatus;
using restraints::RestraintApparatusDpReader;
using topology::AtomGraph;
using topology::NonbondedKit;
using topology::UnitCellType;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Compute the move based on the computed forces.  This is a normalized unit vector along
///        the direction of the forces if the minimization is still in a steepest descent phase,
///        or a conjugate gradient move otherwise.
///
/// \param xfrc       Forces acting on all atoms in the Cartesian X direction (forces acting on
///                   virtual sites must be transmitted to atoms with mass prior to calling this
///                   function)
/// \param yfrc       Forces acting on all atoms in the Cartesian Y direction
/// \param zfrc       Forces acting on all atoms in the Cartesian Z direction
/// \param xmove      The move to apply to each atom in the Cartesian X direction, before scaling
///                   with the step size which will be done by the calling function
/// \param ymove      The move to apply to each atom in the Cartesian Y direction
/// \param zmove      The move to apply to each atom in the Cartesian Z direction
/// \param natom      Number of atoms (trusted length of xcrd, ycrd, ..., ymove, and zmove)
/// \param step       Current step number of the energy minimization
/// \param sd_step    Step number at which steepest descent optimization ends and conjugate
///                   gradient moves begin
void computeGradientMove(const double* xfrc, const double* yfrc, const double* zfrc,
                         double* xmove, double* ymove, double* zmove, double* x_cg_temp,
                         double* y_cg_temp, double* z_cg_temp, int natom, int step, int sd_steps);

/// \brief Move particles based on a direction and a specified distance.  The unit vector spread
///        across xmove, ymove, and zmove provides the direction and dist the length to travel
///        along it.  Replace virtual particles after moving their frame atoms.
///
/// \param xcrd       Cartesian X coordinates of all particles
/// \param ycrd       Cartesian Y coordinates of all particles
/// \param zcrd       Cartesian Z coordinates of all particles
/// \param xmove      Cartesian X moves to apply to each particle
/// \param ymove      Cartesian Y moves to apply to each particle
/// \param zmove      Cartesian Z moves to apply to each particle
/// \param x_cg_temp  Cartesian X moves to apply to each particle
/// \param y_cg_temp  Cartesian Y moves to apply to each particle
/// \param z_cg_temp  Cartesian Z moves to apply to each particle
/// \param umat       Transformation matrix taking coordinates into fractional space
/// \param invu       Transformation matrix taking fractional coordinates back to real space
/// \param unit_cell  Type of simulation cell, to direct re-imaging for virtual site placement
/// \param vsk        Virtual site abstract from the original topology
/// \param natom      Number of atoms (trusted length of xcrd, ycrd, ..., ymove, and zmove)
/// \param dist       Distance along the gradient to travel (a multiplier for the unit vector)
void moveParticles(double* xcrd, double* ycrd, double* zcrd, const double* xmove,
                   const double* ymove, const double* zmove, const double* umat,
                   const double* invu, UnitCellType unit_cell, const VirtualSiteKit<double> &vsk,
                   int natom, double dist);

/// \brief Minimize a set of coordinates governed by a single topology and restraint apparatus.
///        This routine carries out line minimizations.
///
/// Overloaded:
///   - Accept a system in isolated boundary conditions, taking a StaticExclusionMask to describe
///     its permanent non-bonded exclusion list
///   - Accept raw pointers to coordinates, various pre-allocated arrays, and the relevant
///     parameter abstracts
///   - Accept a modifiable PhaseSpace object or one of its abstracts along with the relevant
///     topology and restraint apparatus, or abstracts thereof
///
/// \param xcrd       Cartesian X coordinates of all particles
/// \param ycrd       Cartesian Y coordinates of all particles
/// \param zcrd       Cartesian Z coordinates of all particles
/// \param xfrc       Forces acting on all atoms in the Cartesian X direction (pre-allocated for
///                   speed, probably fed in from a PhaseSapce object or equivalent)
/// \param yfrc       Forces acting on all atoms in the Cartesian Y direction
/// \param zfrc       Forces acting on all atoms in the Cartesian Z direction
/// \param xmove      Cartesian X moves to apply to each particle (pre-allocated for speed and to
///                   replicate the layout of memory available on the GPU)
/// \param ymove      Cartesian Y moves to apply to each particle
/// \param zmove      Cartesian Z moves to apply to each particle
/// \param x_cg_temp  Temporary array reserved for conjugate gradient computations (pre-allocated
///                   for speed and to replicate the layout of memory available on the GPU)
/// \param y_cg_temp  Temporary array reserved for conjugate gradient computations
/// \param z_cg_temp  Temporary array reserved for conjugate gradient computations
/// \param vk         Valence parameters abstract from the original topology
/// \param nbk        Non-bonded parameters abstract from the original topology
/// \param rar        Restraints applicable to the system (supplement to the topology)
/// \param vsk        Virtual site abstract from the original topology
/// \param se         Static exclusion mask for systems with isolated boundary conditions
/// \param mincon     User input from a &minimize namelist
/// \{
void minimize(double* xcrd, double* ycrd, double* zcrd, double* xfrc, double* yfrc, double* zfrc,
              double* xmove, double* ymove, double* zmove, double* x_cg_temp, double* y_cg_temp,
              double* z_cg_temp, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const RestraintApparatusDpReader &rar, const VirtualSiteKit<double> &vsk,
              const StaticExclusionMaskReader &ser, const MinimizeControls &mincon);

void minimize(PhaseSpace *ps, const AtomGraph &ag, const RestraintApparatus &ra,
              const StaticExclusionMask &se, const MinimizeControls &mincon);

void minimize(PhaseSpaceWriter psw, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const RestraintApparatusDpReader &rar, const VirtualSiteKit<double> &vsk,
              const StaticExclusionMaskReader &ser, const MinimizeControls &mincon);
/// \}
  
} // namespace mm
} // namespace omni

#endif
