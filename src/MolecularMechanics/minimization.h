// -*-c++-*-
#ifndef OMNI_MINIMIZATION_H
#define OMNI_MINIMIZATION_H

namespace omni {
namespace mm {

/// \brief Compute the move based on the computed forces.  This is a normalized vector along the
///        direction of the forces if the minimization is still in a steepest descent phase, or a
///        conjugate gradient move otherwise.
///
/// \param xfrc     Forces acting on all atoms in the Cartesian X direction (forces acting on
///                 virtual sites must have been transmitted prior to calling this function)
/// \param yfrc     Forces acting on all atoms in the Cartesian Y direction
/// \param zfrc     Forces acting on all atoms in the Cartesian Z direction
/// \param xmove    The move to apply to each atom in the Cartesian X direction, before scaling
///                 with the step size which will be done by the calling function
/// \param ymove    The move to apply to each atom in the Cartesian Y direction
/// \param zmove    The move to apply to each atom in the Cartesian Z direction
/// \param natom    Number of atoms (trusted length of xcrd, ycrd, ..., ymove, and zmove)
/// \param step     Current step number of the energy minimization
/// \param sd_step  Step number at which steepest descent optimization ends and conjugate gradient
///                 moves begin
void computeGradientMove(const double* xfrc, const double* yfrc, const double* zfrc,
                         double* xmove, double* ymove, double* zmove, int natom, int step,
                         int sd_steps);

/// \brief
void moveParticles(double* xcrd, double* ycrd, double* zcrd, const double* xmove,
                   const double* ymove, const double* zmove, int natom, double dist);

/// \brief
void minimize(double* xcrd, double* ycrd, double* zcrd, double* umat, double* invu,
              UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc, double* xmove,
              double* ymove, double* zmove, double* x_cg_temp, double* y_cg_temp,
              double* z_cg_temp, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const RestraintApparatusDpReader &rar, const VirtualSiteKit<double> &vsk,
              const StaticExclusionMask &se, const MinimizeControls &mincon);
 
} // namespace mm
} // namespace omni

#endif
