#include "minimization.h"
#include "mm_evaluation.h"

namespace omni {
namespace mm {

//-------------------------------------------------------------------------------------------------
void computeGradientMove(const double* xfrc, const double* yfrc, const double* zfrc,
                         double* xmove, double* ymove, double* zmove, const double* umat,
                         const double* invu, const UnitCellType unit_cell,
                         const VirtualSiteKit<double> &vsk, const int natom, const int step,
                         const int sd_steps) {
  if (step < sd_steps) {
    double msum = 0.0;
    for (int i = 0; i < natom; i++) {
      const double fx = xfrc[i];
      const double fy = yfrc[i];
      const double fz = zfrc[i];
      xmove[i] = fx;
      ymove[i] = fy;
      zmove[i] = fz;
      msum += (fx * fx) + (fy * fy) + (fz * fz);
    }
    msum = 1.0 / sqrt(msum);
    for (int i = 0; i < natom; i++) {
      xmove[i] *= msum;
      ymove[i] *= msum;
      zmove[i] *= msum;
    }
  }
  else {

    // Apply the conjugate gradient protocol 
  }
}

//-------------------------------------------------------------------------------------------------
void moveParticles(double* xcrd, double* ycrd, double* zcrd, const double* xmove,
                   const double* ymove, const double* zmove, const double* umat,
                   const double* invu, const UnitCellType unit_cell,
                   const VirtualSiteKit<double> &vsk, const int natom, const double dist) {
  for (int i = 0; i < natom; i++) {
    xcrd[i] += xmove[i] * dist;
    ycrd[i] += ymove[i] * dist;
    zcrd[i] += zmove[i] * dist;
  }
  placeVirtualSites(xcrd, ycrd, zcrd, umat, invu, unit_cell, vsk);
}
  
//-------------------------------------------------------------------------------------------------
void minimize(double* xcrd, double* ycrd, double* zcrd, double* xfrc, double* yfrc, double* zfrc,
              double* xmove, double* ymove, double* zmove, double* x_cg_temp, double* y_cg_temp,
              double* z_cg_temp, const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
              const RestraintApparatusDpReader &rar, const VirtualSiteKit<double> &vsk,
              const StaticExclusionMask &se, const MinimizeControls &mincon) {
  
  // Loop for the requested number of cycles
  ScoreCard sc(1);
  double move_scale = mincon.getInitialStep();
  double evec[4], mvec[4], abcd_coefs[4], d_abcd[3], amat[16], inva[16];
  for (int step = 0; step < mincon.getTotalCycles(); step++) {

    // Compute the forces on all particles
    for (int i = 0; i < vk.natom; i++) {
      xfrc[i] = 0.0;
      yfrc[i] = 0.0;
      zfrc[i] = 0.0;
    }
    evalNonbValeRestMM(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, &sc, vk, nbk,
                       rar, EvaluateForce::YES);
    transmitVirtualSiteForces<double, double>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, umat, invu,
                                              unit_cell, vsk);
    evec[0] = sc.reportTotalEnergy();
    mvec[0] = 0.0;

    // Generate the move    
    computeGradientMove(xfrc, yfrc, zfrc, xmove, ymove, zmove, x_cg_temp, y_cg_temp, z_cg_temp,
                        vk.natom, step, mincon.getSteepestDescentSteps());

    // Implement the move three times, compute the energy, and arrive at a minimum value.
    double move_scale_factor = 1.0;    
    for (int i = 0; i < 3; i++) {
      moveParticles(xcrd, ycrd, zcrd, xmove, ymove, zmove, nullptr, nullptr, UnitCellType::NONE,
                    vk.natom, vsk, move_scale);
      evalNonbValeRestMM(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE, xfrc, yfrc, zfrc,
                         &sc, vk, nbk, rar, EvaluateForce::NO);
      evec[i + 1] = sc.reportTotalEnergy();
      mvec[i + 1] = mvec[i] + move_scale_factor;
      if (evec[i + 1] < evec[i]) {
        move_scale *= 1.05;
        move_scale_factor *= 1.05;
      }
      else {
        move_scale *= 0.90
        move_scale_factor *= 0.90;
      }
    }
    
    // Solve the linear system of equations to arrive at a polynomial Ax^3 + Bx^2 + Cx + D = evec
    // to solve the energy surface for x = 0.0 (gets evec[0]), mvec[1] (gets evec[1]), mvec[2],
    // and mvec[3].  The move vector mvec is kept in units of the initial move, to protect the
    // numerics of solving this system of linear equations as the step size becomes very small.
    for (int i = 0; i < 4; i++) {
      const double x = mvec[i];
      amat[amcon + i     ] = x * x * x;
      amat[amcon + i +  4] = x * x;
      amat[amcon + i +  8] = x;
      amat[amcon + i + 12] = 1.0;
    }
    invertSquareMatrix(amat->data(), inva->data(), 4);
    matrixVectorMultiply(inva->data(), evec->data(), abcd_coefs->data(), 4, 4);

    // The cubic polynomial will have at most one local minimum.  Find the local minimum, and if
    // it is in the range [mvec[0] = 0.0, mvec[3]] (inclusive of the endpoints), then take it by
    // moving all coordinates that distance from their origin at the outset of the move.
    // Otherwise, compare evec[0] and evec[3] to find the new coordinates (leave them where they
    // are after completing the third move thus far, or reset them to where they started).  The
    // step size is adjusting the whole time.
    d_abcd[0] = 3.0 * abcd_coefs[0];
    d_abcd[1] = 2.0 * abcd_coefs[1];
    d_abcd[2] = abcd_coefs[2];
    const double sqrt_arg = (d_abcd[1] * d_abcd[1]) - (4.0 * d_abcd[0] * d_abcd[2]);
    if (sqrt_arg < 0.0) {

      // The cubic equation has no minima or maxima.  Check the extrema of the range to
      // ascertain the correct move.
      if (evec[0] < evec[3]) {
        moveParticles(xcrd, ycrd, zcrd, xmove, ymove, zmove, nullptr, nullptr, UnitCellType::NONE,
                      vk.natom, vsk, -mvec[3] * move_scale);
      }
    }
    else {
      
      // Finish solving the quadratic formula to obtain both solutions and evaluate the cubic
      // polynomial's second derivative at each point.
      const double ext_i  = (-d_abcd[1] + sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double ext_ii = (-d_abcd[1] - sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double min_pos = ((2.0 * d_abcd[0] * ext_i) + d_abcd[1] > 0.0) ? ext_i : ext_ii;
      moveParticles(xcrd, ycrd, zcrd, xmove, ymove, zmove, nullptr, nullptr, UnitCellType::NONE,
                    vk.natom, vsk, (min_pos - mvec[3]) * move_scale);
    }
  }
}

} // namespace mm
} // namespace omni
