#include <cmath>
#include "Math/matrix_ops.h"
#include "Potential/energy_enumerators.h"
#include "Structure/virtual_site_handling.h"
#include "minimization.h"
#include "mm_evaluation.h"

namespace omni {
namespace mm {

using energy::EvaluateForce;
using math::invertSquareMatrix;
using math::matrixVectorMultiply;
using structure::placeVirtualSites;
using structure::transmitVirtualSiteForces;
  
//-------------------------------------------------------------------------------------------------
void computeGradientMove(double* xfrc, double* yfrc, double* zfrc, double* xprv_move,
                         double* yprv_move, double* zprv_move, double* x_cg_temp,
                         double* y_cg_temp, double* z_cg_temp, const int natom, const int step,
                         const int sd_steps) {

  // Apply the conjugate gradient protocol
  if (step >= sd_steps) {
    double gg = 0.0;
    double dgg = 0.0;
    for (int i = 0; i < natom; i++) {
      gg += (xprv_move[i] * xprv_move[i]) + (yprv_move[i] * yprv_move[i]) +
            (zprv_move[i] * zprv_move[i]);
      const double dfx = (xfrc[i] - xprv_move[i]) * xfrc[i];
      const double dfy = (yfrc[i] - yprv_move[i]) * yfrc[i];
      const double dfz = (zfrc[i] - zprv_move[i]) * zfrc[i];
      dgg += dfx + dfy + dfz;
    }
    const double gam = (step == sd_steps) ? 0.0 : dgg / gg;
    for (int i = 0; i < natom; i++) {
      xprv_move[i] = xfrc[i];
      yprv_move[i] = yfrc[i];
      zprv_move[i] = zfrc[i];
      x_cg_temp[i] = xprv_move[i] + (gam * x_cg_temp[i]);
      y_cg_temp[i] = yprv_move[i] + (gam * y_cg_temp[i]);
      z_cg_temp[i] = zprv_move[i] + (gam * z_cg_temp[i]);
      xfrc[i] = x_cg_temp[i];
      yfrc[i] = y_cg_temp[i];
      zfrc[i] = z_cg_temp[i];
    }
  }

  // Normalize the force vector to get the move
  double msum = 0.0;
  for (int i = 0; i < natom; i++) {
    const double fx = xfrc[i];
    const double fy = yfrc[i];
    const double fz = zfrc[i];
    msum += (fx * fx) + (fy * fy) + (fz * fz);
  }
  msum = 1.0 / sqrt(msum);
  for (int i = 0; i < natom; i++) {
    xfrc[i] *= msum;
    yfrc[i] *= msum;
    zfrc[i] *= msum;
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
ScoreCard minimize(double* xcrd, double* ycrd, double* zcrd, double* xfrc, double* yfrc,
                   double* zfrc, double* xprv_move, double* yprv_move, double* zprv_move,
                   double* x_cg_temp, double* y_cg_temp, double* z_cg_temp,
                   const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                   const RestraintApparatusDpReader &rar, const VirtualSiteKit<double> &vsk,
                   const StaticExclusionMaskReader &ser, const MinimizeControls &mincon,
                   const int nrg_scale_bits) {
  
  // Loop for the requested number of cycles
  ScoreCard sc(1, mincon.getTotalCycles() + 1, nrg_scale_bits), sc_temp(1);
  double move_scale = mincon.getInitialStep();
  double evec[4], mvec[4], abcd_coefs[4], d_abcd[3], amat[16], inva[16];
  for (int i = 0; i < vk.natom; i++) {
    xprv_move[i] = 0.0;
    yprv_move[i] = 0.0;
    zprv_move[i] = 0.0;
    x_cg_temp[i] = 0.0;
    y_cg_temp[i] = 0.0;
    z_cg_temp[i] = 0.0;
  }
  for (int step = 0; step < mincon.getTotalCycles(); step++) {

    // Compute the forces on all particles
    for (int i = 0; i < vk.natom; i++) {
      xfrc[i] = 0.0;
      yfrc[i] = 0.0;
      zfrc[i] = 0.0;
    }
    evalNonbValeRestMM(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE, xfrc, yfrc, zfrc,
                       &sc, vk, nbk, ser, rar, EvaluateForce::YES, 0, step);
    transmitVirtualSiteForces<double, double>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                              UnitCellType::NONE, vsk);
    evec[0] = sc.reportTotalEnergy();
    mvec[0] = 0.0;

    // Generate the move    
    computeGradientMove(xfrc, yfrc, zfrc, xprv_move, yprv_move, zprv_move, x_cg_temp, y_cg_temp,
                        z_cg_temp, vk.natom, step, mincon.getSteepestDescentCycles());

    // Implement the move three times, compute the energy, and arrive at a minimum value.
    double move_scale_factor = 1.0;    
    for (int i = 0; i < 3; i++) {
      moveParticles(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr, UnitCellType::NONE,
                    vsk, vk.natom, move_scale);
      evalNonbValeRestMM(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE, xfrc, yfrc, zfrc,
                         &sc_temp, vk, nbk, ser, rar, EvaluateForce::NO, 0, step);
      evec[i + 1] = sc_temp.reportTotalEnergy();
      mvec[i + 1] = mvec[i] + move_scale_factor;
      if (evec[i + 1] < evec[i]) {
        const double idecay = 0.01 * static_cast<double>(i);
        move_scale *= 1.05 - idecay;
        move_scale_factor *= 1.05 - idecay;
      }
      else {
        const double idecay = 0.025 * static_cast<double>(i);        
        move_scale /= 1.05 - idecay;
        move_scale_factor /= 1.05 - idecay;
      }
    }
    
    // Solve the linear system of equations to arrive at a polynomial Ax^3 + Bx^2 + Cx + D = evec
    // to solve the energy surface for x = 0.0 (gets evec[0]), mvec[1] (gets evec[1]), mvec[2],
    // and mvec[3].  The move vector mvec is kept in units of the initial move, to protect the
    // numerics of solving this system of linear equations as the step size becomes very small.
    for (int i = 0; i < 4; i++) {
      const double x = mvec[i];
      amat[i     ] = x * x * x;
      amat[i +  4] = x * x;
      amat[i +  8] = x;
      amat[i + 12] = 1.0;
      abcd_coefs[i] = 0.0;
    }
    invertSquareMatrix(amat, inva, 4);
    matrixVectorMultiply(inva, evec, abcd_coefs, 4, 4);
    
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
        moveParticles(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr, UnitCellType::NONE,
                      vsk, vk.natom, -mvec[3] * move_scale);
      }
    }
    else {
      
      // Finish solving the quadratic formula to obtain both solutions and evaluate the cubic
      // polynomial's second derivative at each point.
      const double ext_i  = (-d_abcd[1] + sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double ext_ii = (-d_abcd[1] - sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double min_pos = ((2.0 * d_abcd[0] * ext_i) + d_abcd[1] > 0.0) ? ext_i : ext_ii;
      if (min_pos <= 0.0) {
        if (evec[3] > evec[0]) {
          moveParticles(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                        UnitCellType::NONE, vsk, vk.natom, -mvec[3] * move_scale);
        }
      }
      else if (min_pos < mvec[3]) {
        moveParticles(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr, UnitCellType::NONE,
                      vsk, vk.natom, (min_pos - mvec[3]) * move_scale);
      }
    }
    sc.incrementSampleCount();
  }
  return sc;
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph &ag, const RestraintApparatus &ra,
                   const StaticExclusionMask &se, const MinimizeControls &mincon,
                   const int nrg_scale_bits) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  PhaseSpaceWriter psw = ps->data();
  return minimize(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc, psw.zfrc, psw.xvel, psw.yvel,
                  psw.zvel, psw.xprv, psw.yprv, psw.zprv, vk, nbk, ra.dpData(), vsk, se.data(),
                  mincon, nrg_scale_bits);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpaceWriter psw, const ValenceKit<double> &vk,
                   const NonbondedKit<double> &nbk, const RestraintApparatusDpReader &rar,
                   const VirtualSiteKit<double> &vsk, const StaticExclusionMaskReader &ser,
                   const MinimizeControls &mincon, int nrg_scale_bits) {
  return minimize(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc, psw.zfrc, psw.xvel, psw.yvel,
                  psw.zvel, psw.xprv, psw.yprv, psw.zprv, vk, nbk, rar, vsk, ser, mincon,
                  nrg_scale_bits);
}
  
} // namespace mm
} // namespace omni
