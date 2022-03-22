#include <cmath>
#include "Math/vector_ops.h"
#include "local_arrangement.h"
#include "virtual_site_handling.h"

namespace omni {
namespace structure {

using math::dot;
using math::project;
using math::crossProduct;
using topology::VirtualSiteKind;

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(double* xcrd, double* ycrd, double* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<double> &vsk) {
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom  = vsk.vs_atoms[i];
    const int parent_atom = vsk.frame1_idx[i];
    const int frame2_atom = vsk.frame2_idx[i];
    const int param_idx = vsk.vs_param_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[param_idx])) {
    case VirtualSiteKind::FLEX_2:
      {
        double dx = xcrd[frame2_atom] - xcrd[parent_atom];
        double dy = ycrd[frame2_atom] - ycrd[parent_atom];
        double dz = zcrd[frame2_atom] - zcrd[parent_atom];
        imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
        const double p_f2_factor = vsk.dim1[param_idx];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz);
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        double dx = xcrd[frame2_atom] - xcrd[parent_atom];
        double dy = ycrd[frame2_atom] - ycrd[parent_atom];
        double dz = zcrd[frame2_atom] - zcrd[parent_atom];
        imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
        const double invr = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const double p_vs_distance = vsk.dim1[param_idx];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dx * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dy * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dz * invr);
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        const int frame3_atom = vsk.frame3_idx[i];
        double dx2 = xcrd[frame2_atom] - xcrd[parent_atom];
        double dy2 = ycrd[frame2_atom] - ycrd[parent_atom];
        double dz2 = zcrd[frame2_atom] - zcrd[parent_atom];
        double dx3 = xcrd[frame3_atom] - xcrd[parent_atom];
        double dy3 = ycrd[frame3_atom] - ycrd[parent_atom];
        double dz3 = zcrd[frame3_atom] - zcrd[parent_atom];
        imageCoordinates(&dx2, &dy2, &dz2, umat, invu, unit_cell, minimum_image);
        imageCoordinates(&dx3, &dy3, &dz3, umat, invu, unit_cell, minimum_image);
        const double p_f2_factor = vsk.dim1[param_idx];
        const double p_f3_factor = vsk.dim2[param_idx];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx2) + (p_f3_factor * dx3);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy2) + (p_f3_factor * dy3);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz2) + (p_f3_factor * dz3);
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        const int frame3_atom = vsk.frame3_idx[i];
        double dx23 = xcrd[frame3_atom] - xcrd[frame2_atom];
        double dy23 = ycrd[frame3_atom] - ycrd[frame2_atom];
        double dz23 = zcrd[frame3_atom] - zcrd[frame2_atom];
        imageCoordinates(&dx23, &dy23, &dz23, umat, invu, unit_cell, minimum_image);
        const double f2_f3_factor  = vsk.dim2[param_idx];
        const double x_midpoint = xcrd[frame2_atom] + (f2_f3_factor * dx23);
        const double y_midpoint = ycrd[frame2_atom] + (f2_f3_factor * dy23);
        const double z_midpoint = zcrd[frame2_atom] + (f2_f3_factor * dz23);
        double dxm = x_midpoint - xcrd[parent_atom];
        double dym = y_midpoint - ycrd[parent_atom];
        double dzm = z_midpoint - zcrd[parent_atom];
        imageCoordinates(&dxm, &dym, &dzm, umat, invu, unit_cell, minimum_image);
        const double invr = 1.0 / sqrt((dxm * dxm) + (dym * dym) + (dzm * dzm));
        const double p_vs_distance = vsk.dim1[param_idx];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dxm * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dym * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dzm * invr);
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        // Use small vectors to make use of a vector operation
        const int frame3_atom  = vsk.frame3_idx[i];
        double p_f2[3], f2_f3[3], f23_t_pf2[3];
        p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
        f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
        f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
        f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
        imageCoordinates( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell, minimum_image);
        const double invr2_p_f2 = 1.0 / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                         (p_f2[2] * p_f2[2]));

        // Compute the projection of f2_f3 (the displacement vector between frame atoms 2 and 3),
        // onto p_f2 (the displacement vector between the parent atom and frame atom 2).  Subtract
        // this from f2_f3 to get the part of f2_f3 that is perpendicular to p_f2.  Use the vector
        // that will ultimately store the result as temporary space to hold the projection.  In
        // more optimized code, the inverse squared displacement p_f2 computed for the projection
        // can be re-used, but it makes things more clear for the reference code to recompute that
        // quantity when needed.
        project(f2_f3, p_f2, f23_t_pf2, 3);
        f23_t_pf2[0] = f2_f3[0] - f23_t_pf2[0];
        f23_t_pf2[1] = f2_f3[1] - f23_t_pf2[1];
        f23_t_pf2[2] = f2_f3[2] - f23_t_pf2[2];
        const double invr_p_f2 = 1.0 / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                            (p_f2[2] * p_f2[2]));
        const double invr_t = 1.0 / sqrt((f23_t_pf2[0] * f23_t_pf2[0]) +
                                         (f23_t_pf2[1] * f23_t_pf2[1]) +
                                         (f23_t_pf2[2] * f23_t_pf2[2]));
        const double p_f2_factor = vsk.dim1[param_idx] * cos(vsk.dim2[param_idx]) * invr_p_f2;
        const double t_factor    = vsk.dim1[param_idx] * sin(vsk.dim2[param_idx]) * invr_t;
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * p_f2[0]) + (t_factor * f23_t_pf2[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * p_f2[1]) + (t_factor * f23_t_pf2[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * p_f2[2]) + (t_factor * f23_t_pf2[2]);
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        double p_f2[3], p_f3[3], f2_f3[3], pf2_x_pf3[3];
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
        imageCoordinates(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell, minimum_image);
        crossProduct(p_f2, p_f3, pf2_x_pf3);
        const double pf2_factor = vsk.dim1[param_idx];
        const double pf3_factor = vsk.dim2[param_idx];
        const double cr_factor  = vsk.dim3[param_idx];
        xcrd[vsite_atom] = xcrd[parent_atom] + (pf2_factor * p_f2[0]) + (pf3_factor * p_f3[0]) +
                           (cr_factor * pf2_x_pf3[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (pf2_factor * p_f2[1]) + (pf3_factor * p_f3[1]) +
                           (cr_factor * pf2_x_pf3[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (pf2_factor * p_f2[2]) + (pf3_factor * p_f3[2]) +
                           (cr_factor * pf2_x_pf3[2]);
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        const int frame4_atom  = vsk.frame4_idx[i];
        double p_f2[3], p_f3[3], p_f4[3], pf3_m_pf2[3], pf4_m_pf2[3], p_vs[3];
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
        p_f4[0] = xcrd[frame4_atom] - xcrd[parent_atom];
        p_f4[1] = ycrd[frame4_atom] - ycrd[parent_atom];
        p_f4[2] = zcrd[frame4_atom] - zcrd[parent_atom];
        imageCoordinates(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&p_f4[0], &p_f4[1], &p_f4[2], umat, invu, unit_cell, minimum_image);
        const double pf3_factor = vsk.dim1[param_idx];
        const double pf4_factor = vsk.dim2[param_idx];
        pf3_m_pf2[0] = (pf3_factor * p_f3[0]) - p_f2[0];
        pf3_m_pf2[1] = (pf3_factor * p_f3[1]) - p_f2[1];
        pf3_m_pf2[2] = (pf3_factor * p_f3[2]) - p_f2[2];
        pf4_m_pf2[0] = (pf4_factor * p_f4[0]) - p_f2[0];
        pf4_m_pf2[1] = (pf4_factor * p_f4[1]) - p_f2[1];
        pf4_m_pf2[2] = (pf4_factor * p_f4[2]) - p_f2[2];
        crossProduct(pf3_m_pf2, pf4_m_pf2, p_vs);
        const double pvs_factor = vsk.dim3[param_idx] /
                                  sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                       (p_vs[2] * p_vs[2]));
        xcrd[vsite_atom] = xcrd[parent_atom] + (pvs_factor * p_vs[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (pvs_factor * p_vs[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (pvs_factor * p_vs[2]);
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag) {
  PhaseSpaceWriter psw = ps->data();
  placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                    ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag) {
  PhaseSpaceWriter psw = ps->data();
  placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell,
                    ag->getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag) {
  CoordinateFrameWriter cfw = cf->data();
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu, cfw.unit_cell,
                    ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag) {
  CoordinateFrameWriter cfw = cf->data();
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu, cfw.unit_cell,
                    ag->getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag) {
  transmitVirtualSiteForces(ps, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag) {
  transmitVirtualSiteForces(ps, ag->getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const VirtualSiteKit<double> &vsk) {
  PhaseSpaceWriter psw = ps->data();
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom  = vsk.frame1_idx[i];
    const int frame2_atom  = vsk.frame2_idx[i];
    const int pidx = vsk.vs_param_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;

    // Copy the force on the virtual site into its own vector.  This will be needed in all
    // subsequent force transmissions.
    const double vs_frc[3] = { psw.xfrc[vsite_atom], psw.yfrc[vsite_atom], psw.zfrc[vsite_atom] };
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pidx])) {
    case VirtualSiteKind::FLEX_2:
      {
        const double p_f2_factor = vsk.dim1[pidx];
        psw.xfrc[parent_atom] += (1.0 - p_f2_factor) * psw.xfrc[vsite_atom];
        psw.yfrc[parent_atom] += (1.0 - p_f2_factor) * psw.yfrc[vsite_atom];
        psw.zfrc[parent_atom] += (1.0 - p_f2_factor) * psw.zfrc[vsite_atom];
        psw.xfrc[frame2_atom] += p_f2_factor * psw.xfrc[vsite_atom];
        psw.yfrc[frame2_atom] += p_f2_factor * psw.yfrc[vsite_atom];
        psw.zfrc[frame2_atom] += p_f2_factor * psw.zfrc[vsite_atom];
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        double p_f2[3], p_vs[3], vs_frc_proj[3], force_partition[3];
        p_f2[0] = psw.xcrd[frame2_atom] - psw.xcrd[parent_atom];
        p_f2[1] = psw.ycrd[frame2_atom] - psw.ycrd[parent_atom];
        p_f2[2] = psw.zcrd[frame2_atom] - psw.zcrd[parent_atom];
        imageCoordinates(&p_f2[0], &p_f2[1], &p_f2[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        p_vs[0] = psw.xcrd[vsite_atom] - psw.xcrd[parent_atom];
        p_vs[1] = psw.ycrd[vsite_atom] - psw.ycrd[parent_atom];
        p_vs[2] = psw.zcrd[vsite_atom] - psw.zcrd[parent_atom];
        imageCoordinates(&p_vs[0], &p_vs[1], &p_vs[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        const double invr_p_f2 = 1.0 / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                            (p_f2[2] * p_f2[2]));

        // Take the projection of the force on the virtual site onto the parent atom -> virtual
        // site displacement vector p_vs.  The force is partitioned between the two frame atoms
        // according to the distance between the parent atom and the virtual site (normalized by
        // the distance between the two frame atoms).
        project(vs_frc, p_vs, vs_frc_proj, 3);
        const double p_vs_distance = vsk.dim1[pidx];
        for (int j = 0; j < 3; j++) {
          force_partition[j] = invr_p_f2 * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
        }
        psw.xfrc[parent_atom] += vs_frc[0] - force_partition[0];
        psw.yfrc[parent_atom] += vs_frc[1] - force_partition[1];
        psw.zfrc[parent_atom] += vs_frc[2] - force_partition[2];
        psw.xfrc[frame2_atom] += force_partition[0];
        psw.yfrc[frame2_atom] += force_partition[1];
        psw.zfrc[frame2_atom] += force_partition[2];
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        const double p_f2_factor = vsk.dim1[pidx];
        const double p_f3_factor = vsk.dim2[pidx];
        psw.xfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * vs_frc[0];
        psw.yfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * vs_frc[1];
        psw.zfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * vs_frc[2];
        psw.xfrc[frame2_atom] += p_f2_factor * vs_frc[0];
        psw.yfrc[frame2_atom] += p_f2_factor * vs_frc[1];
        psw.zfrc[frame2_atom] += p_f2_factor * vs_frc[2];
        psw.xfrc[frame3_atom] += p_f3_factor * vs_frc[0];
        psw.yfrc[frame3_atom] += p_f3_factor * vs_frc[1];
        psw.zfrc[frame3_atom] += p_f3_factor * vs_frc[2];
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        double f2_f3[3], p_vs[3], p_mid[3], vs_frc_proj[3];
        double force_partition[3];
        const int frame3_atom  = vsk.frame3_idx[i];
        f2_f3[0] = psw.xcrd[frame3_atom] - psw.xcrd[frame2_atom];
        f2_f3[1] = psw.ycrd[frame3_atom] - psw.ycrd[frame2_atom];
        f2_f3[2] = psw.zcrd[frame3_atom] - psw.zcrd[frame2_atom];
        imageCoordinates(&f2_f3[0], &f2_f3[1], &f2_f3[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        p_vs[0] = psw.xcrd[vsite_atom] - psw.xcrd[parent_atom];
        p_vs[1] = psw.ycrd[vsite_atom] - psw.ycrd[parent_atom];
        p_vs[2] = psw.zcrd[vsite_atom] - psw.zcrd[parent_atom];
        imageCoordinates(&p_vs[0], &p_vs[1], &p_vs[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        const double f2_f3_factor  = vsk.dim2[pidx];
        p_mid[0] = psw.xcrd[frame2_atom] - psw.xcrd[parent_atom] + (f2_f3_factor * f2_f3[0]);
        p_mid[1] = psw.ycrd[frame2_atom] - psw.ycrd[parent_atom] + (f2_f3_factor * f2_f3[1]);
        p_mid[2] = psw.zcrd[frame2_atom] - psw.zcrd[parent_atom] + (f2_f3_factor * f2_f3[2]);
        imageCoordinates(&p_mid[0], &p_mid[1], &p_mid[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);

        // As with the fixed distance, two-point frame, compute the part of the force on the
        // virtual site that is perpendicular to the parent atom -> virtual site displacement
        // vector and use this to compute the partitioning between the parent atom and some
        // imaginary midpoint on the line between frame atom 2 and frame atom 3.  That's like
        // the FIXED_2 frame type.  The force on the midpoint is subsequently distributed, akin
        // to the FLEX_2 type, between atoms 2 and 3.
        project(vs_frc, p_vs, vs_frc_proj, 3);
        const double p_vs_distance = vsk.dim1[pidx];
        const double invr_p_mid = 1.0 / sqrt((p_mid[0] * p_mid[0]) + (p_mid[1] * p_mid[1]) +
                                             (p_mid[2] * p_mid[2]));
        for (int j = 0; j < 3; j++) {
          force_partition[j] = invr_p_mid * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
        }
        psw.xfrc[parent_atom] += vs_frc[0] - force_partition[0];
        psw.yfrc[parent_atom] += vs_frc[1] - force_partition[1];
        psw.zfrc[parent_atom] += vs_frc[2] - force_partition[2];
        psw.xfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[0];
        psw.yfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[1];
        psw.zfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[2];
        psw.xfrc[frame3_atom] += f2_f3_factor * force_partition[0];
        psw.yfrc[frame3_atom] += f2_f3_factor * force_partition[1];
        psw.zfrc[frame3_atom] += f2_f3_factor * force_partition[2];
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        double p_f2[3], f2_f3[3], f23_t_pf2[3], F1[3], F2[3], F3[3];
        const int frame3_atom  = vsk.frame3_idx[i];
        p_f2[0]  = psw.xcrd[frame2_atom] - psw.xcrd[parent_atom];
        p_f2[1]  = psw.ycrd[frame2_atom] - psw.ycrd[parent_atom];
        p_f2[2]  = psw.zcrd[frame2_atom] - psw.zcrd[parent_atom];
        f2_f3[0] = psw.xcrd[frame3_atom] - psw.xcrd[frame2_atom];
        f2_f3[1] = psw.ycrd[frame3_atom] - psw.ycrd[frame2_atom];
        f2_f3[2] = psw.zcrd[frame3_atom] - psw.zcrd[frame2_atom];
        imageCoordinates( &p_f2[0],  &p_f2[1],  &p_f2[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        imageCoordinates(&f2_f3[0], &f2_f3[1], &f2_f3[2], psw.umat, psw.invu, psw.unit_cell,
                         minimum_image);
        project(f2_f3, p_f2, f23_t_pf2, 3);
        f23_t_pf2[0] = f2_f3[0] - f23_t_pf2[0];
        f23_t_pf2[1] = f2_f3[1] - f23_t_pf2[1];
        f23_t_pf2[2] = f2_f3[2] - f23_t_pf2[2];
        const double invr2_p_f2 = 1.0 / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                         (p_f2[2] * p_f2[2]));
        const double invr_p_f2 = sqrt(invr2_p_f2);
        const double invr2_t = 1.0 / ((f23_t_pf2[0] * f23_t_pf2[0]) +
                                      (f23_t_pf2[1] * f23_t_pf2[1]) +
                                      (f23_t_pf2[2] * f23_t_pf2[2]));
        const double invr_t = sqrt(invr2_t);
        const double f1fac  = dot(p_f2, vs_frc, 3) * invr2_p_f2;
        const double f2fac  = dot(f23_t_pf2, vs_frc, 3) * invr2_t;
        const double p_f2_factor = vsk.dim1[pidx] * cos(vsk.dim2[pidx]) * invr_p_f2;
        const double t_factor    = vsk.dim1[pidx] * sin(vsk.dim2[pidx]) * invr_t;
        const double abbcOabab = dot(p_f2, f2_f3, 3) * invr2_p_f2;
        for (int j = 0; j < 3; j++) {
          F1[j] = vs_frc[j] - (f1fac * p_f2[j]);
          F2[j] = F1[j] - (f2fac * f23_t_pf2[j]);
          F3[j] = f1fac * f23_t_pf2[j];          
        }
        psw.xfrc[parent_atom] += vs_frc[0] - (p_f2_factor * F1[0]) +
                                 (t_factor * ((abbcOabab * F2[0]) + F3[0]));
        psw.yfrc[parent_atom] += vs_frc[1] - (p_f2_factor * F1[1]) +
                                 (t_factor * ((abbcOabab * F2[1]) + F3[1]));
        psw.zfrc[parent_atom] += vs_frc[2] - (p_f2_factor * F1[2]) +
                                 (t_factor * ((abbcOabab * F2[2]) + F3[2]));
        psw.xfrc[frame2_atom] += (p_f2_factor * F1[0]) -
                                 (t_factor * (F2[0] + (abbcOabab * F2[0]) + F3[0]));
        psw.yfrc[frame2_atom] += (p_f2_factor * F1[1]) -
                                 (t_factor * (F2[1] + (abbcOabab * F2[1]) + F3[1]));
        psw.zfrc[frame2_atom] += (p_f2_factor * F1[2]) -
                                 (t_factor * (F2[2] + (abbcOabab * F2[2]) + F3[2]));
        psw.xfrc[frame3_atom] += t_factor * F2[0]; 
        psw.yfrc[frame3_atom] += t_factor * F2[1];
        psw.zfrc[frame3_atom] += t_factor * F2[2];
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        double p_f2[3], p_f3[3], partition_f2[3], partition_f3[3];
        const double mf2_01 = vsk.dim3[pidx] * (psw.zcrd[frame3_atom] - psw.zcrd[parent_atom]);
        const double mf2_02 = vsk.dim3[pidx] * (psw.ycrd[frame3_atom] - psw.ycrd[parent_atom]);
        const double mf2_12 = vsk.dim3[pidx] * (psw.xcrd[frame3_atom] - psw.xcrd[parent_atom]);
        partition_f2[0] = ( vsk.dim1[pidx] * vs_frc[0]) - (mf2_01 * vs_frc[1]) +
                          (mf2_02 * vs_frc[2]);
        partition_f2[1] = ( mf2_01 * vs_frc[0]) + (vsk.dim1[pidx] * vs_frc[1]) -
                          (mf2_12 * vs_frc[2]);
        partition_f2[2] = (-mf2_02 * vs_frc[0]) + (mf2_12 * vs_frc[1]) +
                          (vsk.dim1[pidx] * vs_frc[2]);
        const double mf3_01 = vsk.dim3[pidx] * (psw.zcrd[frame2_atom] - psw.zcrd[parent_atom]);
        const double mf3_02 = vsk.dim3[pidx] * (psw.ycrd[frame2_atom] - psw.ycrd[parent_atom]);
        const double mf3_12 = vsk.dim3[pidx] * (psw.xcrd[frame2_atom] - psw.xcrd[parent_atom]);
        partition_f3[0] = ( vsk.dim2[pidx] * vs_frc[0]) + (mf3_01 * vs_frc[1]) -
                          (mf3_02 * vs_frc[2]);
        partition_f3[1] = (-mf3_01 * vs_frc[0]) + (vsk.dim2[pidx] * vs_frc[1]) +
                          (mf3_12 * vs_frc[2]);
        partition_f3[2] = ( mf3_02 * vs_frc[0]) - (mf3_12 * vs_frc[1]) +
                          (vsk.dim2[pidx] * vs_frc[2]);
        psw.xfrc[parent_atom] += vs_frc[0] - partition_f2[0] - partition_f3[0];
        psw.yfrc[parent_atom] += vs_frc[1] - partition_f2[1] - partition_f3[1];
        psw.zfrc[parent_atom] += vs_frc[2] - partition_f2[2] - partition_f3[2];
        psw.xfrc[frame2_atom] += partition_f2[0];
        psw.yfrc[frame2_atom] += partition_f2[1];
        psw.zfrc[frame2_atom] += partition_f2[2];
        psw.xfrc[frame3_atom] += partition_f3[0];
        psw.yfrc[frame3_atom] += partition_f3[1];
        psw.zfrc[frame3_atom] += partition_f3[2];
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        const int frame4_atom  = vsk.frame4_idx[i];
        double p_f2[3], rj_f3[3], rj_f4[3], rj_f34[3], rm[3], rt[3], fb[3], fc[3], fd[3];
        p_f2[0] = psw.xcrd[frame2_atom] -  psw.xcrd[parent_atom];
        p_f2[1] = psw.ycrd[frame2_atom] -  psw.ycrd[parent_atom];
        p_f2[2] = psw.zcrd[frame2_atom] -  psw.zcrd[parent_atom];
        rj_f3[0] = (vsk.dim1[pidx] * (psw.xcrd[frame3_atom] - psw.xcrd[parent_atom])) - p_f2[0];
        rj_f3[1] = (vsk.dim1[pidx] * (psw.ycrd[frame3_atom] - psw.ycrd[parent_atom])) - p_f2[1];
        rj_f3[2] = (vsk.dim1[pidx] * (psw.zcrd[frame3_atom] - psw.zcrd[parent_atom])) - p_f2[2];
        rj_f4[0] = (vsk.dim2[pidx] * (psw.xcrd[frame4_atom] - psw.xcrd[parent_atom])) - p_f2[0];
        rj_f4[1] = (vsk.dim2[pidx] * (psw.ycrd[frame4_atom] - psw.ycrd[parent_atom])) - p_f2[1];
        rj_f4[2] = (vsk.dim2[pidx] * (psw.zcrd[frame4_atom] - psw.zcrd[parent_atom])) - p_f2[2];
        rj_f34[0] = rj_f4[0] - rj_f3[0];
        rj_f34[1] = rj_f4[1] - rj_f3[1];
        rj_f34[2] = rj_f4[2] - rj_f3[2];
        crossProduct(rj_f3, rj_f4, rm);
        const double invr2_rm = 1.0 / ((rm[0] * rm[0]) + (rm[1] * rm[1]) + (rm[2] * rm[2]));
        const double invr_rm  = sqrt(invr2_rm);
        const double cfx = vsk.dim3[pidx] * invr_rm * vs_frc[0];
        const double cfy = vsk.dim3[pidx] * invr_rm * vs_frc[1];
        const double cfz = vsk.dim3[pidx] * invr_rm * vs_frc[2];
        crossProduct(rm, rj_f34, rt);
        rt[0] *= invr2_rm;
        rt[1] *= invr2_rm;
        rt[2] *= invr2_rm;
        fb[0] = (-rm[0] * rt[0] * cfx) + ((rj_f34[2] - (rm[1] * rt[0])) * cfy) -
                ( (rj_f34[1] + (rm[2] * rt[0])) * cfz);
        fb[1] = (-(rj_f34[2] + (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) +
                ( (rj_f34[0] - (rm[2] * rt[1])) * cfz);
        fb[2] = ((rj_f34[1] - (rm[0] * rt[2])) * cfx) - ((rj_f34[0] + (rm[1] * rt[2])) * cfy) -
                (rm[2] * rt[2] * cfz);
        rt[0] = ((rj_f4[1] * rm[2]) - (rj_f4[2] * rm[1])) * invr2_rm * vsk.dim1[pidx];
        rt[1] = ((rj_f4[2] * rm[0]) - (rj_f4[0] * rm[2])) * invr2_rm * vsk.dim1[pidx];
        rt[2] = ((rj_f4[0] * rm[1]) - (rj_f4[1] * rm[0])) * invr2_rm * vsk.dim1[pidx];
        fc[0] = (-rm[0] * rt[0] * cfx) - (((vsk.dim1[pidx] * rj_f4[2]) + (rm[1] * rt[0])) * cfy) +
                (((vsk.dim1[pidx] * rj_f4[1]) - (rm[2] * rt[0])) * cfz);
        fc[1] = (((vsk.dim1[pidx] * rj_f4[2]) - (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) -
                (((vsk.dim1[pidx] * rj_f4[0]) + (rm[2] * rt[1])) * cfz);
        fc[2] = (-((vsk.dim1[pidx] * rj_f4[1]) + (rm[0] * rt[2])) * cfx) +
                (((vsk.dim1[pidx] * rj_f4[0]) - (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
        rt[0] = ((rm[1] * rj_f3[2]) - (rm[2] * rj_f3[1])) * invr2_rm * vsk.dim2[pidx];
        rt[1] = ((rm[2] * rj_f3[0]) - (rm[0] * rj_f3[2])) * invr2_rm * vsk.dim2[pidx];
        rt[2] = ((rm[0] * rj_f3[1]) - (rm[1] * rj_f3[0])) * invr2_rm * vsk.dim2[pidx];
        fd[0] = (-rm[0] * rt[0] * cfx) + (((vsk.dim2[pidx] * rj_f3[2]) - (rm[1] * rt[0])) * cfy) -
                (((vsk.dim2[pidx] * rj_f3[1]) + (rm[2] * rt[0])) * cfz);
        fd[1] = (-((vsk.dim2[pidx] * rj_f3[2]) + (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) +
                (((vsk.dim2[pidx] * rj_f3[0]) - (rm[2] * rt[1])) * cfz);
        fd[2] = (((vsk.dim2[pidx] * rj_f3[1]) - (rm[0] * rt[2])) * cfx) +
                (-((vsk.dim2[pidx] * rj_f3[0]) + (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
        psw.xfrc[parent_atom] += vs_frc[0] - fb[0] - fc[0] - fd[0];
        psw.yfrc[parent_atom] += vs_frc[1] - fb[1] - fc[1] - fd[1];
        psw.zfrc[parent_atom] += vs_frc[2] - fb[2] - fc[2] - fd[2];
        psw.xfrc[frame2_atom] += fb[0];
        psw.yfrc[frame2_atom] += fb[1];
        psw.zfrc[frame2_atom] += fb[2];
        psw.xfrc[frame3_atom] += fc[0];
        psw.yfrc[frame3_atom] += fc[1];
        psw.zfrc[frame3_atom] += fc[2];
        psw.xfrc[frame4_atom] += fd[0];
        psw.yfrc[frame4_atom] += fd[1];
        psw.zfrc[frame4_atom] += fd[2];
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
    psw.xfrc[vsite_atom] = 0.0;
    psw.yfrc[vsite_atom] = 0.0;
    psw.zfrc[vsite_atom] = 0.0;
  }
}
  
} // namespace structure
} // namespace omni
