#include <cmath>
#include "Math/vector_ops.h"
#include "local_arrangement.h"
#include "virtual_site_placement.h"

namespace omni {
namespace structure {

using math::project;
using math::crossProduct;
using topology::VirtualSiteKind;

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(double* xcrd, double* ycrd, double* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<double> &vsk) {
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom  = vsk.frame1_idx[i];
    const int frame2_atom  = vsk.frame2_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[i])) {
    case VirtualSiteKind::FLEX_2:
      {
        double dx = xcrd[frame2_atom] - xcrd[parent_atom];
        double dy = ycrd[frame2_atom] - ycrd[parent_atom];
        double dz = zcrd[frame2_atom] - zcrd[parent_atom];
        imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
        const double p_f2_factor = vsk.dim1[i];
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
        const double p_vs_distance = vsk.dim1[i];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dx * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dy * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dz * invr);
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        double dx2 = xcrd[frame2_atom] - xcrd[parent_atom];
        double dy2 = ycrd[frame2_atom] - ycrd[parent_atom];
        double dz2 = zcrd[frame2_atom] - zcrd[parent_atom];
        double dx3 = xcrd[frame3_atom] - xcrd[parent_atom];
        double dy3 = ycrd[frame3_atom] - ycrd[parent_atom];
        double dz3 = zcrd[frame3_atom] - zcrd[parent_atom];
        imageCoordinates(&dx2, &dy2, &dz2, umat, invu, unit_cell, minimum_image);
        imageCoordinates(&dx3, &dy3, &dz3, umat, invu, unit_cell, minimum_image);
        const double p_f2_factor = vsk.dim1[i];
        const double p_f3_factor = vsk.dim2[i];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx2) + (p_f3_factor * dx3);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy2) + (p_f3_factor * dy3);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz2) + (p_f3_factor * dz3);
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        double dx23 = xcrd[frame3_atom] - xcrd[frame2_atom];
        double dy23 = ycrd[frame3_atom] - ycrd[frame2_atom];
        double dz23 = zcrd[frame3_atom] - zcrd[frame2_atom];
        imageCoordinates(&dx23, &dy23, &dz23, umat, invu, unit_cell, minimum_image);
        const double f2_f3_factor  = vsk.dim2[i];
        const double x_midpoint = xcrd[frame2_atom] + (f2_f3_factor * dx23);
        const double y_midpoint = ycrd[frame2_atom] + (f2_f3_factor * dx23);
        const double z_midpoint = zcrd[frame2_atom] + (f2_f3_factor * dx23);
        double dxm = x_midpoint - xcrd[parent_atom];
        double dym = y_midpoint - ycrd[parent_atom];
        double dzm = z_midpoint - zcrd[parent_atom];
        imageCoordinates(&dxm, &dym, &dzm, umat, invu, unit_cell, minimum_image);
        const double invr = 1.0 / sqrt((dxm * dxm) + (dym * dym) + (dzm * dzm));
        const double p_vs_distance = vsk.dim1[i];
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dxm * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dym * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dzm * invr);
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        // Use small vectors to make use of a vector operation
        const int frame3_atom  = vsk.frame3_idx[i];
        double p_f2[3], p_f3[3], f2_f3[3], f23_t_pf2[3];
        p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
        p_f3[0]  = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f3[1]  = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f3[2]  = zcrd[frame3_atom] - zcrd[parent_atom];
        f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
        f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
        f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
        imageCoordinates( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates( &p_f3[0],  &p_f3[1],  &p_f3[2], umat, invu, unit_cell, minimum_image);
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
        const double p_f2_factor = vsk.dim1[i] * cos(vsk.dim2[i]) * invr_p_f2;
        const double t_factor    = vsk.dim1[i] * sin(vsk.dim2[i]) * invr_t;
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
        const double pf2_factor = vsk.dim1[i];
        const double pf3_factor = vsk.dim2[i];
        const double cr_factor = vsk.dim3[i];
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
        p_f4[0] = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f4[1] = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f4[2] = zcrd[frame3_atom] - zcrd[parent_atom];
        imageCoordinates(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell, minimum_image);
        imageCoordinates(&p_f4[0], &p_f4[1], &p_f4[2], umat, invu, unit_cell, minimum_image);
        const double pf3_factor = vsk.dim1[i];
        const double pf4_factor = vsk.dim2[i];
        pf3_m_pf2[0] = (pf3_factor * p_f3[0]) - p_f2[0];
        pf3_m_pf2[1] = (pf3_factor * p_f3[1]) - p_f2[1];
        pf3_m_pf2[2] = (pf3_factor * p_f3[2]) - p_f2[2];
        pf4_m_pf2[0] = (pf4_factor * p_f4[0]) - p_f2[0];
        pf4_m_pf2[1] = (pf4_factor * p_f4[1]) - p_f2[1];
        pf4_m_pf2[2] = (pf4_factor * p_f4[2]) - p_f2[2];
        crossProduct(pf3_m_pf2, pf4_m_pf2, p_vs);
        const double pvs_factor = vsk.dim3[i] / sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
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
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag) {
  CoordinateFrameWriter cfw = cf->data();
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu, cfw.unit_cell,
                    ag.getDoublePrecisionVirtualSiteKit());
}

} // namespace structure
} // namespace omni
