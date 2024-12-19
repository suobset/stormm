// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
void drawVirtualSite(const int vsite_atom, const int parent_atom, const int95_t xdisp,
                     const int95_t ydisp, const int95_t zdisp, Tcoord *xcrd, Tcoord *ycrd,
                     Tcoord *zcrd, int* xcrd_ovrf, int* ycrd_ovrf, int* zcrd_ovrf) {
  const int95_t vsite_nx = hostSplitFPSum(xdisp, xcrd[parent_atom], xcrd_ovrf[parent_atom]);
  const int95_t vsite_ny = hostSplitFPSum(ydisp, ycrd[parent_atom], ycrd_ovrf[parent_atom]);
  const int95_t vsite_nz = hostSplitFPSum(zdisp, zcrd[parent_atom], zcrd_ovrf[parent_atom]);
  xcrd[vsite_atom] = vsite_nx.x;
  ycrd[vsite_atom] = vsite_ny.x;
  zcrd[vsite_atom] = vsite_nz.x;
  xcrd_ovrf[vsite_atom] = vsite_nx.y;
  ycrd_ovrf[vsite_atom] = vsite_ny.y;
  zcrd_ovrf[vsite_atom] = vsite_nz.y;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void placeVirtualSite(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                      const double* invu, const UnitCellType unit_cell, const int vsite_atom,
                      const int parent_atom, const int frame2_atom, const int frame3_atom,
                      const int frame4_atom, const VirtualSiteKind frame_type,
                      const Tcalc frame_d1, const Tcalc frame_d2, const Tcalc frame_d3,
                      const Tcalc gpos_scale_factor, int* xcrd_ovrf, int* ycrd_ovrf,
                      int* zcrd_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_factor = value_one / gpos_scale_factor;
  const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
  const bool crd_overflow = (xcrd_ovrf != nullptr || ycrd_ovrf != nullptr || zcrd_ovrf != nullptr);
  switch (frame_type) {
  case VirtualSiteKind::FLEX_2:
    {
      Tcalc dx, dy, dz;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          dx = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          dy = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          dz = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          dx = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        dx = xcrd[frame2_atom] - xcrd[parent_atom];
        dy = ycrd[frame2_atom] - ycrd[parent_atom];
        dz = zcrd[frame2_atom] - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
      const Tcalc p_f2_factor = frame_d1;
      if (tcoord_is_sgnint) {
        const Tcalc disp_mult = p_f2_factor * gpos_scale_factor;
        if (crd_overflow) {
          const int95_t i_dx = hostDoubleToInt95(dx * disp_mult);
          const int95_t i_dy = hostDoubleToInt95(dy * disp_mult);
          const int95_t i_dz = hostDoubleToInt95(dz * disp_mult);
          drawVirtualSite(vsite_atom, parent_atom, i_dx, i_dy, i_dz, xcrd, ycrd, zcrd, xcrd_ovrf,
                          ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(dx * disp_mult);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(dy * disp_mult);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(dz * disp_mult);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz);
      }
    }
    break;
  case VirtualSiteKind::FIXED_2:
    {
      Tcalc dx, dy, dz;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          dx = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          dy = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          dz = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          dx = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        dx = xcrd[frame2_atom] - xcrd[parent_atom];
        dy = ycrd[frame2_atom] - ycrd[parent_atom];
        dz = zcrd[frame2_atom] - zcrd[parent_atom];          
      }
      imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
      const Tcalc invr = value_one / sqrt((dx * dx) + (dy * dy) + (dz * dz));
      const Tcalc p_vs_distance = frame_d1;
      if (tcoord_is_sgnint) {
        const Tcalc disp_mult = p_vs_distance * invr * gpos_scale_factor;
        if (crd_overflow) {
          const int95_t i_dx = hostDoubleToInt95(dx * disp_mult);
          const int95_t i_dy = hostDoubleToInt95(dy * disp_mult);
          const int95_t i_dz = hostDoubleToInt95(dz * disp_mult);
          drawVirtualSite(vsite_atom, parent_atom, i_dx, i_dy, i_dz, xcrd, ycrd, zcrd, xcrd_ovrf,
                          ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * dx);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * dy);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * dz);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dx * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dy * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dz * invr);
      }
    }
    break;
  case VirtualSiteKind::FLEX_3:
    {
      Tcalc dx2, dy2, dz2, dx3, dy3, dz3;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          dx2 = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          dy2 = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          dz2 = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          dx3 = displacement(parent_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          dy3 = displacement(parent_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          dz3 = displacement(parent_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          dx2 = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy2 = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz2 = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          dx3 = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy3 = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz3 = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        dx2 = xcrd[frame2_atom] - xcrd[parent_atom];
        dy2 = ycrd[frame2_atom] - ycrd[parent_atom];
        dz2 = zcrd[frame2_atom] - zcrd[parent_atom];
        dx3 = xcrd[frame3_atom] - xcrd[parent_atom];
        dy3 = ycrd[frame3_atom] - ycrd[parent_atom];
        dz3 = zcrd[frame3_atom] - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&dx2, &dy2, &dz2, umat, invu, unit_cell, minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&dx3, &dy3, &dz3, umat, invu, unit_cell, minimum_image);
      const Tcalc p_f2_factor = frame_d1;
      const Tcalc p_f3_factor = frame_d2;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          const int95_t p_vsx = hostDoubleToInt95(((p_f2_factor * dx2) + (p_f3_factor * dx3)) *
                                                  gpos_scale_factor);
          const int95_t p_vsy = hostDoubleToInt95(((p_f2_factor * dy2) + (p_f3_factor * dy3)) *
                                                  gpos_scale_factor);
          const int95_t p_vsz = hostDoubleToInt95(((p_f2_factor * dz2) + (p_f3_factor * dz3)) *
                                                  gpos_scale_factor);
          drawVirtualSite(vsite_atom, parent_atom, p_vsx, p_vsy, p_vsz, xcrd, ycrd, zcrd,
                          xcrd_ovrf, ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((p_f2_factor * dx2) + (p_f3_factor * dx3)) *
                                     gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((p_f2_factor * dy2) + (p_f3_factor * dy3)) *
                                     gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((p_f2_factor * dz2) + (p_f3_factor * dz3)) *
                                     gpos_scale_factor);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx2) + (p_f3_factor * dx3);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy2) + (p_f3_factor * dy3);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz2) + (p_f3_factor * dz3);
      }
    }
    break;
  case VirtualSiteKind::FIXED_3:
    {
      Tcalc dx23, dy23, dz23;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          dx23 = displacement(frame2_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          dy23 = displacement(frame2_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          dz23 = displacement(frame2_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          dx23 = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          dy23 = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          dz23 = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
      }
      else {
        dx23 = xcrd[frame3_atom] - xcrd[frame2_atom];
        dy23 = ycrd[frame3_atom] - ycrd[frame2_atom];
        dz23 = zcrd[frame3_atom] - zcrd[frame2_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&dx23, &dy23, &dz23, umat, invu, unit_cell, minimum_image);
      const Tcalc f2_f3_factor  = frame_d2;
      Tcoord x_midpoint, y_midpoint, z_midpoint;
      int95_t x95_midpoint, y95_midpoint, z95_midpoint;
      if (tcoord_is_sgnint) {
        const Tcalc disp_mult = f2_f3_factor * gpos_scale_factor;
        if (crd_overflow) {
          x95_midpoint = hostInt95Sum(xcrd[frame2_atom], xcrd_ovrf[frame2_atom], disp_mult * dx23);
          y95_midpoint = hostInt95Sum(ycrd[frame2_atom], ycrd_ovrf[frame2_atom], disp_mult * dy23);
          z95_midpoint = hostInt95Sum(zcrd[frame2_atom], zcrd_ovrf[frame2_atom], disp_mult * dz23);
        }
        else {
          x_midpoint = xcrd[frame2_atom] + llround(disp_mult * dx23);
          y_midpoint = ycrd[frame2_atom] + llround(disp_mult * dy23);
          z_midpoint = zcrd[frame2_atom] + llround(disp_mult * dz23);
        }
      }
      else {
        x_midpoint = xcrd[frame2_atom] + (f2_f3_factor * dx23);
        y_midpoint = ycrd[frame2_atom] + (f2_f3_factor * dy23);
        z_midpoint = zcrd[frame2_atom] + (f2_f3_factor * dz23);
      }
      Tcalc dxm, dym, dzm;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          const int95_t i_dxm = hostSplitFPSubtract(x95_midpoint,
                                                    xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t i_dym = hostSplitFPSubtract(y95_midpoint,
                                                    ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t i_dzm = hostSplitFPSubtract(z95_midpoint,
                                                    zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          dxm = hostSplitFPToReal(i_dxm) * inv_gpos_factor;
          dym = hostSplitFPToReal(i_dym) * inv_gpos_factor;
          dzm = hostSplitFPToReal(i_dzm) * inv_gpos_factor;
        }
        else {
          dxm = static_cast<Tcalc>(x_midpoint - xcrd[parent_atom]) * inv_gpos_factor;
          dym = static_cast<Tcalc>(y_midpoint - ycrd[parent_atom]) * inv_gpos_factor;
          dzm = static_cast<Tcalc>(z_midpoint - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        dxm = x_midpoint - xcrd[parent_atom];
        dym = y_midpoint - ycrd[parent_atom];
        dzm = z_midpoint - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&dxm, &dym, &dzm, umat, invu, unit_cell, minimum_image);
      const Tcalc invr = value_one / sqrt((dxm * dxm) + (dym * dym) + (dzm * dzm));
      const Tcalc p_vs_distance = frame_d1;
      if (tcoord_is_sgnint) {
        const Tcalc disp_mult = p_vs_distance * invr * gpos_scale_factor;
        if (crd_overflow) {
          const int95_t p_vsx = hostDoubleToInt95(disp_mult * dxm);
          const int95_t p_vsy = hostDoubleToInt95(disp_mult * dym);
          const int95_t p_vsz = hostDoubleToInt95(disp_mult * dzm);
          drawVirtualSite(vsite_atom, parent_atom, p_vsx, p_vsy, p_vsz, xcrd, ycrd, zcrd,
                          xcrd_ovrf, ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * dxm);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * dym);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * dzm);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dxm * invr);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dym * invr);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dzm * invr);
      }
    }
    break;
  case VirtualSiteKind::FAD_3:
    {
      // Use small vectors to make use of a vector operation
      Tcalc p_f2[3], f2_f3[3], f23_t_pf2[3];
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          f2_f3[0] = displacement(frame2_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          f2_f3[1] = displacement(frame2_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          f2_f3[2] = displacement(frame2_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          p_f2[0]  = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1]  = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2]  = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
      }
      else {
        p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
        f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
        f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
        f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
      }
      imageCoordinates<Tcalc, Tcalc>( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                     minimum_image);
      const Tcalc invr2_p_f2 = value_one / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
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
      const Tcalc invr_p_f2 = (tcalc_is_double) ?
                              value_one / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                               (p_f2[2] * p_f2[2])) :
                              value_one / sqrtf((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                                (p_f2[2] * p_f2[2]));
          
      const Tcalc invr_t = (tcalc_is_double) ?
                           value_one / sqrt((f23_t_pf2[0] * f23_t_pf2[0]) +
                                            (f23_t_pf2[1] * f23_t_pf2[1]) +
                                            (f23_t_pf2[2] * f23_t_pf2[2])) :
                           value_one / sqrtf((f23_t_pf2[0] * f23_t_pf2[0]) +
                                             (f23_t_pf2[1] * f23_t_pf2[1]) +
                                             (f23_t_pf2[2] * f23_t_pf2[2]));          
      const Tcalc p_f2_factor = frame_d1 * cos(frame_d2) * invr_p_f2;
      const Tcalc t_factor    = frame_d1 * sin(frame_d2) * invr_t;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          const int95_t p_vsx = hostDoubleToInt95(((p_f2_factor * p_f2[0]) +
                                                   (t_factor * f23_t_pf2[0])) *
                                                  gpos_scale_factor);
          const int95_t p_vsy = hostDoubleToInt95(((p_f2_factor * p_f2[1]) +
                                                   (t_factor * f23_t_pf2[1])) *
                                                  gpos_scale_factor);
          const int95_t p_vsz = hostDoubleToInt95(((p_f2_factor * p_f2[2]) +
                                                   (t_factor * f23_t_pf2[2])) *
                                                  gpos_scale_factor);
          drawVirtualSite(vsite_atom, parent_atom, p_vsx, p_vsy, p_vsz, xcrd, ycrd, zcrd,
                          xcrd_ovrf, ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[0]) + (t_factor * f23_t_pf2[0])) *
                                     gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[1]) + (t_factor * f23_t_pf2[1])) *
                                     gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[2]) + (t_factor * f23_t_pf2[2])) *
                                     gpos_scale_factor);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * p_f2[0]) + (t_factor * f23_t_pf2[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * p_f2[1]) + (t_factor * f23_t_pf2[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * p_f2[2]) + (t_factor * f23_t_pf2[2]);
      }
    }
    break;
  case VirtualSiteKind::OUT_3:
    {
      Tcalc p_f2[3], p_f3[3], f2_f3[3], pf2_x_pf3[3];
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          p_f3[0] = displacement(parent_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f3[1] = displacement(parent_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f3[2] = displacement(parent_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell,
                                     minimum_image);
      crossProduct(p_f2, p_f3, pf2_x_pf3);
      const Tcalc pf2_factor = frame_d1;
      const Tcalc pf3_factor = frame_d2;
      const Tcalc cr_factor  = frame_d3;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          const int95_t p_vsx = hostDoubleToInt95(((pf2_factor * p_f2[0]) +
                                                   (pf3_factor * p_f3[0]) +
                                                   (cr_factor * pf2_x_pf3[0])) *
                                                  gpos_scale_factor);
          const int95_t p_vsy = hostDoubleToInt95(((pf2_factor * p_f2[1]) +
                                                   (pf3_factor * p_f3[1]) +
                                                   (cr_factor * pf2_x_pf3[1])) *
                                                  gpos_scale_factor);
          const int95_t p_vsz = hostDoubleToInt95(((pf2_factor * p_f2[2]) +
                                                   (pf3_factor * p_f3[2]) +
                                                   (cr_factor * pf2_x_pf3[2])) *
                                                  gpos_scale_factor);
          drawVirtualSite(vsite_atom, parent_atom, p_vsx, p_vsy, p_vsz, xcrd, ycrd, zcrd,
                          xcrd_ovrf, ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((pf2_factor * p_f2[0]) + (pf3_factor * p_f3[0]) +
                                     (cr_factor * pf2_x_pf3[0])) * gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((pf2_factor * p_f2[1]) + (pf3_factor * p_f3[1]) +
                                     (cr_factor * pf2_x_pf3[1])) * gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((pf2_factor * p_f2[2]) + (pf3_factor * p_f3[2]) +
                                     (cr_factor * pf2_x_pf3[2])) * gpos_scale_factor);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (pf2_factor * p_f2[0]) + (pf3_factor * p_f3[0]) +
                           (cr_factor * pf2_x_pf3[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (pf2_factor * p_f2[1]) + (pf3_factor * p_f3[1]) +
                           (cr_factor * pf2_x_pf3[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (pf2_factor * p_f2[2]) + (pf3_factor * p_f3[2]) +
                           (cr_factor * pf2_x_pf3[2]);
      }
    }
    break;
  case VirtualSiteKind::FIXED_4:
    {
      Tcalc p_f2[3], p_f3[3], p_f4[3], pf3_m_pf2[3], pf4_m_pf2[3], p_vs[3];
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          p_f3[0] = displacement(parent_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f3[1] = displacement(parent_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f3[2] = displacement(parent_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          p_f4[0] = displacement(parent_atom, frame4_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f4[1] = displacement(parent_atom, frame4_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f4[2] = displacement(parent_atom, frame4_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f4[0] = static_cast<Tcalc>(xcrd[frame4_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f4[1] = static_cast<Tcalc>(ycrd[frame4_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f4[2] = static_cast<Tcalc>(zcrd[frame4_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
        p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
        p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
        p_f4[0] = xcrd[frame4_atom] - xcrd[parent_atom];
        p_f4[1] = ycrd[frame4_atom] - ycrd[parent_atom];
        p_f4[2] = zcrd[frame4_atom] - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_f4[0], &p_f4[1], &p_f4[2], umat, invu, unit_cell,
                                     minimum_image);
      const Tcalc pf3_factor = frame_d1;
      const Tcalc pf4_factor = frame_d2;
      pf3_m_pf2[0] = (pf3_factor * p_f3[0]) - p_f2[0];
      pf3_m_pf2[1] = (pf3_factor * p_f3[1]) - p_f2[1];
      pf3_m_pf2[2] = (pf3_factor * p_f3[2]) - p_f2[2];
      pf4_m_pf2[0] = (pf4_factor * p_f4[0]) - p_f2[0];
      pf4_m_pf2[1] = (pf4_factor * p_f4[1]) - p_f2[1];
      pf4_m_pf2[2] = (pf4_factor * p_f4[2]) - p_f2[2];
      crossProduct(pf3_m_pf2, pf4_m_pf2, p_vs);
      const Tcalc pvs_factor = (tcalc_is_double) ?
                               frame_d3 / sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                               (p_vs[2] * p_vs[2])) :
                               frame_d3 / sqrtf((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                                (p_vs[2] * p_vs[2]));
      if (tcoord_is_sgnint) {
        const Tcalc disp_mult = pvs_factor * gpos_scale_factor;
        if (crd_overflow) {
          const int95_t p_vsx = hostDoubleToInt95(disp_mult * p_vs[0]);
          const int95_t p_vsy = hostDoubleToInt95(disp_mult * p_vs[1]);
          const int95_t p_vsz = hostDoubleToInt95(disp_mult * p_vs[2]);
          drawVirtualSite(vsite_atom, parent_atom, p_vsx, p_vsy, p_vsz, xcrd, ycrd, zcrd,
                          xcrd_ovrf, ycrd_ovrf, zcrd_ovrf);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * p_vs[0]);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * p_vs[1]);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * p_vs[2]);
        }
      }
      else {
        xcrd[vsite_atom] = xcrd[parent_atom] + (pvs_factor * p_vs[0]);
        ycrd[vsite_atom] = ycrd[parent_atom] + (pvs_factor * p_vs[1]);
        zcrd[vsite_atom] = zcrd[parent_atom] + (pvs_factor * p_vs[2]);
      }
    }
    break;
  case VirtualSiteKind::NONE:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void placeVirtualSites(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc gpos_scale_factor) {
  for (int i = 0; i < vsk.nsite; i++) {
    const int param_idx = vsk.vs_param_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
    placeVirtualSite(xcrd, ycrd, zcrd, umat, invu, unit_cell, vsk.vs_atoms[i], vsk.frame1_idx[i],
                     vsk.frame2_idx[i], vsk.frame3_idx[i], vsk.frame4_idx[i],
                     static_cast<VirtualSiteKind>(vsk.vs_types[param_idx]), vsk.dim1[param_idx],
                     vsk.dim2[param_idx], vsk.dim3[param_idx], gpos_scale_factor);
  }
}
    
//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void placeVirtualSites(PhaseSpaceWriter *psw, const VirtualSiteKit<Tcalc> vsk,
                       const bool edit_alternate_image) {
  if (edit_alternate_image) {
    placeVirtualSites(psw->xalt, psw->yalt, psw->zalt, psw->umat_alt, psw->invu_alt,
                      psw->unit_cell, vsk);
  }
  else {
    placeVirtualSites(psw->xcrd, psw->ycrd, psw->zcrd, psw->umat, psw->invu, psw->unit_cell, vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void placeVirtualSites(CoordinateFrameWriter *cfw, const VirtualSiteKit<Tcalc> vsk) {
  placeVirtualSites(cfw->xcrd, cfw->ycrd, cfw->zcrd, cfw->umat, cfw->invu, cfw->unit_cell, vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void placeVirtualSites(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                       const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk) {

  // Allocate space to mock the GPU's __shared__ memory buffers
  std::vector<int2> vwu_map(vwu_abstract_length);
  std::vector<llint> sh_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xcrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_ycrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zcrd_ovrf(maximum_valence_work_unit_atoms);

  // Set pointers with regard to the detail present in the synthesis
  llint* xcrd_ptr = sh_xcrd.data();
  llint* ycrd_ptr = sh_ycrd.data();
  llint* zcrd_ptr = sh_zcrd.data();
  int *xcrd_ovrf_ptr, *ycrd_ovrf_ptr, *zcrd_ovrf_ptr;
  if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
    xcrd_ovrf_ptr = sh_xcrd_ovrf.data();
    ycrd_ovrf_ptr = sh_ycrd_ovrf.data();
    zcrd_ovrf_ptr = sh_zcrd_ovrf.data();
  }
  else {
    xcrd_ovrf_ptr = nullptr;
    ycrd_ovrf_ptr = nullptr;
    zcrd_ovrf_ptr = nullptr;
  }

  // Determine the transform stride.
  const int xfrm_stride = roundUp(9, warp_size_int);

  // Loop over all valence work units.
  for (int i = 0; i < poly_vk.nvwu; i++) {

    // Extract the valence work unit's abstract.
    for (int j = 0; j < vwu_abstract_length; j++) {
      vwu_map[j] = poly_vk.vwu_abstracts[(i * vwu_abstract_length) + j];
    }

    // Determine the system index.
    const int sys_idx = vwu_map[static_cast<size_t>(VwuAbstractMap::SYSTEM_ID)].x;

    // Import the atoms of the work unit.
    const int2 import_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::IMPORT)];
    for (int j = import_limits.x; j < import_limits.y; j++) {
      const size_t synth_atom = poly_vk.vwu_imports[j];
      const size_t local_atom = j - import_limits.x;
      sh_xcrd[local_atom] = poly_psw->xalt[synth_atom];
      sh_ycrd[local_atom] = poly_psw->yalt[synth_atom];
      sh_zcrd[local_atom] = poly_psw->zalt[synth_atom];
      if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        sh_xcrd_ovrf[local_atom] = poly_psw->xalt_ovrf[synth_atom];
        sh_ycrd_ovrf[local_atom] = poly_psw->yalt_ovrf[synth_atom];
        sh_zcrd_ovrf[local_atom] = poly_psw->zalt_ovrf[synth_atom];
      }
    }
    double* umat = &poly_psw->umat[xfrm_stride * sys_idx];
    double* invu = &poly_psw->invu[xfrm_stride * sys_idx];

    // Perform virtual site placement
    const int2 vs_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::VSITE)];
    for (int j = vs_limits.x; j < vs_limits.y; j++) {
      const uint2 t_insr = poly_auk.vste_insr[j];
      const int vsite_atom  = (t_insr.x & 0x3ff);
      const int parent_atom = ((t_insr.x >> 10) & 0x3ff);
      const int frame2_atom = ((t_insr.x >> 20) & 0x3ff);
      const int frame3_atom = (t_insr.y & 0x3ff);
      const int frame4_atom = ((t_insr.y >> 10) & 0x3ff);
      const uint param_idx = (t_insr.y >> 20);
      const Tcalc4 vs_frame = poly_auk.vs_params[param_idx];
      const VirtualSiteKind frame_type = static_cast<VirtualSiteKind>(vs_frame.w);
      placeVirtualSite<llint, Tcalc>(xcrd_ptr, ycrd_ptr, zcrd_ptr, umat, invu, poly_psw->unit_cell,
                                     vsite_atom, parent_atom, frame2_atom, frame3_atom,
                                     frame4_atom, frame_type, vs_frame.x, vs_frame.y, vs_frame.z,
                                     poly_psw->gpos_scale, xcrd_ovrf_ptr, ycrd_ovrf_ptr,
                                     zcrd_ovrf_ptr);

      // The updated virtual site can be written back to the global position arrays immediately.
      const int manip_segment = (vsite_atom >> 5);
      const int manip_bitpos = (vsite_atom & 0x1f);
      const int2 manip_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::MANIPULATE)];
      const uint2 manip_mask = poly_auk.vwu_manip[manip_limits.x + manip_segment];
      if ((manip_mask.y >> manip_bitpos) & 0x1) {
        const size_t synth_atom = poly_vk.vwu_imports[import_limits.x + vsite_atom];
        poly_psw->xalt[synth_atom] = sh_xcrd[vsite_atom];
        poly_psw->yalt[synth_atom] = sh_ycrd[vsite_atom];
        poly_psw->zalt[synth_atom] = sh_zcrd[vsite_atom];
        if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
          poly_psw->xalt_ovrf[synth_atom] = sh_xcrd_ovrf[vsite_atom];
          poly_psw->yalt_ovrf[synth_atom] = sh_ycrd_ovrf[vsite_atom];
          poly_psw->zalt_ovrf[synth_atom] = sh_zcrd_ovrf[vsite_atom];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, const UnitCellType unit_cell,
                               const int vsite_atom, const int parent_atom, const int frame2_atom,
                               const int frame3_atom, const int frame4_atom,
                               const VirtualSiteKind frame_type, const Tcalc frame_d1,
                               const Tcalc frame_d2, const Tcalc frame_d3,
                               const Tcalc gpos_scale_factor, const Tcalc force_scale_factor,
                               const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                               int* xfrc_ovrf, int* yfrc_ovrf, int* zfrc_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  const bool crd_overflow = (xcrd_ovrf != nullptr || ycrd_ovrf != nullptr || zcrd_ovrf != nullptr);
  const bool frc_overflow = (xfrc_ovrf != nullptr || yfrc_ovrf != nullptr || zfrc_ovrf != nullptr);
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_factor = value_one / gpos_scale_factor;
  const Tcalc inv_force_factor = value_one / force_scale_factor;
  const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
  
  // Copy the force on the virtual site into its own vector.  This will be needed in most
  // subsequent force transmissions as a real-valued representation of the forces (scaling
  // factor for fixed precision removed).
  Tcalc vs_frc[3];
  if (tforce_is_sgnint) {
    if (frc_overflow) {
      vs_frc[0] = hostInt95ToDouble(xfrc[vsite_atom], xfrc_ovrf[vsite_atom]) * inv_force_factor;
      vs_frc[1] = hostInt95ToDouble(yfrc[vsite_atom], yfrc_ovrf[vsite_atom]) * inv_force_factor;
      vs_frc[2] = hostInt95ToDouble(zfrc[vsite_atom], zfrc_ovrf[vsite_atom]) * inv_force_factor;
    }
    else {
      vs_frc[0] = static_cast<Tcalc>(xfrc[vsite_atom]) * inv_force_factor;
      vs_frc[1] = static_cast<Tcalc>(yfrc[vsite_atom]) * inv_force_factor;
      vs_frc[2] = static_cast<Tcalc>(zfrc[vsite_atom]) * inv_force_factor;
    }
  }
  else {
    vs_frc[0] = xfrc[vsite_atom];
    vs_frc[1] = yfrc[vsite_atom];
    vs_frc[2] = zfrc[vsite_atom];
  }
  switch (frame_type) {
  case VirtualSiteKind::FLEX_2:
    {
      const Tcalc p_f2_factor = frame_d1;
      if (tforce_is_sgnint) {
        const Tcalc part_mult = (value_one - p_f2_factor);
        if (frc_overflow) {
          const Tcalc vs_fx = hostInt95ToDouble(xfrc[vsite_atom], xfrc_ovrf[vsite_atom]);
          const Tcalc vs_fy = hostInt95ToDouble(yfrc[vsite_atom], yfrc_ovrf[vsite_atom]);
          const Tcalc vs_fz = hostInt95ToDouble(zfrc[vsite_atom], zfrc_ovrf[vsite_atom]);
          const int95_t xpart = hostDoubleToInt95(part_mult * vs_fx);
          const int95_t ypart = hostDoubleToInt95(part_mult * vs_fy);
          const int95_t zpart = hostDoubleToInt95(part_mult * vs_fz);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], xpart);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], ypart);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], zpart);
          const int95_t ivsr_fx = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                                    xpart.x, xpart.y);
          const int95_t ivsr_fy = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                                    ypart.x, ypart.y);
          const int95_t ivsr_fz = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                                    zpart.x, zpart.y);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], ivsr_fx);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], ivsr_fy);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], ivsr_fz);
        }
        else {
          const Tforce xpart = llround(part_mult * static_cast<Tcalc>(xfrc[vsite_atom]));
          const Tforce ypart = llround(part_mult * static_cast<Tcalc>(yfrc[vsite_atom]));
          const Tforce zpart = llround(part_mult * static_cast<Tcalc>(zfrc[vsite_atom]));
          xfrc[parent_atom] += xpart;
          yfrc[parent_atom] += ypart;
          zfrc[parent_atom] += zpart;
          xfrc[frame2_atom] += xfrc[vsite_atom] - xpart;
          yfrc[frame2_atom] += yfrc[vsite_atom] - ypart;
          zfrc[frame2_atom] += zfrc[vsite_atom] - zpart;
        }
      }
      else {
        xfrc[parent_atom] += (1.0 - p_f2_factor) * xfrc[vsite_atom];
        yfrc[parent_atom] += (1.0 - p_f2_factor) * yfrc[vsite_atom];
        zfrc[parent_atom] += (1.0 - p_f2_factor) * zfrc[vsite_atom];
        xfrc[frame2_atom] += p_f2_factor * xfrc[vsite_atom];
        yfrc[frame2_atom] += p_f2_factor * yfrc[vsite_atom];
        zfrc[frame2_atom] += p_f2_factor * zfrc[vsite_atom];
      }
    }
    break;
  case VirtualSiteKind::FIXED_2:
    {
      Tcalc p_f2[3], p_vs[3], vs_frc_proj[3], force_partition[3];
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          p_vs[0] = displacement(parent_atom, vsite_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_vs[1] = displacement(parent_atom, vsite_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_vs[2] = displacement(parent_atom, vsite_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_vs[0] = static_cast<Tcalc>(xcrd[vsite_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_vs[1] = static_cast<Tcalc>(ycrd[vsite_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_vs[2] = static_cast<Tcalc>(zcrd[vsite_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
      }
      else {
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        p_vs[0] = xcrd[vsite_atom] - xcrd[parent_atom];
        p_vs[1] = ycrd[vsite_atom] - ycrd[parent_atom];
        p_vs[2] = zcrd[vsite_atom] - zcrd[parent_atom];
      }
      imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_vs[0], &p_vs[1], &p_vs[2], umat, invu, unit_cell,
                                     minimum_image);
      const Tcalc invr_p_f2 = (tcalc_is_double) ?
                              value_one / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                               (p_f2[2] * p_f2[2])) :
                              value_one / sqrtf((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                                (p_f2[2] * p_f2[2]));

      // Take the projection of the force on the virtual site onto the parent atom -> virtual
      // site displacement vector p_vs.  The force is partitioned between the two frame atoms
      // according to the distance between the parent atom and the virtual site (normalized by
      // the distance between the two frame atoms).
      project(vs_frc, p_vs, vs_frc_proj, 3);
      const Tcalc p_vs_distance = frame_d1;
      for (int j = 0; j < 3; j++) {
        force_partition[j] = invr_p_f2 * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
      }
      if (tforce_is_sgnint) {
        if (frc_overflow) {
          const int95_t xpart = hostDoubleToInt95(force_partition[0] * force_scale_factor);
          const int95_t ypart = hostDoubleToInt95(force_partition[1] * force_scale_factor);
          const int95_t zpart = hostDoubleToInt95(force_partition[2] * force_scale_factor);
          const int95_t ivsr_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                                   xpart.x, xpart.y);
          const int95_t ivsr_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                                   ypart.x, ypart.y);
          const int95_t ivsr_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                                   zpart.x, zpart.y);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], ivsr_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], ivsr_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], ivsr_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], xpart);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], ypart);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], zpart);
        }
        else {
          const Tforce xpart = llround(force_partition[0] * force_scale_factor);
          const Tforce ypart = llround(force_partition[1] * force_scale_factor);
          const Tforce zpart = llround(force_partition[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - zpart;
          xfrc[frame2_atom] += xpart;
          yfrc[frame2_atom] += ypart;
          zfrc[frame2_atom] += zpart;
        }
      }
      else {
        xfrc[parent_atom] += vs_frc[0] - force_partition[0];
        yfrc[parent_atom] += vs_frc[1] - force_partition[1];
        zfrc[parent_atom] += vs_frc[2] - force_partition[2];
        xfrc[frame2_atom] += force_partition[0];
        yfrc[frame2_atom] += force_partition[1];
        zfrc[frame2_atom] += force_partition[2];
      }
    }
    break;
  case VirtualSiteKind::FLEX_3:
    {
      if (tforce_is_sgnint) {
        const Tcalc p_f2_factor = frame_d1 * force_scale_factor;
        const Tcalc p_f3_factor = frame_d2 * force_scale_factor;
        if (frc_overflow) {
          const int95_t f2x_part = hostDoubleToInt95(p_f2_factor * vs_frc[0]);
          const int95_t f2y_part = hostDoubleToInt95(p_f2_factor * vs_frc[1]);
          const int95_t f2z_part = hostDoubleToInt95(p_f2_factor * vs_frc[2]);
          const int95_t f3x_part = hostDoubleToInt95(p_f3_factor * vs_frc[0]);
          const int95_t f3y_part = hostDoubleToInt95(p_f3_factor * vs_frc[1]);
          const int95_t f3z_part = hostDoubleToInt95(p_f3_factor * vs_frc[2]);
          const int95_t vs_dfx = hostSplitFPSum(f2x_part, f3x_part);
          const int95_t vs_dfy = hostSplitFPSum(f2y_part, f3y_part);
          const int95_t vs_dfz = hostSplitFPSum(f2z_part, f3z_part);
          const int95_t ivsr_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                                   vs_dfx.x, vs_dfx.y);
          const int95_t ivsr_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                                   vs_dfy.x, vs_dfy.y);
          const int95_t ivsr_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                                   vs_dfz.x, vs_dfz.y);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], ivsr_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], ivsr_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], ivsr_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], f2x_part);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], f2y_part);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], f2z_part);
          hostSplitFPSum(&xfrc[frame3_atom], &xfrc_ovrf[frame3_atom], f3x_part);
          hostSplitFPSum(&yfrc[frame3_atom], &yfrc_ovrf[frame3_atom], f3y_part);
          hostSplitFPSum(&zfrc[frame3_atom], &zfrc_ovrf[frame3_atom], f3z_part);
        }
        else {
          const Tforce f2x_part = llround(p_f2_factor * vs_frc[0]);
          const Tforce f2y_part = llround(p_f2_factor * vs_frc[1]);
          const Tforce f2z_part = llround(p_f2_factor * vs_frc[2]);
          const Tforce f3x_part = llround(p_f3_factor * vs_frc[0]);
          const Tforce f3y_part = llround(p_f3_factor * vs_frc[1]);
          const Tforce f3z_part = llround(p_f3_factor * vs_frc[2]);
          xfrc[parent_atom] += xfrc[vsite_atom] - f2x_part - f3x_part;
          yfrc[parent_atom] += yfrc[vsite_atom] - f2y_part - f3y_part;
          zfrc[parent_atom] += zfrc[vsite_atom] - f2z_part - f3z_part;
          xfrc[frame2_atom] += f2x_part;
          yfrc[frame2_atom] += f2y_part;
          zfrc[frame2_atom] += f2z_part;
          xfrc[frame3_atom] += f3x_part;
          yfrc[frame3_atom] += f3y_part;
          zfrc[frame3_atom] += f3z_part;
        }
      }
      else {
        const Tcalc p_f2_factor = frame_d1;
        const Tcalc p_f3_factor = frame_d2;
        xfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * xfrc[vsite_atom];
        yfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * yfrc[vsite_atom];
        zfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * zfrc[vsite_atom];
        xfrc[frame2_atom] += p_f2_factor * xfrc[vsite_atom];
        yfrc[frame2_atom] += p_f2_factor * yfrc[vsite_atom];
        zfrc[frame2_atom] += p_f2_factor * zfrc[vsite_atom];
        xfrc[frame3_atom] += p_f3_factor * xfrc[vsite_atom];
        yfrc[frame3_atom] += p_f3_factor * yfrc[vsite_atom];
        zfrc[frame3_atom] += p_f3_factor * zfrc[vsite_atom];
      }
    }
    break;
  case VirtualSiteKind::FIXED_3:
    {
      Tcalc f2_f3[3], p_vs[3], p_mid[3], vs_frc_proj[3], force_partition[3];
      const Tcalc f2_f3_factor  = frame_d2;
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          f2_f3[0] = displacement(frame2_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          f2_f3[1] = displacement(frame2_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          f2_f3[2] = displacement(frame2_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          p_vs[0] = displacement(parent_atom, vsite_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_vs[1] = displacement(parent_atom, vsite_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_vs[2] = displacement(parent_atom, vsite_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          const int95_t dpf2x = hostInt95Subtract(xcrd[frame2_atom], xcrd_ovrf[frame2_atom],
                                                  xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t dpf2y = hostInt95Subtract(ycrd[frame2_atom], ycrd_ovrf[frame2_atom],
                                                  ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t dpf2z = hostInt95Subtract(zcrd[frame2_atom], zcrd_ovrf[frame2_atom],
                                                  zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          p_mid[0] = static_cast<Tcalc>(hostSplitFPToReal(dpf2x) * inv_gpos_factor) +
                     (f2_f3_factor * f2_f3[0]);
          p_mid[1] = static_cast<Tcalc>(hostSplitFPToReal(dpf2y) * inv_gpos_factor) +
                     (f2_f3_factor * f2_f3[1]);
          p_mid[2] = static_cast<Tcalc>(hostSplitFPToReal(dpf2z) * inv_gpos_factor) +
                     (f2_f3_factor * f2_f3[2]);
        }
        else {
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
          p_vs[0]  = static_cast<Tcalc>(xcrd[vsite_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_vs[1]  = static_cast<Tcalc>(ycrd[vsite_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_vs[2]  = static_cast<Tcalc>(zcrd[vsite_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_mid[0] = (static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[0]);
          p_mid[1] = (static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[1]);
          p_mid[2] = (static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[2]);
        }
      }
      else {
        f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
        f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
        f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
        p_vs[0] = xcrd[vsite_atom] - xcrd[parent_atom];
        p_vs[1] = ycrd[vsite_atom] - ycrd[parent_atom];
        p_vs[2] = zcrd[vsite_atom] - zcrd[parent_atom];
        p_mid[0] = xcrd[frame2_atom] - xcrd[parent_atom] + (f2_f3_factor * f2_f3[0]);
        p_mid[1] = ycrd[frame2_atom] - ycrd[parent_atom] + (f2_f3_factor * f2_f3[1]);
        p_mid[2] = zcrd[frame2_atom] - zcrd[parent_atom] + (f2_f3_factor * f2_f3[2]);
      }
      imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_vs[0], &p_vs[1], &p_vs[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&p_mid[0], &p_mid[1], &p_mid[2], umat, invu, unit_cell,
                                     minimum_image);

      // As with the fixed distance, two-point frame, compute the part of the force on the
      // virtual site that is perpendicular to the parent atom -> virtual site displacement
      // vector and use this to compute the partitioning between the parent atom and some
      // imaginary midpoint on the line between frame atom 2 and frame atom 3.  That's like
      // the FIXED_2 frame type.  The force on the midpoint is subsequently distributed, akin
      // to the FLEX_2 type, between atoms 2 and 3.
      project(vs_frc, p_vs, vs_frc_proj, 3);
      const Tcalc p_vs_distance = frame_d1;
      const Tcalc invr_p_mid = 1.0 / sqrt((p_mid[0] * p_mid[0]) + (p_mid[1] * p_mid[1]) +
                                          (p_mid[2] * p_mid[2]));
      for (int j = 0; j < 3; j++) {
        force_partition[j] = invr_p_mid * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
      }
      if (tforce_is_sgnint) {
        if (frc_overflow) {
          const int95_t xpart = hostDoubleToInt95(force_partition[0] * force_scale_factor);
          const int95_t ypart = hostDoubleToInt95(force_partition[1] * force_scale_factor);
          const int95_t zpart = hostDoubleToInt95(force_partition[2] * force_scale_factor);
          const int95_t f2f3_xprt = hostDoubleToInt95(f2_f3_factor * force_partition[0] *
                                                      force_scale_factor);
          const int95_t f2f3_yprt = hostDoubleToInt95(f2_f3_factor * force_partition[1] *
                                                      force_scale_factor);
          const int95_t f2f3_zprt = hostDoubleToInt95(f2_f3_factor * force_partition[2] *
                                                      force_scale_factor);
          const int95_t dparent_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                                      xpart.x, xpart.y);
          const int95_t dparent_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                                      ypart.x, ypart.y);
          const int95_t dparent_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                                      zpart.x, zpart.y);
          const int95_t dframe2_x = hostSplitFPSubtract(xpart, f2f3_xprt);
          const int95_t dframe2_y = hostSplitFPSubtract(ypart, f2f3_yprt);
          const int95_t dframe2_z = hostSplitFPSubtract(zpart, f2f3_zprt);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], dparent_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], dparent_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], dparent_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], dframe2_x);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], dframe2_y);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], dframe2_z);
          hostSplitFPSum(&xfrc[frame3_atom], &xfrc_ovrf[frame3_atom], f2f3_xprt);
          hostSplitFPSum(&yfrc[frame3_atom], &yfrc_ovrf[frame3_atom], f2f3_yprt);
          hostSplitFPSum(&zfrc[frame3_atom], &zfrc_ovrf[frame3_atom], f2f3_zprt);
        }
        else {
          const Tforce xpart = llround(force_partition[0] * force_scale_factor);
          const Tforce ypart = llround(force_partition[1] * force_scale_factor);
          const Tforce zpart = llround(force_partition[2] * force_scale_factor);
          const Tforce f2f3_xprt = llround(f2_f3_factor * force_partition[0] * force_scale_factor);
          const Tforce f2f3_yprt = llround(f2_f3_factor * force_partition[1] * force_scale_factor);
          const Tforce f2f3_zprt = llround(f2_f3_factor * force_partition[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - zpart;
          xfrc[frame2_atom] += xpart - f2f3_xprt;
          yfrc[frame2_atom] += ypart - f2f3_yprt;
          zfrc[frame2_atom] += zpart - f2f3_zprt;
          xfrc[frame3_atom] += f2f3_xprt;
          yfrc[frame3_atom] += f2f3_yprt;
          zfrc[frame3_atom] += f2f3_zprt;
        }
      }
      else {
        xfrc[parent_atom] += vs_frc[0] - force_partition[0];
        yfrc[parent_atom] += vs_frc[1] - force_partition[1];
        zfrc[parent_atom] += vs_frc[2] - force_partition[2];
        xfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[0];
        yfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[1];
        zfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[2];
        xfrc[frame3_atom] += f2_f3_factor * force_partition[0];
        yfrc[frame3_atom] += f2_f3_factor * force_partition[1];
        zfrc[frame3_atom] += f2_f3_factor * force_partition[2];
      }
    }
    break;
  case VirtualSiteKind::FAD_3:
    {
      Tcalc p_f2[3], f2_f3[3], f23_t_pf2[3], F1[3], F2[3], F3[3];
      if (tcoord_is_sgnint) {
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          f2_f3[0] = displacement(frame2_atom, frame3_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          f2_f3[1] = displacement(frame2_atom, frame3_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          f2_f3[2] = displacement(frame2_atom, frame3_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
        }
        else {
          p_f2[0]  = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1]  = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2]  = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
      }
      else {
        p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
        f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
        f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
        f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
      }
      imageCoordinates<Tcalc, Tcalc>( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell,
                                     minimum_image);
      imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                     minimum_image);
      project(f2_f3, p_f2, f23_t_pf2, 3);
      f23_t_pf2[0] = f2_f3[0] - f23_t_pf2[0];
      f23_t_pf2[1] = f2_f3[1] - f23_t_pf2[1];
      f23_t_pf2[2] = f2_f3[2] - f23_t_pf2[2];
      const Tcalc invr2_p_f2 = value_one / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                            (p_f2[2] * p_f2[2]));
      const Tcalc invr2_t = value_one / ((f23_t_pf2[0] * f23_t_pf2[0]) +
                                         (f23_t_pf2[1] * f23_t_pf2[1]) +
                                         (f23_t_pf2[2] * f23_t_pf2[2]));
      Tcalc invr_p_f2, invr_t, p_f2_factor, t_factor;
      if (tcalc_is_double) {
        invr_p_f2 = sqrt(invr2_p_f2);
        invr_t = sqrt(invr2_t);
        p_f2_factor = frame_d1 * cos(frame_d2) * invr_p_f2;
        t_factor    = frame_d1 * sin(frame_d2) * invr_t;
      }
      else {
        invr_p_f2 = sqrtf(invr2_p_f2);
        invr_t = sqrtf(invr2_t);
        p_f2_factor = frame_d1 * cosf(frame_d2) * invr_p_f2;
        t_factor    = frame_d1 * sinf(frame_d2) * invr_t;
      }
      const Tcalc f1fac  = dot(p_f2, vs_frc, 3) * invr2_p_f2;
      const Tcalc f2fac  = dot(f23_t_pf2, vs_frc, 3) * invr2_t;
      const Tcalc abbcOabab = dot(p_f2, f2_f3, 3) * invr2_p_f2;
      for (int j = 0; j < 3; j++) {
        F1[j] = vs_frc[j] - (f1fac * p_f2[j]);
        F2[j] = F1[j] - (f2fac * f23_t_pf2[j]);
        F3[j] = f1fac * f23_t_pf2[j];          
      }
      if (tforce_is_sgnint) {
        if (frc_overflow) {
          const int95_t f1_xpart = hostDoubleToInt95(((p_f2_factor * F1[0]) -
                                                      (t_factor * ((abbcOabab * F2[0]) + F3[0]))) *
                                                     force_scale_factor);
          const int95_t f1_ypart = hostDoubleToInt95(((p_f2_factor * F1[1]) -
                                                      (t_factor * ((abbcOabab * F2[1]) + F3[1]))) *
                                                     force_scale_factor);
          const int95_t f1_zpart = hostDoubleToInt95(((p_f2_factor * F1[2]) -
                                                      (t_factor * ((abbcOabab * F2[2]) + F3[2]))) *
                                                     force_scale_factor);
          const int95_t f3_xpart = hostDoubleToInt95(t_factor * F2[0] * force_scale_factor);
          const int95_t f3_ypart = hostDoubleToInt95(t_factor * F2[1] * force_scale_factor);
          const int95_t f3_zpart = hostDoubleToInt95(t_factor * F2[2] * force_scale_factor);
          const int95_t dparent_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                                      f1_xpart.x, f1_xpart.y);
          const int95_t dparent_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                                      f1_ypart.x, f1_ypart.y);
          const int95_t dparent_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                                      f1_zpart.x, f1_zpart.y);
          const int95_t dframe2_x = hostSplitFPSubtract(f1_xpart, f3_xpart);
          const int95_t dframe2_y = hostSplitFPSubtract(f1_ypart, f3_ypart);
          const int95_t dframe2_z = hostSplitFPSubtract(f1_zpart, f3_zpart);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], dparent_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], dparent_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], dparent_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], dframe2_x);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], dframe2_y);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], dframe2_z);
          hostSplitFPSum(&xfrc[frame3_atom], &xfrc_ovrf[frame3_atom], f3_xpart);
          hostSplitFPSum(&yfrc[frame3_atom], &yfrc_ovrf[frame3_atom], f3_ypart);
          hostSplitFPSum(&zfrc[frame3_atom], &zfrc_ovrf[frame3_atom], f3_zpart);
        }
        else {
          const Tforce f1_xpart = llround(((p_f2_factor * F1[0]) -
                                           (t_factor * ((abbcOabab * F2[0]) + F3[0]))) *
                                          force_scale_factor);
          const Tforce f1_ypart = llround(((p_f2_factor * F1[1]) -
                                           (t_factor * ((abbcOabab * F2[1]) + F3[1]))) *
                                          force_scale_factor);
          const Tforce f1_zpart = llround(((p_f2_factor * F1[2]) -
                                           (t_factor * ((abbcOabab * F2[2]) + F3[2]))) *
                                          force_scale_factor);
          const Tforce f3_xpart = llround(t_factor * F2[0] * force_scale_factor);
          const Tforce f3_ypart = llround(t_factor * F2[1] * force_scale_factor);
          const Tforce f3_zpart = llround(t_factor * F2[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - f1_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - f1_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - f1_zpart;
          xfrc[frame2_atom] += f1_xpart - f3_xpart;
          yfrc[frame2_atom] += f1_ypart - f3_ypart;
          zfrc[frame2_atom] += f1_zpart - f3_zpart;
          xfrc[frame3_atom] += f3_xpart;
          yfrc[frame3_atom] += f3_ypart;
          zfrc[frame3_atom] += f3_zpart;
        }
      }
      else {
        xfrc[parent_atom] += vs_frc[0] - (p_f2_factor * F1[0]) +
                             (t_factor * ((abbcOabab * F2[0]) + F3[0]));
        yfrc[parent_atom] += vs_frc[1] - (p_f2_factor * F1[1]) +
                             (t_factor * ((abbcOabab * F2[1]) + F3[1]));
        zfrc[parent_atom] += vs_frc[2] - (p_f2_factor * F1[2]) +
                             (t_factor * ((abbcOabab * F2[2]) + F3[2]));
        xfrc[frame2_atom] += (p_f2_factor * F1[0]) -
                             (t_factor * (F2[0] + (abbcOabab * F2[0]) + F3[0]));
        yfrc[frame2_atom] += (p_f2_factor * F1[1]) -
                             (t_factor * (F2[1] + (abbcOabab * F2[1]) + F3[1]));
        zfrc[frame2_atom] += (p_f2_factor * F1[2]) -
                             (t_factor * (F2[2] + (abbcOabab * F2[2]) + F3[2]));
        xfrc[frame3_atom] += t_factor * F2[0]; 
        yfrc[frame3_atom] += t_factor * F2[1];
        zfrc[frame3_atom] += t_factor * F2[2];
      }
    }
    break;
  case VirtualSiteKind::OUT_3:
    {
      Tcalc p_f2[3], p_f3[3], partition_f2[3], partition_f3[3];
      Tcalc mf2_01, mf2_02, mf2_12, mf3_01, mf3_02, mf3_12;
      if (tcoord_is_sgnint) {
        const Tcalc d3_factor = frame_d3 * inv_gpos_factor;
        if (crd_overflow) {
          const int95_t dpf2x = hostInt95Subtract(xcrd[frame2_atom], xcrd_ovrf[frame2_atom],
                                                  xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t dpf2y = hostInt95Subtract(ycrd[frame2_atom], ycrd_ovrf[frame2_atom],
                                                  ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t dpf2z = hostInt95Subtract(zcrd[frame2_atom], zcrd_ovrf[frame2_atom],
                                                  zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          const int95_t dpf3x = hostInt95Subtract(xcrd[frame3_atom], xcrd_ovrf[frame3_atom],
                                                  xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t dpf3y = hostInt95Subtract(ycrd[frame3_atom], ycrd_ovrf[frame3_atom],
                                                  ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t dpf3z = hostInt95Subtract(zcrd[frame3_atom], zcrd_ovrf[frame3_atom],
                                                  zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          mf2_01 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3z));
          mf2_02 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3y));
          mf2_12 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3x));
          mf3_01 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf2z));
          mf3_02 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf2y));
          mf3_12 = d3_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf2x));
        }
        else {
          mf2_01 = d3_factor * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]);
          mf2_02 = d3_factor * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]);
          mf2_12 = d3_factor * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]);
          mf3_01 = d3_factor * static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]);
          mf3_02 = d3_factor * static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]);
          mf3_12 = d3_factor * static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]);
        }
      }
      else {
        mf2_01 = frame_d3 * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]);
        mf2_02 = frame_d3 * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]);
        mf2_12 = frame_d3 * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]);
        mf3_01 = frame_d3 * static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]);
        mf3_02 = frame_d3 * static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]);
        mf3_12 = frame_d3 * static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]);
      }
      partition_f2[0] = ( frame_d1 * vs_frc[0]) - (mf2_01 * vs_frc[1]) + (mf2_02 * vs_frc[2]);
      partition_f2[1] = ( mf2_01 * vs_frc[0]) + (frame_d1 * vs_frc[1]) - (mf2_12 * vs_frc[2]);
      partition_f2[2] = (-mf2_02 * vs_frc[0]) + (mf2_12 * vs_frc[1]) + (frame_d1 * vs_frc[2]);
      partition_f3[0] = ( frame_d2 * vs_frc[0]) + (mf3_01 * vs_frc[1]) - (mf3_02 * vs_frc[2]);
      partition_f3[1] = (-mf3_01 * vs_frc[0]) + (frame_d2 * vs_frc[1]) + (mf3_12 * vs_frc[2]);
      partition_f3[2] = ( mf3_02 * vs_frc[0]) - (mf3_12 * vs_frc[1]) + (frame_d2 * vs_frc[2]);
      if (tforce_is_sgnint) {
        if (frc_overflow) {
          const int95_t f2_xpart = hostDoubleToInt95(partition_f2[0] * force_scale_factor);
          const int95_t f2_ypart = hostDoubleToInt95(partition_f2[1] * force_scale_factor);
          const int95_t f2_zpart = hostDoubleToInt95(partition_f2[2] * force_scale_factor);
          const int95_t f3_xpart = hostDoubleToInt95(partition_f3[0] * force_scale_factor);
          const int95_t f3_ypart = hostDoubleToInt95(partition_f3[1] * force_scale_factor);
          const int95_t f3_zpart = hostDoubleToInt95(partition_f3[2] * force_scale_factor);
          int95_t ivsr_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom],
                                             f2_xpart.x, f2_xpart.y);
          int95_t ivsr_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom],
                                             f2_ypart.x, f2_ypart.y);
          int95_t ivsr_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom],
                                             f2_zpart.x, f2_zpart.y);
          ivsr_x = hostSplitFPSubtract(ivsr_x, f3_xpart);
          ivsr_y = hostSplitFPSubtract(ivsr_y, f3_ypart);
          ivsr_z = hostSplitFPSubtract(ivsr_z, f3_zpart);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], ivsr_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], ivsr_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], ivsr_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], f2_xpart);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], f2_ypart);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], f2_zpart);
          hostSplitFPSum(&xfrc[frame3_atom], &xfrc_ovrf[frame3_atom], f3_xpart);
          hostSplitFPSum(&yfrc[frame3_atom], &yfrc_ovrf[frame3_atom], f3_ypart);
          hostSplitFPSum(&zfrc[frame3_atom], &zfrc_ovrf[frame3_atom], f3_zpart);
        }
        else {
          const Tforce f2_xpart = llround(partition_f2[0] * force_scale_factor);
          const Tforce f2_ypart = llround(partition_f2[1] * force_scale_factor);
          const Tforce f2_zpart = llround(partition_f2[2] * force_scale_factor);
          const Tforce f3_xpart = llround(partition_f3[0] * force_scale_factor);
          const Tforce f3_ypart = llround(partition_f3[1] * force_scale_factor);
          const Tforce f3_zpart = llround(partition_f3[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - f2_xpart - f3_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - f2_ypart - f3_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - f2_zpart - f3_zpart;
          xfrc[frame2_atom] += f2_xpart;
          yfrc[frame2_atom] += f2_ypart;
          zfrc[frame2_atom] += f2_zpart;
          xfrc[frame3_atom] += f3_xpart;
          yfrc[frame3_atom] += f3_ypart;
          zfrc[frame3_atom] += f3_zpart;
        }
      }
      else {
        xfrc[parent_atom] += vs_frc[0] - partition_f2[0] - partition_f3[0];
        yfrc[parent_atom] += vs_frc[1] - partition_f2[1] - partition_f3[1];
        zfrc[parent_atom] += vs_frc[2] - partition_f2[2] - partition_f3[2];
        xfrc[frame2_atom] += partition_f2[0];
        yfrc[frame2_atom] += partition_f2[1];
        zfrc[frame2_atom] += partition_f2[2];
        xfrc[frame3_atom] += partition_f3[0];
        yfrc[frame3_atom] += partition_f3[1];
        zfrc[frame3_atom] += partition_f3[2];
      }
    }
    break;
  case VirtualSiteKind::FIXED_4:
    {
      Tcalc p_f2[3], rj_f3[3], rj_f4[3], rj_f34[3], rm[3], rt[3], fb[3], fc[3], fd[3];
      if (tcoord_is_sgnint) {
        const Tcalc d1_factor = frame_d1 * inv_gpos_factor;
        const Tcalc d2_factor = frame_d2 * inv_gpos_factor;
        if (crd_overflow) {
          p_f2[0] = displacement(parent_atom, frame2_atom, xcrd, xcrd_ovrf, inv_gpos_factor);
          p_f2[1] = displacement(parent_atom, frame2_atom, ycrd, ycrd_ovrf, inv_gpos_factor);
          p_f2[2] = displacement(parent_atom, frame2_atom, zcrd, zcrd_ovrf, inv_gpos_factor);
          const int95_t dpf3x = hostInt95Subtract(xcrd[frame3_atom], xcrd_ovrf[frame3_atom],
                                                  xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t dpf3y = hostInt95Subtract(ycrd[frame3_atom], ycrd_ovrf[frame3_atom],
                                                  ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t dpf3z = hostInt95Subtract(zcrd[frame3_atom], zcrd_ovrf[frame3_atom],
                                                  zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          rj_f3[0] = (d1_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3x))) - p_f2[0];
          rj_f3[1] = (d1_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3y))) - p_f2[1];
          rj_f3[2] = (d1_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf3z))) - p_f2[2];
          const int95_t dpf4x = hostInt95Subtract(xcrd[frame4_atom], xcrd_ovrf[frame4_atom],
                                                  xcrd[parent_atom], xcrd_ovrf[parent_atom]);
          const int95_t dpf4y = hostInt95Subtract(ycrd[frame4_atom], ycrd_ovrf[frame4_atom],
                                                  ycrd[parent_atom], ycrd_ovrf[parent_atom]);
          const int95_t dpf4z = hostInt95Subtract(zcrd[frame4_atom], zcrd_ovrf[frame4_atom],
                                                  zcrd[parent_atom], zcrd_ovrf[parent_atom]);
          rj_f4[0] = (d2_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf4x))) - p_f2[0];
          rj_f4[1] = (d2_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf4y))) - p_f2[1];
          rj_f4[2] = (d2_factor * static_cast<Tcalc>(hostSplitFPToReal(dpf4z))) - p_f2[2];
        }
        else {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          rj_f3[0] = (d1_factor * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom])) -
                     p_f2[0];
          rj_f3[1] = (d1_factor * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom])) -
                     p_f2[1];
          rj_f3[2] = (d1_factor * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom])) -
                     p_f2[2];
          rj_f4[0] = (d2_factor * static_cast<Tcalc>(xcrd[frame4_atom] - xcrd[parent_atom])) -
                     p_f2[0];
          rj_f4[1] = (d2_factor * static_cast<Tcalc>(ycrd[frame4_atom] - ycrd[parent_atom])) -
                     p_f2[1];
          rj_f4[2] = (d2_factor * static_cast<Tcalc>(zcrd[frame4_atom] - zcrd[parent_atom])) -
                     p_f2[2];
        }
      }
      else {
        p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
        p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
        p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
        rj_f3[0] = (frame_d1 * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom])) -
                   p_f2[0];
        rj_f3[1] = (frame_d1 * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom])) -
                   p_f2[1];
        rj_f3[2] = (frame_d1 * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom])) -
                   p_f2[2];
        rj_f4[0] = (frame_d2 * static_cast<Tcalc>(xcrd[frame4_atom] - xcrd[parent_atom])) -
                   p_f2[0];
        rj_f4[1] = (frame_d2 * static_cast<Tcalc>(ycrd[frame4_atom] - ycrd[parent_atom])) -
                   p_f2[1];
        rj_f4[2] = (frame_d2 * static_cast<Tcalc>(zcrd[frame4_atom] - zcrd[parent_atom])) -
                   p_f2[2];
      }
      rj_f34[0] = rj_f4[0] - rj_f3[0];
      rj_f34[1] = rj_f4[1] - rj_f3[1];
      rj_f34[2] = rj_f4[2] - rj_f3[2];
      crossProduct(rj_f3, rj_f4, rm);
      const Tcalc invr2_rm = (tcalc_is_double) ?
                             value_one / ((rm[0] * rm[0]) + (rm[1] * rm[1]) + (rm[2] * rm[2])) :
                             value_one / ((rm[0] * rm[0]) + (rm[1] * rm[1]) + (rm[2] * rm[2]));
      const Tcalc invr_rm  = (tcalc_is_double) ? sqrt(invr2_rm) : sqrtf(invr2_rm);
      const Tcalc cfx = frame_d3 * invr_rm * vs_frc[0];
      const Tcalc cfy = frame_d3 * invr_rm * vs_frc[1];
      const Tcalc cfz = frame_d3 * invr_rm * vs_frc[2];
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
      rt[0] = ((rj_f4[1] * rm[2]) - (rj_f4[2] * rm[1])) * invr2_rm * frame_d1;
      rt[1] = ((rj_f4[2] * rm[0]) - (rj_f4[0] * rm[2])) * invr2_rm * frame_d1;
      rt[2] = ((rj_f4[0] * rm[1]) - (rj_f4[1] * rm[0])) * invr2_rm * frame_d1;
      fc[0] = (-rm[0] * rt[0] * cfx) - (((frame_d1 * rj_f4[2]) + (rm[1] * rt[0])) * cfy) +
              (((frame_d1 * rj_f4[1]) - (rm[2] * rt[0])) * cfz);
      fc[1] = (((frame_d1 * rj_f4[2]) - (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) -
              (((frame_d1 * rj_f4[0]) + (rm[2] * rt[1])) * cfz);
      fc[2] = (-((frame_d1 * rj_f4[1]) + (rm[0] * rt[2])) * cfx) +
              (((frame_d1 * rj_f4[0]) - (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
      rt[0] = ((rm[1] * rj_f3[2]) - (rm[2] * rj_f3[1])) * invr2_rm * frame_d2;
      rt[1] = ((rm[2] * rj_f3[0]) - (rm[0] * rj_f3[2])) * invr2_rm * frame_d2;
      rt[2] = ((rm[0] * rj_f3[1]) - (rm[1] * rj_f3[0])) * invr2_rm * frame_d2;
      fd[0] = (-rm[0] * rt[0] * cfx) + (((frame_d2 * rj_f3[2]) - (rm[1] * rt[0])) * cfy) -
              (((frame_d2 * rj_f3[1]) + (rm[2] * rt[0])) * cfz);
      fd[1] = (-((frame_d2 * rj_f3[2]) + (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) +
              (((frame_d2 * rj_f3[0]) - (rm[2] * rt[1])) * cfz);
      fd[2] = (((frame_d2 * rj_f3[1]) - (rm[0] * rt[2])) * cfx) +
              (-((frame_d2 * rj_f3[0]) + (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
      if (tforce_is_sgnint) {
        if (frc_overflow) {
          const int95_t fb_xpart = hostDoubleToInt95(fb[0] * force_scale_factor);
          const int95_t fb_ypart = hostDoubleToInt95(fb[1] * force_scale_factor);
          const int95_t fb_zpart = hostDoubleToInt95(fb[2] * force_scale_factor);
          const int95_t fc_xpart = hostDoubleToInt95(fc[0] * force_scale_factor);
          const int95_t fc_ypart = hostDoubleToInt95(fc[1] * force_scale_factor);
          const int95_t fc_zpart = hostDoubleToInt95(fc[2] * force_scale_factor);
          const int95_t fd_xpart = hostDoubleToInt95(fd[0] * force_scale_factor);
          const int95_t fd_ypart = hostDoubleToInt95(fd[1] * force_scale_factor);
          const int95_t fd_zpart = hostDoubleToInt95(fd[2] * force_scale_factor);
          int95_t ivsr_x = fb_xpart + fc_xpart + fd_xpart;
          int95_t ivsr_y = fb_ypart + fc_ypart + fd_ypart;
          int95_t ivsr_z = fb_zpart + fc_zpart + fd_zpart;
          ivsr_x = hostInt95Subtract(xfrc[vsite_atom], xfrc_ovrf[vsite_atom], ivsr_x.x, ivsr_x.y);
          ivsr_y = hostInt95Subtract(yfrc[vsite_atom], yfrc_ovrf[vsite_atom], ivsr_y.x, ivsr_y.y);
          ivsr_z = hostInt95Subtract(zfrc[vsite_atom], zfrc_ovrf[vsite_atom], ivsr_z.x, ivsr_z.y);
          hostSplitFPSum(&xfrc[parent_atom], &xfrc_ovrf[parent_atom], ivsr_x);
          hostSplitFPSum(&yfrc[parent_atom], &yfrc_ovrf[parent_atom], ivsr_y);
          hostSplitFPSum(&zfrc[parent_atom], &zfrc_ovrf[parent_atom], ivsr_z);
          hostSplitFPSum(&xfrc[frame2_atom], &xfrc_ovrf[frame2_atom], fb_xpart);
          hostSplitFPSum(&yfrc[frame2_atom], &yfrc_ovrf[frame2_atom], fb_ypart);
          hostSplitFPSum(&zfrc[frame2_atom], &zfrc_ovrf[frame2_atom], fb_zpart);
          hostSplitFPSum(&xfrc[frame3_atom], &xfrc_ovrf[frame3_atom], fc_xpart);
          hostSplitFPSum(&yfrc[frame3_atom], &yfrc_ovrf[frame3_atom], fc_ypart);
          hostSplitFPSum(&zfrc[frame3_atom], &zfrc_ovrf[frame3_atom], fc_zpart);
          hostSplitFPSum(&xfrc[frame4_atom], &xfrc_ovrf[frame4_atom], fd_xpart);
          hostSplitFPSum(&yfrc[frame4_atom], &yfrc_ovrf[frame4_atom], fd_ypart);
          hostSplitFPSum(&zfrc[frame4_atom], &zfrc_ovrf[frame4_atom], fd_zpart);
        }
        else {
          const Tforce fb_xpart = llround(fb[0] * force_scale_factor);
          const Tforce fb_ypart = llround(fb[1] * force_scale_factor);
          const Tforce fb_zpart = llround(fb[2] * force_scale_factor);
          const Tforce fc_xpart = llround(fc[0] * force_scale_factor);
          const Tforce fc_ypart = llround(fc[1] * force_scale_factor);
          const Tforce fc_zpart = llround(fc[2] * force_scale_factor);
          const Tforce fd_xpart = llround(fd[0] * force_scale_factor);
          const Tforce fd_ypart = llround(fd[1] * force_scale_factor);
          const Tforce fd_zpart = llround(fd[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - fb_xpart - fc_xpart - fd_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - fb_ypart - fc_ypart - fd_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - fb_zpart - fc_zpart - fd_zpart;
          xfrc[frame2_atom] += fb_xpart;
          yfrc[frame2_atom] += fb_ypart;
          zfrc[frame2_atom] += fb_zpart;
          xfrc[frame3_atom] += fc_xpart;
          yfrc[frame3_atom] += fc_ypart;
          zfrc[frame3_atom] += fc_zpart;
          xfrc[frame4_atom] += fd_xpart;
          yfrc[frame4_atom] += fd_ypart;
          zfrc[frame4_atom] += fd_zpart;
        }
      }
      else {
        xfrc[parent_atom] += vs_frc[0] - fb[0] - fc[0] - fd[0];
        yfrc[parent_atom] += vs_frc[1] - fb[1] - fc[1] - fd[1];
        zfrc[parent_atom] += vs_frc[2] - fb[2] - fc[2] - fd[2];
        xfrc[frame2_atom] += fb[0];
        yfrc[frame2_atom] += fb[1];
        zfrc[frame2_atom] += fb[2];
        xfrc[frame3_atom] += fc[0];
        yfrc[frame3_atom] += fc[1];
        zfrc[frame3_atom] += fc[2];
        xfrc[frame4_atom] += fd[0];
        yfrc[frame4_atom] += fd[1];
        zfrc[frame4_atom] += fd[2];
      }
    }
    break;
  case VirtualSiteKind::NONE:
    break;
  }

  // Eliminate any force on the virtual site--it has been transferred to the frame atoms.
  if (tforce_is_sgnint) {
    xfrc[vsite_atom] = 0LL;
    yfrc[vsite_atom] = 0LL;
    zfrc[vsite_atom] = 0LL;
    if (frc_overflow) {
      xfrc_ovrf[vsite_atom] = 0;
      yfrc_ovrf[vsite_atom] = 0;
      zfrc_ovrf[vsite_atom] = 0;
    }
  }
  else {
    xfrc[vsite_atom] = 0.0;
    yfrc[vsite_atom] = 0.0;
    zfrc[vsite_atom] = 0.0;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, const UnitCellType unit_cell,
                               const VirtualSiteKit<Tcalc> &vsk,
                               const Tcalc gpos_scale_factor, const Tcalc force_scale_factor) {
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom  = vsk.frame1_idx[i];
    const int frame2_atom  = vsk.frame2_idx[i];
    const int param_idx = vsk.vs_param_idx[i];
    transmitVirtualSiteForces(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, umat, invu, unit_cell,
                              vsk.vs_atoms[i], vsk.frame1_idx[i], vsk.frame2_idx[i],
                              vsk.frame3_idx[i], vsk.frame4_idx[i],
                              static_cast<VirtualSiteKind>(vsk.vs_types[param_idx]),
                              vsk.dim1[param_idx], vsk.dim2[param_idx], vsk.dim3[param_idx],
                              gpos_scale_factor, force_scale_factor);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void transmitVirtualSiteForces(PhaseSpaceWriter *psw, const VirtualSiteKit<Tcalc> vsk) {
  transmitVirtualSiteForces(psw->xcrd, psw->ycrd, psw->zcrd, psw->xfrc, psw->yfrc, psw->zfrc,
                            psw->umat, psw->invu, psw->unit_cell, vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void transmitVirtualSiteForces(PsSynthesisWriter *poly_psw, const SyValenceKit<Tcalc> &poly_vk,
                               const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk) {

  // Allocate space to mock the GPU's __shared__ memory buffers
  std::vector<int2> vwu_map(vwu_abstract_length);
  std::vector<llint> sh_xfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zfrc(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_yfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xcrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_ycrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zcrd_ovrf(maximum_valence_work_unit_atoms);

  // Set pointers with regard to the detail present in the synthesis
  llint* xfrc_ptr = sh_xfrc.data();
  llint* yfrc_ptr = sh_yfrc.data();
  llint* zfrc_ptr = sh_zfrc.data();
  int *xfrc_ovrf_ptr, *yfrc_ovrf_ptr, *zfrc_ovrf_ptr;
  if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
    xfrc_ovrf_ptr = sh_xfrc_ovrf.data();
    yfrc_ovrf_ptr = sh_yfrc_ovrf.data();
    zfrc_ovrf_ptr = sh_zfrc_ovrf.data();
  }
  else {
    xfrc_ovrf_ptr = nullptr;
    yfrc_ovrf_ptr = nullptr;
    zfrc_ovrf_ptr = nullptr;
  }
  llint* xcrd_ptr = sh_xcrd.data();
  llint* ycrd_ptr = sh_ycrd.data();
  llint* zcrd_ptr = sh_zcrd.data();
  int *xcrd_ovrf_ptr, *ycrd_ovrf_ptr, *zcrd_ovrf_ptr;
  if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
    xcrd_ovrf_ptr = sh_xcrd_ovrf.data();
    ycrd_ovrf_ptr = sh_ycrd_ovrf.data();
    zcrd_ovrf_ptr = sh_zcrd_ovrf.data();
  }
  else {
    xcrd_ovrf_ptr = nullptr;
    ycrd_ovrf_ptr = nullptr;
    zcrd_ovrf_ptr = nullptr;
  }

  // Determine the transform stride.
  const int xfrm_stride = roundUp(9, warp_size_int);

  // Loop over all valence work units.
  for (int i = 0; i < poly_vk.nvwu; i++) {

    // Extract the valence work unit's abstract.
    for (int j = 0; j < vwu_abstract_length; j++) {
      vwu_map[j] = poly_vk.vwu_abstracts[(i * vwu_abstract_length) + j];
    }

    // Determine the system index.
    const int sys_idx = vwu_map[static_cast<size_t>(VwuAbstractMap::SYSTEM_ID)].x;

    // Import the atoms of the work unit.
    const int2 import_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::IMPORT)];
    for (int j = import_limits.x; j < import_limits.y; j++) {
      const size_t synth_atom = poly_vk.vwu_imports[j];
      const size_t local_atom = j - import_limits.x;
      sh_xfrc[local_atom] = poly_psw->xfrc[synth_atom];
      sh_yfrc[local_atom] = poly_psw->yfrc[synth_atom];
      sh_zfrc[local_atom] = poly_psw->zfrc[synth_atom];
      if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
        sh_xfrc_ovrf[local_atom] = poly_psw->xfrc_ovrf[synth_atom];
        sh_yfrc_ovrf[local_atom] = poly_psw->yfrc_ovrf[synth_atom];
        sh_zfrc_ovrf[local_atom] = poly_psw->zfrc_ovrf[synth_atom];
      }
      sh_xcrd[local_atom] = poly_psw->xcrd[synth_atom];
      sh_ycrd[local_atom] = poly_psw->ycrd[synth_atom];
      sh_zcrd[local_atom] = poly_psw->zcrd[synth_atom];
      if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        sh_xcrd_ovrf[local_atom] = poly_psw->xcrd_ovrf[synth_atom];
        sh_ycrd_ovrf[local_atom] = poly_psw->ycrd_ovrf[synth_atom];
        sh_zcrd_ovrf[local_atom] = poly_psw->zcrd_ovrf[synth_atom];
      }
    }
    double* umat = &poly_psw->umat[xfrm_stride * sys_idx];
    double* invu = &poly_psw->invu[xfrm_stride * sys_idx];

    // Perform virtual site force transmission
    const int2 vs_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::VSITE)];
    for (int j = vs_limits.x; j < vs_limits.y; j++) {
      const uint2 t_insr = poly_auk.vste_insr[j];
      const int vsite_atom  = (t_insr.x & 0x3ff);
      const int parent_atom = ((t_insr.x >> 10) & 0x3ff);
      const int frame2_atom = ((t_insr.x >> 20) & 0x3ff);
      const int frame3_atom = (t_insr.y & 0x3ff);
      const int frame4_atom = ((t_insr.y >> 10) & 0x3ff);
      const uint param_idx = (t_insr.y >> 20);
      const Tcalc4 vs_frame = poly_auk.vs_params[param_idx];
      const VirtualSiteKind frame_type = static_cast<VirtualSiteKind>(vs_frame.w);
      transmitVirtualSiteForces<llint,
                                llint, Tcalc>(xcrd_ptr, ycrd_ptr, zcrd_ptr, xfrc_ptr, yfrc_ptr,
                                              zfrc_ptr, umat, invu, poly_psw->unit_cell,
                                              vsite_atom, parent_atom, frame2_atom, frame3_atom,
                                              frame4_atom, frame_type, vs_frame.x, vs_frame.y,
                                              vs_frame.z, poly_psw->gpos_scale,
                                              poly_psw->frc_scale, xcrd_ovrf_ptr, ycrd_ovrf_ptr,
                                              zcrd_ovrf_ptr, xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                              zfrc_ovrf_ptr);
    }
    
    // Write results back to the main force arrays.  In serial processing mode all forces can be
    // copied back without risk of race conditions.
    for (int j = import_limits.x; j < import_limits.y; j++) {
      const size_t synth_atom = poly_vk.vwu_imports[j];
      const size_t local_atom = j - import_limits.x;
      poly_psw->xfrc[synth_atom] = sh_xfrc[local_atom];
      poly_psw->yfrc[synth_atom] = sh_yfrc[local_atom];
      poly_psw->zfrc[synth_atom] = sh_zfrc[local_atom];
      if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
        poly_psw->xfrc_ovrf[synth_atom] = sh_xfrc_ovrf[local_atom];
        poly_psw->yfrc_ovrf[synth_atom] = sh_yfrc_ovrf[local_atom];
        poly_psw->zfrc_ovrf[synth_atom] = sh_zfrc_ovrf[local_atom];
      }
    }
  }
}

} // namespace structure
} // namespace stormm
