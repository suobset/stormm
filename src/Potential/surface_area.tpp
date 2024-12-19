// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc3> std::vector<double>
surfaceArea(PsSynthesisWriter *poly_psw, const SyNonbondedKit<Tcalc, Tcalc2> &synbk,
            const std::vector<Tcalc3> &sphere_pts, const double probe_radius, const double weight,
            const EvaluateForce eval_frc, const SasaReference radii_source) {
  std::vector<double> result(poly_psw->system_count, 0.0);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc3> std::vector<double>
surfaceArea(PsSynthesisWriter *poly_psw, const SyNonbondedKit<Tcalc, Tcalc2> &synbk,
            const int n_pts, const double probe_radius, const double weight,
            const EvaluateForce eval_frc, const SasaReference radii_source) {
  std::vector<double> result(poly_psw->system_count, 0.0);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const ImplicitSolventKit<Tcalc> &isk,
                   const std::vector<Tcalc3> &sphere_pts, const double probe_radius,
                   const double weight, const EvaluateForce eval_frc) {
  double result = 0.0;
  const size_t nsph_pts = sphere_pts.size();
  for (int i = 0; i < psw->natom; i++) {
    const Tcalc atom_x = psw->xcrd[i];
    const Tcalc atom_y = psw->ycrd[i];
    const Tcalc atom_z = psw->zcrd[i];
    for (size_t j = 0; j < nsph_pts; j++) {
      const Tcalc pt_x = (sphere_pts[j].x * isk.pb_radii[j]) + atom_x;
      const Tcalc pt_y = (sphere_pts[j].y * isk.pb_radii[j]) + atom_y;
      const Tcalc pt_z = (sphere_pts[j].z * isk.pb_radii[j]) + atom_z;

    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const NonbondedKit<Tcalc> &nbk,
                   const std::vector<Tcalc3> &sphere_pts, const double probe_radius,
                   const double weight, const EvaluateForce eval_frc) {
  double result = 0.0;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc3>
double surfaceArea(PhaseSpace *ps, const AtomGraph *ag, const std::vector<Tcalc3> &sphere_pts,
                   const double probe_radius, const double weight, const EvaluateForce eval_frc,
                   const SasaReference radii_source, const PrecisionModel prec) {
  double result = 0.0;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc3>
double surfaceArea(PhaseSpace *ps, const AtomGraph &ag, const std::vector<Tcalc3> &sphere_pts,
                   const double probe_radius, const double weight, const EvaluateForce eval_frc,
                   const SasaReference radii_source, const PrecisionModel prec) {
  return surfaceArea(ps, ag.getSelfPointer(), sphere_pts, probe_radius, weight, radii_source,
                     prec);
}

} // namespace energy
} // namespace stormm
