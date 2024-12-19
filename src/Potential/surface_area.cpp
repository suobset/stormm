#include "copyright.h"
#include "surface_area.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
double surfaceArea(PhaseSpace *ps, const AtomGraph *ag, const int n_pts, const double probe_radius,
                   const double weight, const EvaluateForce eval_frc,
                   const SasaReference radii_source, const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  switch (radii_source) {
  case SasaReference::LENNARD_JONES_SIGMA:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const std::vector<double3> sphere_pts = surfaceDistribution<double, double3>(n_pts);
        const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
        return surfaceArea(&psw, nbk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const std::vector<float3> sphere_pts = surfaceDistribution<float, float3>(n_pts);
        const NonbondedKit<float> nbk = ag->getSinglePrecisionNonbondedKit();
        return surfaceArea(&psw, nbk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    }
    break;
  case SasaReference::BORN_RADII:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const std::vector<double3> sphere_pts = surfaceDistribution<double, double3>(n_pts);
        const ImplicitSolventKit<double> isk = ag->getDoublePrecisionImplicitSolventKit();
        return surfaceArea(&psw, isk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const std::vector<float3> sphere_pts = surfaceDistribution<float, float3>(n_pts);
        const ImplicitSolventKit<float> isk = ag->getSinglePrecisionImplicitSolventKit();
        return surfaceArea(&psw, isk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double surfaceArea(PhaseSpace *ps, const AtomGraph &ag, const int n_pts, const double probe_radius,
                   const double weight, const EvaluateForce eval_frc,
                   const SasaReference radii_source, const PrecisionModel prec) {
  return surfaceArea(ps, ag.getSelfPointer(), n_pts, probe_radius, weight, eval_frc, radii_source,
                     prec);
}
  
} // namespace energy
} // namespace stormm
