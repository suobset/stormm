#include "copyright.h"
#include "Math/rounding.h"
#include "surface_area.h"

namespace stormm {
namespace energy {

using stmath::roundUp;
  
//-------------------------------------------------------------------------------------------------
void surfaceArea(PhaseSpaceSynthesis *poly_ps, ScoreCard *sc, const AtomGraphSynthesis *poly_ag,
                 const std::vector<int> &molecule_idx, const int n_pts, const double probe_radius,
                 const double weight, const EvaluateForce eval_frc,
                 const SasaReference radii_source, const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  const size_t last_sys = poly_psw.system_count - 1;
  const int padded_atom_count = roundUp(poly_psw.atom_starts[last_sys] +
                                        poly_psw.atom_counts[last_sys], 32);
  const int nmol = molecule_idx.size();
  std::vector<uint> sasa_mask((padded_atom_count >> 5), (nmol == 0) * 0xffffffffU);
  if (molecule_idx.size() > 0) {
    for (int i = 0; i < nmol; i++) {
    }
  }
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double, double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit();
      const std::vector<double3> sphere_pts = surfaceDistribution<double, double3>(n_pts);
      surfaceArea<double, double2, double3>(&poly_psw, sc, poly_nbk, sasa_mask.data(), sphere_pts,
                                            probe_radius, weight, eval_frc, radii_source);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float, float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit();
      const std::vector<float3> sphere_pts = surfaceDistribution<float, float3>(n_pts);
      surfaceArea<float, float2, float3>(&poly_psw, sc, poly_nbk, sasa_mask.data(), sphere_pts,
                                         probe_radius, weight, eval_frc, radii_source);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
double surfaceArea(PhaseSpace *ps, const AtomGraph *ag, const std::vector<int> &molecule_idx,
                   const int n_pts, const double probe_radius, const double weight,
                   const EvaluateForce eval_frc, const SasaReference radii_source,
                   const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  const size_t nmol = molecule_idx.size();
  std::vector<bool> sasa_mask(psw.natom, (nmol == 0));
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  for (size_t i = 0; i < nmol; i++) {
    for (int j = cdk.mol_limits[molecule_idx[i]]; j < cdk.mol_limits[molecule_idx[i] + 1]; j++) {
      sasa_mask[cdk.mol_contents[j]] = true;
    }
  }
  switch (radii_source) {
  case SasaReference::LENNARD_JONES_SIGMA:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const std::vector<double3> sphere_pts = surfaceDistribution<double, double3>(n_pts);
        const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
        return surfaceArea(&psw, sasa_mask, nbk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const std::vector<float3> sphere_pts = surfaceDistribution<float, float3>(n_pts);
        const NonbondedKit<float> nbk = ag->getSinglePrecisionNonbondedKit();
        return surfaceArea(&psw, sasa_mask, nbk, sphere_pts, probe_radius, weight, eval_frc);
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
        return surfaceArea(&psw, sasa_mask, isk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const std::vector<float3> sphere_pts = surfaceDistribution<float, float3>(n_pts);
        const ImplicitSolventKit<float> isk = ag->getSinglePrecisionImplicitSolventKit();
        return surfaceArea(&psw, sasa_mask, isk, sphere_pts, probe_radius, weight, eval_frc);
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double surfaceArea(PhaseSpace *ps, const AtomGraph &ag, const std::vector<int> &molecule_idx,
                   const int n_pts, const double probe_radius, const double weight,
                   const EvaluateForce eval_frc, const SasaReference radii_source,
                   const PrecisionModel prec) {
  return surfaceArea(ps, ag.getSelfPointer(), molecule_idx, n_pts, probe_radius, weight, eval_frc,
                     radii_source, prec);
}

} // namespace energy
} // namespace stormm
