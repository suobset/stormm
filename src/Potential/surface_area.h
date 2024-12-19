// -*-c++-*-
#ifndef STORMM_SURFACE_AREA_H
#define STORMM_SURFACE_AREA_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Math/geodesic.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;
using stmath::surfaceDistribution;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using trajectory::PhaseSpaceWriter;
using trajectory::PhaseSpace;
  
/// \brief The default amount of energy enclosed within a given solvent-accessible surface area.
constexpr double default_sasa_energy = 1.0;

/// \brief The default probe radius for solvent accessible surface area
constexpr double default_sasa_probe_radius = 1.4;

/// \brief The default number of points to distribute over each sphere surface when computing
///        solvent-accessible surface area.  This is used only for computing the overall value of
///        the surface area, and will have no effect on the accuracy of the derivative.
constexpr int default_sasa_point_count = 272;

/// \brief THe default source of SASA radii is the system's Lennard-Jones sigma parameters.
constexpr SasaReference default_sasa_radii = SasaReference::LENNARD_JONES_SIGMA;

/// \brief Single precision is appropriate for a very rough approximation like the
///        solvent-accessible surface area.
constexpr PrecisionModel default_sasa_precision = PrecisionModel::SINGLE;
  
/// \brief Compute the surface area of a molecule based on the area covered by the center of a
///        sphere rolling across the surface.
///
/// Overloaded:
///   - Provide a single system or a collection of systems
///   - Provide abstracts of the relevant objects or the objects themselves
///   - Provide a pre-determined unit sphere with equidistributed points, or indicate the number of
///     of points to use in covering the surface of any given atom
///   - Use different sources to define the reference atomic radii
///
/// \param psw           Abstract of the coordinate object for a single system
/// \param ps            Coordinates of the molecular system of interest, also contains arrays of
///                      forces for each atom
/// \param isk           Implicit solvent parameters for each atom in the molecular system,
///                      indicating that the atomic radii are to be taken from the atoms' inherent
///                      Generalized Born radii
/// \param nbk           Non-bonded parameters for each atom in the molecular system, indicating
///                      that the radii should be taken from the Lennard-Jones sigma values
/// \param ag            Topology from which to obtain atomic radii
/// \param sphere_pts    An array of points distributed over the surface of the unit sphere, which
///                      can be scaled as needed for any given atom with the intrinsic radius and
///                      solvent probe radius
/// \param probe_radius  The radius of the solvent probe used in computing the solvent-accessible
///                      surface area
/// \param weight        The energy implied by a given amount of surface area, expressed in
///                      kcal/mol-A^2
/// \param radii_source  Array of atomic radii from which to calculate the surface area
/// \{
template <typename Tcalc, typename Tcalc2, typename Tcalc3> std::vector<double>
surfaceArea(PsSynthesisWriter *poly_psw, const SyNonbondedKit<Tcalc, Tcalc2> &synbk,
            const std::vector<Tcalc3> &sphere_pts, double probe_radius = default_sasa_probe_radius,
            double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
            SasaReference radii_source = default_sasa_radii);

template <typename Tcalc, typename Tcalc2> std::vector<double>
surfaceArea(PsSynthesisWriter *poly_psw, const SyNonbondedKit<Tcalc, Tcalc2> &synbk,
            int n_pts = default_sasa_point_count, double probe_radius = default_sasa_probe_radius,
            double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
            SasaReference radii_source = default_sasa_radii);

template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const ImplicitSolventKit<Tcalc> &isk,
                   const std::vector<Tcalc3> &sphere_pts,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy,
                   EvaluateForce eval_frc = EvaluateForce::NO);

template <typename Tcalc, typename Tcalc3>
double surfaceArea(PhaseSpaceWriter *psw, const NonbondedKit<Tcalc> &nbk,
                   const std::vector<Tcalc3> &sphere_pts,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy,
                   EvaluateForce eval_frc = EvaluateForce::NO);

template <typename Tcalc3>
double surfaceArea(PhaseSpace *ps, const AtomGraph *ag, const std::vector<Tcalc3> &sphere_pts,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
                   SasaReference radii_source = default_sasa_radii,
                   PrecisionModel prec = default_sasa_precision);

template <typename Tcalc3>
double surfaceArea(PhaseSpace *ps, const AtomGraph &ag, const std::vector<Tcalc3> &sphere_pts,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
                   SasaReference radii_source = default_sasa_radii,
                   PrecisionModel prec = default_sasa_precision);

double surfaceArea(PhaseSpace *ps, const AtomGraph *ag, int n_pts = default_sasa_point_count,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
                   SasaReference radii_source = default_sasa_radii,
                   PrecisionModel prec = default_sasa_precision);
                   
double surfaceArea(PhaseSpace *ps, const AtomGraph &ag, int n_pts = default_sasa_point_count,
                   double probe_radius = default_sasa_probe_radius,
                   double weight = default_sasa_energy, EvaluateForce eval_frc = EvaluateForce::NO,
                   SasaReference radii_source = default_sasa_radii,
                   PrecisionModel prec = default_sasa_precision);
/// \}
  
} // namespace energy
} // namespace energy

#include "surface_area.tpp"

#endif
