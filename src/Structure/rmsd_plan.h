// -*-c++-*-
#ifndef STORMM_RMSD_PLAN_H
#define STORMM_RMSD_PLAN_H

#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Chemistry/atom_equivalence.h"
#include "Chemistry/chemical_features.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::AtomEquivalence;
using chemistry::ChemicalFeatures;
using synthesis::PhaseSpaceSynthesis;
using topology::AtomGraph;

/// \brief Assume that, if 75% of the mass (or coordinates) of a molecule are asymmetric or
///        otherwise placed in definie positions, the remainder of the molecule's symmetric atom
///        groups can be tested for best fit based on the pose that aligning the rest of the
///        molecule leaves them in.
constexpr double default_required_mass_fraction = 0.75;

/// \brief Read-only abstract for the RmsdPlan object, containing pointers to various buckets of
///        atom indices.
struct RmsdPlanReader {

  /// \brief The constructor takes inputs for all arguments.
  RmsdPlanReader(int plan_count_in, RmsdMethod strategy_in, double mass_fraction_in,
                 const int* alignment_steps_in, const int* core_atoms_in,
                 const int* core_counts_in, const int* core_starts_in, const int* symm_atoms_in,
                 const int4* symm_bounds_in, const int* symm_ranges_in);
  
  /// \brief The presence of const members will implicitly delete the copy and move assignment
  ///        operators, but the default copy and move constructors will apply.
  /// \{
  RmsdPlanReader(const RmsdPlanReader &original);
  RmsdPlanReader(RmsdPlanReader &&original);
  /// \}

  const int plan_count;        ///< The number of individual system RMSD calculation plans
  const RmsdMethod strategy;   ///< RMSD calculations are mass- or coordinate-weighted, aligned or
                               ///<   not
  const double mass_fraction;  ///< Total mass fraction needed to provide a critical alignment for
                               ///<   subsequent determination of symmetry-related atom groups
  const int* alignment_steps;  ///< Protocols for perfomring RMSD calculations under each plan
  const int* core_atoms;       ///< Concatenated lists of asymmetric atoms in each system
  const int* core_counts;      ///< Exact numbers of asymmetric atoms in each system (each system's
                               ///<   core atom list is padded by the warp size)
  const int* core_starts;      ///< Bounds array for core_atoms, above (use this for offsets, and
                               ///<   core counts for exact counts
  const int* symm_atoms;       ///< Concatenated lists of each system's 
  const int4* symm_bounds;     ///< Bounds for symmetry-related atom group lists in the symm_atoms
                               ///<   array.  The "x" and "y" members list the lower and upper
                               ///<   limits of all symmetry-related domains in each group, while
                               ///<   the "z" members list the domain size and the "w" member gives
                               ///<   the total number of domains.
  const int* symm_ranges;      ///< Bounds array on the symmetric atom bounds array--indicating the
                               ///<   lower and upper limits of symmetric atom numbered groups in
                               ///<   each system
};
  
/// \brief Collect instructions for one or more systems (intend to work with any coordinate object,
///        including the PhaseSpaceSynthesis)
class RmsdPlan {
public:

  /// \brief The constructor can take any number of AtomEquivalence objects and build plans from
  ///        them.  Topology pointers will be harvested from these objects and, if necessary,
  ///        matched to topologies in a synthesis object in order to connect each system to the
  ///        appropriate instruction set.
  ///
  /// \param strategy_in  The type of RMSD calculations that this plan will guide
  /// \param rmf_in       The required mass fraction that constitutes enough of the molecule to
  ///                     perform placement of small symmetric atom groups without further
  ///                     alignment calculations
  /// \param ag_in        Topology for which to draw a positional RMSD computation plan
  /// \param cf_in        Input coordinates (used for making the necessary AtomEquivalence objects
  ///                     if they are not supplied explicitly)
  /// \param eq_in        Pre-computed breakdown of all symmetry-related atoms
  /// \param poly_ps_in   PhaseSpaceSynthesis object for which to draw plans covering all systems.
  ///                     Topology pointers will be harvested from the object.
  /// \param eq_list_in   List of symmetry-related atom breakdowns for all systems in poly_ps_in
  /// \param gpu          Details of the GPU that will be used in computing RMSD values
  /// \{
  RmsdPlan(RmsdMethod strategy_in = RmsdMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction);

  RmsdPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
           RmsdMethod strategy_in = RmsdMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);

  RmsdPlan(const AtomEquivalence &eq_in, RmsdMethod strategy_in = RmsdMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);
  
  RmsdPlan(const PhaseSpaceSynthesis &poly_ps_in, RmsdMethod strategy_in = RmsdMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);

  RmsdPlan(const PhaseSpaceSynthesis &poly_ps_in, const std::vector<AtomEquivalence> &eq_list_in,
           RmsdMethod strategy_in = RmsdMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief With no POINTER-kind Hybrid or pointers to its own member variables to repair, and no
  ///        const members, the RmsdPlan can make use of default copy and move constructors as well
  ///        as copy and move assignment operators.
  ///
  /// \param original  Another object to copy or move in constructing this one
  /// \param other     The right-hand side object of an assignment
  /// \{
  RmsdPlan(const RmsdPlan &original) = default;
  RmsdPlan(RmsdPlan &&original) = default;
  RmsdPlan& operator=(const RmsdPlan &original) = default;
  RmsdPlan& operator=(RmsdPlan &&original) = default;
  /// \}
  
  /// \brief Get the number of plans kept within this object.
  int getPlanCount() const;

  /// \brief Get the general, prescribed strategy for computing RMSDs.
  RmsdMethod getGeneralStrategy() const;

  /// \brief Get the mass fraction required to test symmetry-related atom arrangements without
  ///        specific alignments.
  double getRequiredMassFraction() const;

  /// \brief Get the read-only abstract of the system
  ///
  /// \param tier  Get pointers at the level of the CPU host or GPU device
  const RmsdPlanReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
private:
  int plan_count;                      ///< The number of distinct plans, one per system that this
                                       ///<   object is to serve
  RmsdMethod general_strategy;         ///< Align systems by mass (ALIGN_MASS) or by coordinates of
                                       ///<   all particles (ALIGN_GEOM), no do not align and take
                                       ///<   the mass-weighted (NO_ALIGN_MASS) or the
                                       ///<   coordinate-averaged (NO_ALIGN_GEOM) RMSD
  double required_mass_fraction;       ///< The mass fraction (or number of atoms, if mass is not a
                                       ///<   a consideration)
  Hybrid<int> alignment_steps;         ///< Integer cast of RmsdAlignmentProtocol enumerations for
                                       ///<   each system
  Hybrid<int> asymmetric_core_atoms;   ///< Atoms found in each system's asymmetric core, given in
                                       ///<   terms of the system's own topological indices.
  Hybrid<int> asymmetric_core_counts;  ///< Atom counts in each 
  Hybrid<int> asymmetric_core_starts;  ///< Starting indices in asymmetric_core_atoms for each
                                       ///<   system's list of asymmetric core atoms
  Hybrid<int> symmetry_group_atoms;    ///< Atoms involved in each collection of interchangeable
                                       ///<   groups needed to compute a best-fit RMSD.  These
                                       ///<   indicies work in terms of each system's local
                                       ///<   topology numbering.  Unlike the energy calculation
                                       ///<   work units, which focus on one system at a time and
                                       ///<   contain instructions for pulling out specific atoms
                                       ///<   from the overall list, symmetry group plans will
                                       ///<   describe how to compare one snapshot of a system to
                                       ///<   another in terms of the one system's structure.
  Hybrid<int4> symmetry_group_bounds;  ///< Bounds array for symmetry_group_atoms, indicating the
                                       ///<   limits of atom indices involved in each collection of
                                       ///<   equivalent atoms.  Because each system's
                                       ///<   symmetry-related groups are defined back-to-back,
                                       ///<   without padding, each element of this array will give
                                       ///<   the starting index of symmetry_group_atoms in the "x"
                                       ///<   member and the upper limit in the "y" member.
                                       ///<   Between systems, however, the symmetry-related atom
                                       ///<   lists will be padded by the warp size.  The "z"
                                       ///<   member of each element gives the specific size of
                                       ///<   domains in each symmetry-related group (the number of
                                       ///<   symmetry-related atoms involved in each swap), while
                                       ///<   the "w" member gives the number of domains (the
                                       ///<   number of interchangeable sub-structures).
  Hybrid<int> symmetry_group_ranges;   ///< A bounds array on symmetry_group_bounds, indicating
                                       ///<   the lower and upper limits of symmetry groups
                                       ///<   relevant to each system.

  /// List of all topologies described by this plan
  std::vector<AtomGraph*> ag_pointers; 
};
  
} // namespace structure
} // namespace stormm

#endif
