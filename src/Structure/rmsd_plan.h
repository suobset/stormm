// -*-c++-*-
#ifndef STORMM_RMSD_PLAN_H
#define STORMM_RMSD_PLAN_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Chemistry/atom_equivalence.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::AtomEquivalence;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;

/// \brief Assume that, if 75% of the mass (or coordinates) of a molecule are asymmetric or
///        otherwise placed in definie positions, the remainder of the molecule's symmetric atom
///        groups can be tested for best fit based on the pose that aligning the rest of the
///        molecule leaves them in.
constexpr double default_required_mass_fraction = 0.75;

/// \brief Read-only abstract for the RMSDPlan object, containing pointers to various buckets of
///        atom indices.
template <typename T> struct RMSDPlanReader {

  /// \brief The constructor takes inputs for all arguments.
  RMSDPlanReader(int plan_count_in, RMSDMethod strategy_in, double mass_fraction_in,
                 const T* masses_in, const int* atom_counts_in, const int* atom_starts_in,
                 const int* alignment_steps_in, const int* core_atoms_in,
                 const int* core_counts_in, const int* core_starts_in, const int* symm_atoms_in,
                 const int4* symm_bounds_in, const int2* symm_ranges_in,
                 const int* frame_counts_in, const int* atr_starts_in,
                 const size_t* ata_starts_in);
  
  /// \brief The presence of const members will implicitly delete the copy and move assignment
  ///        operators, but the default copy and move constructors will apply.
  /// \{
  RMSDPlanReader(const RMSDPlanReader &original);
  RMSDPlanReader(RMSDPlanReader &&original);
  /// \}

  const int plan_count;        ///< The number of individual system RMSD calculation plans
  const RMSDMethod strategy;   ///< RMSD calculations are mass- or coordinate-weighted, aligned or
                               ///<   not
  const double mass_fraction;  ///< Total mass fraction needed to provide a critical alignment for
                               ///<   subsequent determination of symmetry-related atom groups
  const T* masses;             ///< Masses of particles in all systems
  const int* atom_counts;      ///< Numbers of atoms involved in the systems served by each plan
  const int* atom_starts;      ///< Starting positions of the atom series associated with each plan
                               ///<   (each plan pertains to a unique topology)
  const int* alignment_steps;  ///< Protocols for performing RMSD calculations under each plan
  const int* core_atoms;       ///< Concatenated lists of asymmetric atoms in each system
  const int* core_counts;      ///< Exact numbers of asymmetric atoms in each system (each system's
                               ///<   core atom list is padded by the warp size)
  const int* core_starts;      ///< Bounds array for core_atoms, above (use this for offsets, and
                               ///<   core counts for exact counts
  const int* symm_atoms;       ///< Concatenated lists of each system's symmetry-related atoms
  const int4* symm_bounds;     ///< Bounds for symmetry-related atom group lists in the symm_atoms
                               ///<   array.  The "x" and "y" members list the lower and upper
                               ///<   limits of all symmetry-related domains in each group, while
                               ///<   the "z" members list the domain size and the "w" member gives
                               ///<   the total number of domains.
  const int2* symm_ranges;     ///< Bounds array on the symmetric atom bounds array--indicating the
                               ///<   lower and upper limits of symmetric atom numbered groups in
                               ///<   each system in the "x" and "y" members of each tuple.
  const int* frame_counts;     ///< The number of frames referencing each unique topology handled
                               ///<   by this plan
  const int* atr_starts;       ///< Starting indices of the results array for storing reference
                               ///<   RMSD calculation results based on each topology
  const size_t* ata_starts;    ///< Starting indices of the results array for storing all-to-all
                               ///<   matrix RMSD calculation results based on each topology
};
  
/// \brief Collect instructions for one or more systems (intend to work with any coordinate object,
///        including the PhaseSpaceSynthesis)
class RMSDPlan {
public:

  /// \brief The constructor can take any number of AtomEquivalence objects and build plans from
  ///        them.  Topology pointers will be harvested from these objects and, if necessary,
  ///        matched to topologies in a synthesis object in order to connect each system to the
  ///        appropriate instruction set.
  ///
  /// \param strategy_in   The type of RMSD calculations that this plan will guide
  /// \param rmf_in        The required mass fraction that constitutes enough of the molecule to
  ///                      perform placement of small symmetric atom groups without further
  ///                      alignment calculations
  /// \param ag_in         Topology for which to draw a positional RMSD computation plan
  /// \param cf_in         Input coordinates (used for making the necessary AtomEquivalence objects
  ///                      if they are not supplied explicitly)
  /// \param eq_in         Pre-computed breakdown of all symmetry-related atoms
  /// \param poly_ps_in    PhaseSpaceSynthesis object for which to draw plans covering all systems.
  ///                      Topology pointers will be harvested from the object.
  /// \param eq_list_in    List of symmetry-related atom breakdowns for all systems in poly_ps_in
  /// \param gpu           Details of the GPU that will be used in computing RMSD values
  /// \param low_mol_idx   Lower limit of molecules from the topology to assess for symmetry groups
  ///                      (this, and high_mol_idx, are only used if a single topology is provided
  ///                      for automatic deduction of symmetry groups)
  /// \param high_mol_idx  Upper limit of molecules from the topology to assess for symmetry
  ///                      groups.  The default of -1 indicates that only one molecule, indicated
  ///                      by the lower index, should be assessed.
  /// \{
  RMSDPlan(RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction);

  RMSDPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
           RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           int low_mol_idx = 0, int high_mol_idx = -1);

  RMSDPlan(const AtomEquivalence &eq_in, RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);
  
  RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           int low_mol_idx = 0, int high_mol_idx = -1);

  RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const std::vector<AtomEquivalence> &eq_list_in,
           RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief With no POINTER-kind Hybrid or pointers to its own member variables to repair, and no
  ///        const members, the RMSDPlan can make use of default copy and move constructors as well
  ///        as copy and move assignment operators.
  ///
  /// \param original  Another object to copy or move in constructing this one
  /// \param other     The right-hand side object of an assignment
  /// \{
  RMSDPlan(const RMSDPlan &original) = default;
  RMSDPlan(RMSDPlan &&original) = default;
  RMSDPlan& operator=(const RMSDPlan &original) = default;
  RMSDPlan& operator=(RMSDPlan &&original) = default;
  /// \}
  
  /// \brief Get the number of plans kept within this object.
  int getPlanCount() const;

  /// \brief Get the size of the vector of RMSD results for all frames to a single reference, for
  ///        size considerations on an array of results.
  size_t getReferenceRMSDSize() const;
  
  /// \brief Get the size of the array filled with RMSD matrix results for all frames to a single
  ///        reference, for size considerations on an array of results.
  size_t getRMSDMatrixSize() const;

  /// \brief Get the starting index for catalogging RMSD calculations to a reference structure.
  ///
  /// \param plan_index  The system / plan of interest
  size_t getReferenceRMSDStart(const int plan_index) const;
  
  /// \brief Get the starting index for catalogging matrix RMSD calculations.
  ///
  /// \param plan_index  The system / plan of interest
  size_t getRMSDMatrixStart(const int plan_index) const;

  /// \brief Get the number of frames making use of a particular plan, matching the reference RMSD
  ///        bounds set forth in the plan.
  ///
  /// \param plan_index  The system / plan of interest
  int getFrameCount(const int plan_index) const;
  
  /// \brief Get the general, prescribed strategy for computing RMSDs.
  RMSDMethod getGeneralStrategy() const;

  /// \brief Get the mass fraction required to test symmetry-related atom arrangements without
  ///        specific alignments.
  double getRequiredMassFraction() const;

  /// \brief Get one of the topology pointers used by the RMSD plan.
  ///
  /// \param plan_index  The system / plan of interest
  const AtomGraph* getTopologyPointer(int plan_index) const;

  /// \brief Get one of the alignment protocols.
  ///
  /// \param plan_index  The system / plan of interest
  RMSDAlignmentProtocol getAlignmentProtocol(int plan_index) const;

  /// \brief Get the read-only abstract of the system in double precision.  The masses of particles
  ///        are the only templated type.
  ///
  /// \param tier  Get pointers at the level of the CPU host or GPU device
  const RMSDPlanReader<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the read-only abstract of the system in single precision.  The masses of particles
  ///        are the only templated type.
  ///
  /// \param tier  Get pointers at the level of the CPU host or GPU device
  const RMSDPlanReader<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Resize the results arrays to accommodate a given number of RMSD calculations, or an
  ///        assortment of frames.
  ///
  /// \param frame_count  A trusted number of frames to allocate for
  /// \param cs           A set of cooordinate frames (the number of frames in the series will be
  ///                     taken as the resizing value)
  /// \param poly_ps      A synthesis of many systems
  /// \{
  void formatResults(int frame_count);

  template <typename T>
  void formatResults(const CoordinateSeries<T> &cs);
  
  void formatResults(const PhaseSpaceSynthesis &poly_ps);
  /// \}

#ifdef STORMM_USE_HPC
  /// \brief Upload the object's contents to the GPU device.
  void upload();

  /// \brief Download the object's contents from the GPU device.
  void download();
#endif
  
private:
  int plan_count;                      ///< The number of distinct plans, one per system that this
                                       ///<   object is to serve
  int total_rmsd_scores;               ///< The number of RMSD scores that the plan anticipates to
                                       ///<   be computed.
  size_t total_rmsd_pair_scores;       ///< The number of pairwise RMSD scores between frames of
                                       ///<   the same system that the plan anticipates to be
                                       ///<   computed.
  RMSDMethod general_strategy;         ///< Align systems by mass (ALIGN_MASS) or by coordinates of
                                       ///<   all particles (ALIGN_GEOM), no do not align and take
                                       ///<   the mass-weighted (NO_ALIGN_MASS) or the
                                       ///<   coordinate-averaged (NO_ALIGN_GEOM) RMSD
  double required_mass_fraction;       ///< The mass fraction (or number of atoms, if mass is not a
                                       ///<   a consideration)
  Hybrid<double> masses;               ///< Masses of each particle represented in double precision
  Hybrid<float> sp_masses;             ///< Masses of each particle represented in single precision
  Hybrid<int> atom_counts;             ///< Overall numbers of atoms in each unique topology
  Hybrid<int> atom_starts;             ///< Starting indices of each unique system in the masses
                                       ///<   and sp_masses arrays.  The masses are padded up to
                                       ///<   the warp size in each system.
  Hybrid<int> asymmetric_core_atoms;   ///< Atoms found in ach system's asymmetric core, given in
                                       ///<   terms of the system's own topological indices.
  Hybrid<int> asymmetric_core_counts;  ///< Atom counts in each system's core, not part of any of
                                       ///<   its symmetry-related groups
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
  Hybrid<int> alignment_steps;         ///< Integer cast of EquivalenceSwap enumerations for
                                       ///<   the groups of each system.  This is indexes in step
                                       ///<   with symmetry_group_bounds above and each entry can
                                       ///<   be thought of as a fifth member of the tuples in that
                                       ///<   array.
  Hybrid<int2> symmetry_group_ranges;  ///< A bounds array on symmetry_group_bounds, indicating
                                       ///<   the lower and upper limits of symmetry groups
                                       ///<   relevant to each system in its "x" and "y" members,
                                       ///<   respectively.

  /// The number of symmetry-related groups that must be placed in order to ensure that the
  /// required mass fraction of each system is satisfied before subsequent groups are placed.
  /// Because the symmetry-related groups are arranged in order of decreasing size, these numbers
  /// indicate the number of groups to take from the front of the list.  When running the
  /// combinatorial permutations of all such groups, the progress will be stored in an unsigned
  /// 32-bit integers devoting four bits to each symmetry-related group.  This will break if any
  /// of the groups have greater than 16-fold symmetry, which is unheard of in chemistry, and at
  /// most eight distinct groups can be sampled in combinatorial fashion.  At the low end, this
  /// implies at most 2^8 = 256 graph isomorphisms to test for the core alignment, and if just
  /// half of the groups have three-fold, as opposed to two-fold, symmetry, the number of
  /// isomorphism tests rises to 3^4 * 2^4 = 1296, an enormous number of computations just to lay
  /// the foundation for a single RMSD calculation.  But, it implies that a single warp running
  /// each RMSD calculation standas a good chance of working a compute-bound problem.
  Hybrid<int> alignment_prerequisites;

  /// Indicate the approach that must be taken to compute the overall RMSD.  This is an integer
  /// re-cast of the RMSDAlignmentProtocol markers for each system.
  Hybrid<int> alignment_protocols;

  /// The number of symmetry levels indicates the depth of rearrangement that any given atom could
  /// undergo while evaluating different settings for symmmetry-equivalent groups and their
  /// dependents.  The symmetry-equivalent groups to which each atom belongs at any given level of
  /// symmetry are recorded in the array symmetry_group_membership.
  Hybrid<int> symmetry_depth;

  /// The membership of atoms in each symmetry-equivalent group for a given level of evaluation.
  /// The symmetry group is given in the "x" member and the domain of the group in the "y" member.
  /// Entries of -1 in the "x" member indicate that, at that level of symmetry, the atom is not
  /// involved in any symmetry group.  There are entries of this array for every atom of the system
  /// and all possible layers of the system.  Any given atom may be a member of as many distinct
  /// symmetry groups as the system has symmetry layers, but all of those groups will be
  /// dependents, in some tree, upon one symmetry group in particular.  Each system's entries are
  /// padded by the warp size.
  Hybrid<int2> symmetry_group_membership;

  /// Bounds array for symmetry_group_membership above.  The multiplicative effect that the number
  /// of symmetry layers may have on array sizes, plus the need to store additional information
  /// about the upper bound of each atom's applicable groups, forces the type to be ullint rather
  /// than int, even though the total number of atoms in the synthesis cannot exceed INT_MAX.  Each
  /// entry stores the lower limit of the particular atom's symmetry group indices in the low 48
  /// bits.  The number of groups that the atom is included in makes up the high 16 bits.
  /// If SGMB[i] is the ith element of symmetry_group_membership_bounds, then for most atoms,
  ///
  ///   (SGMB[i] & 0xffffffffffffULL) + ((SGMB[i] >> 48) & 0xffffULL) =
  ///   (SGMB[i + 1] & 0xffffffffffffULL).
  ///
  /// In the padded spaces between systems, that relationship will not hold, but it is simply
  /// illustrative.
  Hybrid<ullint> symmetry_group_membership_bounds;

  /// The number of frames based on each plan (each plan pertains to a particular topology).  While
  /// the frames for multiple plans are not guaranteed to exist back-to-back in a
  /// PhaseSpaceSynthesis object, the object will contain lists of its systems that enumerate the
  /// indices corresponding to each topology.
  Hybrid<int> frame_counts;

  /// Offsets for results from each frame for each unique topology, if computing RMSDs of all
  /// frames to a reference structure.
  Hybrid<int> all_to_ref_offsets;
  
  /// When computing all-to-all RMSD matrices for each system, the results are organized in
  /// segments of (Ni * (Ni - 1) / 2) numbers each, the lower triangles of matrices whee Ni is the
  /// number of frames using each unique topology or plan.  The position of the RMSD result for
  /// frames i and j (for i < j) is given by ((j - 1) * (j - 2) / 2) + i.
  Hybrid<size_t> all_to_all_offsets;
  
  /// List of all topologies described by this plan
  std::vector<AtomGraph*> ag_pointers; 

  /// \brief Validate the plan index.
  ///
  /// \param plan_index  Index of the plan that will be requested
  /// \param caller      Name of the calling function
  void validatePlanIndex(int plan_index, const char* caller = nullptr) const;
  
  /// \brief Resize the storage arrays and add groups of symmetry-related atoms to the object's
  ///        Hybrid arrays.
  ///
  /// \param eq_tables  List of all relevant AtomEquivalence objects
  void chartSymmetryGroups(const std::vector<AtomEquivalence> &eq_tables);
  
  /// \brief Compute codes for each atom which, given a series of numbers indicating the
  ///        configurations that each symmetry-related group is to take on, will guide the
  ///        creation of a molecular isomorphism suitable for a direct RMSD calculation.  Greater
  ///        description of the method is found in the implementation.
  void writeSymmetryGroupCodes();
};

} // namespace structure
} // namespace stormm

#include "rmsd_plan.tpp"

#endif
