// -*-c++-*-
#ifndef STORMM_SYNTHESIS_PERMUTOR_H
#define STORMM_SYNTHESIS_PERMUTOR_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Chemistry/chemical_features.h"
#include "Math/tickcounter.h"
#include "Namelists/nml_conformer.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace synthesis {

using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::ChemicalFeatures;
using chemistry::ConformationEdit;
using chemistry::CoupledEdit;
using chemistry::IsomerPlan;
using math::TickCounter;
using namelist::ConformerControls;
using topology::AtomGraph;

/// \brief Read-only abstract for the SynthesisPermutor.  Access to all of the atom groups,
///        critical measurements, and settings for mutable components is contained herein.
template <typename T> struct SyPermutorKit {
public:

  /// \brief As with other abstracts, the constructor takes input arguments for every member
  ///        variable.
  explicit SyPermutorKit(int system_count_in, int perm_map_count_in,
                         const int* perm_map_idx_in, const int2* perm_elements_in,
                         const int* perm_element_bounds_in, const int* system_settings_in,
                         const int* rot_grp_atoms_in, const int* rot_grp_bounds_in,
                         const int* prm_rot_grp_bounds_in, const int* ctx_grp_atoms_in,
                         const int* ctx_grp_bounds_in, const int* prm_ctx_grp_bounds_in,
                         const int* inv_grp_atoms_in, const int* inv_grp_bounds_in,
                         const int* prm_inv_grp_bounds_in, const int4* rot_bond_markers_in,
                         const int4* ctx_bond_markers_in, const int4* chiral_markers_in,
                         const T* rot_bond_settings_in, const T* ctx_bond_settings_in,
                         const int* rot_bond_settings_bounds_in,
                         const int* ctx_bond_settings_bounds_in,
                         const int* chiral_settings_in, const int* chiral_settings_bounds_in);

  // The member variables store size constants for arrays (the numbers of systems in the synthesis
  // and the number of unique prmologies (translating to the number of unique permutor maps)
  // involved in manipulating each structure of the synthesis.
  int system_count;                     ///< Number of systems in the synthesis
  int perm_map_count;                   ///< Number of permutation maps (unique topologies in the
                                        ///<   synthesis)
  const int* perm_map_idx;              ///< Permutation map indices for every system in the
                                        ///<   synthesis
  const int2* perm_elements;            ///< Array identities and array indices of each element in
                                        ///<   each permutor map.  This information directs threads
                                        ///<   on where to search in rot_bond_settings,
                                        ///<   ctx_bond_settings, or chiral_settings after reading
                                        ///<   limits set forth in those arrays' bounds.
  const int* perm_element_bounds;       ///< Bounds array for perm_elements
  const int* system_settings;           ///< Array of settings indices for all elements of all
                                        ///<   systems in the synthesis.
  const int* rot_grp_atoms;             ///< Rotatable bond group associated atoms
  const int* rot_grp_bounds;            ///< Bounds array for rot_grp_atoms
  const int* prm_rot_grp_bounds;        ///< Bounds array for rot_grp_bounds--giving the limits of
                                        ///<   rot_grp_atoms limits within each permutor map
  const int* ctx_grp_atoms;             ///< Cis-trans isomeric bond group atoms
  const int* ctx_grp_bounds;            ///< Bounds array for ctx_grp_atoms
  const int* prm_ctx_grp_bounds;        ///< Bounds array for ctx_grp_bounds--giving the limits of
                                        ///<   ctx_grp_atoms limits within each permutor map
  const int* inv_grp_atoms;             ///< Chiral center inversion group atoms
  const int* inv_grp_bounds;            ///< Bounds array for inv_grp_atoms
  const int* prm_inv_grp_bounds;        ///< Bounds array for inv_grp_bounds--giving the limits of
                                        ///<   inv_grp_atoms limits within each permutor map
  const int4* rot_bond_markers;         ///< Marker atoms for each rotatable bond, used to define
                                        ///<   the angle about that bond
  const int4* ctx_bond_markers;         ///< Marker atoms for each cis-trans isomeric bond, used to
                                        ///<   define the angle about that bond
  const int4* chiral_markers;           ///< Marker atoms for each chiral center, giving chiral arm
                                        ///<   priorities in the order .x = 0, .y = 3, .z = 2, and
                                        ///<   .w = 1.  See the ChemicalFeatures object
                                        ///<   documentation for more details.
  const T* rot_bond_settings;           ///< Specific values that each rotatable bond can take on
  const T* ctx_bond_settings;           ///< Specific values that each cis-trans isomeric bond can
                                        ///<   take on
  const int* chiral_settings;           ///< Specific values that each chiral center can take on
  const int* rot_bond_settings_bounds;  ///< Bounds array for rot_bond_settings
  const int* ctx_bond_settings_bounds;  ///< Bounds array for ctx_bond_settings
  const int* chiral_settings_bounds;    ///< Bounds array for chiral_settings
};
  
/// \brief An object for tracking the states and The object itself stores a series of permutor
///        maps, detailing the atoms that move as a consequence of rotating about some bond or
///        inverting a chiral center.  One such map is kept for each unique topology, and the
///        atom indices of a particular system are then obtained by adding the atom offset from
///        the PhaseSpaceSynthesis, Condensate, or AtomGraphSynthesis.
class SynthesisPermutor {
public:

  /// \brief The constructor accepts an array of pointers to existing chemical features objects,
  ///        in which case a pointer to the original object for each system will be retained.  The
  ///        constructor also accepts an array of topologies, in which case the chemical features
  ///        object will be temporarily created for the purpose of extracting critical arrays of
  ///        atoms used in manipulations while the pointer to the chemical features object itself
  ///        will go on to read null.
  ///
  /// \param poly_ag          Topology synthesis (sufficient for generating the object, but will
  ///                         calculate and discard chemical features for all individual
  ///                         topologies)
  /// \param poly_ps          Coordinate synthesis (also sufficient for generating the object,
  ///                         but will again imply creating and discarding chemical features)
  /// \param chemfe_in        A series of pointers to chemical features objects for each system
  /// \param retain_pointers  Indicate whether a pointer to the provided chemical features object
  ///                         should be retained.  This has a default value of true to make it
  ///                         transparent to future development, with any SynthesisPermutor created
  ///                         based on a topology passing a value of false.
  /// \param timer            Wall time tracking object, passed down to ChemicalFeatures
  ///                         computations
  /// \{
  SynthesisPermutor();
  SynthesisPermutor(const AtomGraphSynthesis *poly_ag, StopWatch *timer = nullptr);
  SynthesisPermutor(const AtomGraphSynthesis &poly_ag, StopWatch *timer = nullptr);
  SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps, StopWatch *timer = nullptr);
  SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps, StopWatch *timer = nullptr);
  SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in, bool retain_pointers = true,
                    StopWatch *timer = nullptr);
  SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in, bool retain_pointers = true,
                    StopWatch *timer = nullptr);
  /// \}

  /// \brief The presence of POINTER-kind Hybrid objects and multiple ARRAY-kind members means
  ///        that the copy and move constructors as well as assignment operators must be written
  ///        out explicitly.
  ///
  /// \param original  The object to copy or move
  /// \param other     Another object sitting on the right hand side of the assignment operator
  /// \{
  SynthesisPermutor(const SynthesisPermutor &original);
  SynthesisPermutor(SynthesisPermutor &&original);
  SynthesisPermutor& operator=(const SynthesisPermutor &other);
  SynthesisPermutor& operator=(SynthesisPermutor &&other);
  /// \}
  
  /// \brief Get the total number of systems tracked in the synthesis.
  int getSystemCount() const;

  /// \brief Get the number of unique permutor sets (the number of unique topologies served by the
  ///        object).
  int getPermutorSetCount() const;

  /// \brief Get the default number of rotatable bond angle samples.
  int getRotatableBondSampleCount() const;
  
  /// \brief Get the default number of cis-trans isomeric bond angle samples.
  int getCisTransBondSampleCount() const;

  /// \brief Get the numbers of options for each mutable element in a specific system.
  ///
  /// Overloaded:
  ///   - Search by system index
  ///   - Search by topology
  ///
  /// \param system_index  Index of the system of interest
  /// \param query_ag      Topology of interest (must match one of the topology pointers guiding a
  ///                      map, or the function will return an error)
  /// \{
  std::vector<int> getElementSampleCounts(int system_index) const;
  std::vector<int> getElementSampleCounts(const AtomGraph *query_ag) const;
  std::vector<int> getElementSampleCounts(const AtomGraph &query_ag) const;
  /// \}

  /// \brief Get the total number of variable elements in one of the systems: rotatable bonds,
  ///        cis-trans bonds, and chiral centers (even centers that cannot invert are counted).
  ///
  /// Overloaded:
  ///   - Search by system index
  ///   - Search by topology
  ///
  /// \param system_index  Index of the system of interest
  int getVariableElementCount(int system_index) const;

  /// \brief Get the number of rotatable groups.
  ///
  /// \param system_index  Index of the system of interest
  int getRotatableBondCount(int system_index) const;

  /// \brief Get the number of cis-trans isomeric groups.
  ///
  /// \param system_index  Index of the system of interest
  int getCisTransBondCount(int system_index) const;

  /// \brief Get the number of chiral centers.
  ///
  /// \param system_index  Index of the system of interest
  int getChiralCenterCount(int system_index) const;

  /// \brief Get one of the rotatable groups.
  ///
  /// \param system_index  Index of the system of interest
  /// \param group_indx    The rotatable group of interest
  const IsomerPlan& getRotatableGroup(int system_index, int group_index);

  /// \brief Get the additional degrees of freedom coupled to one in particular.
  std::vector<CoupledEdit> getCoupledEdits(int system_index, ConformationEdit motion,
                                           int element_index) const;
  
  /// \brief Get a const reference to the state tracker for the molecule.  Also useable to get a
  ///        copy of the state tracker which is then mutable.
  ///
  /// \param system_index  Index of the system of interest
  const TickCounter<double>& getStateTracker(int system_index) const;

  /// \brief Load specific settings into the available states of rotatable bonds for a particular
  ///        system.  If the vectors differ from the originally stated rotatable bond sample
  ///        counts, the specific system will be updated.
  ///
  /// \param system_index             Index of the system of interest
  /// \param rotatable_bond_settings  A vector of vectors containing the acceptable
  void defineRotatableBondSettings(int system_index,
                                   const std::vector<std::vector<double>> &settings);

  /// \brief Load specific settings into the available states of cis-trans bonds for a particular
  ///        system.  If the vectors are not of length 2 (the assumption for the number of choices
  ///        in a cis-trans rotatable bond), the number of states will be updated to accommodate
  ///        the provided choices.
  ///
  /// \param system_index             Index of the system of interest
  /// \param rotatable_bond_settings  A vector of vectors containing the acceptable
  void defineCisTransBondSettings(int system_index,
                                  const std::vector<std::vector<double>> &settings);

  /// \brief Load specific settings into the available states of an invertible chiral center for a
  ///        particular system.
  ///
  /// \param system_index             Index of the system of interest
  /// \param rotatable_bond_settings  A vector of vectors containing the acceptable settings as
  ///                                 0 or 1
  void defineInvertibleCenterSettings(int system_index,
                                      const std::vector<std::vector<int>> &settings);

  /// \brief Set the general number of rotation samples per rotatable bond
  ///
  /// \param samples  The number of samples to make for each rotatable bond in each system
  void setRotatableBondSamples(int samples_in);
  
  /// \brief Get the pointer to the chemical features object, which can in turn reveal the original
  ///        topology governing the molecule.
  const ChemicalFeatures* getChemicalFeaturesPointer() const;

  /// \brief Get the double-precision read-only abstract for the object.
  ///
  /// \param tier  Get data pointers relevant to the CPU, or to the GPU.
  const SyPermutorKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision read-only abstract for the object.
  ///
  /// \param tier  Get data pointers relevant to the CPU, or to the GPU.
  const SyPermutorKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU device.
  void upload();

  /// \brief Download data from the GPU device.
  void download();
#endif

private:
  int system_count;                              ///< Total number of systems, and the length of
                                                 ///<   the permutor_indices array
  int permutor_map_count;                        ///< The number of permutor maps for different
                                                 ///<   systems in use by this object
  int rotatable_bond_samples;                    ///< General number of times to sample each
                                                 ///<   rotatable bond (default 3)
  int cis_trans_bond_samples;                    ///< General number of times to sample each
                                                 ///<   cis-trans isomeric bond (default 2)
  bool features_objects_accessible;              ///< Indicates whether the ChemicalFeatures
                                                 ///<   objects upon which the permutor maps are
                                                 ///<   based can be accessed by other functions
                                                 ///<   once the object is created.  Some of the
                                                 ///<   constructors create the objects temporarily
                                                 ///<   and the SynthesisPermutors thus created
                                                 ///<   cannot return valid pointers to the
                                                 ///<   features afterwards.
  Hybrid<int> permutor_map_indices;              ///< Index of the permutation map for each system
                                                 ///<   in the synthesis
  Hybrid<int2> permutor_elements;                ///< CoupledEdit objects, translated to int2,
                                                 ///<   indicating the appropriate arrays, and
                                                 ///<   indices therein, where the settings for
                                                 ///<   each mutable element of a given permtuor
                                                 ///<   map can be found.
  Hybrid<int> permutor_element_bounds;           ///< Bounds array for permutor_elements.  The
                                                 ///<   total number of variable elements in the
                                                 ///<   ith permutor map is given by
                                                 ///<   permutor_element_bounds[i + 1] -
                                                 ///<   permutor_element_bounds[i]
  Hybrid<int> system_settings;                   ///< Settings for each system in the synthesis
                                                 ///<   referencing the appropriate TickCounter
                                                 ///<   states.  One TickCounter, a CPU-only
                                                 ///<   object, is stored for each unique permutor
                                                 ///<   map, each unique topology in the synthesis.
                                                 ///<   Feed the settings in this array into the
                                                 ///<   TickCounter to see the actual angle values,
                                                 ///<   or reference the permutor_elements array,
                                                 ///<   then find the appropriate value series in
                                                 ///<   rotatable_bond_settings,
                                                 ///<   cis_trans_bond_settings, or
                                                 ///<   chiral_settings to get the values in memory
                                                 ///<   available to the GPU.
  Hybrid<int> synthesis_data;                    ///< An ARRAY-kind Hybrid targeted by the above
                                                 ///<   Hybrid<int> objects supporting the
                                                 ///<   synthesis.
  Hybrid<int> rotatable_group_atoms;             ///< Atoms present in each rotatable group of
                                                 ///<   each molecular system, starting with the
                                                 ///<   root and pivot atoms defining the rotatable
                                                 ///<   bond / axis of rotation and continuing to
                                                 ///<   list all that move as a consequence of
                                                 ///<   rotating about the bond.  All such groups
                                                 ///<   are listed back-to-back, with no padding
                                                 ///<   between them.
  Hybrid<int> rotatable_group_bounds;            ///< Bounds array for rotatable_group_atoms.  This
                                                 ///<   array knows where each list of topological
                                                 ///<   atom indices in a rotatable group starts
                                                 ///<   and stops, but not which topology that list
                                                 ///<   pertains to.  For that information, one
                                                 ///<   must follow the array below,
                                                 ///<   system_rotatable_group_bounds, stepping
                                                 ///<   over the range of indices for a particular
                                                 ///<   system.
  Hybrid<int> permutor_rotatable_group_bounds;   ///< Bounds array for rotatable_group_atoms.  This
                                                 ///<   array defines the range of rotatable groups
                                                 ///<   for the ith system in its ith and (i + 1)th
                                                 ///<   elements.
  Hybrid<int> cis_trans_group_atoms;             ///< Atoms present in each isomeric group of each
                                                 ///<   molecular system, beginning with the root
                                                 ///<   and pivot atoms of the bond itself.  All
                                                 ///<   such systems are listed back-to-back, with
                                                 ///<   no padding between them.
  Hybrid<int> cis_trans_group_bounds;            ///< Bounds array for cis_trans_group_atoms, laid
                                                 ///<   out in an analogous manner to
                                                 ///<   rotatable_group_bounds
  Hybrid<int> permutor_cis_trans_group_bounds;   ///< Bounds array for cis_trans_group_bounds, laid
                                                 ///<   out in an analogous manner to
                                                 ///<   system_rotatable_group_bounds
  Hybrid<int> invertible_group_atoms;            ///< Atoms present in each chiral inversion group.
                                                 ///<   The list of moving atoms is prefaced by
                                                 ///<   the indices of the anchor atoms for the two
                                                 ///<   smallest chains.
  Hybrid<int> invertible_group_bounds;           ///< Bounds array for invertible_group_atoms
  Hybrid<int> permutor_invertible_group_bounds;  ///< Bounds array for invertible_group_bounds
  Hybrid<int4> rotatable_bond_markers;           ///< Atoms that define the angle made by a
                                                 ///<   particular rotatable bond
  Hybrid<int4> cis_trans_bond_markers;           ///< Atoms that define the angle made by a
                                                 ///<   particular cis-trans isomeric bond
  Hybrid<int4> chiral_markers;                   ///< Atoms that define chirality of a particular
                                                 ///<   center, with the priority order going
                                                 ///<   y > z > w > x (see the chiral_arm_atoms
                                                 ///<   member variable in the ChemicalFeatures
                                                 ///<   object)
  Hybrid<double> rotatable_bond_settings;        ///< Angle settings for each rotatable bond, as
                                                 ///<   defined by atoms in rotatable_bond_markers
  Hybrid<double> cis_trans_bond_settings;        ///< Angle settings for each cis-trans bond, as
                                                 ///<   defined by atoms in cis_trans_bond_markers
  Hybrid<float> sp_rotatable_bond_settings;      ///< Angle settings for each rotatable bond, as
                                                 ///<   defined by atoms in rotatable_bond_markers
                                                 ///<   (single-precision)
  Hybrid<float> sp_cis_trans_bond_settings;      ///< Angle settings for each cis-trans bond, as
                                                 ///<   defined by atoms in cis_trans_bond_markers
                                                 ///<   (single-precision)
  Hybrid<int> rotatable_bond_settings_bounds;    ///< Bounds array for rotatable_bond_settings
  Hybrid<int> cis_trans_bond_settings_bounds;    ///< Bounds array for cis_trans_bond_settings
  Hybrid<int> group_data;                        ///< ARRAY-kind Hybrid object targeted by the
                                                 ///<   preceding POINTER-kind Hybrid<int> objects
  Hybrid<int4> marker_data;                      ///< ARRAY-kind Hybrid object targeted by the
                                                 ///<   preceding POINTER-kind Hybrid<int4> objects
  Hybrid<int> chiral_settings;                   ///< ChiralOrientation enumerations translated to
                                                 ///<   int values for each possible setting of
                                                 ///<   each chiral center.  If chiral sampling is
                                                 ///<   active for a specific topology, each center
                                                 ///<   will have two options.  If chiral sampling
                                                 ///<   is not active for a specific topology, each
                                                 ///<   center will have one option for every
                                                 ///<   orientation (D- or L-) found in structures
                                                 ///<   of the original synthesis.  This is an
                                                 ///<   ARRAY-kind Hybrid object that does not
                                                 ///<   target group_data.
  Hybrid<int> chiral_settings_bounds;            ///< Bounds array for chiral_settings.  This is an
                                                 ///<   ARRAY-kind Hybrid object that does not
                                                 ///<   target group_data.

  // The following CPU-based data arrays help to stage the Hybrid data.
  std::vector<double> general_rotation_settings;  ///< Settings for rotatable bonds, if they are
                                                  ///<   not given for a specific system
  std::vector<double> general_cis_trans_settings; ///< Settings for cis-trans rotatable bonds, if
                                                  ///<   they are not given for a specific system
  std::vector<IsomerPlan> rotatable_groups;       ///< A list of rotatable groups, pulled from
                                                  ///<   each system's chemical features.  The
                                                  ///<   permutor_rotatable_group_bounds array can
                                                  ///<   also be thought of as bounding this.
  std::vector<IsomerPlan> cis_trans_groups;       ///< A list of cis-trans isomeric groups, pulled
                                                  ///<   from each system's chemical features.  The
                                                  ///<   permutor_cis_trans_group_bounds array can
                                                  ///<   also be thought of as bounding this.
  std::vector<IsomerPlan> invertible_groups;      ///< A list of invertible groups for chiral
                                                  ///<   centers in each system, pulled from their
                                                  ///<   respective chemical features objects.  The
                                                  ///<   permutor_invertible_group_bounds array can
                                                  ///<   also be thought of as bounding this.

  /// Exemplary counter wheels for looping over each unique system's mutable variables.
  std::vector<TickCounter<double>> state_trackers;

  /// A list of pointers to the unique topologies guiding each permutor map
  std::vector<AtomGraph*> topologies;
  
  /// Chemical features for each unique topology
  std::vector<ChemicalFeatures*> features;       

  /// \brief Repair (rebase) POINTER-kind HYBRID objects to the object's present ARRAY-kind Hybrid
  ///        targets.
  void rebasePointers();
  
  /// \brief Create a temporary vector of ChemicalFeatures objects based on the unique topologies
  ///        at hand, for use in constructing the object.
  ///
  /// \param timer  Wall time tracker
  std::vector<ChemicalFeatures> temporaryFeatures(StopWatch *timer = nullptr) const;

  /// \brief Create pointers to each member of a list of ChemicalFeatures objects.
  ///
  /// \param feat_in  The list of ChemicalFeatures objects
  const std::vector<ChemicalFeatures*>
  temporaryFeaturesPointers(const std::vector<ChemicalFeatures> &feat_in) const;

  /// \brief Fill the permutor maps based on a defined collection of chemical features.
  ///
  /// \param chemfe_in  The chemical features for each unique system
  /// \param timer      Wall time tracker
  void fillPermutorMaps(const std::vector<ChemicalFeatures*> &chemfe_in, StopWatch *timer);

  /// \brief Set the ranges of variables (rotatable bonds, cis-trans bonds, chiral centers) for
  ///        each permutor map.
  ///
  /// \param confcon  User input derived from a &conformer namelist
  void setVariableRanges(const ConformerControls &confcon);
  
  /// \brief Validate a request for a specific system index.
  ///
  /// \param system_index  The index of interest
  void validateSystemIndex(int system_index) const;
};

} // namespace synthesis
} // namespace stormm

#include "synthesis_permutor.tpp"

#endif

