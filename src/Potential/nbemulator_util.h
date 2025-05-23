// -*-c++-*-
#ifndef STORMM_NBEMULATOR_UTIL_H
#define STORMM_NBEMULATOR_UTIL_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Topology/atomgraph.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using chemistry::ChemicalFeatures;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SynthesisCacheMap;
using synthesis::SystemCache;

/// \brief Default settings for non-bonded emulation potentials.
/// \{
const double default_emul_support_width = 2.0;
const double default_emul_support_start = 4.0;
const double default_emul_basis_zero_hold = 0.0001;
const double default_emul_target_weight = 1.0;
const double default_emul_restraint_weight = 1.0;
const double default_force_zero_hold = 0.001;
const char default_emul_subset_mask[] = ":*";
/// \}

/// \brief This class serves the non-bonded potential emulation features.  Name a structure (by
///        index in an NBEmulator's set of input structures) and provide a list of atoms within
///        that structure which should have no net force on them due to the emulated pairwise
///        potentials.
class ForceZero {
public:
  
  /// \brief The constructor takes initializations for all member variables.
  ForceZero(const std::string &structure_label_in, int structure_index_in,
            const std::string &mask_in, double zero_hold = default_force_zero_hold);

  /// \brief Get the structure label, referencing a label group from a &files namelist.
  const std::string& getStructureLabel() const;
  
  /// \brief Get the sturcture index for which atomic forces are to be near zero.
  int getStructureIndex() const;
  
  /// \brief Get a reference to mask defining which atoms should have their net forces due to
  ///        emulation pair potentials restrained towards zero.  The mask must be interpreted
  ///        against the correct structure and topology using a ChemicalFeatures object.
  const std::string& getMask() const;

  /// \brief Get the regularization holding forces on the selected atoms towards zero.
  double getRegularization() const;
  
private:
  std::string structure_label;  ///< The label group of the structure containing atoms to balance,
                                ///<   as found in a SystemCache obtained from a &files namelist
                                ///<   control block
  int structure_index;          ///< The index of the structure among a larger label group
  std::string mask;             ///< Atom selection within the structure of interest
  double zero_hold;             ///< Weight by which the forces on atoms structure and selected by
                                ///<   the mask will be restrained towards zero.  This weight is
                                ///<   scaled by the number of structures named by the structure
                                ///<   label and index, if the index is negative and therefore
                                ///<   implies that all members of the named label group are to be
                                ///<   included.
};
  
/// \brief The atom types making up sources for each potential must be defined in terms of their
///        distinguishing chemical features.  This class provides options for defining each atom
///        type and then finding its footprint in any given topology.
class NBAtomSource {
public:

  /// \brief The constructor accepts values for all member variables.
  NBAtomSource(const std::string &label_in, const std::vector<int> &atomic_number_in,
               const std::vector<int> &bonds_to_other_atoms_in,
               const std::vector<bool> &in_aromatic_group_in,
               const std::vector<bool> &has_formal_charge_in);

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable.
  /// \{
  NBAtomSource(const NBAtomSource &original) = default;
  NBAtomSource(NBAtomSource &&original) = default;
  NBAtomSource& operator=(const NBAtomSource &original) = default;
  NBAtomSource& operator=(NBAtomSource &&original) = default;
  /// \}

  /// \brief Get the label for this potential source.
  const std::string& getLabel() const;
  
  /// \brief Get one or all options for the atomic number of the source atoms.
  ///
  /// Overloaded:
  ///   - Get the list of all atomic number options
  ///   - Get a specific atomic number specification by index
  ///
  /// \param index  Index of the option of interest
  /// \{
  const std::vector<int>& getAtomicNumber() const;
  int getAtomicNumber(int index) const;
  /// \}

  /// \brief Get the number of bonds to which the defined source atoms must make with other atoms.
  ///        Different indicies indicate different options, e.g. a source may be carbon atoms with
  ///        2, 3, or 4 bonds to other atoms.
  ///
  /// Overloaded:
  ///   - Get the list of all bond count options
  ///   - Get a specific bond count option by index
  ///
  /// \param index  Index of the option of interest
  /// \{
  const std::vector<int>& getBondCount() const;
  int getBondCount(int index) const;
  /// \}

  /// \brief Get the matching conditions related to aromaticity.
  const std::vector<bool>& getAromaticConditions() const;

  /// \brief Get the matching conditions related to formal charge.
  const std::vector<bool>& getFormalChargeConditions() const;

  /// \brief Merge the current object with another, taking unions of all member variable arrays.
  ///
  /// \param other  Another object to merge into this one
  void merge(const NBAtomSource &other);

  /// \brief Check whether all details of this object's source definition are present in another.
  bool subsumedBy(const NBAtomSource &other) const;
 
  /// \brief Determine whether an atom in a topology is a match for this source.
  bool atomMatches(int atom_index, const ChemicalFeatures &chemfe,
                   const std::vector<bool> &aromatic_mask = {}) const;
  
private:
  std::string label;                      ///< The name of this atom source type
  std::vector<int> atomic_number;         ///< List of matching atomic numbers
  std::vector<int> bonds_to_other_atoms;  ///< List of matching numbers of bonds to other atoms
                                          ///<   (virtual sites will not be counted)
  std::vector<bool> in_aromatic_group;    ///< List of aromaticity statuses
  std::vector<bool> has_formal_charge;    ///< Adding an entry of TRUE will include atoms which
                                          ///<   have nonzero formal charge states.  This includes
                                          ///<   any partial formal charges due to resonance.  An
                                          ///<   entry of FALSE in this array will include atoms
                                          ///<   which have no formal charge state.
};

/// \brief The typical target defines a single row of the fitting matrix.
class EmulationTarget {
public:

  /// \brief The constructor again takes intializations for all member variables.
  ///
  /// Overloaded:
  ///   - Indicate a complex, ligand, and receptor, implying that the difference in energy between
  ///     the complex and the sum of ligand and receptor is the target.
  ///   - Indicate a single structure.
  ///   - Indicate a single structure and an alternate structure, implying that the difference in
  ///     energy between the structures is the target.
  /// \{
  EmulationTarget(const std::string &label_in, double target_energy_in, double weight_in,
                  const std::string &first_label_in, const std::string &second_label_in,
                  const std::string &third_label_in, int first_frame_in = 0,
                  int second_frame_in = 0, int third_frame_in = 0,
                  const std::string &first_mask_in = std::string(default_emul_subset_mask),
                  const std::string &second_mask_in = std::string(default_emul_subset_mask),
                  const std::string &third_mask_in = std::string(default_emul_subset_mask),
                  const bool floating_baseline_in = false);

  EmulationTarget(const std::string &label_in, double target_energy_in, double weight_in,
                  const std::string &first_label_in, const std::string &second_label_in,
                  int first_frame_in = 0, int second_frame_in = 0,
                  const std::string &first_mask_in = std::string(default_emul_subset_mask),
                  const std::string &second_mask_in = std::string(default_emul_subset_mask),
                  const bool floating_baseline_in = false);

  EmulationTarget(const std::string &label_in, double target_energy_in, double weight_in,
                  const std::string &first_label_in, int first_frame_in = 0,
                  const std::string &first_mask_in = std::string(default_emul_subset_mask),
                  const bool floating_baseline_in = false);
  /// \}

  /// \brief Get the label for this target.
  const std::string& getLabel() const;

  /// \brief Get the target energy, the right hand side of the fitting equation.
  double getTargetEnergy() const;

  /// \brief Get the weight assigned to this target.
  double getWeight() const;

  /// \brief Get the first system's label in the systems cache (as carried over from a &files
  ///        namelist).  This first structure could be the complex of ligand and receptor, the
  ///        first structure of a pair for which the difference first - second is the target, or
  ///        the only named structure.
  const std::string& getFirstStructureLabel() const;

  /// \brief Get the second system's label in the systems cache.  This second structure could be
  ///        the ligand from a complex, ligand, and receptor triad, or the alternate structure in
  ///        a target where the energy of the first structure minus the alternate (second)
  ///        structure is the objective.
  const std::string& getSecondStructureLabel() const;

  /// \brief Get the third system's label in the systems cache.  This last structure is the
  ///        receptor of a complex, ligand, and receptor triad where the energy difference of
  ///        forming the complex is the objective.
  const std::string& getThirdStructureLabel() const;
  
  /// \brief Get the first system referenced by the target's fitting equation.
  int getFirstStructureFrame() const;

  /// \brief Get the second system referenced by the target's fitting equation.
  int getSecondStructureFrame() const;

  /// \brief Get the third system referenced by the target's fitting equation.
  int getThirdStructureFrame() const;

  /// \brief Get the first system referenced by the target's fitting equation.
  const std::string& getFirstStructureMask() const;

  /// \brief Get the second system referenced by the target's fitting equation.
  const std::string& getSecondStructureMask() const;

  /// \brief Get the third system referenced by the target's fitting equation.
  const std::string& getThirdStructureMask() const;

  /// \brief Determine whether this target involves a relative or absolute energy value.
  bool hasAbsoluteEnergy() const;

  /// \brief Indicate whether the target incorporates a floating baseline energy.  Such baselines
  ///        are incorporated into the data fitting to relieve the fitted parameters which will
  ///        become the new force field from having to "lift" or "lower" the total energy of their
  ///        respective structures, instead putting that effort off onto an extra parameter which
  ///        is then discarded.  This can be helpful when specifying multiple conformations of the
  ///        same structure with absolute energies, even though a force field with accuracy in
  ///        relative energies is desired.  Targets which do not have absolute energies cannot
  ///        incorporate a floating baseline.  A general parameter in the overarching
  ///        EmulatorControls object can be used to specify that all targets with absolute energies
  ///        incorporate a floating baseline energy.  In the fitting matrix, unique baselines will
  ///        be allocated their own columns and cover all similar atom mask selections in
  ///        structures with the same topology.
  bool useFloatingBaseline() const;
  
private:
  std::string label;         ///< The label for this target.  If unspecified a unique label
                             ///<   will be generated automatically.
  double target_energy;      ///< The target energy of the first structure or difference in
                             ///<   energies of the first and second or third structures
  double weight;             ///< Scaling factor applied to all elements of the fitting matrix and
                             ///<   the right-hand side vector related to this target.  This is the
                             ///<   means by which some data points can be given more weight than
                             ///<   others.
  std::string first_label;   ///< Label for the obligatory structure from which the fitted
                             ///<   non-bonded potentials are to derive the target energy.  If
                             ///<   given without second_label and / or third_label, the target
                             ///<   energy is considered to be an absolute energy.  Otherwise, the
                             ///<   target energy is considered to be a relative energy.
  std::string second_label;  ///< Label for the second (optional) structure to consider in the
                             ///<   target energy computation.  If second_label is present and no
                             ///<   third structure is given, the target energy will be assumed to
                             ///<   be the difference in results from first minus second.
  std::string third_label;   ///< Label for the third (optional) structure to consider in the
                             ///<   target energy computation.  If specified, second_label
                             ///<   second must also name a valid label group and perhaps a frame
                             ///<   number within that group.  The target energy will be assumed
                             ///<   emerge from the energy of the first structure, minus the sum of
                             ///<   energies of the second and third combined.
  int first_frame;           ///< Optional frame number within the first structure's label group
  int second_frame;          ///< Optional frame number within the second structure's label group
  int third_frame;           ///< Optional frame number within the third structure's label group
  std::string first_mask;    ///< Mask of atoms in the first structure (by default, the whole
                             ///<   structure), which shall contribute to the fitted energy.
                             ///<   Non-bonded potentials between these atoms and only these atoms
                             ///<   will populate data in the target's row of the fitting matrix.
  std::string second_mask;   ///< Mask of atoms in the second structure that shall contribute to
                             ///<   the fitted energy
  std::string third_mask;    ///< Mask of atoms in the third structure that shall contribute to the
                             ///<   fitted energy
  bool floating_baseline;    ///< Flag to incorporate a floating baseline into the fitted target
                             ///<   energy.  See the accessor function useFloatingBaseline() for a
                             ///<   more detailed description.
};

} // namespace energy
} // namespace stormm

#endif
