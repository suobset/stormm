#include "copyright.h"
#include "Topology/atomgraph_abstracts.h"
#include "nbemulator.h"

namespace stormm {
namespace energy {

using chemistry::ChemicalFeaturesReader;
using stmath::foundIn;
using stmath::locateValue;
using stmath::reduceUniqueValues;
using synthesis::SynthesisMapReader;
using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
ForceZero::ForceZero(const std::string &structure_label_in, const int structure_index_in,
                     const std::string &mask_in, const double zero_hold_in) :
    structure_label{structure_label_in}, structure_index{structure_index_in}, mask{mask_in},
    zero_hold{zero_hold_in}
{}

//-------------------------------------------------------------------------------------------------
const std::string& ForceZero::getStructureLabel() const {
  return structure_label;
}

//-------------------------------------------------------------------------------------------------
int ForceZero::getStructureIndex() const {
  return structure_index;
}

//-------------------------------------------------------------------------------------------------
const std::string& ForceZero::getMask() const {
  return mask;
}

//-------------------------------------------------------------------------------------------------
double ForceZero::getRegularization() const {
  return zero_hold;
}

//-------------------------------------------------------------------------------------------------
NBAtomSource::NBAtomSource(const std::string &label_in, const std::vector<int> &atomic_number_in,
                           const std::vector<int> &bonds_to_other_atoms_in,
                           const std::vector<bool> &in_aromatic_group_in,
                           const std::vector<bool> &has_formal_charge_in) :
    label{label_in}, atomic_number{atomic_number_in},
    bonds_to_other_atoms{bonds_to_other_atoms_in}, in_aromatic_group{in_aromatic_group_in},
    has_formal_charge{has_formal_charge_in}
{}

//-------------------------------------------------------------------------------------------------
const std::string& NBAtomSource::getLabel() const {
  return label;
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<int>& NBAtomSource::getAtomicNumber() const {
  return atomic_number;
}

//-------------------------------------------------------------------------------------------------
int NBAtomSource::getAtomicNumber(int index) const {
  if (index < 0 || index >= atomic_number.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a list of " +
          std::to_string(atomic_number.size()) + " selected atomic numbers.", "NBAtomSource",
          "getAtomicNumber");
  }
  return atomic_number[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& NBAtomSource::getBondCount() const {
  return bonds_to_other_atoms;
}

//-------------------------------------------------------------------------------------------------
int NBAtomSource::getBondCount(int index) const {
  if (index < 0 || index >= bonds_to_other_atoms.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a list of " +
          std::to_string(bonds_to_other_atoms.size()) + " selected bond counts.", "NBAtomSource",
          "getBondCount");
  }
  return bonds_to_other_atoms[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<bool>& NBAtomSource::getAromaticConditions() const {
  return in_aromatic_group;
}

//-------------------------------------------------------------------------------------------------
const std::vector<bool>& NBAtomSource::getFormalChargeConditions() const {
  return has_formal_charge;
}

//-------------------------------------------------------------------------------------------------
void NBAtomSource::merge(const NBAtomSource &other) {
  atomic_number.insert(atomic_number.end(), other.atomic_number.begin(),
                       other.atomic_number.end());
  reduceUniqueValues(&atomic_number);
  bonds_to_other_atoms.insert(bonds_to_other_atoms.end(), other.bonds_to_other_atoms.begin(),
                       other.bonds_to_other_atoms.end());
  reduceUniqueValues(&bonds_to_other_atoms);
  in_aromatic_group = reduceUniqueValues(in_aromatic_group, other.in_aromatic_group);
  has_formal_charge = reduceUniqueValues(has_formal_charge, other.has_formal_charge);
}

//-------------------------------------------------------------------------------------------------
bool NBAtomSource::subsumedBy(const NBAtomSource &other) const {
  return (foundIn(atomic_number, other.atomic_number) &&
          foundIn(bonds_to_other_atoms, other.bonds_to_other_atoms) &&
          foundIn(in_aromatic_group, other.in_aromatic_group) &&
          foundIn(has_formal_charge, other.has_formal_charge));
}

//-------------------------------------------------------------------------------------------------
bool NBAtomSource::atomMatches(const int atom_index, const ChemicalFeatures &chemfe,
                               const std::vector<bool> &aromatic_mask) const {
  const AtomGraph *ag = chemfe.getTopologyPointer();
  const int natom = ag->getAtomCount();
  if (atom_index < 0 || atom_index >= natom) {
    rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a topology of " +
          std::to_string(natom) + " atoms.", "NBAtomSource", "atomMatches");
  }
  const int znum = ag->getAtomicNumber(atom_index);
  if (locateValue(atomic_number, znum) < atomic_number.size()) {
    int nreal_bond = 0;
    const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
    for (int i = nbk.nb12_bounds[atom_index]; i < nbk.nb12_bounds[atom_index + 1]; i++) {
      nreal_bond += (ag->getAtomicNumber(nbk.nb12x[i]) > 0);
    }
    const size_t nbond_options = bonds_to_other_atoms.size();
    if (nbond_options == 0 || locateValue(bonds_to_other_atoms, nreal_bond) < nbond_options) {
      const ChemicalFeaturesReader chemfer = chemfe.data();
      const size_t nq_state = has_formal_charge.size();
      const double qval = chemfer.formal_charges[atom_index];
      bool chrg_match = false;
      if (nq_state >= 2) {
        chrg_match = true;
      }
      else {
        for (size_t i = 0; i < nq_state; i++) {
          chrg_match = (chrg_match || (has_formal_charge[i] && fabs(qval) >= 1.0e-6) ||
                        (has_formal_charge[i] == false && fabs(qval) < 1.0e-6));
        }
      }
      const size_t narom_state = in_aromatic_group.size();
      bool arom_match;
      if (narom_state >= 2) {
        arom_match = true;
      }
      else if (in_aromatic_group[0]) {
        arom_match = false;
        if (chemfer.ring_inclusion[atom_index] > 0) {
          if (aromatic_mask.size() == 0) {
            size_t i = 0;
            while (arom_match == false &&
                   i < chemfer.aromatic_group_bounds[chemfer.aromatic_group_count]) {
              arom_match = (atom_index == chemfer.aromatic_groups[i]);
              i++;
            }
          }
          else {
            if (aromatic_mask.size() != natom) {
              rtErr("The size of the aromatic mask (" + std::to_string(aromatic_mask.size()) +
                    ") must match the number of particles in the system (" +
                    std::to_string(natom) + ").", "NBAtomSource", "atomMatches");
            }
            arom_match = aromatic_mask[atom_index];
          }
        }
      }
      else {
        arom_match = true;
        if (chemfer.ring_inclusion[atom_index] > 0) {
          if (aromatic_mask.size() == 0) {
            size_t i = 0;
            while (arom_match && i < chemfer.aromatic_group_bounds[chemfer.aromatic_group_count]) {
              arom_match = (arom_match && (atom_index != chemfer.aromatic_groups[i]));
              i++;
            }
          }
          else {
            if (aromatic_mask.size() != natom) {
              rtErr("The size of the aromatic mask (" + std::to_string(aromatic_mask.size()) +
                    ") must match the number of particles in the system (" +
                    std::to_string(natom) + ").", "NBAtomSource", "atomMatches");
            }
            arom_match = (aromatic_mask[atom_index] == false);
          }
        }
        else {
          arom_match = true;
        }
      }
      return (chrg_match && arom_match);
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
EmulationTarget::EmulationTarget(const std::string &label_in, const double target_energy_in,
                                 const double weight_in, const std::string &first_label_in,
                                 const std::string &second_label_in,
                                 const std::string &third_label_in, const int first_frame_in,
                                 const int second_frame_in, const int third_frame_in,
                                 const std::string &first_mask_in,
                                 const std::string &second_mask_in,
                                 const std::string &third_mask_in,
                                 const bool floating_baseline_in) :
    label{label_in}, target_energy{target_energy_in}, weight{weight_in},
    first_label{first_label_in}, second_label{second_label_in}, third_label{third_label_in},
    first_frame{first_frame_in}, second_frame{second_frame_in}, third_frame{third_frame_in},
    first_mask{first_mask_in}, second_mask{second_mask_in}, third_mask{third_mask_in},
    floating_baseline{floating_baseline_in}
{}

//-------------------------------------------------------------------------------------------------
EmulationTarget::EmulationTarget(const std::string &label_in, const double target_energy_in,
                                 const double weight_in, const std::string &first_label_in,
                                 const std::string &second_label_in, const int first_frame_in,
                                 const int second_frame_in, const std::string &first_mask_in,
                                 const std::string &second_mask_in,
                                 const bool floating_baseline_in) :
    EmulationTarget(label_in, target_energy_in, weight_in, first_label_in, second_label_in,
                    std::string(""), first_frame_in, second_frame_in, 0, first_mask_in,
                    second_mask_in, std::string(default_emul_subset_mask), floating_baseline_in)
{}

//-------------------------------------------------------------------------------------------------
EmulationTarget::EmulationTarget(const std::string &label_in, const double target_energy_in,
                                 const double weight_in, const std::string &first_label_in,
                                 const int first_frame_in, const std::string &first_mask_in,
                                 const bool floating_baseline_in) :
    EmulationTarget(label_in, target_energy_in, weight_in, first_label_in, std::string(""),
                    std::string(""), first_frame_in, 0, 0, first_mask_in,
                    std::string(default_emul_subset_mask), std::string(default_emul_subset_mask),
                    floating_baseline_in)
{}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getLabel() const {
  return label;
}

//-------------------------------------------------------------------------------------------------
double EmulationTarget::getTargetEnergy() const {
  return target_energy;
}

//-------------------------------------------------------------------------------------------------
double EmulationTarget::getWeight() const {
  return weight;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getFirstStructureLabel() const {
  return first_label;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getSecondStructureLabel() const {
  return second_label;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getThirdStructureLabel() const {
  return third_label;
}

//-------------------------------------------------------------------------------------------------
int EmulationTarget::getFirstStructureFrame() const {
  return first_frame;
}

//-------------------------------------------------------------------------------------------------
int EmulationTarget::getSecondStructureFrame() const {
  return second_frame;
}

//-------------------------------------------------------------------------------------------------
int EmulationTarget::getThirdStructureFrame() const {
  return third_frame;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getFirstStructureMask() const {
  return first_mask;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getSecondStructureMask() const {
  return second_mask;
}

//-------------------------------------------------------------------------------------------------
const std::string& EmulationTarget::getThirdStructureMask() const {
  return third_mask;
}

//-------------------------------------------------------------------------------------------------
bool EmulationTarget::hasAbsoluteEnergy() const {
  return (second_label.size() > 0);
}

//-------------------------------------------------------------------------------------------------
bool EmulationTarget::useFloatingBaseline() const {
  return floating_baseline;
}

} // namespace energy
} // namespace stormm
