#include <algorithm>
#include <cmath>
#include <limits.h>
#include <string>
#include <vector>
#include "Constants/symbol_values.h"
#include "DataTypes/omni_vector_types.h"
#include "Math/rounding.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "indigo.h"
#include "periodic_table.h"

namespace omni {
namespace chemistry {

using card::HybridKind;
using math::minValue;
using math::maxValue;
using math::prefixSumInPlace;
using math::PrefixSumType;
using math::roundUp;
using math::sum;
using parse::char4ToString;
using topology::ChemicalDetailsKit;
using topology::ValenceKit;
using topology::NonbondedKit;
  
//-------------------------------------------------------------------------------------------------
int IndigoFormalChargeKey(const int atomic_number, const int formal_charge) {
  if (atomic_number < 0 || atomic_number > maximum_indigo_atomic_number ||
      formal_charge < minimum_formal_charge || formal_charge > maximum_formal_charge) {
    rtErr("The atomic number must be positive and the formal charge between " +
          std::to_string(minimum_formal_charge) + " and +" +
          std::to_string(maximum_formal_charge) + ".  Atomic number " +
          std::to_string(atomic_number) + " and / or formal charge " +
          std::to_string(formal_charge) + " are invalid.", "IndigoFormalChargeKey");
  }
  return (formal_charge_range * atomic_number) + formal_charge - minimum_formal_charge;
}

//-------------------------------------------------------------------------------------------------
int IndigoBondOrderKey(const int atomic_number_i, const int atomic_number_j,
                       const int bond_order) {
  if (atomic_number_i < 0 || atomic_number_j < 0 ||
      atomic_number_i > maximum_indigo_atomic_number ||
      atomic_number_j > maximum_indigo_atomic_number || bond_order < 1 ||
      bond_order > maximum_bond_order) {
    rtErr("Atomic numbers must be greater than zero and bond orders within {1, 2, 3}.  "
          "Atomic numbers " + std::to_string(atomic_number_i) + " and " +
          std::to_string(atomic_number_j) + ", with bond order " + std::to_string(bond_order) +
          " are invalid.", "IndigoBondOrderKey");
  }
  return (((atomic_number_i * indigo_atomic_number_range) + atomic_number_j) * bond_order_range) +
         bond_order;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> indigoFormalChargeScores() {

  // Pre-allocate the array with very high values, which will be detected later as signifying that
  // no value for a particular charge state is available.
  std::vector<int> result(indigo_atomic_number_range * formal_charge_range, INT_MAX);

  // Hydrogen
  result[IndigoFormalChargeKey( 1, -1)] =    -264;
  result[IndigoFormalChargeKey( 1,  0)] =       0;
  result[IndigoFormalChargeKey( 1,  1)] =   49981;

  // Carbon
  result[IndigoFormalChargeKey( 6, -4)] =  110560;
  result[IndigoFormalChargeKey( 6, -3)] =   48600;
  result[IndigoFormalChargeKey( 6, -2)] =   12642;
  result[IndigoFormalChargeKey( 6, -1)] =   -9020;
  result[IndigoFormalChargeKey( 6,  0)] =   -5385;
  result[IndigoFormalChargeKey( 6,  1)] =   35759;
  result[IndigoFormalChargeKey( 6,  2)] =  124978;
  result[IndigoFormalChargeKey( 6,  3)] =  300539;
  result[IndigoFormalChargeKey( 6,  4)] =  537111;

  // Nitrogen
  result[IndigoFormalChargeKey( 7, -3)] =   70853;
  result[IndigoFormalChargeKey( 7, -2)] =   20738;
  result[IndigoFormalChargeKey( 7, -1)] =   -6710;
  result[IndigoFormalChargeKey( 7,  0)] =  -10055;
  result[IndigoFormalChargeKey( 7,  1)] =   43125;
  result[IndigoFormalChargeKey( 7,  2)] =  151559;
  result[IndigoFormalChargeKey( 7,  3)] =  325361;
  result[IndigoFormalChargeKey( 7,  4)] =  609490;
  result[IndigoFormalChargeKey( 7,  5)] =  968634;

  // Oxygen
  result[IndigoFormalChargeKey( 8, -2)] =    9397;
  result[IndigoFormalChargeKey( 8, -1)] =  -12915;
  result[IndigoFormalChargeKey( 8,  0)] =   -8232;
  result[IndigoFormalChargeKey( 8,  1)] =   40942;
  result[IndigoFormalChargeKey( 8,  2)] =  169748;
  result[IndigoFormalChargeKey( 8,  3)] =  371155;
  result[IndigoFormalChargeKey( 8,  4)] =  654831;
  result[IndigoFormalChargeKey( 8,  5)] = 1072593;
  result[IndigoFormalChargeKey( 8,  6)] = 1579333;

  // Fluorine
  result[IndigoFormalChargeKey( 9, -1)] =  -11934;
  result[IndigoFormalChargeKey( 9,  0)] =       0;
  result[IndigoFormalChargeKey( 9,  1)] =   63230;
  result[IndigoFormalChargeKey( 9,  2)] =  190703;
  result[IndigoFormalChargeKey( 9,  3)] =  420886;
  result[IndigoFormalChargeKey( 9,  4)] =  740644;
  result[IndigoFormalChargeKey( 9,  5)] = 1159413;
  result[IndigoFormalChargeKey( 9,  6)] = 1735880;
  result[IndigoFormalChargeKey( 9,  7)] = 2415246;

  // Phosphorus
  result[IndigoFormalChargeKey(15, -3)] =   47698;
  result[IndigoFormalChargeKey(15, -2)] =   11858;
  result[IndigoFormalChargeKey(15, -1)] =   -7367;
  result[IndigoFormalChargeKey(15,  0)] =   -6395;
  result[IndigoFormalChargeKey(15,  1)] =   32033;
  result[IndigoFormalChargeKey(15,  2)] =  104397;
  result[IndigoFormalChargeKey(15,  3)] =  214834;
  result[IndigoFormalChargeKey(15,  4)] =  403146;
  result[IndigoFormalChargeKey(15,  5)] =  641082;
  
  // Sulfur
  result[IndigoFormalChargeKey(16, -2)] =    2470;
  result[IndigoFormalChargeKey(16, -1)] =  -12008;
  result[IndigoFormalChargeKey(16,  0)] =   -5080;
  result[IndigoFormalChargeKey(16,  1)] =   32100;
  result[IndigoFormalChargeKey(16,  2)] =  117760;
  result[IndigoFormalChargeKey(16,  3)] =  245359;
  result[IndigoFormalChargeKey(16,  4)] =  418644;
  result[IndigoFormalChargeKey(16,  5)] =  684375;
  result[IndigoFormalChargeKey(16,  6)] = 1006529;

  // Chlorine
  result[IndigoFormalChargeKey(17, -1)] =  -12662;
  result[IndigoFormalChargeKey(17,  0)] =       0;
  result[IndigoFormalChargeKey(17,  1)] =   46802;
  result[IndigoFormalChargeKey(17,  2)] =  133168;
  result[IndigoFormalChargeKey(17,  3)] =  278608;
  result[IndigoFormalChargeKey(17,  4)] =  473544;
  result[IndigoFormalChargeKey(17,  5)] =  721469;
  result[IndigoFormalChargeKey(17,  6)] = 1076125;
  result[IndigoFormalChargeKey(17,  7)] = 1493864;

  // Bromine
  result[IndigoFormalChargeKey(35, -1)] =  -12263;
  result[IndigoFormalChargeKey(35,  0)] =       0;
  result[IndigoFormalChargeKey(35,  1)] =   42725;
  result[IndigoFormalChargeKey(35,  2)] =  120276;
  result[IndigoFormalChargeKey(35,  3)] =  249427;
  result[IndigoFormalChargeKey(35,  4)] =  420793;
  result[IndigoFormalChargeKey(35,  5)] =  636814;
  result[IndigoFormalChargeKey(35,  6)] =  951383;
  result[IndigoFormalChargeKey(35,  7)] = 1321067;

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> indigoBondOrderScores() {

  // Again, pre-allocate the array with very high values, which will be detected later as
  // signifying that no value for a particular charge state is available.
  std::vector<int> result(indigo_atomic_number_range * indigo_atomic_number_range *
                          bond_order_range, INT_MAX);

  // Bromine, single bonds
  result[IndigoBondOrderKey(35, 35, 1)] =   -7253;
  result[IndigoBondOrderKey(35, 17, 1)] =   -7827;
  result[IndigoBondOrderKey(17, 35, 1)] =   -7827;
  result[IndigoBondOrderKey(35,  9, 1)] =   -9073;
  result[IndigoBondOrderKey( 9, 35, 1)] =   -9073;
  result[IndigoBondOrderKey(35,  1, 1)] =  -14372;
  result[IndigoBondOrderKey( 1, 35, 1)] =  -14372;
  result[IndigoBondOrderKey(35,  7, 1)] =   -8401;
  result[IndigoBondOrderKey( 7, 35, 1)] =   -8401;
  result[IndigoBondOrderKey(35,  8, 1)] =   -7888;
  result[IndigoBondOrderKey( 8, 35, 1)] =   -7888;
  result[IndigoBondOrderKey(35, 15, 1)] =  -10189;
  result[IndigoBondOrderKey(15, 35, 1)] =  -10189;
  result[IndigoBondOrderKey(35, 16, 1)] =   -8366;
  result[IndigoBondOrderKey(16, 35, 1)] =   -8366;
  result[IndigoBondOrderKey(35,  6, 1)] =  -11532;
  result[IndigoBondOrderKey( 6, 35, 1)] =  -11532;

  // Carbon, single bonds
  result[IndigoBondOrderKey( 6,  6, 1)] =  -15022;
  result[IndigoBondOrderKey( 6, 17, 1)] =  -13227;
  result[IndigoBondOrderKey(17,  6, 1)] =  -13227;
  result[IndigoBondOrderKey( 6,  9, 1)] =  -17572;
  result[IndigoBondOrderKey( 9,  6, 1)] =  -17572;
  result[IndigoBondOrderKey( 6,  1, 1)] =  -17646;
  result[IndigoBondOrderKey( 1,  6, 1)] =  -17646;
  result[IndigoBondOrderKey( 6,  7, 1)] =  -14014;
  result[IndigoBondOrderKey( 7,  6, 1)] =  -14014;
  result[IndigoBondOrderKey( 6,  8, 1)] =  -14967;
  result[IndigoBondOrderKey( 8,  6, 1)] =  -14967;
  result[IndigoBondOrderKey( 6, 15, 1)] =  -11640;
  result[IndigoBondOrderKey(15,  6, 1)] =  -11640;
  result[IndigoBondOrderKey( 6, 16, 1)] =  -11951;
  result[IndigoBondOrderKey(16,  6, 1)] =  -11951;

  // Chlorine, single bonds
  result[IndigoBondOrderKey(17, 17, 1)] =   -8290;
  result[IndigoBondOrderKey(17,  9, 1)] =   -8916;
  result[IndigoBondOrderKey( 9, 17, 1)] =   -8916;
  result[IndigoBondOrderKey(17,  1, 1)] =  -16569;
  result[IndigoBondOrderKey( 1, 17, 1)] =  -16569;
  result[IndigoBondOrderKey(17,  7, 1)] =   -9453;
  result[IndigoBondOrderKey( 7, 17, 1)] =   -9453;
  result[IndigoBondOrderKey(17,  8, 1)] =   -8340;
  result[IndigoBondOrderKey( 8, 17, 1)] =   -8340;
  result[IndigoBondOrderKey(17, 15, 1)] =  -11867;
  result[IndigoBondOrderKey(15, 17, 1)] =  -11867;
  result[IndigoBondOrderKey(17, 16, 1)] =   -9491;
  result[IndigoBondOrderKey(16, 17, 1)] =   -9491;

  // Fluorine, single bonds
  result[IndigoBondOrderKey( 9,  9, 1)] =   -5322;
  result[IndigoBondOrderKey( 9,  1, 1)] =  -21939;
  result[IndigoBondOrderKey( 1,  9, 1)] =  -21939;
  result[IndigoBondOrderKey( 9,  7, 1)] =  -10913;
  result[IndigoBondOrderKey( 7,  9, 1)] =  -10913;
  result[IndigoBondOrderKey( 9,  8, 1)] =   -7280;
  result[IndigoBondOrderKey( 8,  9, 1)] =   -7280;
  result[IndigoBondOrderKey( 9, 15, 1)] =  -17200;
  result[IndigoBondOrderKey(15,  9, 1)] =  -17200;
  result[IndigoBondOrderKey( 9, 16, 1)] =  -12469;
  result[IndigoBondOrderKey(16,  9, 1)] =  -12469;

  // Hydrogen, single bonds
  result[IndigoBondOrderKey( 1,  1, 1)] =  -17291;
  result[IndigoBondOrderKey( 1,  7, 1)] =  -17960;
  result[IndigoBondOrderKey( 7,  1, 1)] =  -17960;
  result[IndigoBondOrderKey( 1,  8, 1)] =  -19492;
  result[IndigoBondOrderKey( 8,  1, 1)] =  -19492;
  result[IndigoBondOrderKey( 1, 15, 1)] =  -13608;
  result[IndigoBondOrderKey(15,  1, 1)] =  -13608;
  result[IndigoBondOrderKey( 1, 16, 1)] =  -14783;
  result[IndigoBondOrderKey(16,  1, 1)] =  -14783;

  // Nitrogen, single bonds
  result[IndigoBondOrderKey( 7,  7, 1)] =  -10899;
  result[IndigoBondOrderKey( 7,  8, 1)] =  -10376;
  result[IndigoBondOrderKey( 8,  7, 1)] =  -10376;
  result[IndigoBondOrderKey( 7, 15, 1)] =  -12181;
  result[IndigoBondOrderKey(15,  7, 1)] =  -12181;
  result[IndigoBondOrderKey( 7, 16, 1)] =  -10518;
  result[IndigoBondOrderKey(16,  7, 1)] =  -10518;

  // Oxygen, single bonds
  result[IndigoBondOrderKey( 8,  8, 1)] =   -7877;
  result[IndigoBondOrderKey( 8, 15, 1)] =  -14082;
  result[IndigoBondOrderKey(15,  8, 1)] =  -14082;
  result[IndigoBondOrderKey( 8, 16, 1)] =  -10816;
  result[IndigoBondOrderKey(16,  8, 1)] =  -10816;

  // Phosphorus, single bonds
  result[IndigoBondOrderKey(15, 15, 1)] =   -9084;
  result[IndigoBondOrderKey(15, 16, 1)] =  -10308;
  result[IndigoBondOrderKey(16, 15, 1)] =  -10308;

  // Sulfur, single bonds
  result[IndigoBondOrderKey(16, 16, 1)] =   -9677;

  // Carbon, double bonds
  result[IndigoBondOrderKey( 6,  6, 2)] =  -28182;
  result[IndigoBondOrderKey( 6,  7, 2)] =  -25265;
  result[IndigoBondOrderKey( 7,  6, 2)] =  -25265;
  result[IndigoBondOrderKey( 6,  8, 2)] =  -28007;
  result[IndigoBondOrderKey( 8,  6, 2)] =  -28007;
  result[IndigoBondOrderKey( 6, 15, 2)] =  -19184;
  result[IndigoBondOrderKey(15,  6, 2)] =  -19184;
  result[IndigoBondOrderKey( 6, 16, 2)] =  -20381;
  result[IndigoBondOrderKey(16,  6, 2)] =  -20381;

  // Nitrogen, double bonds
  result[IndigoBondOrderKey( 7,  7, 2)] =  -19511;
  result[IndigoBondOrderKey( 7,  8, 2)] =  -18116;
  result[IndigoBondOrderKey( 8,  7, 2)] =  -18116;
  result[IndigoBondOrderKey( 7, 15, 2)] =  -16692;
  result[IndigoBondOrderKey(15,  7, 2)] =  -16692;
  result[IndigoBondOrderKey( 7, 16, 2)] =  -14236;
  result[IndigoBondOrderKey(16,  7, 2)] =  -14236;

  // Oxygen, double bonds
  result[IndigoBondOrderKey( 8,  8, 2)] =  -17967;
  result[IndigoBondOrderKey( 8, 15, 2)] =  -20583;
  result[IndigoBondOrderKey(15,  8, 2)] =  -20583;
  result[IndigoBondOrderKey( 8, 16, 2)] =  -18580;
  result[IndigoBondOrderKey(16,  8, 2)] =  -18580;

  // Phosphorus, double bonds
  result[IndigoBondOrderKey(15, 15, 2)] =  -12738;
  result[IndigoBondOrderKey(15, 16, 2)] =  -14290;
  result[IndigoBondOrderKey(16, 15, 2)] =  -14290;

  // Sulfur, double bonds
  result[IndigoBondOrderKey(16, 16, 2)] =  -15083;

  // Carbon, triple bonds
  result[IndigoBondOrderKey( 6,  6, 2)] =  -41714;
  result[IndigoBondOrderKey( 6,  7, 2)] =  -37482;
  result[IndigoBondOrderKey( 7,  6, 2)] =  -37482;
  result[IndigoBondOrderKey( 6, 15, 2)] =  -27955;
  result[IndigoBondOrderKey(15,  6, 2)] =  -27955;

  // Nitrogen, triple bonds
  result[IndigoBondOrderKey( 7,  7, 2)] =  -34502;
  result[IndigoBondOrderKey( 7, 15, 2)] =  -21302;
  result[IndigoBondOrderKey(15,  7, 2)] =  -21302;

  // Phosphorus, triple bonds
  result[IndigoBondOrderKey(15, 15, 2)] =  -16750;

  return result;
}

//-------------------------------------------------------------------------------------------------
uint IndigoAtomCenter::computeState(const int fc_value, const std::vector<uint> &bond_orders) {
  const uint charge_comp = static_cast<uint>(fc_value - minimum_formal_charge);
  uint bond_comp = 0U;
  for (int i = 0; i < bond_count; i++) {
    bond_comp |= (bond_orders[i] << (2*i + formal_charge_bits));
  }
  return (charge_comp | bond_comp);
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::computeScore(const uint state, const std::vector<int> &atom_scores,
                                   const std::vector<int> &bond_scores) {

  // Leave the formal charge in the [0, 16) range.  Here, it is an index into a table, not a
  // quantity to manipulate in a scientific equation.
  const int fc_value = (state & 0xf);
  int result = 2 * atom_scores[(table_index * formal_charge_range) + fc_value];
  for (int i = 0; i < bond_count; i++) {
    const int ord = ((state >> (2*i + formal_charge_bits)) & 0x3);
    result += bond_scores[(relevant_bonds[i] * bond_order_range) + ord];
  }
  return result;
}
                                   
//-------------------------------------------------------------------------------------------------
IndigoAtomCenter::IndigoAtomCenter(const int table_index_in, const int z_number_in,
                                   const int bond_count_in, const int valence_electrons_in,
                                   const std::vector<int> &relevant_bonds_in,
                                   const std::vector<int> &partner_atoms_in,
                                   const std::vector<int> &atom_scores,
                                   const std::vector<int> &bond_scores) :
  table_index{table_index_in},
  z_number{z_number_in},
  bond_count{bond_count_in},
  valence_electrons{valence_electrons_in},
  possible_states{0},
  relevant_bonds{relevant_bonds_in},
  partner_atoms{partner_atoms_in},
  states{},
  scores{}
{
  // Check that the number of bonds will not break the 32-bit unsigned int format
  if ((2 * bond_count) + 8 > sizeof(uint) * 8) {
    rtErr("An atom center with " + std::to_string(bond_count) + " bonds cannot be encoded in " +
          std::to_string(sizeof(uint) * 8) + "-bit unsigned int format.", "IndigoAtomCenter");
  }
  
  // Each bond accounts for at least two electrons out of the valence target.
  const int unclaimed_valence_e = valence_electrons - (2 * bond_count);
  
  // If there are no unclaimed valence electrons simply after accounting for the number of
  // bonds, then every bond must be single order and the formal charge of the atom is fixed.
  // There is only one state.
  possible_states = 1;
  std::vector<uint> bond_orders(bond_count, 0x1);
  if (unclaimed_valence_e == 0) {
    int fc_value;
    switch(z_number) {
    case 1:
      fc_value = 0;
      break;
    case 6:
    case 15:
    case 16:

      // Carbon, phosphorus, and sulfur will all have zero net charge if they fulfill their
      // valence shells in this manner.  If phosphorus or sulfur to have three or fewer bonds,
      // their valence shell targets will shift to eight electrons but the bonds alone cannot
      // be guaranteed to fulfill the valence shells in those cases.
      fc_value = 0;
      break;
    case 7:
      fc_value = 1;
      break;
    case 8:
      fc_value = 2;
      break;
    case 9:
    case 17:
    case 35:
      fc_value = 3;
      break;
    default:
      {
        const char2 esymb = elemental_symbols[z_number];
        rtErr("Element " + std::to_string(esymb.x) + std::to_string(esymb.y) + " is not coevered "
              "by the Indigo scoring function.", "IndigoAtomCenter");
      }
    }

    // Shift the formal charge state to a number in the range [0, 16) with 8 being neutral
    states.push_back(computeState(fc_value, bond_orders));
    scores.push_back(computeScore(states[0], atom_scores, bond_scores));
  }
  else {

    // Deal with possible new states by incrementing each and every bond up to order 3.  Many of
    // these attempts will fail right away because they make too many increments to the bond order,
    // but they can be culled quickly.  This will take (3^bond_count) attempts, but that number is
    // well bounded by chemical features of any atom.
    std::vector<uint> trial_bond_orders(bond_count, 1);
    const uint trial_limit = maximum_bond_order * bond_count;
    uint trial_sum = sum<uint>(trial_bond_orders);
    do {
      
      // Test the new state
      const int free_pairs = (valence_electrons - (2 * static_cast<int>(trial_sum))) / 2;
      if (free_pairs >= 0) {
        int fc_value;
        switch(z_number) {
        case 1:
          rtErr("Element H cannot have a non-neutral formal charge state, or anything other than "
                "one single bond connecting it to other atoms.", "IndigoAtomCenter");
          break;
        case 6:
          fc_value = -free_pairs;
          break;
        case 7:
          fc_value = 1 - free_pairs;
          break;
        case 8:
          fc_value = 2 - free_pairs;
          break;
        case 9:
        case 17:
        case 35:
          fc_value = 3 - free_pairs;
          break;
        case 15:
          if (bond_count >= 3) {

            // The target valence electron population is ten
            fc_value = -free_pairs;
          }
          else {

            // The target valence electron population is eight
            fc_value = 1 - free_pairs;
          }
          break;
        case 16:
          if (bond_count >= 3) {

            // The target valence electron population is twelve
            fc_value = -free_pairs;
          }
          else {

            // The target valence electron population is eight
            fc_value = 2 - free_pairs;
          }
          break;
        default:
          {
            const char2 esymb = elemental_symbols[z_number];
            rtErr("Element " + std::to_string(esymb.x) + std::to_string(esymb.y) + " is not "
                  "covered by the Indigo scoring function.", "IndigoAtomCenter");
          }
        }

        // Compute the new state and add it (with its score) to the list if it is unique
        const uint new_state = computeState(fc_value, trial_bond_orders);
        bool not_found = true;
        const size_t state_count = states.size();
        for (size_t i = 0; i < state_count; i++) {
          not_found = (not_found && new_state != states[i]);
        }
        if (not_found) {
          states.push_back(new_state);
          scores.push_back(computeScore(new_state, atom_scores, bond_scores));
        }
      }
      
      // Increment the vector of bond upgrades
      trial_bond_orders[0] += 1;
      for (int i = 0; i < bond_count; i++) {
        if (trial_bond_orders[i] > maximum_bond_order) {
          trial_bond_orders[i] = 1;
          if (i < bond_count - 1) {
            trial_bond_orders[i + 1] += 1;
          }
        }
      }
      trial_sum = sum<uint>(trial_bond_orders);
    } while (trial_sum < trial_limit);
    possible_states = states.size();
  }
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getIndigoTableIndex() const {
  return table_index;
}
  
//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getStateCount() const {
  return possible_states;
}

//-------------------------------------------------------------------------------------------------
uint IndigoAtomCenter::getState(const int state_index) const {
  return states[state_index];
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getCharge(const int state_index) const {
  return (static_cast<int>(states[state_index] & 0xf) + minimum_formal_charge);
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getScore(const int state_index) const {
  return scores[state_index];
}

//-------------------------------------------------------------------------------------------------
void IndigoAtomCenter::sortStates() {

  // Bubble sort for very small things
  if (possible_states == 1) {
    return;
  }
  else if (possible_states <= 4) {
    for (int i = 0; i < possible_states - 1; i++) {
      for (int j = i + 1; j < possible_states; j++) {
        if (scores[i] > scores[j]) {
          std::swap(scores[i], scores[j]);
          std::swap(states[i], states[j]);
        }
      }
    }
  }
  else {

    // Use a long long int to ensure that the packed bits in the
    // unsigned int state do not overflow the format.
    std::vector<longlong2> at_data(possible_states);
    for (int i = 0; i < possible_states; i++) {
      at_data[i].x = scores[i];
      at_data[i].y = states[i];
    }
    std::sort(at_data.begin(), at_data.end(), [](longlong2 a, longlong2 b) { return a.x < b.x; });
    for (int i = 0; i < possible_states; i++) {
      scores[i] = at_data[i].x;
      states[i] = at_data[i].y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getPartnerCount() const {
  return static_cast<int>(partner_atoms.size());
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getRelevantBond(const int index) const {
  return relevant_bonds[index];
}
  
//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getPartner(const int index) const {
  return partner_atoms[index];
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::findPartnerIndex(const int partner_atom_number) const {
  
  for (int i = 0; i < bond_count; i++) {
    if (partner_atoms[i] == partner_atom_number) {
      return i;
    }
  }
  rtErr("Partner atom " + std::to_string(partner_atom_number) + " was not found in atom center " +
        std::to_string(table_index) + ".", "IndigoAtomCenter", "findPartnerIndex");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::getBondOrderOfState(const int bond_index, const int state_index) const {
  return static_cast<int>((states[state_index] >> ((2 * bond_index) + formal_charge_bits)) & 0x3);
}
  
//-------------------------------------------------------------------------------------------------
bool IndigoAtomCenter::queryBondToPartner(const int partner_atom_number,
                                          const int partner_bo) const {
  int local_index = findPartnerIndex(partner_atom_number);

  // With the partner atom's local index found, scan the available states and determine whether
  // any of them contain the right bond order reaching out to the partner atom.
  bool matched = false;
  for (int i = 0; i < possible_states; i++) {
    const int local_bo = ((states[i] >> ((2 * local_index) + formal_charge_bits)) & 0x3);
    matched = (matched || local_bo == partner_bo);
  }
  return matched;
}
  
//-------------------------------------------------------------------------------------------------
int IndigoAtomCenter::cullByPartners(const std::vector<IndigoAtomCenter> &acen) {
  const int partner_count = partner_atoms.size();
  std::vector<bool> retain_states(possible_states, true);
  for (int i = 0; i < partner_count; i++) {
    for (int j = 0; j < possible_states; j++) {

      // Isolate the bond going to the partner in this state--does the partner have a state with
      // the same bond order pointing back?
      const int tbo = ((states[j] >> (2*i + formal_charge_bits)) & 0x3);
      if (! acen[partner_atoms[i]].queryBondToPartner(table_index, tbo)) {
        retain_states[j] = false;
      }
    }
  }

  // Cull states that are incompatible with one or more partner atoms
  int k = 0;
  int culled_states = 0;
  for (int i = 0; i < possible_states; i++) {
    if (retain_states[i]) {
      states[k] = states[i];
      scores[k] = scores[i];
      k++;
    }
    else {
      culled_states++;
    }
  }
  possible_states = k;
  states.resize(k);
  scores.resize(k);

  return culled_states;
}

//-------------------------------------------------------------------------------------------------
IndigoFragment::IndigoFragment(const std::vector<int> &centers_list_in,
                               const std::vector<IndigoAtomCenter> &all_centers) :
  center_count{static_cast<int>(centers_list_in.size())},
  possible_states{0},
  centers_list{centers_list_in},
  states{},
  states_data{},
  net_charges{},
  scores{}
{
  // Map the interconnections among the centers.  Joined pairs of atom centers are stored in the
  // vector relevant_pairs, which stores the global atom center numbers based on the list of
  // atom centers specific to this fragment and avoids double-counting the interactions.  Use
  // const representations even of local variables to avoid dereferencing &this.
  std::vector<int> settings(center_count, 0);
  std::vector<int> max_settings(center_count);
  std::vector<int4> tmp_relevant_pairs;
  std::vector<int2> tmp_local_pair_idx;
  const int ccenter_count = center_count;
  for (int i = 0; i < ccenter_count; i++) {
    const int tc_index = centers_list_in[i];
    max_settings[i] = all_centers[tc_index].getStateCount();
    const int n_partners = all_centers[tc_index].getPartnerCount();
    for (int j = 0; j < n_partners; j++) {
      const int partner_index = all_centers[tc_index].getPartner(j);

      // Avoid double-counting pairs by requiring that the second member of the pair come later
      // in the IndigoTable.
      if (partner_index < tc_index) {
        continue;
      }
      bool partner_in_fragment = false;
      int local_partner_index;
      for (int k = 0; k < ccenter_count; k++) {
        if (centers_list_in[k] == partner_index) {
          partner_in_fragment = true;
          local_partner_index = k;
          break;
        }
      }
      if (partner_in_fragment) {

        // Find the bond of the partner atom that links back to atom center tc_index
        const int reciprocal_j = all_centers[partner_index].findPartnerIndex(tc_index);

        // The relevant pairs vector stores the indices of two atom centers in its x and y members,
        // providing their numbers in the overarching IndigoTable.  Its z and w members provide the
        // local bond indices of the aotm centers identified in x and y, respectively, which form
        // the bond between atom centers x and y or between y and x.
        tmp_relevant_pairs.push_back({ tc_index, partner_index, j, reciprocal_j });

        // What is still needed are the local indices of each atom center, within this fragment.
        // Another array stores that information.  Together, these arrays define the potential
        // function for this fragment.
        tmp_local_pair_idx.push_back({ i, local_partner_index });
      }
    }
  }
  const int total_options = sum<int>(max_settings) - ccenter_count;
  const int total_pairs = tmp_relevant_pairs.size();

  // Create a key for all pairs involving a particular local atom index.
  std::vector<int> relevant_pair_bounds(ccenter_count + 1, 0);
  for (int i = 0; i < total_pairs; i++) {
    if (tmp_local_pair_idx[i].x > tmp_local_pair_idx[i].y) {
      relevant_pair_bounds[tmp_local_pair_idx[i].x] += 1;
    }
    else {
      relevant_pair_bounds[tmp_local_pair_idx[i].y] += 1;
    }
  }
  prefixSumInPlace<int>(&relevant_pair_bounds, PrefixSumType::EXCLUSIVE, "IndigoFragment");
  std::vector<int4> relevant_pairs(total_pairs);
  std::vector<int2> local_pair_idx(total_pairs);
  for (int i = 0; i < total_pairs; i++) {
    if (tmp_local_pair_idx[i].x > tmp_local_pair_idx[i].y) {
      const int inc_idx = relevant_pair_bounds[tmp_local_pair_idx[i].x];
      relevant_pairs[inc_idx] = tmp_relevant_pairs[i];
      local_pair_idx[inc_idx] = tmp_local_pair_idx[i];
      relevant_pair_bounds[tmp_local_pair_idx[i].x] = inc_idx + 1;
    }
    else {
      const int inc_idx = relevant_pair_bounds[tmp_local_pair_idx[i].y];
      relevant_pairs[inc_idx] = tmp_relevant_pairs[i];
      local_pair_idx[inc_idx] = tmp_local_pair_idx[i];
      relevant_pair_bounds[tmp_local_pair_idx[i].y] = inc_idx + 1;
    }
  }
  for (int i = ccenter_count; i > 0; i--) {
    relevant_pair_bounds[i] = relevant_pair_bounds[i - 1];
  }
  relevant_pair_bounds[0] = 0;
  
  // Create a vector to run through all states of the fragment.  Accumulate results in
  // preliminary arrays before commiting them to the actual object.
  std::vector<int> prelim_net_charges;
  std::vector<int> prelim_scores;
  int centers_participating = 0;
  do {

    // Test whether this combination of settings is viable.  For each pair, make sure that both
    // atoms, whatever state they may be in, are participating in the current fragment structure.
    bool viable = true;
    if (centers_participating >= 1) {
      for (int pos = relevant_pair_bounds[centers_participating - 1];
           pos < relevant_pair_bounds[centers_participating]; pos++) {
        const int local_atom_i = local_pair_idx[pos].x;
        const int local_atom_j = local_pair_idx[pos].y;
        const int atom_i = relevant_pairs[pos].x;
        const int atom_j = relevant_pairs[pos].y;
        const int bond_out_of_i = relevant_pairs[pos].z;
        const int bond_out_of_j = relevant_pairs[pos].w;
        if (all_centers[atom_i].getBondOrderOfState(bond_out_of_i, settings[local_atom_i]) !=
            all_centers[atom_j].getBondOrderOfState(bond_out_of_j, settings[local_atom_j])) {
          viable = false;
          break;
        }
      }
    }
    if (viable) {

      // If not all atom centers are participating, increment the number of centers that are
      // included in the fragment and start the search for the state of the next atom at its
      // first position.  Otherwise, compute the score and net charge of this state, then push
      // those results and the state itself onto the stack for this fragment.
      if (centers_participating < ccenter_count) {
        settings[centers_participating] = 0;
        centers_participating++;
        continue;
      }
      else {
        int st_q = 0;
        int st_score = 0;
        for (int i = 0; i < ccenter_count; i++) {
          const int c_idx = centers_list_in[i];
          const int s_idx = settings[i];
          st_q += all_centers[c_idx].getCharge(s_idx);
          st_score += all_centers[c_idx].getScore(s_idx);
          states_data.push_back(s_idx);
        }
        prelim_net_charges.push_back(st_q);
        prelim_scores.push_back(st_score);
      }
    }

    // Increment the state of the final participating atom, until all of its states are
    // exhausted.  If the final atom has been incremented past its avaiable states, reduce the
    // number of participating atoms by one.  Increment the state of the final participating atom
    // (which was once the next-to-last participating atom).  The next round will determine
    // whether the fragment state is viable, and roll forward again with a new attempt to add
    // atoms.
    int last_participant = centers_participating - 1;
    settings[last_participant] += 1;
    while (settings[last_participant] >= max_settings[last_participant] &&
           centers_participating > 1) {
      centers_participating--;
      last_participant--;
      settings[last_participant] += 1;
    }    
  } while (settings[0] < max_settings[0]);

  // Sort states by energy, then by charge.  Let the maximum charge of any fragment be +/-8192.
  possible_states = prelim_net_charges.size();
  std::vector<longlong2> frag_pkg(possible_states);
  for (int i = 0; i < possible_states; i++) {
    frag_pkg[i].x = (static_cast<llint>(prelim_scores[i]) * 8192LL) +
                    static_cast<llint>(prelim_net_charges[i]);
    frag_pkg[i].y = i;
  }
  std::sort(frag_pkg.begin(), frag_pkg.end(), [](longlong2 a, longlong2 b) { return a.x < b.x; });

  // Commit results to the object.  Set pointers to the state data.
  net_charges.resize(possible_states);
  scores.resize(possible_states);
  states.resize(possible_states);
  int* sdd_ptr = states_data.data();
  for (int i = 0; i < possible_states; i++) {
    const size_t fidx = frag_pkg[i].y;
    net_charges[i] = prelim_net_charges[fidx];
    scores[i] = prelim_scores[fidx];
    states[i] = &sdd_ptr[fidx * ccenter_count];
  }
}

//-------------------------------------------------------------------------------------------------
int IndigoFragment::getStateCount() const {
  return possible_states;
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> IndigoFragment::getState(const int state_index) const {
  std::vector<int2> result(center_count);
  for (int i = 0; i < center_count; i++) {
    result[i] = { centers_list[i], states[state_index][i] };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int IndigoFragment::getCharge(const int state_index) const {
  return net_charges[state_index];
}

//-------------------------------------------------------------------------------------------------
int IndigoFragment::getScore(const int state_index) const {
  return scores[state_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> IndigoFragment::getChargeStates() const {
  std::vector<bool> covered(possible_states, false);
  std::vector<int> result;
  for (int i = 0; i < possible_states; i++) {
    if (covered[i]) {
      continue;
    }
    const int ichg = net_charges[i];
    result.push_back(ichg);
    for (int j = i; j < possible_states; j++) {
      covered[j] = (covered[j] || net_charges[j] == ichg);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int3> IndigoFragment::getChargesAndBestEnergies() const {

  // Find the minimum and maximum possible charges.
  const int min_q = minValue(net_charges);
  const int max_q = maxValue(net_charges);

  // Accumulate results
  std::vector<int3> result(max_q - min_q + 1);
  std::vector<bool> marked(max_q - min_q + 1, false);
  for (int i = 0; i < possible_states; i++) {
    const int pos = net_charges[i] - min_q;
    if (marked[pos] == false || result[pos].y > scores[i]) {
      marked[pos] = true;
      result[pos].x = net_charges[i];
      result[pos].y = scores[i];
      result[pos].z = i;
    }
  }

  // The initial pass at charge and energy states may have left a lot of holes in the result.
  // Tighten it up by removing unrepresented charge states.
  int j = 0;
  for (int i = 0; i <= max_q - min_q; i++) {
    if (marked[i]) {
      result[j] = result[i];
      j++;
    }
  }
  result.resize(j);
  
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> IndigoFragment::getStatesBearingCharge(const int charge_value) const {
  std::vector<int2> result;
  for (int i = 0; i < possible_states; i++) {
    if (net_charges[i] == charge_value) {
      result.push_back({scores[i], i});
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool IndigoFragment::testEquivalence(const IndigoFragment &other,
                                     const std::vector<int> &real_atom_map,
                                     const AtomGraph *ag_pointer,
                                     const std::vector<IndigoAtomCenter> &atom_centers) const {
  if (other.center_count != center_count || other.possible_states != possible_states) {
    return false;
  }
  for (int i = 0; i < possible_states; i++) {
    if (scores[i] != other.scores[i]) {
      return false;
    }
    if (net_charges[i] != other.net_charges[i]) {
      return false;
    }
  }
  for (int i = 0; i < center_count; i++) {
    if (ag_pointer->getAtomicNumber(real_atom_map[centers_list[i]]) !=
        ag_pointer->getAtomicNumber(real_atom_map[other.centers_list[i]])) {
      return false;
    }
    const int my_state_count = atom_centers[centers_list[i]].getStateCount();
    if (my_state_count != atom_centers[other.centers_list[i]].getStateCount()) {
      return false;
    }
    for (int j = 0; j < my_state_count; j++) {
      if (atom_centers[centers_list[i]].getState(j) !=
          atom_centers[other.centers_list[i]].getState(j)) {
        return false;
      }
    }
  }
  
  return true;
}

//-------------------------------------------------------------------------------------------------
int IndigoFragment::cullStatesBearingCharge(const int charge_value) {
  std::vector<bool> keep_state(possible_states);
  for (int i = 0; i < possible_states; i++) {
    keep_state[i] = (net_charges[i] != charge_value);
  }
  int j = 0;
  for (int i = 0; i < possible_states; i++) {
    if (keep_state[i]) {
      states[j] = states[i];
      net_charges[j] = net_charges[i];
      scores[j] = scores[i];
      j++;
    }
  }
  const int n_culled = possible_states - j;
  possible_states = j;
  return n_culled;
}

//-------------------------------------------------------------------------------------------------
int IndigoFragment::cullHigherEnergyStatesByCharge(const int charge_value, const int score_delta) {
  int best_score;
  bool best_score_set = false;
  for (int i = 0; i < possible_states; i++) {
    if (net_charges[i] == charge_value) {
      best_score = (best_score_set) ? std::min(best_score, scores[i]) : scores[i];
      best_score_set = true;
    }
  }
  std::vector<bool> keep_state(possible_states);
  for (int i = 0; i < possible_states; i++) {
    keep_state[i] = (net_charges[i] != charge_value || scores[i] <= best_score + score_delta);
  }
  int j = 0;
  for (int i = 0; i < possible_states; i++) {
    if (keep_state[i]) {
      states[j] = states[i];
      net_charges[j] = net_charges[i];
      scores[j] = scores[i];
      j++;
    }
  }
  const int n_culled = possible_states - j;
  possible_states = j;
  return n_culled;
}

//-------------------------------------------------------------------------------------------------
IndigoTable::IndigoTable(const AtomGraph *ag_in, const int molecule_index,
                         const double temperature_in) :
    atom_count{0}, bond_count{0}, net_charge{0}, fragment_count{0}, temperature{temperature_in},
    ground_state_energy{0.0}, real_atom_map{}, real_bond_map{}, bond_i_atoms{}, bond_j_atoms{},
    atom_connect_bounds{}, atom_relevant_bonds{}, atom_bond_partners{}, atom_scores{},
    bond_scores{}, valence_electrons{}, atom_centers{}, mutable_fragments{},
    ground_state_formal_charges{}, ground_state_bond_orders{}, ground_state_free_electrons{},
    ag_pointer{ag_in}
{
  // Get the number of real atoms in the molecule.  Filter out virtual sites and raise an error
  // if atoms that the scoring function does not cover are encountered.
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const int2 mol_limits = ag_pointer->getMoleculeLimits(molecule_index);
  std::vector<int> ag_atom_idx_map(cdk.natom, -1);
  int j = 0;
  for (int i = mol_limits.x; i < mol_limits.y; i++) {
    const int iznum = cdk.z_numbers[i];
    if (iznum == 0) {
      continue;
    }
    if (iznum == 1  || iznum == 6  || iznum == 7  || iznum == 8  || iznum == 9  || iznum == 15  ||
        iznum == 16 || iznum == 17 || iznum == 35) {
      ag_atom_idx_map[i] = j;
      real_atom_map.push_back(i);

      // Push a guess as to the valence electron count onto the array for now.  Phosphorus and
      // sulfur atoms with sufficient coordination will have to be updated once the bonding
      // patterns are ascertained.
      const int val_e = (iznum == 1) ? 2 : 8;
      valence_electrons.push_back(val_e);
      j++;
    }
    else {
      const char2 esymb = elemental_symbols[iznum];
      const int res_index = ag_pointer->getResidueIndex(i);
      const int res_number = ag_pointer->getResidueNumber(i);
      rtErr("Element " + std::to_string(esymb.x) + std::to_string(esymb.y) + " is not coevered "
            "by the Indigo scoring function.  Atom " + char4ToString(ag_pointer->getAtomName(i)) +
            " in residue " + char4ToString(ag_pointer->getResidueName(res_index)) + " " +
            std::to_string(res_number) + " of topology " + ag_pointer->getFileName() +
            " cannot be treated with Indigo.", "IndigoAtomCenter");
    }
  }
  atom_count = j;
  
  // Allocate the bond partners bounds with the known number of atoms in this Lewis structure
  atom_connect_bounds.resize(atom_count + 1, 0);

  // Get the net charge of the molecule in question
  net_charge = static_cast<int>(round(sum<double>(ag_in->getPartialCharge<double>(mol_limits.x,
                                                                                  mol_limits.y))));

  // Make a map of the bonds.  Filter down to bonds that are within the relevant atoms.  Record,
  // for each atom in this system, how many bonds it makes to other atoms.
  const ValenceKit<double> vk = ag_pointer->getDoublePrecisionValenceKit();
  std::vector<int> ag_bond_idx_map(vk.nbond, -1);
  j = 0;
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int idg_i_atom = ag_atom_idx_map[vk.bond_i_atoms[pos]];
    const int idg_j_atom = ag_atom_idx_map[vk.bond_j_atoms[pos]];
    if (idg_i_atom >= 0 && idg_j_atom >= 0 &&
        (cdk.z_numbers[vk.bond_i_atoms[pos]] > 1 || cdk.z_numbers[vk.bond_j_atoms[pos]] > 1 ||
         atom_count == 2)) {
      ag_bond_idx_map[pos] = j;
      real_bond_map.push_back(pos);
      bond_i_atoms.push_back(idg_i_atom);
      bond_j_atoms.push_back(idg_j_atom);
      atom_connect_bounds[idg_i_atom] += 1;
      atom_connect_bounds[idg_j_atom] += 1;
      j++;
    }
  }
  bond_count = j;
  
  // Take this opportunity to identify phosphorus and sulfur with sufficient coordination
  for (int i = mol_limits.x; i < mol_limits.y; i++) {
    const int iznum = cdk.z_numbers[i];
    const int it_index = ag_atom_idx_map[i];
    if (it_index == -1) {
      continue;
    }
    if (iznum == 15 && atom_connect_bounds[it_index] >= 3) {
      valence_electrons[it_index] = 10;
    }
    if (iznum == 16 && atom_connect_bounds[it_index] >= 3) {
      valence_electrons[it_index] = 12;
    }
  }

  // Compute the prefix sum over atom 1:2 neighbors to get the bounds array
  prefixSumInPlace<int>(&atom_connect_bounds, PrefixSumType::EXCLUSIVE, "IndigoTable");  

  // Make lists of all bond partners for each atom
  atom_relevant_bonds.resize(atom_connect_bounds[atom_count]);
  atom_bond_partners.resize(atom_connect_bounds[atom_count]);
  std::vector<int> partner_counters(atom_count, 0);
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int idg_i_atom = ag_atom_idx_map[vk.bond_i_atoms[pos]];
    const int idg_j_atom = ag_atom_idx_map[vk.bond_j_atoms[pos]];
    if (idg_i_atom >= 0 && idg_j_atom >= 0 &&
        (cdk.z_numbers[vk.bond_i_atoms[pos]] > 1 || cdk.z_numbers[vk.bond_j_atoms[pos]] > 1 ||
         atom_count == 2)) {
      const int idg_i_offset = atom_connect_bounds[idg_i_atom] + partner_counters[idg_i_atom];
      const int idg_j_offset = atom_connect_bounds[idg_j_atom] + partner_counters[idg_j_atom];
      atom_relevant_bonds[idg_i_offset] = ag_bond_idx_map[pos];
      atom_relevant_bonds[idg_j_offset] = ag_bond_idx_map[pos];
      atom_bond_partners[idg_i_offset] = idg_j_atom;
      atom_bond_partners[idg_j_offset] = idg_i_atom;
      partner_counters[idg_i_atom] += 1;
      partner_counters[idg_j_atom] += 1;
    }
  }
  
  // Just get the maps here, rather than make more global instances
  const std::vector<int> fc_scores = indigoFormalChargeScores();
  const std::vector<int> bo_scores = indigoBondOrderScores();

  // Determine the score options for each atom
  std::vector<int> atom_scores(atom_count * formal_charge_range);
  for (int i = 0; i < cdk.natom; i++) {
    const int it_index = ag_atom_idx_map[i];
    if (it_index == -1) {
      continue;
    }
    const int iznum = cdk.z_numbers[i];
    const int ij_base_key = IndigoFormalChargeKey(iznum, minimum_formal_charge);
    for (int j = 0; j < formal_charge_range; j++) {
      atom_scores[(it_index * formal_charge_range) + j] = fc_scores[ij_base_key + j];
    }
  }

  // Determine the score options for each bond
  std::vector<int> bond_scores(bond_count * bond_order_range, INT_MAX);
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int it_index = ag_bond_idx_map[pos];
    if (it_index == -1) {
      continue;
    }
    const int iznum = cdk.z_numbers[vk.bond_i_atoms[pos]];
    const int jznum = cdk.z_numbers[vk.bond_j_atoms[pos]];
    for (int ord = 1; ord <= maximum_bond_order; ord++) {
      const int pos_key = IndigoBondOrderKey(iznum, jznum, ord);
      bond_scores[(it_index * bond_order_range) + ord] = bo_scores[pos_key];
    }
  }
  
  // Create the atom centers array
  for (int i = mol_limits.x; i < mol_limits.y; i++) {
    const int it_index = ag_atom_idx_map[i];
    if (it_index == -1) {
      continue;
    }
    const int llim = atom_connect_bounds[it_index];
    const int hlim = atom_connect_bounds[it_index + 1];
    atom_centers.push_back(IndigoAtomCenter(it_index, cdk.z_numbers[i], hlim - llim,
                                            valence_electrons[it_index],
                                            std::vector<int>(atom_relevant_bonds.begin() + llim,
                                                             atom_relevant_bonds.begin() + hlim),
                                            std::vector<int>(atom_bond_partners.begin() + llim,
                                                             atom_bond_partners.begin() + hlim),
                                            atom_scores, bond_scores));

  }
  
  // Cull states that are incompatible with their neighbors
  int n_culled = 1;
  while (n_culled > 0) {
    n_culled = 0;
    for (int i = 0; i < atom_count; i++) {
      n_culled += atom_centers[i].cullByPartners(atom_centers);
    }
  }

  // Check that all atoms have at least one viable state--if not, abort the construction.
  bool all_atoms_viable = true;
  for (int i = 0; i < atom_count; i++) {
    all_atoms_viable = (all_atoms_viable && atom_centers[i].getStateCount() > 0);
  }
  if (all_atoms_viable == false) {
    ground_state_energy = 0.0;
    ground_state_formal_charges.resize(atom_count, 0.0);
    ground_state_bond_orders.resize(bond_count, 0.0);
    ground_state_free_electrons.resize(atom_count, 0.0);
    return;
  }

  // Sort the states within each atom in increasing order of energy
  for (int i = 0; i < atom_count; i++) {
    atom_centers[i].sortStates();
  }
  
  // Identify fragments with atoms containing more than one possible state
  std::vector<bool> fragment_sweep(atom_count, false);
  std::vector<bool> layer_sweep(atom_count, false);
  for (int i = 0; i < atom_count; i++) {
    if (fragment_sweep[i]) {
      continue;
    }
    fragment_sweep[i] = true;
    if (atom_centers[i].getStateCount() > 1) {
      size_t previous_fragment_size = 0;
      std::vector<int> fragment(1, i);
      while (fragment.size() > previous_fragment_size) {
        previous_fragment_size = fragment.size();
        for (size_t j = 0; j < previous_fragment_size; j++) {
          const int j_atom = fragment[j];
          for (int k = atom_connect_bounds[j_atom]; k < atom_connect_bounds[j_atom + 1]; k++) {
            const int k_atom = atom_bond_partners[k];
            if ((! fragment_sweep[k_atom]) && atom_centers[k_atom].getStateCount() == 1) {

              // Mark this atom right away to prevent further effort evaluating it 
              fragment_sweep[k_atom] = true;
            }
            else if (! (fragment_sweep[k_atom] || layer_sweep[k_atom])) {
              fragment.push_back(k_atom);
              layer_sweep[k_atom] = true;
            }
          }
        }

        // Mark all atoms added as part of this layer.  No cleaning in the layer_sweep array is
        // needed as those atoms will now be excluded due to their markings in fragment_sweep.
        const size_t current_fragment_size = fragment.size();
        for (size_t j = previous_fragment_size; j < current_fragment_size; j++) {
          fragment_sweep[fragment[j]] = true;
        }
      }

      // Add this fragment to the growing list
      mutable_fragments.push_back(IndigoFragment(fragment, atom_centers));
    }
  }
  fragment_count = mutable_fragments.size();
  
  // Compute the charge of all fixed atom centers
  int target_charge = round(sum<double>(ag_pointer->getPartialCharge<double>(mol_limits.x,
                                                                             mol_limits.y)));
  int baseline_charge = 0;
  for (int i = 0; i < atom_count; i++) {
    if (atom_centers[i].getStateCount() == 1) {
      baseline_charge += atom_centers[i].getCharge(0);
    }
  }

  // CHECK
  if (atom_count == 304) {
    printf("Target / baseline charge = %2d / %2d\n", target_charge, baseline_charge);
  }
  // END CHECK
  
  // Make an initial pass to cull fragment states that exceed the ground state energy for any
  // given charge state.
  std::vector<int> q_options;
  std::vector<int> options_bounds(fragment_count + 1, 0);
  for (int i = 0; i < fragment_count; i++) {
    const std::vector<int> frag_q_states = mutable_fragments[i].getChargeStates();
    const size_t nqs = frag_q_states.size();
    q_options.insert(q_options.end(), frag_q_states.begin(), frag_q_states.end());
    options_bounds[i + 1] = options_bounds[i] + nqs;
    for (size_t j = 0; j < nqs; j++) {
      mutable_fragments[i].cullHigherEnergyStatesByCharge(frag_q_states[j]);
    }
  }

  // Prepare a guardrail, a pair of bounds indicating the maximum and minimum current charge
  // that the previously included fragments must have in order to reach the correct total charge
  // on the molecule.
  std::vector<int2> q_guardrail(fragment_count);
  int ncull;
  do {
    int cumulative_min = baseline_charge - target_charge;
    int cumulative_max = baseline_charge - target_charge;
    for (int i = fragment_count - 1; i >= 0; i--) {
      int fragment_minq = q_options[options_bounds[i]];
      int fragment_maxq = q_options[options_bounds[i]];
      for (int j = options_bounds[i] + 1; j < options_bounds[i + 1]; j++) {
        fragment_minq = std::min(fragment_minq, q_options[j]);
        fragment_maxq = std::max(fragment_maxq, q_options[j]);
      }
      cumulative_min += fragment_minq;
      cumulative_max += fragment_maxq;
      q_guardrail[i].x = cumulative_min;
      q_guardrail[i].y = cumulative_max;
    }

    if (atom_count == 304) {
      printf("Guardrail = [\n");
      for (int i = 0; i < fragment_count; i++) {
        printf("  %4d %4d    ", q_guardrail[i].x, q_guardrail[i].y);
        for (int j = options_bounds[i]; j < options_bounds[i + 1]; j++) {
          printf("%2d ", q_options[j]);
        }
        printf("\n");
      }
      printf("];\n");
    }
    
    // Loop over all fragments
    ncull = 0;
    for (int i = 0; i < fragment_count; i++) {
      for (int j = options_bounds[i]; j < options_bounds[i + 1]; j++) {
        if (q_options[j] < q_guardrail[i].x || q_options[j] > q_guardrail[i].y) {
          ncull += mutable_fragments[i].cullStatesBearingCharge(q_options[j]);

          // CHECK
          if (atom_count == 304) {
            printf("Cull states in fragment %2d bearing charge %2d\n", i, q_options[j]);
          }
          // END CHECK
        }
      }
    }
    if (ncull > 0) {
      q_options.resize(0);
      options_bounds[0] = 0;
      for (int i = 0; i < fragment_count; i++) {
        const std::vector<int> frag_q_states = mutable_fragments[i].getChargeStates();
        const size_t nqs = frag_q_states.size();
        q_options.insert(q_options.end(), frag_q_states.begin(), frag_q_states.end());
        options_bounds[i + 1] = options_bounds[i] + nqs;
      }
    }
  } while (ncull > 0);

  // CHECK
  if (atom_count == 304) {
    printf("Charge states of fragments:\n");
    for (int i = 0; i < fragment_count; i++) {
      printf("    %3d : ", i);
      for (int j = options_bounds[i]; j < options_bounds[i + 1]; j++) {
        printf(" %2d", q_options[j]);
      }
      printf("\n");
    }
    printf("\n");
  }
  // END CHECK
  
  // Begin to form the state based on atom centers that are not part of any mutable fragment.
  ground_state_formal_charges.resize(atom_count, 0.0);
  ground_state_bond_orders.resize(bond_count, 0.0);
  ground_state_free_electrons.resize(atom_count, 0.0);
  double total_score = 0.0;
  for (int i = 0; i < atom_count; i++) {
    if (atom_centers[i].getStateCount() > 1) {
      continue;
    }
    total_score += addToGroundState(atom_centers[i], 0);
  }

  // Obtain Boltzmann-weighted average energies for each fragment's charge state with degenerate
  // solutions.
  std::vector<double> e_options(options_bounds[fragment_count]);
  for (int i = 0; i < fragment_count; i++) {
    for (int j = options_bounds[i]; j < options_bounds[i + 1]; j++) {
      const int ifrag_q = q_options[j];
      const std::vector<int2> frag_states = mutable_fragments[i].getStatesBearingCharge(ifrag_q);
      const int n_degen_state = frag_states.size();
      std::vector<double> frag_state_prob(n_degen_state);
      double frag_state_sum = 0.0;
      double min_nrg;
      for (int k = 0; k < n_degen_state; k++) {

        // Convert to Hartrees, then develop the Boltzmann distribution.  The factor of 200,000
        // is based on the conversion of Indigo scores to integers from their original values in
        // Hartrees (the scores are given to five decimal places, thus multipled by 100,000) and
        // the fact that all bonds are double-counted, thus the atom formal charge scores are also
        // doubled when assembling the scores for atom centers and fragment states.
        const double nrg = static_cast<double>(frag_states[k].x) / 200000.0 *
                           symbols::hartree_to_kcal;
        frag_state_prob[k] = nrg;
        if (k == 0 || nrg < min_nrg) {
          min_nrg = nrg;
        }
      }

      // Work in the relative energy of each state to keep the numerics stable
      for (int k = 0; k < n_degen_state; k++) {
        const double rel_nrg = frag_state_prob[k] - min_nrg;
        frag_state_prob[k] = exp(-rel_nrg / (symbols::boltzmann_constant * temperature));
        frag_state_sum += frag_state_prob[k];
      }
      frag_state_sum = 1.0 / frag_state_sum;
      double weighted_energy = 0.0;
      for (int k = 0; k < n_degen_state; k++) {
        weighted_energy += frag_states[k].x * frag_state_prob[k] * frag_state_sum;
      }
      e_options[j] = weighted_energy;
    }
  }
  
  // Incorporate fragments into the solution, grinding through each viable charge state of each
  // fragment.  As in other parts of this workflow, this is a do... while loop to ensure that it
  // executes at least once.  If there are many degenerate states of the system with equal
  // energies, this will select the final such state.  This is about the final situation where the
  // code can still run into a combinatorial explosion, but even then calculation of each energy
  // is rapid.
  double best_score = total_score;
  for (int i = 0; i < fragment_count; i++) {

    // CHECK
    if (i >= static_cast<int>(options_bounds.size()) || i < 0) {
      printf("Trying to access element %4d of options_bounds, when there are only %4zu "
             "elements.\n", i, options_bounds.size());
    }
    if (options_bounds[i] >= static_cast<int>(e_options.size()) || options_bounds[i] < 0) {
      printf("Trying to access element %4d of e_options, when there are only %4zu elements.\n",
             options_bounds[i], e_options.size());
    }
    // END CHECK

    double fragment_max = (options_bounds[i] < options_bounds[i + 1]) ?
                          e_options[options_bounds[i]] : 0.0;
    for (int j = options_bounds[i] + 1; j < options_bounds[i + 1]; j++) {
      fragment_max = std::max(fragment_max, e_options[j]);
    }
    best_score += fragment_max;
  }
  std::vector<int> max_frag_settings(fragment_count);
  for (int i = 0; i < fragment_count; i++) {
    max_frag_settings[i] = options_bounds[i + 1] - options_bounds[i];
  }
  const int max_frag_setval = options_bounds[fragment_count] - fragment_count;
  std::vector<int> frag_settings(fragment_count, 0);
  std::vector<int> best_settings(fragment_count, 0);

  // CHECK
  llint n_tested = 0LL;
  // END CHECK
  
  if (fragment_count > 0) {
    do {

      // Create the state based on these fragment charge settings
      double trial_score = total_score;
      for (int i = 0; i < fragment_count; i++) {
        trial_score += e_options[options_bounds[i] + frag_settings[i]];
      }
      if (trial_score < best_score) {
        best_score = trial_score;
        for (int i = 0; i < fragment_count; i++) {
          best_settings[i] = frag_settings[i];
        }
      }
    
      // Increment the charge settings for various fragments
      frag_settings[0] += 1;
      int pos = 0;
      while (pos < fragment_count && frag_settings[pos] >= max_frag_settings[pos]) {
        frag_settings[pos] = 0;
        if (pos < fragment_count - 1) {
          frag_settings[pos + 1] += 1;
        }
      }

      // CHECK
      if (atom_count == 304) {
        printf("State = [ ");
        for (int i = 0; i < fragment_count; i++) {
          printf("%d ", frag_settings[i]); 
        }
        printf("];\n");
      }
      // END CHECK
      
    } while (sum<int>(frag_settings) < max_frag_setval);
  }

  // CHECK
  if (atom_count == 304) {
    printf("%lld states were tested.\n", n_tested);
  }
  // END CHECK
  
  // Incorporate fragments into the solution with their best detected charge states.  This
  // requires re-working the energetics for the best state.
  for (int i = 0; i < fragment_count; i++) {
    const int ibest_q = q_options[options_bounds[i] + best_settings[i]];
    const std::vector<int2> best_states = mutable_fragments[i].getStatesBearingCharge(ibest_q);
    const int n_degen_state = best_states.size();
    std::vector<double> best_state_prob(n_degen_state);
    double best_state_sum = 0.0;
    double min_nrg;
    for (int j = 0; j < n_degen_state; j++) {

      // Convert to Hartrees, then develop the Boltzmann distribution.  The factor of 200,000
      // is based on the conversion of Indigo scores to integers from their original values in
      // Hartrees (the scores are given to five decimal places, thus multipled by 100,000) and
      // the fact that all bonds are double-counted, thus the atom formal charge scores are also
      // doubled when assembling the scores for atom centers and fragment states.
      const double nrg = static_cast<double>(best_states[j].x) / 200000.0 *
                         symbols::hartree_to_kcal;
      best_state_prob[j] = nrg;
      if (j == 0 || nrg < min_nrg) {
        min_nrg = nrg;
      }
    }

    // Work in the relative energy of each state to keep the numerics stable
    for (int j = 0; j < n_degen_state; j++) {
      const double rel_nrg = best_state_prob[j] - min_nrg;
      best_state_prob[j] = exp(-rel_nrg / (symbols::boltzmann_constant * temperature));
      best_state_sum += best_state_prob[j];
    }
    best_state_sum = 1.0 / best_state_sum;
    for (int j = 0; j < n_degen_state; j++) {
      const double fs_prob = best_state_prob[j] * best_state_sum;

      // With the probabilities of each state known, evaluate the atomic formal charges and bond
      // orders as a sum of all available states.
      const std::vector<int2> t_stt = mutable_fragments[i].getState(j);
      const int n_center = t_stt.size();
      for (int k = 0; k < n_center; k++) {
        total_score += addToGroundState(atom_centers[t_stt[k].x], t_stt[k].y, j, fs_prob);
      }
    }
  }

  // Divide the total score by two for the ground state energy in Hartrees, then convert to
  // kcal/mol.  See above for the origin of the factor of 200,000.
  total_score /= 200000.0;
  total_score *= 0.5 * symbols::hartree_to_kcal;
  ground_state_energy = total_score;
}

//-------------------------------------------------------------------------------------------------
int IndigoTable::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int IndigoTable::getBondCount() const {
  return bond_count;
}

//-------------------------------------------------------------------------------------------------
int IndigoTable::getNetCharge() const {
  return net_charge;
}

//-------------------------------------------------------------------------------------------------
std::vector<CombineIDp> IndigoTable::getGroundStateFormalCharges() const {
  std::vector<CombineIDp> result(atom_count);
  for (int i = 0; i < atom_count; i++) {
    result[i] = { real_atom_map[i], ground_state_formal_charges[i] };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<CombineIDp> IndigoTable::getGroundStateBondOrders() const {
  std::vector<CombineIDp> result(bond_count);
  for (int i = 0; i < bond_count; i++) {
    result[i] = { real_bond_map[i], ground_state_bond_orders[i] };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<CombineIDp> IndigoTable::getGroundStateFreeElectrons() const {
  std::vector<CombineIDp> result(atom_count);
  for (int i = 0; i < atom_count; i++) {
    result[i] = { real_atom_map[i], ground_state_free_electrons[i] };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double IndigoTable::getGroundStateEnergy() const {
  return ground_state_energy;
}

//-------------------------------------------------------------------------------------------------
double IndigoTable::addToGroundState(const IndigoAtomCenter &ac, const int state_index,
                                     const int accumulation, const double probability) {
  const uint t_stt = ac.getState(state_index);

  // Adding an int to an unsigned int would "promote" the int to unsigned and create a HUGE
  // charge if the formal charge is intended to be negative.  Re-cast the unsigned int from the
  // state mask explicitly to int to avoid this problem.
  const double q_inc = static_cast<double>(static_cast<int>(t_stt & 0xf) + minimum_formal_charge) *
                       probability;
  const int ac_idx = ac.getIndigoTableIndex();
  if (accumulation == 0) {
    ground_state_formal_charges[ac_idx] = q_inc;
  }
  else {
    ground_state_formal_charges[ac_idx] += q_inc;
  }
  const int nbonds = ac.getPartnerCount();
  double valence_e = 0.0;
  for (int i = 0; i < nbonds; i++) {
    const int local_bond_idx = ac.getRelevantBond(i);

    // Re-cast the unsigned integer direct to double here, as no other int value is being added
    const double tbo = static_cast<double>((t_stt >> (2*i + formal_charge_bits)) & 0x3);
    valence_e += tbo * 2.0;
    
    // Check whether this atom center should be the one to dictate the bond order accumulation in
    // the ground state, by being the i atom of the bond.
    if (bond_i_atoms[local_bond_idx] != ac_idx) {
      continue;
    }
    
    const double b_inc = tbo * probability;
    if (accumulation == 0) {
      ground_state_bond_orders[local_bond_idx] = b_inc;
    }
    else {
      ground_state_bond_orders[local_bond_idx] += b_inc;
    }
  }
  const int znum = ag_pointer->getAtomicNumber(real_atom_map[ac_idx]);
  switch (znum) {
  case 1:
    ground_state_free_electrons[ac_idx] += (2.0 - valence_e) * probability;
    break;
  case 6:
  case 7:
  case 8:
  case 9:
  case 17:
  case 35:
    ground_state_free_electrons[ac_idx] += (8.0 - valence_e) * probability;
    break;
  case 15:
    {
      const double target_valence_e = (nbonds >= 3) ? 10.0 : 8.0;
      ground_state_free_electrons[ac_idx] += (target_valence_e - valence_e) * probability;
    }
    break;
  }
  return (static_cast<double>(ac.getScore(state_index)) * probability);
}

} // namespace chemistry
} // namespace omni
