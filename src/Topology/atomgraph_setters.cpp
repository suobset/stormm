#include <cmath>
#include <cstdio>
#include <climits>
#include "Constants/generalized_born.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"

namespace omni {
namespace topology {

using parse::CaseSensitivity;
using parse::char4ToString;
using parse::stringToChar4;
using parse::strncmpCased;
using namespace generalized_born_defaults;

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const MobilitySetting movement) {

  // Just feed into the next, most general case.  The cost of the unneccessary bounds check is
  // trivial compared to the simplicity this affords in the code.
  modifyAtomMobility(0, atom_count, movement);
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int low_index, const int high_index,
                                   const MobilitySetting movement) {

  // Range check as this will use the pointer
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot change atom mobility "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "modifyAtomMobility");
  }
  const int int_bits = sizeof(int) * 8;
  int* m_ptr = mobile_atoms.data();
  for (int i = low_index; i < high_index; i++) {
    const int access_index = i / int_bits;
    const int bshift = i - (access_index * int_bits);
    const uint orig_mask = static_cast<uint>(m_ptr[access_index]);
    switch (movement) {
    case MobilitySetting::OFF:
      m_ptr[access_index] = (orig_mask & (~(0x1 << bshift)));
      break;
    case MobilitySetting::ON:
      m_ptr[access_index] = (orig_mask | (0x1 << bshift));
      break;
    case MobilitySetting::TOGGLE:
      m_ptr[access_index] = (orig_mask ^ (0x1 << bshift));
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int index, const MobilitySetting movement) {
  const int int_bits = sizeof(int) * 8;
  const int access_index = index / int_bits;
  const int m_val = mobile_atoms.readHost(access_index);
  const int bshift = index - (access_index * int_bits);
  switch (movement) {
  case MobilitySetting::OFF:
    mobile_atoms.putHost(access_index, m_val & (~(0x1 << bshift)));
    break;
  case MobilitySetting::ON:
    mobile_atoms.putHost(access_index, m_val | (0x1 << bshift));
    break;
  case MobilitySetting::TOGGLE:
    mobile_atoms.putHost(access_index, m_val ^ (0x1 << bshift));
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setSource(const std::string &new_source) {
  source = new_source;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                        const double dielectric_in, const double saltcon_in,
                                        const AtomicRadiusSet radii_set,
                                        const ExceptionResponse policy) {
  gb_style = igb_in;

  // Trap GB use cases with Parse radii (these are much older, for Poisson-Boltzmann calculations)
  switch (radii_set) {
  case AtomicRadiusSet::BONDI:
  case AtomicRadiusSet::AMBER6:
  case AtomicRadiusSet::MBONDI:
  case AtomicRadiusSet::MBONDI2:
  case AtomicRadiusSet::MBONDI3:
    break;
  case AtomicRadiusSet::PARSE:

    // The only use of Parse radii is PB calculations, which are currently not supported in OMNI.
    // Nonetheless, the radius set is included for completeness.  Just make sure that no one is
    // using it for Generalized Born calculations.
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Parse radii are not intended for Generalized Born calculations and cannot "
              "therefore be applied to topology " + source + ".", "AtomGraph",
              "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("Parse radii are not intended for Generalized Born calculations and should not "
               "therefore be applied to topology " + source + ".", "AtomGraph",
               "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  case AtomicRadiusSet::NONE:
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      {
        // Radii are not being applied, so check that nonzero radii are at least present
        bool nonzero_radii_found = false;
        const double* radii_ptr = atomic_pb_radii.data();
        for (int i = 0; i < atom_count; i++) {
          nonzero_radii_found = (nonzero_radii_found || radii_ptr[i] > constants::tiny);
        }
        if (nonzero_radii_found == false) {

          // This is probably always an error, but use the policy switch in case someone is
          // trying something clever.
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("No nonzero radii were found in topology file " + source + ", or applied "
                  "when setting a Generalized Born implicit solvent model.", "AtomGraph",
                  "setImplicitSolventModel");
          case ExceptionResponse::WARN:
            rtWarn("No nonzero radii were found in topology file " + source + ", or applied "
                   "when setting a Generalized Born implicit solvent model.  This is likely to "
                   "cause problems later.", "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
      }
    }
  }

  // Check for offending radii, specifically those that might not fit within the bounds of the
  // "neck" GB models.  If such radii are found, recursively call this function and set the
  // appropriate (Bondi or mBondi3) radii and screening parameters.
  if ((gb_style == ImplicitSolventModel::NECK_GB    && radii_set != AtomicRadiusSet::BONDI) ||
      (gb_style == ImplicitSolventModel::NECK_GB_II && radii_set != AtomicRadiusSet::MBONDI3)) {
    const double* radii_ptr = atomic_pb_radii.data();
    bool bad_radii_found = false;
    for (int i = 0; i < atom_count; i++) {
      bad_radii_found = (bad_radii_found || radii_ptr [i] < 1.0 || radii_ptr[i] > 2.0);
    }
    if (bad_radii_found) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
              "Anstroms or larger than 2.0 Angstroms, which were found in the topology.",
              "AtomGraph", "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
               "Anstroms or larger than 2.0 Angstroms, which were found in the topology.",
               "AtomGraph", "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, AtomicRadiusSet::BONDI, policy);
      return;
    }
  }

  double* alpha_ptr = gb_alpha_parameters.data();
  double* beta_ptr  = gb_beta_parameters.data();
  double* gamma_ptr = gb_gamma_parameters.data();
  double* screen_ptr = gb_screening_factors.data();
  const int* znum_ptr = atomic_numbers.data();

  // Set the dielectric constant and salt concentration, if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    dielectric_constant = dielectric_in;
    salt_concentration = saltcon_in;
    break;
  }

  // Set the radius set and screening parameters (a radius set of NONE will leave these
  // parameters unchanged in the topology)
  const int* nb12_bounds_ptr = nb12_exclusion_bounds.data();
  const int* nb12_excl_ptr = nb12_exclusion_list.data();
  const int* nb13_bounds_ptr = nb13_exclusion_bounds.data();
  const int* nb13_excl_ptr = nb13_exclusion_list.data();
  const char4* atom_type_ptr = atom_types.data();
  const char4* atom_name_ptr = atom_names.data();
  const char4* res_name_ptr = residue_names.data();
  const double* mass_ptr = atomic_masses.data();
  const std::vector<int> atom_residue_idx = getResidueIndex();
  for (int i = 0; i < atom_count; i++) {

    // Impart the radius set 
    double atom_rad = atomic_pb_radii.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.2;
        if (nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] > 0) {
          const int bonded_atom_idx = nb12_excl_ptr[nb12_bounds_ptr[i]];
          if (radii_set == AtomicRadiusSet::AMBER6 || radii_set == AtomicRadiusSet::MBONDI) {
            switch (znum_ptr[bonded_atom_idx]) {
            case 1:
              if (strncmpCased(char4ToString(atom_type_ptr[bonded_atom_idx]).c_str(), "HW", 2,
                               CaseSensitivity::NO)) {
                atom_rad = 0.8;
              }
              break;
            case 6:
              atom_rad = 1.3;
              break;
            case 7:
              if (radii_set == AtomicRadiusSet::MBONDI) {
                atom_rad = 1.3;
              }
              break;
            case 8:
            case 16:
              atom_rad = 0.8;
              break;
            default:
              break;
            }
          }
          else if (radii_set == AtomicRadiusSet::MBONDI2 ||
                   radii_set == AtomicRadiusSet::MBONDI3) {
            if (znum_ptr[bonded_atom_idx] == 7) {
              atom_rad = 1.3;
              const char4 atmc4 = atom_name_ptr[i];
              if (radii_set == AtomicRadiusSet::MBONDI3 &&
                  char4ToString(res_name_ptr[atom_residue_idx[bonded_atom_idx]]) == "ARG " &&
                  (atmc4.x == 'H' && (atmc4.y == 'H' || atmc4.y == 'E'))) {
                atom_rad = 1.17;
              }
            }
          }
        }
        else {
          switch (policy) {
          case ExceptionResponse::DIE:
            break;
          case ExceptionResponse::WARN:
            rtErr("Unbonded hydrogen atom " + char4ToString(atom_name_ptr[i]) + " detected in "
                  "residue " + char4ToString(res_name_ptr[atom_residue_idx[i]]) + ".  This "
                  "hydrogen's radius will keep the default Bondi value of 1.2 Angstroms.\n",
                  "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
        break;
      case 6:
        {
          // Detect repartitioning in hydrogen masses when calculating the atom's true mass
          double atomi_mass = mass_ptr[i];
          for (int j = nb12_bounds_ptr[i]; j < nb12_bounds_ptr[i + 1]; j++) {
            const int neighbor_atom_idx = nb12_excl_ptr[j];
            if (znum_ptr[neighbor_atom_idx] == 1 && mass_ptr[neighbor_atom_idx] > 1.008) {
              atomi_mass += mass_ptr[neighbor_atom_idx] - 1.008;
            }
          }

          // Assuming all masses are based on the natural abundances, identify unified carbon
          // atoms as having (mass of C + mass of one or more H).
          const char4 attc4 = atom_type_ptr[i];
          if ((attc4.x == 'C' && attc4.y == '1' && atomi_mass > 13.018) ||
              (attc4.x == 'C' && attc4.y == '2' && atomi_mass > 14.026) ||
              (attc4.x == 'C' && attc4.y == '3' && atomi_mass > 15.034)) {

            // United atom carbon radius
            atom_rad = 2.2;
          }
          else {

            // Standard carbon radius
            atom_rad = 1.7;
          }
        }
        break;
      case 7:
        atom_rad = 1.55;
        break;
      case 8:
        atom_rad = 1.5;
        if (radii_set == AtomicRadiusSet::MBONDI3) {

          // Adjust carboxylic oxygens on proteins (side-chains as well as C-termini)
          const std::string resc4 = char4ToString(res_name_ptr[atom_residue_idx[i]]);
          const std::string atmc4 = char4ToString(atom_name_ptr[i]);
          if (((resc4 == "ASP " || resc4 == "AS4 ") && (atmc4 == "OD1 " || atmc4 == "OD2 ")) ||
              ((resc4 == "GLU " || resc4 == "GL4 ") && (atmc4 == "OE1 " || atmc4 == "OE2 ")) ||
              atmc4 == "OXT ") {
            atom_rad = 1.4;
          }
          if (znum_ptr[i] == 8 && nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] == 1 &&
              znum_ptr[nb12_excl_ptr[nb12_bounds_ptr[i]]] == 6) {
            for (int j = nb13_bounds_ptr[i]; j < nb13_bounds_ptr[i + 1]; j++) {
              if (znum_ptr[nb13_excl_ptr[j]] == 8) {
                atom_rad = 1.4;
              }
            }
          }
        }
        break;
      case 9:
        atom_rad = 1.5;
        break;
      case 14:
        atom_rad = 2.1;
        break;
      case 15:
        atom_rad = 1.85;
        break;
      case 16:
        atom_rad = 1.8;
        break;
      case 17:
        atom_rad = 1.7;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::PARSE:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.00;
        break;
      case 6:
        atom_rad = 1.7;
        break;
      case 7:
        atom_rad = 1.5;
        break;
      case 8:
        atom_rad = 1.4;
        break;
      case 16:
        atom_rad = 1.85;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::NONE:
      break;
    }
    atomic_pb_radii.putHost(atom_rad, i);

    // Impart the GB screening parameters
    double atom_screen = gb_screening_factors.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_screen = 0.85;
        break;
      case 6:
        atom_screen = 0.72;
        break;
      case 7:
        atom_screen = 0.79;
        break;
      case 8:
        atom_screen = 0.85;
        break;
      case 9:
        atom_screen = 0.88;
        break;
      case 15:
        atom_screen = 0.86;
        break;
      case 16:
        atom_screen = 0.96;
        break;
      default:
        atom_screen = 0.8;
        break;
      }
    case AtomicRadiusSet::PARSE:
    case AtomicRadiusSet::NONE:
      break;
    }
    gb_screening_factors.putHost(atom_screen, i);
  }

  // Set alpha, beta, gamma, and screening factors if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = 0.0;
      beta_ptr[i]  = 0.0;
      gamma_ptr[i] = 0.0;
    }
    break;
  case ImplicitSolventModel::OBC_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_i_alpha;
      beta_ptr[i]  = gb_obc_i_beta;
      gamma_ptr[i] = gb_obc_i_gamma;
    }
    break;
  case ImplicitSolventModel::OBC_GB_II:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_ii_alpha;
      beta_ptr[i]  = gb_obc_ii_beta;
      gamma_ptr[i] = gb_obc_ii_gamma;
    }
    break;
  case ImplicitSolventModel::NECK_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_neck_i_alpha;
      beta_ptr[i]  = gb_neck_i_beta;
      gamma_ptr[i] = gb_neck_i_gamma;
      switch(znum_ptr[i]) {
      case 1:
        screen_ptr[i] = gb_neck_i_screen_h;
        break;
      case 6:
        screen_ptr[i] = gb_neck_i_screen_c;
        break;
      case 7:
        screen_ptr[i] = gb_neck_i_screen_n;
        break;
      case 8:
        screen_ptr[i] = gb_neck_i_screen_o;
        break;
      case 16:
        screen_ptr[i] = gb_neck_i_screen_s;
        break;
      default:
        screen_ptr[i] = gb_neck_i_screen_default;
        break;
      }
    }
    break;
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      switch(znum_ptr[i]) {
      case 1:
        alpha_ptr[i]  = gb_neck_ii_alpha_h;
        beta_ptr[i]   = gb_neck_ii_beta_h;
        gamma_ptr[i]  = gb_neck_ii_gamma_h;
        screen_ptr[i] = gb_neck_ii_screen_h;
        break;
      case 6:
        alpha_ptr[i]  = gb_neck_ii_alpha_c;
        beta_ptr[i]   = gb_neck_ii_beta_c;
        gamma_ptr[i]  = gb_neck_ii_gamma_c;
        screen_ptr[i] = gb_neck_ii_screen_c;
        break;
      case 7:
        alpha_ptr[i]  = gb_neck_ii_alpha_n;
        beta_ptr[i]   = gb_neck_ii_beta_n;
        gamma_ptr[i]  = gb_neck_ii_gamma_n;
        screen_ptr[i] = gb_neck_ii_screen_n;
        break;
      case 8:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_o;
        break;
      case 15:
        alpha_ptr[i]  = gb_neck_ii_alpha_p;
        beta_ptr[i]   = gb_neck_ii_beta_p;
        gamma_ptr[i]  = gb_neck_ii_gamma_p;
        screen_ptr[i] = gb_neck_ii_screen_p;
        break;
      case 16:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_s;
        break;
      default:
        alpha_ptr[i]  = 1.0;
        beta_ptr[i]   = 0.8;
        gamma_ptr[i]  = 4.85;
        screen_ptr[i] = 0.5;
        break;
      }
    }
    break;
  }

  // Note the radius set in the topology, if it has indeed changed
  if (pb_radii_set.size() == 0 || radii_set != AtomicRadiusSet::NONE) {
    pb_radii_set = getAtomicRadiusSetName(radii_set);
  }

  // Compute the neck GB indices based on the baseline atomic PB radii.  These values must later
  // be checked against the available table size.
  int* neck_idx_ptr = neck_gb_indices.data();
  double* radii_ptr = atomic_pb_radii.data();
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      neck_idx_ptr[i] = static_cast<int>(((radii_ptr[i] - 1.0) * 20.0) + 0.5);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setWaterResidueName(const char4 new_name) {
  water_residue_name = new_name;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setWaterResidueName(const std::string &new_name) {
  if (new_name.size() > 4) {
    rtErr("The proposed water model name (" + new_name + ") cannot contain more than four "
          "characters.", "AtomGraph", "setWaterResidueName");
  }
  water_residue_name = stringToChar4(new_name);
}

} // namespace topology
} // namespace omni
