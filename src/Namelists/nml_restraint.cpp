#include <cmath>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Structure/local_arrangement.h"
#include "namelist_element.h"
#include "nml_restraint.h"

namespace stormm {
namespace namelist {

using structure::imageValue;
using structure::ImagingMethod;
  
//-------------------------------------------------------------------------------------------------
RestraintControls::RestraintControls(const ExceptionResponse policy_in) :
    policy{policy_in}, mask_i{std::string("")}, mask_j{std::string("")}, mask_k{std::string("")},
    mask_l{std::string("")}, atom_i{-1}, atom_j{-1}, atom_k{-1}, atom_l{-1}, order{0},
    initiation_step{0}, maturation_step{0}, initial_k2{0.0}, initial_k3{0.0}, initial_r1{0.0},
    initial_r2{0.0}, initial_r3{0.0}, initial_r4{0.0}, mature_k2{0.0}, mature_k3{0.0},
    mature_r1{0.0}, mature_r2{0.0}, mature_r3{0.0}, mature_r4{0.0}, initial_crd{0.0, 0.0, 0.0},
    mature_crd{0.0, 0.0, 0.0}
{}

//-------------------------------------------------------------------------------------------------
RestraintControls::RestraintControls(const TextFile &tf, int *start_line,
                                     const ExceptionResponse policy_in) :
    RestraintControls()
{
  const int starting_line = *start_line + 1;
  NamelistEmulator t_nml = restraintInput(tf, start_line, policy);

  // Take in atom selections
  if (t_nml.getKeywordStatus("mask1") != InputStatus::MISSING) {
    mask_i = t_nml.getStringValue("mask1");
  }
  if (t_nml.getKeywordStatus("mask2") != InputStatus::MISSING) {
    mask_j = t_nml.getStringValue("mask2");
  }
  if (t_nml.getKeywordStatus("mask3") != InputStatus::MISSING) {
    mask_k = t_nml.getStringValue("mask3");
  }
  if (t_nml.getKeywordStatus("mask4") != InputStatus::MISSING) {
    mask_l = t_nml.getStringValue("mask4");
  }
  if (t_nml.getKeywordStatus("iat1") != InputStatus::MISSING) {
    atom_i = t_nml.getIntValue("iat1");
  }
  if (t_nml.getKeywordStatus("iat2") != InputStatus::MISSING) {
    atom_j = t_nml.getIntValue("iat2");
  }
  if (t_nml.getKeywordStatus("iat3") != InputStatus::MISSING) {
    atom_k = t_nml.getIntValue("iat3");
  }
  if (t_nml.getKeywordStatus("iat4") != InputStatus::MISSING) {
    atom_l = t_nml.getIntValue("iat4");
  }

  // Check the validity of the atom selections
  if (mask_i.size() == 0LLU && atom_i < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A restraint must have at least a valid first atom, but no such atom was specified "
            "for the &restraint namelist beginning on line " + std::to_string(starting_line) +
            "of " + tf.getFileName() + ".", "RestraintControls", "validateRestraint");
    case ExceptionResponse::WARN:
      rtWarn("A restraint must have at least a valid first atom, but no such atom was specified.  "
             "This may create problems downstream, or the restraint input beginning on line " +
             std::to_string(starting_line) + " of " + tf.getFileName() + " will be ignored.",
             "RestraintControls", "validateRestraint");
      restraint_is_valid = false;
      return;
    case ExceptionResponse::SILENT:
      restraint_is_valid = false;
      return;
    }
  }  

  // Compute the order of the restraint based on the number of atom specifications.  This will not
  // validate the atom specifications (in particular, that at least the first atom must be
  // specified in some way), leaving that work for the validation private member function.  The
  // actual identities of the atoms are not discernible at this point.  Rather, the order will be
  // inferred based on attempts to specify them.
  order = (atom_i >= 0 || mask_i.size() >= 0LLU) + (atom_j >= 0 || mask_j.size() >= 0LLU) +
          (atom_k >= 0 || mask_k.size() >= 0LLU) + (atom_l >= 0 || mask_l.size() >= 0LLU);

  // Enforce consistency in the means of atom specification
  if (getAtomSpecification() == RestraintAnchoring::MIXED) {
    enforceSpecification();    
  }
    
  // Get time-dependencies
  bool nstep1_found = false;
  bool nstep2_found = false;
  if (t_nml.getKeywordStatus("nstep1") != InputStatus::MISSING) {
    nstep1_found = true;
    initiation_step = t_nml.getIntValue("nstep1");
  }
  if (t_nml.getKeywordStatus("nstep2") != InputStatus::MISSING) {
    nstep2_found = true;
    maturation_step = t_nml.getIntValue("nstep2");
  }
  if (nstep1_found != nstep2_found) {
    rtErr("Both nstep1 and nstep2 must be provided in a &restraint namelist for time-dependent "
          "restraints, or neither may be provided for a time-independent restraint.  To apply a "
          "restraint that takes immediate effect at some time during the simulation, set both "
          "nstep1 and nstep2 to the same non-zero value at line " + std::to_string(starting_line) +
          " of file " + tf.getFileName() + ".", "RestraintControls");
  }

  // Get restraint parameters
  bool r1_found = false;
  bool r2_found = false;
  bool r3_found = false;
  bool r4_found = false;
  bool x0_found = false;
  bool y0_found = false;
  bool z0_found = false;
  bool k2a_found = false;
  bool k3a_found = false;
  bool r1a_found = false;
  bool r2a_found = false;
  bool r3a_found = false;
  bool r4a_found = false;
  bool x0a_found = false;
  bool y0a_found = false;
  bool z0a_found = false;
  if (t_nml.getKeywordStatus("rk2") == InputStatus::USER_SPECIFIED) {
    initial_k2 = t_nml.getRealValue("rk2");
  }
  if (t_nml.getKeywordStatus("rk3") == InputStatus::USER_SPECIFIED) {
    initial_k3 = t_nml.getRealValue("rk3");
  }
  if (t_nml.getKeywordStatus("r1") == InputStatus::USER_SPECIFIED) {
    initial_r1 = t_nml.getRealValue("r1");
    r1_found = true;
  }
  if (t_nml.getKeywordStatus("r2") == InputStatus::USER_SPECIFIED) {
    initial_r2 = t_nml.getRealValue("r2");
    r2_found = true;
  }
  if (t_nml.getKeywordStatus("r3") == InputStatus::USER_SPECIFIED) {
    initial_r3 = t_nml.getRealValue("r3");
    r3_found = true;
  }
  if (t_nml.getKeywordStatus("r4") == InputStatus::USER_SPECIFIED) {
    initial_r4 = t_nml.getRealValue("r4");
    r4_found = true;
  }
  if (t_nml.getKeywordStatus("rx0") == InputStatus::USER_SPECIFIED) {
    initial_crd.x = t_nml.getRealValue("rx0");
    x0_found = true;
  }
  if (t_nml.getKeywordStatus("ry0") == InputStatus::USER_SPECIFIED) {
    initial_crd.y = t_nml.getRealValue("ry0");
    y0_found = true;
  }
  if (t_nml.getKeywordStatus("rz0") == InputStatus::USER_SPECIFIED) {
    initial_crd.z = t_nml.getRealValue("rz0");
    z0_found = true;
  }
  bool time_dependent_parameters_found = false;
  if (t_nml.getKeywordStatus("rk2a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_k2 = t_nml.getRealValue("rk2a");
    k2a_found = true;
  }
  if (t_nml.getKeywordStatus("rk3a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_k3 = t_nml.getRealValue("rk3a");
    k3a_found = true;
  }
  if (t_nml.getKeywordStatus("r1a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_r1 = t_nml.getRealValue("r1a");
    r1a_found = true;
  }
  if (t_nml.getKeywordStatus("r2a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_r2 = t_nml.getRealValue("r2a");
    r2a_found = true;
  }
  if (t_nml.getKeywordStatus("r3a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_r3 = t_nml.getRealValue("r3a");
    r3a_found = true;
  }
  if (t_nml.getKeywordStatus("r4a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_r4 = t_nml.getRealValue("r4a");
    r4a_found = true;
  }
  if (t_nml.getKeywordStatus("rx0a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_crd.x = t_nml.getRealValue("rx0a");
    x0a_found = true;
  }
  if (t_nml.getKeywordStatus("ry0a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_crd.y = t_nml.getRealValue("ry0a");
    y0a_found = true;
  }
  if (t_nml.getKeywordStatus("rz0a") == InputStatus::USER_SPECIFIED) {
    time_dependent_parameters_found = true;
    mature_crd.z = t_nml.getRealValue("rz0a");
    z0a_found = true;
  }
  if (nstep1_found == false && time_dependent_parameters_found) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Time-dependent parameters were found in a &restraint namelist beginning on line " +
            std::to_string(starting_line) + " of " + tf.getFileName() + ".", "RestraintControls");
    case ExceptionResponse::WARN:
      rtWarn("Time-dependent parameters were found in a &restraint namelist beginning on line " +
             std::to_string(starting_line) + " of " + tf.getFileName() + ".  The restraint will "
             "be made static based on its initial values.", "RestraintControls");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (order == 1) {
    if (x0_found == false || y0_found == false || z0_found == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("A positional restraint specified at line " + std::to_string(starting_line) +
              " of " + tf.getFileName() + " must have an x0, y0, and z0 for the target position.",
              "RestraintControls");
      case ExceptionResponse::WARN:
        rtErr("A positional restraint specified at line " + std::to_string(starting_line) +
              " of " + tf.getFileName() + " must have an x0, y0, and z0 for the target position.  "
              "The restraint will be ignored.", "RestraintControls");
        restraint_is_valid = false;
        break;
      case ExceptionResponse::SILENT:
        restraint_is_valid = false;
        break;
      }
    }
    if (nstep1_found && (x0a_found == false || y0a_found == false || z0a_found == false)) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("A time-dependent positional restraint specified at line " +
              std::to_string(starting_line) + " of " + tf.getFileName() + " must have a final "
              "x0a, y0a, and z0a for the target position.", "RestraintControls");
      case ExceptionResponse::WARN:
        rtWarn("A time-dependent positional restraint specified at line " +
               std::to_string(starting_line) + " of " + tf.getFileName() + " must have a final "
               "x0a, y0a, and z0a for the target position.  The restraint will be ignored.",
               "RestraintControls");
        restraint_is_valid = false;
        break;
      case ExceptionResponse::SILENT:
        restraint_is_valid = false;
        break;
      }
    }
  }
  if (order == 2) {
    checkFinalRestraintSettings(nstep1_found, r2a_found, r3a_found, k2a_found, k3a_found,
                                starting_line, tf.getFileName());
    if (r1_found == false) {
      initial_r1 = initial_r2 - 1000.0;
    }
    if (r4_found == false) {
      initial_r4 = initial_r3 + 1000.0;
    }
    if (nstep1_found) {
      if (r1a_found == false) {
        mature_r1 = mature_r2 - 1000.0;
      }
      if (r4_found == false) {
        mature_r4 = mature_r3 + 1000.0;
      }
    }
  }
  if (order == 3) {
    checkFinalRestraintSettings(nstep1_found, r2a_found, r3a_found, k2a_found, k3a_found,
                                starting_line, tf.getFileName());
    initial_r2 = fabs(imageValue(initial_r2, symbols::twopi, ImagingMethod::MINIMUM_IMAGE));
    initial_r3 = fabs(imageValue(initial_r3, symbols::twopi, ImagingMethod::MINIMUM_IMAGE));
    mature_r2  = fabs(imageValue(mature_r2, symbols::twopi, ImagingMethod::MINIMUM_IMAGE));
    mature_r3  = fabs(imageValue(mature_r3, symbols::twopi, ImagingMethod::MINIMUM_IMAGE));
    if (r1_found == false) {
      initial_r1 = 0.0;
    }
    if (r4_found == false) {
      initial_r4 = symbols::pi;
    }
    if (nstep1_found) {
      if (r1a_found == false) {
        mature_r1 = 0.0;
      }
      if (r4_found == false) {
        mature_r4 = symbols::pi;
      }
    }
  }
  if (order == 4) {
    checkFinalRestraintSettings(nstep1_found, r2a_found, r3a_found, k2a_found, k3a_found,
                                starting_line, tf.getFileName());
    initial_r2 = imageValue(initial_r2, symbols::twopi, ImagingMethod::PRIMARY_UNIT_CELL);
    initial_r3 = imageValue(initial_r3, symbols::twopi, ImagingMethod::PRIMARY_UNIT_CELL);
    mature_r2  = imageValue(mature_r2, symbols::twopi, ImagingMethod::PRIMARY_UNIT_CELL);
    mature_r3  = imageValue(mature_r3, symbols::twopi, ImagingMethod::PRIMARY_UNIT_CELL);
    if (r1_found == false) {
      initial_r1 = 0.0;
    }
    if (r4_found == false) {
      initial_r4 = symbols::twopi;
    }
    if (nstep1_found) {
      if (r1a_found == false) {
        mature_r1 = 0.0;
      }
      if (r4_found == false) {
        mature_r4 = symbols::twopi;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
RestraintAnchoring RestraintControls::getAtomSpecification() const {
  const int indexed_order = (atom_i >= 0) + (atom_j >= 0) + (atom_k >= 0) + (atom_l >= 0);
  const int masked_order = (mask_i.size() >= 0LLU) + (mask_j.size() >= 0LLU) +
                           (mask_k.size() >= 0LLU) + (mask_l.size() >= 0LLU);
  const int mixed_order = (atom_i >= 0 || mask_i.size() >= 0LLU) +
                          (atom_j >= 0 || mask_j.size() >= 0LLU) +
                          (atom_k >= 0 || mask_k.size() >= 0LLU) +
                          (atom_l >= 0 || mask_l.size() >= 0LLU);
    
  if (indexed_order == order) {
    return RestraintAnchoring::INDICES;
  }
  else if (masked_order == order) {
    return RestraintAnchoring::ATOMMASK;
  }
  else if (mixed_order == order) {
    return RestraintAnchoring::MIXED;
  }
  else {
    rtErr("The order of this restraint cannot be explained by either atom mask or atom index "
          "specifications.", "RestraintControls", "getAtomSpecification");
  }
  return RestraintAnchoring::UNKNOWN;
}

//-------------------------------------------------------------------------------------------------
void RestraintControls::enforceSpecification() {
  if (getAtomSpecification() == RestraintAnchoring::MIXED) {
    if (atom_i >= 0 && mask_i.size() == 0) {
      mask_i = std::string("@") + std::to_string(atom_i + 1);
      atom_i = -1;
    }
    if (atom_j >= 0 && mask_j.size() == 0) {
      mask_j = std::string("@") + std::to_string(atom_j + 1);
      atom_j = -1;
    }
    if (atom_k >= 0 && mask_k.size() == 0) {
      mask_k = std::string("@") + std::to_string(atom_k + 1);
      atom_k = -1;
    }
    if (atom_l >= 0 && mask_l.size() == 0) {
      mask_l = std::string("@") + std::to_string(atom_l + 1);
      atom_l = -1;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void RestraintControls::checkFinalRestraintSettings(const bool nstep1_found, const bool r2a_found,
                                                    const bool r3a_found, const bool k2a_found,
                                                    const bool k3a_found, const int starting_line,
                                                    const std::string &filename) {
  if (nstep1_found && (r2a_found == false || r3a_found == false ||
                       k2a_found == false || k3a_found == false)) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("A time-dependent restraint specified at line " +
             std::to_string(starting_line) + " of " + filename + " must have final r2a and r3a, "
             "rk2a and rk3a values.  Missing values will be set based on the initial "
             "values--beware that those may have their default values.");
      if (r2a_found == false) {
        mature_r2 = initial_r2;
      }
      if (r3a_found == false) {
        mature_r3 = initial_r3;
      }
      if (k2a_found == false) {
        mature_k2 = initial_k2;
      }
      if (k3a_found == false) {
        mature_k3 = initial_k3;
      }
      break;
    case ExceptionResponse::SILENT:
      if (r2a_found == false) {
        mature_r2 = initial_r2;
      }
      if (r3a_found == false) {
        mature_r3 = initial_r3;
      }
      if (k2a_found == false) {
        mature_k2 = initial_k2;
      }
      if (k3a_found == false) {
        mature_k3 = initial_k3;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint RestraintControls::getRestraint(const AtomGraph *ag,
                                                 const ChemicalFeatures *chemfe,
                                                 const CoordinateFrameReader &cfr) const {
  switch (getAtomSpecification()) {
  case RestraintAnchoring::INDICES:
    return BoundedRestraint(atom_i, atom_j, atom_k, atom_l, ag, initiation_step,
                            maturation_step, initial_k2, initial_k3, initial_r1, initial_r2,
                            initial_r3, initial_r4, mature_k2, mature_k3, mature_r1, mature_r2,
                            mature_r3, mature_r4, initial_crd, mature_crd);
  case RestraintAnchoring::ATOMMASK:
    return BoundedRestraint(mask_i, mask_j, mask_k, mask_l, ag, chemfe, cfr, initiation_step,
                            maturation_step, initial_k2, initial_k3, initial_r1, initial_r2,
                            initial_r3, initial_r4, mature_k2, mature_k3, mature_r1, mature_r2,
                            mature_r3, mature_r4, initial_crd, mature_crd);
  case RestraintAnchoring::MIXED:
  case RestraintAnchoring::UNKNOWN:
    rtErr("The restraint atom encoding could not be determined.", "RestraintControls",
          "getRestraint");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint RestraintControls::getRestraint(const AtomGraph *ag) const {
  switch (getAtomSpecification()) {
  case RestraintAnchoring::INDICES:
    return BoundedRestraint(atom_i, atom_j, atom_k, atom_l, ag, initiation_step,
                            maturation_step, initial_k2, initial_k3, initial_r1, initial_r2,
                            initial_r3, initial_r4, mature_k2, mature_k3, mature_r1, mature_r2,
                            mature_r3, mature_r4, initial_crd, mature_crd);
  case RestraintAnchoring::ATOMMASK:
    rtErr("A restraint cannot be specified with atom masks unless its chemical features and "
          "coordinates are also provided.  An alternative overload of getRestraint() must be "
          "called instead.", "RestraintControls", "getRestraint");
  case RestraintAnchoring::MIXED:
  case RestraintAnchoring::UNKNOWN:
    rtErr("The restraint atom encoding could not be determined.", "RestraintControls",
          "getRestraint");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator restraintInput(const TextFile &tf, int *start_line,
                                const ExceptionResponse policy) {
  NamelistEmulator t_nml("restraint", CaseSensitivity::AUTOMATIC, policy,
                         "Replicates the Amber NMR restraint namelist within STORMM.");
  t_nml.addKeyword(NamelistElement("iat1", NamelistType::INTEGER, "MISSING"));
  t_nml.addKeyword(NamelistElement("iat2", NamelistType::INTEGER, "MISSING"));
  t_nml.addKeyword(NamelistElement("iat3", NamelistType::INTEGER, "MISSING"));
  t_nml.addKeyword(NamelistElement("iat4", NamelistType::INTEGER, "MISSING"));
  t_nml.addKeyword(NamelistElement("nstep1", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("nstep2", NamelistType::INTEGER, "0"));
  t_nml.addKeyword(NamelistElement("r1", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("r2", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r3", NamelistType::REAL, "0.0"));
  t_nml.addKeyword(NamelistElement("r4", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("r1a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("r2a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("r3a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("r4a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rk2", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rk3", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rk2a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rk3a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("mask1", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("mask2", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("mask3", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("mask4", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("rx0", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("ry0", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rz0", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rx0a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("ry0a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("rz0a", NamelistType::REAL, "MISSING"));
  t_nml.addKeyword(NamelistElement("system", NamelistType::STRING, "MISSING"));
  t_nml.addHelp("iat1", "The first atom in the restraint (this or mask_i is required)");
  t_nml.addHelp("iat2", "The second atom in the restraint (this or mask_j is required)");
  t_nml.addHelp("iat3", "The third atom in the restraint (optional, will convert a distance "
                "restraint into an angle restraint if provided)");
  t_nml.addHelp("iat4", "The fourth atom in the restraint (optional, will convert an angle "
                "bending restraint into a dihedral angle restraint if provided)");
  t_nml.addHelp("nstep1", "Step number at which to begin applying the restraint, with "
                "displacement parameters r1, r2, r3, and r4, plus stiffness parameter rk2 and "
                "rk3.");
  t_nml.addHelp("nstep2", "Step number by which to finish applying the restraint, with "
                "displacement parameters r1a, r2a, r3a, and r4a, plus stiffness parameter rk2a "
                "and rk3a (if r[1-4]a and rk[2,3]a are unspecified, they are taken to be the same "
                "as their initial counterparts and the restraint does not scale between nstep1 "
                "and nstep2).");
  t_nml.addHelp("r1", "Leftmost distance at which the left-hand harmonic restraint linearizes");
  t_nml.addHelp("r2", "Middle distance parameter, to the left of which the left-hand harmonic "
                "restraint applies.  The left-hand harmonic restraint evaluates as "
                "k2 * (R - r2)^2 for r1 <= R <= r2.");
  t_nml.addHelp("r3", "Middle distance parameter, to the left of which the potential is flat and "
                "to the right of which the right-hand harmonic restraint is applied.  The "
                "potential is flat, force 0, for a distance R such that r2 <= R <= r3.");
  t_nml.addHelp("r4", "Rightmost distance parameter, defining the right-hand harmonic restraint "
                "on the interval r3 <= R <= r4.  To the right of r4 the right-hand harmonic "
                "restraint linearizes.");
  t_nml.addHelp("rk2", "Stiffness constant for the left-hand harmonic restraint.");
  t_nml.addHelp("rk2", "Stiffness constant for the right-hand harmonic restraint.");
  t_nml.addHelp("r1a", "If specified along with distinct nstep1 and nstep2, this specifies "
                "the final value of the left-most restraint distance parameter.");
  t_nml.addHelp("r2a", "Final value of the right-most limit of the left-hand harmonic restraint, "
                "given appropriate nstep values.");
  t_nml.addHelp("r3a", "Final value of the left-most limit of the right-hand harmonic restraint, "
                "given appropriate nstep values.");
  t_nml.addHelp("r4a", "Final value of the right-most limit of the right-hand harmonic restraint "
                "(beyond which it linearizes), given appropriate nstep values.");
  t_nml.addHelp("rk2a", "Final value of the left-hand harmonic restraint stiffness, given "
                "appropriate nstep values.");
  t_nml.addHelp("rk3a", "Final value of the right-hand harmonic restraint stiffness, given "
                "appropriate nstep values.");
  t_nml.addHelp("mask1", "Ambmask string evaluating to the first atom in the restraint (may be "
                "given in place of iatm1, but one of these is required)");
  t_nml.addHelp("mask2", "Ambmask string evaluating to the second atom in the restraint (may be "
                "given in place of iatm2, but one of these is required)");
  t_nml.addHelp("mask3", "Ambmask string evaluating to the third atom in the restraint (may be "
                "given in place of iatm3, and either will convert a distance restraint into an "
                "angle bending restraint)");
  t_nml.addHelp("mask4", "Ambmask string evaluating to the fourth atom in the restraint (may be "
                "given in place of iatm4, and either will convert a distance restraint into an "
                "angle bending restraint)");
  t_nml.addHelp("rx0", "Cartesian X coordinate of the anchor point for a positional restraint.");
  t_nml.addHelp("ry0", "Cartesian Y coordinate of the anchor point for a positional restraint.");
  t_nml.addHelp("rz0", "Cartesian Z coordinate of the anchor point for a positional restraint.");
  t_nml.addHelp("rx0a", "Final value of the Cartesian X coordinate of the anchor point for a "
                "positional restraint, given appropriate nstep values.");
  t_nml.addHelp("ry0a", "Final value of the Cartesian Y coordinate of the anchor point for a "
                "positional restraint, given appropriate nstep values.");
  t_nml.addHelp("rz0a", "Final value of the Cartesian Z coordinate of the anchor point for a "
                "positional restraint, given appropriate nstep values.");
  t_nml.addHelp("system", "The system to which this restraint shall apply.  The value of this "
                "keyword should match one of the labels given to a system with the -sys keyword "
                "of the &files namelist, or use one of the reserved values 'all' or "
                "'all_possible' (both case-insensitive) to apply the restraint to all systems, or "
                "all systems meeting the atom requirements.");
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount());

  return t_nml;
}

} // namespace namelist
} // namespace stormm
