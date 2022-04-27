#include <cmath>
#include "Parsing/parse.h"
#include "forcefield_element.h"

namespace omni {
namespace modeling {

using parse::uppercase;

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in) :
    kind{kind_in}, atom_name_i{' ', ' ', ' ', ' '}, atom_name_j{' ', ' ', ' ', ' '},
    atom_name_k{' ', ' ', ' ', ' '}, atom_name_l{' ', ' ', ' ', ' '},
    atom_name_m{' ', ' ', ' ', ' '}, residue_name_i{' ', ' ', ' ', ' '},
    residue_name_j{' ', ' ', ' ', ' '}, residue_name_k{' ', ' ', ' ', ' '},
    residue_name_l{' ', ' ', ' ', ' '}, residue_name_m{' ', ' ', ' ', ' '}, property_a{0.0},
    property_b{0.0}, property_c{0.0}, surface{}, surface_dimension{0},
    frame_type{VirtualSiteKind::NONE}
{}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const double prop_a_in, const double prop_b_in,
                                     const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const double prop_a_in, const double prop_b_in,
                                     const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 residue_i_in,
                                     const char4 residue_j_in, const char4 residue_k_in,
                                     const char4 residue_l_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in, const double prop_a_in,
                                     const double prop_b_in, const double prop_c_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  property_a = prop_a_in;
  property_b = prop_b_in;
  property_c = prop_c_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in,
                                     const std::vector<double> &surface_in) :
  ForceFieldElement(kind_in)
{
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  const size_t nsurf = surface_in.size();
  surface.resize(nsurf);
  for (size_t i = 0; i < nsurf; i++) {
    surface[i] = surface_in[i];
  }
  surface_dimension = round(sqrt(static_cast<double>(nsurf))) + 0.0001;
  if (surface_dimension * surface_dimension != nsurf) {
    rtErr("The provided surface of " + std::to_string(nsurf) + " elements cannot be formed into a "
          "square array.", "ForceFieldElement");
  }
}

//-------------------------------------------------------------------------------------------------
ParameterKind ForceFieldElement::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getNameOfAtom(const char atom_rank) const {
  switch (kind) {
  case ParameterKind::CMAP:
    switch (uppercase(atom_rank)) {
    case 'I':
      return atom_name_i;
    case 'J':
      return atom_name_j;
    case 'K':
      return atom_name_k;
    case 'L':
      return atom_name_l;
    case 'M':
      return atom_name_m;
    }
    break;
  case ParameterKind::CHARGE:
    switch(uppercase(atom_rank)) {
    case 'I':
      return atom_name_i;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a charge parameter.  To get "
            "partial charge information, call getNameOfAtom() with no argument.",
            "ForceFieldElement", "getNameOfAtom");
    }
    break;
  case ParameterKind::VIRTUAL_SITE_FRAME:
    int number_rank;
    switch (uppercase(atom_rank)) {
    case 'I':
      number_rank = 0;
    case 'J':
      number_rank = 1;
      break;
    case 'K':
      number_rank = 2;
      break;
    case 'L':
      number_rank = 3;
      break;
    case 'M':
      number_rank = 4;
      break;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a virtual site frame.",
            "ForceFieldElement", "getNameOfAtom");
    }
    if (number_rank == 0) {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by "
             "getVirtualSiteAtom().", "ForceFieldElement", "getNameOfAtom");
    }
    else {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by getFrameAtom(" +
             std::to_string(number_rank) + ").", "ForceFieldElement", "getNameOfAtom");
    }
    switch (number_rank) {
    case 0:
      return atom_name_i;
    case 1:
      return atom_name_j;
    case 2:
      return atom_name_k;
    case 3:
      return atom_name_l;
    case 4:
      return atom_name_m;
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getParameterKindName(kind) + "\".", "ForceFieldElement", "getNameOfAtom");
  }
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getTypeOfAtom(const char atom_rank) const {
  bool problem = false;
  switch (uppercase(atom_rank)) {
  case 'I':
    switch (kind) {
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
      return atom_name_i;
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'J':
    switch (kind) {
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_j;
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'K':
    switch (kind) {
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_k;
    case ParameterKind::BOND:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'L':
    switch (kind) {
    case ParameterKind::DIHEDRAL:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_l;
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  default:
    problem = true;
    break;
  }
  if (problem) {
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getParameterKindName(kind) + "\".", "ForceFieldElement", "getTypeOfAtom");
  }
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getNameOfResidue(const char atom_rank) const {
  switch (kind) {
  case ParameterKind::CMAP:
    switch (uppercase(atom_rank)) {
    case 'I':
      return residue_name_i;
    case 'J':
      return residue_name_j;
    case 'K':
      return residue_name_k;
    case 'L':
      return residue_name_l;
    case 'M':
      return residue_name_m;
    }
    break;
  case ParameterKind::CHARGE:
    switch(uppercase(atom_rank)) {
    case 'I':
      return residue_name_i;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a charge parameter.  To get "
            "partial charge information, call getNameOfAtom() with no argument.",
            "ForceFieldElement", "getNameOfAtom");
    }
    break;
  case ParameterKind::VIRTUAL_SITE_FRAME:
    int number_rank;
    switch (uppercase(atom_rank)) {
    case 'I':
      number_rank = 0;
    case 'J':
      number_rank = 1;
      break;
    case 'K':
      number_rank = 2;
      break;
    case 'L':
      number_rank = 3;
      break;
    case 'M':
      number_rank = 4;
      break;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a virtual site frame.",
            "ForceFieldElement", "getNameOfResidue");
    }
    if (number_rank == 0) {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by "
             "getVirtualSiteResidue().", "ForceFieldElement", "getNameOfResidue");
    }
    else {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by getFrameAtom(" +
             std::to_string(number_rank) + ").", "ForceFieldElement", "getNameOfResidue");
    }
    switch (number_rank) {
    case 0:
      return residue_name_i;
    case 1:
      return residue_name_j;
    case 2:
      return residue_name_k;
    case 3:
      return residue_name_l;
    case 4:
      return residue_name_m;
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getParameterKindName(kind) + "\".", "ForceFieldElement", "getNameOfResidue");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getCharge() const {
  switch (kind) {
  case ParameterKind::CHARGE:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no charge property.",
          "ForceFieldElement", "getCharge");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getSigma() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no sigma property.",
          "ForceFieldElement", "getSigma");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getEpsilon() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no epsilon property.",
          "ForceFieldElement", "getEpsilon");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getRho() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no rho property.",
          "ForceFieldElement", "getRho");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getStiffnessConstant() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return property_a;
  case ParameterKind::DIHEDRAL:
    rtWarn("A " + getParameterKindName(kind) + "\" has an amplitude, which is more correctly "
           "accessed with getAmplitude(), but the amplitude is equivalent to a stiffness.",
           "ForceFieldElement", "getStiffnessConstant");
    return property_a;
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no stiffness constant.",
          "ForceFieldElement", "getStiffnessConstant");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getEquilibriumConstant() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return property_b;
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no equilibrium "
          "constant.", "ForceFieldElement", "getEquilibriumConstant");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getAmplitude() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no amplitude.",
          "ForceFieldElement", "getAmplitude");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getPhaseAngle() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no phase angle.",
          "ForceFieldElement", "getPhaseAngle");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getPeriodicity() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no periodicity.",
          "ForceFieldElement", "getPeriodicity");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getElectrostaticScaling() const {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no 1:4 scaling factor.",
          "ForceFieldElement", "getElectrostaticScaling");
  }
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getVanDerWaalsScaling() const {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" has no 1:4 scaling factor.",
          "ForceFieldElement", "getVanDerWaalsScaling");
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ForceFieldElement::getSurface() const {
  switch (kind) {
  case ParameterKind::CMAP:
    return surface;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" does not have a "
          "two-dimensional potential surface.", "ForceFieldElement", "getSurface");
  }
}

//-------------------------------------------------------------------------------------------------
int ForceFieldElement::getSurfaceDimension() const {
  switch (kind) {
  case ParameterKind::CMAP:
    return surface_dimension;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" does not have a "
          "two-dimensional potential surface.", "ForceFieldElement", "getSurface");
  }
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKind ForceFieldElement::getVirtualSiteFrameType() const {
  switch (kind) {
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return frame_type;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getParameterKindName(kind) + "\" does not have a "
          "virtual site frame type.", "ForceFieldElement", "getVirtualSiteFrameType");
  }
}

} // namespace modeling
} // namespace omni
