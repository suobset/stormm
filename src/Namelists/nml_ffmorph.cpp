#include "Parsing/parse.h"
#include "namelist_element.h"
#include "nml_ffmorph.h"

namespace omni {
namespace namelist {

using parse::stringToChar4;

//-------------------------------------------------------------------------------------------------
FFMorphControls::FFMorphControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    harmonic_bonds{}, harmonic_angles{}, cosine_dihedrals{}, urey_bradley_angles{},
    charmm_impropers{}, cmap_surfaces{}, attn14_scalings{}, charge_properties{},
    van_der_waals_properties{}, virtual_sites{}
{}

//-------------------------------------------------------------------------------------------------
FFMorphControls::FFMorphControls(const TextFile &tf, int *start_line, bool *found_nml,
                                 const ExceptionResponse policy_in) :
    FFMorphControls(policy_in)
{
  NamelistEmulator t_nml = ffmorphInput(tf, start_line, found_nml, policy);

  // Load each kind of parameter edit
  const int nbond    = t_nml.getKeywordEntries("bond");
  const int nangl    = t_nml.getKeywordEntries("angle");
  const int ndihe    = t_nml.getKeywordEntries("dihedral");
  const int nubrd    = t_nml.getKeywordEntries("urey_bradley");
  const int ncimp    = t_nml.getKeywordEntries("charmm_improper");
  const int ncmap    = t_nml.getKeywordEntries("cmap");
  const int nattn_14 = t_nml.getKeywordEntries("attenuation");
  const int ncharge  = t_nml.getKeywordEntries("charge");
  const int nljparm  = t_nml.getKeywordEntries("lennard_jones");
  const int nvsite   = t_nml.getKeywordEntries("virtual_site");
  harmonic_bonds.reserve(nbond);
  harmonic_angles.reserve(nangl);
  cosine_dihedrals.reserve(ndihe);
  urey_bradley_angles.reserve(nubrd);
  charmm_impropers.reserve(ncimp);
  cmap_surfaces.reserve(ncmap);
  attn14_scalings.reserve(nattn_14);
  charge_properties.reserve(ncharge);
  van_der_waals_properties.reserve(nljparm);
  virtual_sites.reserve(nvsite);
  const InputStatus stt_missing = InputStatus::MISSING;
  for (int i = 0; i < nbond; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("bond", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("bond", "-tj", i));
    const bool keq_provided = (t_nml.getKeywordStatus("bond", "-k") != stt_missing);
    const bool leq_provided = (t_nml.getKeywordStatus("bond", "-l0") != stt_missing);
    const double keq = (keq_provided) ? t_nml.getRealValue("bond", "-k") : 0.0;
    const double leq = (leq_provided) ? t_nml.getRealValue("bond", "-l0") : 0.0;
    harmonic_bonds.emplace_back(ParameterKind::BOND, ti, tj, keq, leq);
  }
}

//-------------------------------------------------------------------------------------------------
int FFMorphControls::getEditCount(const ParameterKind kind) const {
  switch (kind) {
  case ParameterKind::BOND:
    return harmonic_bonds.size();
  case ParameterKind::ANGLE:
    return harmonic_angles.size();
  case ParameterKind::DIHEDRAL:
    return cosine_dihedrals.size();
  case ParameterKind::UREY_BRADLEY:
    return urey_bradley_angles.size();
  case ParameterKind::CHARMM_IMPROPER:
    return charmm_impropers.size();
  case ParameterKind::CMAP:
    return cmap_surfaces.size();
  case ParameterKind::ATTN_14_SCALE:
    return attn14_scalings.size();
  case ParameterKind::CHARGE:
    return charge_properties.size();
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return van_der_waals_properties.size();
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return virtual_sites.size();
  case ParameterKind::NONE:
    rtErr("A valid parameter kind must be specified.", "FFMorphControls", "getEditCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement FFMorphControls::getModelEdit(const ParameterKind kind, const int index) const {
  const int nedits = getEditCount(kind);
  if (index >= nedits) {
    rtErr("Index " + std::to_string(index) + " was requested from an array of " +
          std::to_string(nedits) + " stated edits for " + getParameterKindName(kind) +
          " parameters.", "FFMorphControls", "getModelEdit");
  }
  switch (kind) {
  case ParameterKind::BOND:
    return harmonic_bonds[index];
  case ParameterKind::ANGLE:
    return harmonic_angles[index];
  case ParameterKind::DIHEDRAL:
    return cosine_dihedrals[index];
  case ParameterKind::UREY_BRADLEY:
    return urey_bradley_angles[index];
  case ParameterKind::CHARMM_IMPROPER:
    return charmm_impropers[index];
  case ParameterKind::CMAP:
    return cmap_surfaces[index];
  case ParameterKind::ATTN_14_SCALE:
    return attn14_scalings[index];
  case ParameterKind::CHARGE:
    return charge_properties[index];
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return van_der_waals_properties[index];
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return virtual_sites[index];
  case ParameterKind::NONE:

    // This case is trapped in the call to getEditCount() above.
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator ffmorphInput(const TextFile &tf, int *start_line, bool *found,
                              const ExceptionResponse policy) {
  NamelistEmulator t_nml("ffmorph", CaseSensitivity::AUTOMATIC, policy, "Permits user control of "
                         "specific parameters within the topologies at hand.  The control extends "
                         "as far as changing individual parameters' values, but not the atoms to "
                         "which the parameters apply or the insertion of new potential terms.  "
                         "This namelist drives a metamorphosis of the available topologies, but "
                         "not a change in the non-bonded pair list or valence parameter array "
                         "sizes.  The modified topologies are not written back to disk.");
  const std::string bond_help("Alter the parameters of a bond between two atom types.");
  const std::vector<std::string> bond_keys_help = {
    "Atom type for the I atom in the bond (required)",
    "Atom type for the J atom in the bond (required)",
    "Bond stiffness constant (if unspecified, the existing parameter will not change)",
    "Bond equilibrium constant (if unspecified, the existing parameter will not change)" };
  t_nml.addKeyword(NamelistElement("bond", { "-ti", "-tj", "-k", "-l0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string("") }, DefaultIsObligatory::NO, InputRepeats::YES,
                                   bond_help, bond_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL }));
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount(),
                             found);
  return t_nml;
}
  
} // namespace namelist 
} // namespace omni
