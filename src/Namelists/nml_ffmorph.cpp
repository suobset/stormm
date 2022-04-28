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
  //const int ncmap    = t_nml.getKeywordEntries("cmap");
  const int nattn_14 = t_nml.getKeywordEntries("attenuation");
  const int ncharge  = t_nml.getKeywordEntries("charge");
  const int nljparm  = t_nml.getKeywordEntries("vdw");
  const int nvsite   = t_nml.getKeywordEntries("virtual_site");
  harmonic_bonds.reserve(nbond);
  harmonic_angles.reserve(nangl);
  cosine_dihedrals.reserve(ndihe);
  urey_bradley_angles.reserve(nubrd);
  charmm_impropers.reserve(ncimp);
  //cmap_surfaces.reserve(ncmap);
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
    harmonic_bonds.emplace_back(ParameterKind::BOND, ti, tj);
    if (keq_provided) {
      harmonic_bonds[i].setStiffness(keq);
    }
    if (leq_provided) {
      harmonic_bonds[i].setEquilibrium(leq);
    }
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
                         "sizes.  In all contextes, unspecified parameters will leave any "
                         "existing settings unchanged.  The modified topologies are not written "
                         "back to disk.");
  const std::string bond_help("Modify the parameters of a bond between two atom types.");
  const std::vector<std::string> bond_keys_help = {
    "Atom type for the I atom in the bond (required)",
    "Atom type for the J atom in the bond (required)",
    "Bond stiffness constant, units of kcal/mol-Angstrom",
    "Bond equilibrium constant, units of Angstroms" };
  t_nml.addKeyword(NamelistElement("bond", { "-ti", "-tj", "-k", "-l0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string("") }, DefaultIsObligatory::NO, InputRepeats::YES,
                                   bond_help, bond_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL }));
  const std::string angle_help("Modify the parameters of an angle between three atom types.");
  const std::vector<std::string> angle_keys_help = {
    "Atom type for the I atom in the angle (required)",
    "Atom type for the J atom in the angle (required)",
    "Atom type for the K atom in the angle (required)",
    "Angle stiffness constant, units of kcal/mol-radian",
    "Angle equilibrium constant, units of degrees" };
  t_nml.addKeyword(NamelistElement("angle", { "-ti", "-tj", "-tk", "-k", "-theta0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("") }, DefaultIsObligatory::NO,
                                   InputRepeats::YES, angle_help, angle_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL }));
  const std::string dihedral_help("Modify the parameters of a dihedral interaction between four "
                                  "atom types.");
  const std::vector<std::string> dihedral_keys_help = {
    "Atom type for the I atom in the dihedral angle (required)",
    "Atom type for the J atom in the dihedral angle (required)",
    "Atom type for the K atom in the dihedral angle (required)",
    "Atom type for the L atom in the dihedral angle (required)",
    "Cosine function periodicity", "Dihedral cosine function amplitude (kcal/mol)",
    "Dihedral phase angle (degrees)", "Indicator of whether the parameter set pertains to a "
    "proper or improper torsion" };
  t_nml.addKeyword(NamelistElement("dihedral",
                                   { "-ti", "-tj", "-tk", "-tl", "-n", "-amp", "-phi", "-kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::STRING },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("PROPER") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, dihedral_help,
                                   dihedral_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::REQUIRED }));
  const std::string ubrd_help("Modify the parameters of a Urey-Bradley interaction between three "
                              "atom types (the spring constant applies only between the first and "
                              "third atom types, but the central atom type is essential for "
                              "determining where the parameter should be applied).");
  const std::vector<std::string> ubrd_keys_help = {
    "Atom type for the I atom in the Urey-Bradley interaction (required)",
    "Atom type for the J atom in the Urey-Bradley interaction (required)",
    "Atom type for the K atom in the Urey-Bradley interaction (required)",
    "Urey-Bradley stretching constant (kcal/mol-Angstrom)",
    "Urey-Bradley equilibrium constant (Angstrom)" };
  t_nml.addKeyword(NamelistElement("urey_bradley", { "-ti", "-tj", "-tk", "-k", "-theta0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("") }, DefaultIsObligatory::NO,
                                   InputRepeats::YES, ubrd_help, ubrd_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL }));
  const std::string cimp_help("Modify the parameters of a CHARMM improper dihedral between four "
                              "atom types.");
  const std::vector<std::string> cimp_keys_help = {
    "Atom type for the I atom in the CHARMM improper dihedral (required)",
    "Atom type for the J atom in the CHARMM improper dihedral (required)",
    "Atom type for the K atom in the CHARMM improper dihedral (required)",
    "Atom type for the L atom in the CHARMM improper dihedral (required)",
    "Harmonic stiffness of the restoring force (kcal/mol-radian)",
    "Phase angle for the equilibrium of the angle between two planes (degrees)" };
  t_nml.addKeyword(NamelistElement("charmm_improper", { "-ti", "-tj", "-tk", "-tl", "-k", "-phi" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, cimp_help,
                                   cimp_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL }));
  const std::string attn14_help("Modify the non-bonded scaling parameters of an attenuated 1:4 "
                                "interaction.");
  const std::vector<std::string> attn14_keys_help = {
    "Atom type for the I atom in the CHARMM improper dihedral (required)",
    "Atom type for the J atom in the CHARMM improper dihedral (required)",
    "Atom type for the K atom in the CHARMM improper dihedral (required)",
    "Atom type for the L atom in the CHARMM improper dihedral (required)",
    "Electrostatic non-bonded inverse scaling factor (a factor of two sets electrostatic "
    "interactions to occur at half strength", "Inverse scaling factor for van-der Waals "
    "interactions" };
  t_nml.addKeyword(NamelistElement("attenuation", { "-ti", "-tj", "-tk", "-tl", "-qq", "-vdw" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, attn14_help,
                                   attn14_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL }));
  const std::string charge_help("Modify the properties of a charged atom.");
  const std::vector<std::string> charge_keys_help = {
    "Name of the atom to alter (this is not an atom type, and is required)", "Residue name of the "
    "atom to alter (required)", "Charge value, atomic units (e, the proton charge)" };
  t_nml.addKeyword(NamelistElement("charge", { "-atom", "-resi", "-q" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, charge_help,
                                   charge_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL }));
  const std::string vdw_help("Modify the van-der Waals properties of an atom.");
  const std::vector<std::string> vdw_keys_help = {
    "Atom type to alter (the van-der Waals atom types are one and the same with those used in "
    "valence parameters)", "Sigma (Lennard-Jones radius) value (Angstrom)",
    "Epsilon (Lennard-Jones well depth) value, units of kcal/mol", "Rho (tertiary parameter for "
    "some Lennard-Jones models or Buckingham potentials--the definition will depend on the "
    "context of each topology)"};
  t_nml.addKeyword(NamelistElement("vdw", { "-t", "-sig", "-eps", "-rho" },
                                   { NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, vdw_help,
                                   vdw_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL }));
  const std::string vsite_help("Modify the dimensions of a virtual site frame (the non-bonded "
                               "properties of the virtual site can be altered with the 'charge' "
                               "and 'vdw' keywords--this keyword pertains to the geometry of the "
                               "virtual site amidst its frame atoms)");
  const std::vector<std::string> vsite_keys_help = {
    "Virtual site atom name", "Parent atom name.  This is the real atom to which the virtual site "
    "likely transfers most of the forces it accumulates, although in principle the distinction "
    "between the parent atom and other frame atoms is the way they factor into various chain "
    "rules.  This and other atom names are case sensitive.", "Frame atom 2 name", "Frame atom 3 "
    "name (if applicable to the frame type)", "Frame atom 4 name (if applicable to the frame "
    "type)", "The virtual site frame type, i.e. Flex-2, case-insensitive)", "First frame "
    "dimension (the meaning and units of this and other dimensions depend on the context of the "
    "frame type)", "Second frame dimension", "Third frame dimension" };
  t_nml.addKeyword(NamelistElement("virtual_site", { "-vsatom", "-parent", "-frame2", "-frame3",
                                                     "-frame4", "-ft", "-dim1", "-dim2", "-dim3" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, vsite_help,
                                   vsite_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::REQUIRED, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL }));

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::YES, tf.getLineCount(),
                             found);
  return t_nml;
}
  
} // namespace namelist 
} // namespace omni
