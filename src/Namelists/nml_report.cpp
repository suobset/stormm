#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#  include <windows.h>
#  include <Lmcons.h>
#else
#  include <unistd.h>
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "namelist_enumerators.h"
#include "nml_report.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
ReportControls::ReportControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  policy{policy_in}, report_layout{OutputSyntax::STANDALONE}, report_scope{OutputScope::AVERAGES},
  username{std::string("")}, start_date{}, print_walltime_data{true}, reported_quantities{},
  sdf_addons{}
{
  // Detect the time of day
  gettimeofday(&start_date, nullptr);

  // Get the username
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  std::string buffer(UNLEN+1, ' ');
  DWORD buffer_len = UNLEN+1;
  GetUserName(buffer.data(), &buffer_len);
#else
  std::string buffer(64, ' ');
  getlogin_r(buffer.data(), buffer.size());
  buffer.resize(strlen(buffer.data()));
#endif
  username = buffer;
}

//-------------------------------------------------------------------------------------------------
ReportControls::ReportControls(const TextFile &tf, int *start_line, bool *found_nml,
                               const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  ReportControls(policy_in, wrap)
{
  NamelistEmulator t_nml = reportInput(tf, start_line, found_nml, policy, wrap);
  if (t_nml.getKeywordStatus("syntax") != InputStatus::MISSING) {
    setOutputSyntax(t_nml.getStringValue("syntax"));
  }
  if (t_nml.getKeywordStatus("scope") != InputStatus::MISSING) {
    setOutputScope(t_nml.getStringValue("scope"));
  }
  if (t_nml.getKeywordStatus("username") != InputStatus::MISSING) {
    setUsername(t_nml.getStringValue("username"));
  }
  if (t_nml.getKeywordStatus("timings") != InputStatus::MISSING) {
    setWallTimeData(t_nml.getStringValue("timings"));
  }
  const int ndetail = t_nml.getKeywordEntries("detail");
  for (int i = 0; i < ndetail; i++) {
    setReportedQuantities(t_nml.getStringValue("detail", i));
  }
  const int ndata = t_nml.getKeywordEntries("sdf_item");
  for (int i = 0; i < ndata; i++) {
    addDataItem(translateSdfKeywordInput(t_nml, i));
  }

  // Construct the list of state variables that will be presented in the output file (the actual
  // reporting will be contingent on whether these quantities are indeed relevant to the
  // calculation).
  const int n_reported_sv = reported_quantities.size();
  bool specific_mm   = false;
  bool specific_temp = false;
  bool specific_vir  = false;
  for (int i = 0; i < n_reported_sv; i++) {
    switch (reported_quantities[i]) {
    case StateVariable::BOND:
    case StateVariable::ANGLE:
    case StateVariable::PROPER_DIHEDRAL:
    case StateVariable::IMPROPER_DIHEDRAL:
    case StateVariable::UREY_BRADLEY:
    case StateVariable::CHARMM_IMPROPER:
    case StateVariable::CMAP:
    case StateVariable::VDW:
    case StateVariable::VDW_ONE_FOUR:
    case StateVariable::ELECTROSTATIC:
    case StateVariable::ELECTROSTATIC_ONE_FOUR:
    case StateVariable::GENERALIZED_BORN:
    case StateVariable::RESTRAINT:
      specific_mm   = true;
      break;
    case StateVariable::TEMPERATURE_PROTEIN:
    case StateVariable::TEMPERATURE_LIGAND:
    case StateVariable::TEMPERATURE_SOLVENT:
      specific_temp = true;
      break;
    case StateVariable::VIRIAL_11:
    case StateVariable::VIRIAL_12:
    case StateVariable::VIRIAL_22:
    case StateVariable::VIRIAL_13:
    case StateVariable::VIRIAL_23:
    case StateVariable::VIRIAL_33:
      specific_vir  = true;
      break;
    case StateVariable::KINETIC:
    case StateVariable::PRESSURE:
    case StateVariable::VOLUME:
    case StateVariable::TEMPERATURE_ALL:
    case StateVariable::DU_DLAMBDA:
    case StateVariable::POTENTIAL_ENERGY:
    case StateVariable::TOTAL_ENERGY:
    case StateVariable::ALL_STATES:
      break;
    }
  }
  const std::vector<StateVariable> tmp_rq = reported_quantities;
  reported_quantities.resize(0);
  reported_quantities.push_back(StateVariable::TOTAL_ENERGY);
  reported_quantities.push_back(StateVariable::POTENTIAL_ENERGY);
  reported_quantities.push_back(StateVariable::KINETIC);
  reported_quantities.push_back(StateVariable::TEMPERATURE_ALL);
  reported_quantities.push_back(StateVariable::PRESSURE);
  reported_quantities.push_back(StateVariable::VOLUME);
  reported_quantities.push_back(StateVariable::DU_DLAMBDA);
  if (specific_mm == false) {
    reported_quantities.push_back(StateVariable::BOND);
    reported_quantities.push_back(StateVariable::ANGLE);
    reported_quantities.push_back(StateVariable::PROPER_DIHEDRAL);
    reported_quantities.push_back(StateVariable::IMPROPER_DIHEDRAL);
    reported_quantities.push_back(StateVariable::UREY_BRADLEY);
    reported_quantities.push_back(StateVariable::CHARMM_IMPROPER);
    reported_quantities.push_back(StateVariable::CMAP);
    reported_quantities.push_back(StateVariable::VDW);
    reported_quantities.push_back(StateVariable::VDW_ONE_FOUR);
    reported_quantities.push_back(StateVariable::ELECTROSTATIC);
    reported_quantities.push_back(StateVariable::ELECTROSTATIC_ONE_FOUR);
    reported_quantities.push_back(StateVariable::GENERALIZED_BORN);
    reported_quantities.push_back(StateVariable::RESTRAINT);
  }
  reported_quantities.insert(reported_quantities.end(), tmp_rq.begin(), tmp_rq.end());
}

//-------------------------------------------------------------------------------------------------
OutputSyntax ReportControls::getOutputSyntax() const {
  return report_layout;
}

//-------------------------------------------------------------------------------------------------
OutputScope ReportControls::getOutputScope() const {
  return report_scope;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReportControls::getUsername() const {
  return username;
}

//-------------------------------------------------------------------------------------------------
timeval ReportControls::getStartDate() const {
  return start_date;
}

//-------------------------------------------------------------------------------------------------
bool ReportControls::printWallTimeData() const {
  return print_walltime_data;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getReportedQuantityCount() const {
  return reported_quantities.size();
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<StateVariable>& ReportControls::getReportedQuantities() const {
  return reported_quantities;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getSDFileDataRequestCount() const {
  return sdf_addons.size();
}

//-------------------------------------------------------------------------------------------------
const std::vector<MolObjDataRequest>& ReportControls::getSDFileDataRequests() const {
  return sdf_addons;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputSyntax(const OutputSyntax report_layout_in) {
  report_layout = report_layout_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputSyntax(const std::string &report_layout_in) {
  if (strcmpCased(report_layout_in, "matplotlib", CaseSensitivity::NO) ||
      strcmpCased(report_layout_in, "mat_plot_lib", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::MATPLOTLIB;
  }
  else if (strcmpCased(report_layout_in, "matrix_prg", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "matrix_program", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "matrix_algebra", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::MATRIX_PKG;
  }
  else if (strcmpCased(report_layout_in, "standalone", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "stormm", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::STANDALONE;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(report_layout_in + " was not recognized as a valid output syntax.", "ReportControls",
            "setOutputSyntax");
    case ExceptionResponse::WARN:
      rtWarn(report_layout_in + " was not recognized as a valid output syntax.  The " +
             getEnumerationName(OutputSyntax::STANDALONE) + " format will be taken instead.",
             "ReportControls", "setOutputSyntax");
      report_layout = OutputSyntax::STANDALONE;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputScope(const OutputScope report_scope_in) {
  report_scope = report_scope_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputScope(const std::string &report_scope_in) {
  if (strcmpCased(report_scope_in, "average", CaseSensitivity::NO) ||
      strcmpCased(report_scope_in, "averages", CaseSensitivity::NO) ||
      strcmpCased(report_scope_in, "mean", CaseSensitivity::NO) ||
      strcmpCased(report_scope_in, "means", CaseSensitivity::NO)) {
    report_scope = OutputScope::AVERAGES;
  }
  else if (strcmpCased(report_scope_in, "outlier", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "outliers", CaseSensitivity::NO)) {
    report_scope = OutputScope::OUTLIERS;
  }
  else if (strcmpCased(report_scope_in, "cluster_average", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "cluster_averages", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "clusteraverage", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "clusteraverages", CaseSensitivity::NO)) {
    report_scope = OutputScope::CLUSTER_AVERAGES;
  }
  else if (strcmpCased(report_scope_in, "cluster_outlier", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "cluster_outliers", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "clusteroutlier", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "clusteroutliers", CaseSensitivity::NO)) {
    report_scope = OutputScope::CLUSTER_OUTLIERS;
  }
  else if (strcmpCased(report_scope_in, "all", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "full", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "entire", CaseSensitivity::NO) ||
           strcmpCased(report_scope_in, "complete", CaseSensitivity::NO)) {
    report_scope = OutputScope::FULL;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(report_scope_in + " was not recognized as a valid scope for reporting energies and "
            "other system diagnostics.", "ReportControls", "setOutputSyntax");
    case ExceptionResponse::WARN:
      rtWarn(report_scope_in + " was not recognized as a valid output syntax.  The " +
             getEnumerationName(OutputScope::AVERAGES) + " format will be taken instead.",
             "ReportControls", "setOutputScope");
      report_scope = OutputScope::AVERAGES;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setUsername(const std::string &username_in) {
  username = username_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setWallTimeData(const bool preference) {
  print_walltime_data = preference;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setWallTimeData(const std::string &preference) {
  if (strcmpCased(preference, "on") || strcmpCased(preference, "active")) {
    print_walltime_data = true;
  }
  else if (strcmpCased(preference, "off")) {
    print_walltime_data = false;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(preference + " is not a recognized setting for the \"timings\" keyword.  Use ON / "
            "ACTIVE or OFF.", "ReportControls", "setWallTimeData");
    case ExceptionResponse::WARN:
      rtWarn(preference + " is not a recognized setting for the \"timings\" keyword.  Timings "
             "display will remain ON / ACTIVE unless set to OFF.", "ReportControls",
             "setWallTimeData");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const std::string &quantity_in) {
  std::vector<StateVariable> sv_deciphered;
  if (strcmpCased(quantity_in, "bond") || strcmpCased(quantity_in, "bonds")) {
    sv_deciphered.push_back(StateVariable::BOND);
  }
  else if (strcmpCased(quantity_in, "angle") || strcmpCased(quantity_in, "angles")) {
    sv_deciphered.push_back(StateVariable::ANGLE);
    sv_deciphered.push_back(StateVariable::UREY_BRADLEY);
  }
  else if (strcmpCased(quantity_in, "dihedral") || strcmpCased(quantity_in, "torsion") ||
           strcmpCased(quantity_in, "dihedrals") || strcmpCased(quantity_in, "torsions")) {
    sv_deciphered.push_back(StateVariable::PROPER_DIHEDRAL);
    sv_deciphered.push_back(StateVariable::IMPROPER_DIHEDRAL);
    sv_deciphered.push_back(StateVariable::CHARMM_IMPROPER);
  }
  else if (strcmpCased(quantity_in, "proper_dihedral") ||
           strcmpCased(quantity_in, "proper_dihedrals") ||
           strcmpCased(quantity_in, "properdihedral") ||
           strcmpCased(quantity_in, "properdihedrals")) {
    sv_deciphered.push_back(StateVariable::PROPER_DIHEDRAL);    
  }
  else if (strcmpCased(quantity_in, "improper_dihedral") ||
           strcmpCased(quantity_in, "improper_dihedrals") ||
           strcmpCased(quantity_in, "improperdihedral") ||
           strcmpCased(quantity_in, "improperdihedrals")) {
    sv_deciphered.push_back(StateVariable::IMPROPER_DIHEDRAL);
    sv_deciphered.push_back(StateVariable::CHARMM_IMPROPER);
  }
  else if (strcmpCased(quantity_in, "cmap") || strcmpCased(quantity_in, "cmaps")) {
    sv_deciphered.push_back(StateVariable::CMAP);
  }
  else if (strcmpCased(quantity_in, "vdw") || strcmpCased(quantity_in, "van_der_waals") ||
           strcmpCased(quantity_in, "vanderwaals") || strcmpCased(quantity_in, "lennardjones") ||
           strcmpCased(quantity_in, "lj")) {
    sv_deciphered.push_back(StateVariable::VDW);
    sv_deciphered.push_back(StateVariable::VDW_ONE_FOUR);
  }
  else if (strcmpCased(quantity_in, "elec") || strcmpCased(quantity_in, "electrostatic")) {
    sv_deciphered.push_back(StateVariable::ELECTROSTATIC);
    sv_deciphered.push_back(StateVariable::ELECTROSTATIC_ONE_FOUR);
  }
  else if (strcmpCased(quantity_in, "gb") || strcmpCased(quantity_in, "generalized_born") ||
           strcmpCased(quantity_in, "solvent") || strcmpCased(quantity_in, "generalizedborn")) {
    sv_deciphered.push_back(StateVariable::GENERALIZED_BORN);
  }
  else if (strcmpCased(quantity_in, "restraint") || strcmpCased(quantity_in, "nmr")) {
    sv_deciphered.push_back(StateVariable::RESTRAINT);
  }
  else if (strcmpCased(quantity_in, "local_temperature")) {
    sv_deciphered.push_back(StateVariable::TEMPERATURE_PROTEIN);
    sv_deciphered.push_back(StateVariable::TEMPERATURE_LIGAND);
    sv_deciphered.push_back(StateVariable::TEMPERATURE_SOLVENT);
  }
  else if (strcmpCased(quantity_in, "virial_components")) {
    sv_deciphered.push_back(StateVariable::VIRIAL_11);
    sv_deciphered.push_back(StateVariable::VIRIAL_12);
    sv_deciphered.push_back(StateVariable::VIRIAL_22);
    sv_deciphered.push_back(StateVariable::VIRIAL_13);
    sv_deciphered.push_back(StateVariable::VIRIAL_23);
    sv_deciphered.push_back(StateVariable::VIRIAL_33);
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(quantity_in + " was not recognized as a molecular mechanics energy component valid "
            "for focused reporting.  Providing no input in this section will have all molecular "
            "mechanics terms be reported.", "ReportControls", "setReportedQuantities");
    case ExceptionResponse::WARN:
      rtWarn(quantity_in + " was not recognized as a molecular mechanics energy component valid "
             "for focused reporting.  Providing no input in this section will have all molecular "
             "mechanics terms be reported.  The present input will be ignored.", "ReportControls",
             "setReportedQuantities");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  setReportedQuantities(sv_deciphered);
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const StateVariable quantities_in) {
  setReportedQuantities(std::vector<StateVariable>(1, quantities_in));
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const std::vector<StateVariable> &quantities_in) {
  const int ns = static_cast<int>(StateVariable::ALL_STATES);
  std::vector<bool> activated(ns, false);
  const int nq = quantities_in.size();
  int n_monitored = 0;
  for (int i = 0LLU; i < nq; i++) {
    const int qno = static_cast<int>(quantities_in[i]);
    if (quantities_in[i] != StateVariable::ALL_STATES && activated[qno] == false) {
      activated[qno] = true;
      n_monitored++;
    }
  }
  reported_quantities.resize(n_monitored);
  int j = 0;
  for (int i = 0LLU; i < ns; i++) {
    if (activated[i]) {
      reported_quantities[j] = static_cast<StateVariable>(i);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::addDataItem(const MolObjDataRequest &ask) {

  // Check that there is no other data item with a conflicting title.
  bool problem = false;
  const int n_items = sdf_addons.size();
  for (int i = 0; i < n_items; i++) {
    problem = (problem || (sdf_addons[i].getTitle() == ask.getTitle()));
  }
  if (problem == false) {
    sdf_addons.push_back(ask);
  }
}

//-------------------------------------------------------------------------------------------------
MolObjDataRequest ReportControls::translateSdfKeywordInput(const NamelistEmulator &t_nml,
                                                           const int index) {
  MolObjDataRequest result;

  return result;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator reportInput(const TextFile &tf, int *start_line, bool *found,
                             const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("report", CaseSensitivity::AUTOMATIC, policy, "Collects directives "
                         "pertaining to the content and layout of the diagnostics output file.  "
                         "This degree of control, while optional, can be important for managing "
                         "the voluminous output that can come from a program which runs multiple "
                         "simulations in a single runtime process.");
  const std::string sdf_help("");
  const std::vector<std::string> sdf_keys_help = {
  };
  t_nml.addKeyword(NamelistElement("syntax", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("scope", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("username", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("timings", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("detail", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("sdf_item", { "-title", "-label", "-detail", "-message",
                                                 "-mask", "-parameter", "-type1", "-type2",
                                                 "-type3", "-type4", "-type5" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sdf_help,
                                   sdf_keys_help,
                                   { SubkeyRequirement::REQUIRED, SubkeyRequirement::REQUIRED,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL, SubkeyRequirement::OPTIONAL,
                                     SubkeyRequirement::OPTIONAL }));
  t_nml.addHelp("syntax", "Layout of the diagnostics report file, intended to make it amenable "
                "to one of a variety of plotting programs for further analysis.");
  t_nml.addHelp("scope", "The extent of reporting that shall take place for the energies and "
                "other properties of individual systems.");
  t_nml.addHelp("username", "Name of the user driving the run (if different from that which would "
                "be detected automatically).");
  t_nml.addHelp("timings", "By default, the wall time devoted to various aspects of a calculation "
                "will be displayed at the end of the run.  Set to ON or ACTIVE to ensure this "
                "behavior, or OFF to decline printed timings.");
  t_nml.addHelp("detail", "In addition to certain obligatory outputs, a user can select that only "
                "specific molecular mechanics energy components be printed, to help focus the "
                "output and reduce file sizes.  In all cases, the total, total potential, and "
                "total kinetic energies will be reported, as will thermodynamic integration "
                "derivatives (if applicable), pressure (if applicable), overall temperature, and "
                "volume.  By default, all molecular mechanics energy components are also "
                "reported, but a user may request specific components and thereby disable the "
                "printing of unrequested terms.  Acceptable arguments include: BOND, ANGLE, "
                "TORSION / DIHEDRAL, PROPER_DIHEDRAL, IMPROPER_DIHEDRAL, CMAP, VDW / LJ, ELEC / "
                "ELECTROSTATIC, GB / SOLVENT, and NMR / RESTRAINT.  Various plurals or omissions "
                "of underscores may also be recognized, and the arguments are not case "
                "sensitive.  A user may also opt to print temperatures of specific regions of the "
                "simulation by supplying the LOCAL_TEMPERATURE argument, and components of the "
                "virial using the VIRIAL_COMPONENTS argument.");
  t_nml.addHelp("sdf_item", "");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  
  return t_nml;
}

} // namespace namelist
} // namespace namelist
