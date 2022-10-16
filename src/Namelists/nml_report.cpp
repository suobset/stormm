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
  gettimeofday(&start_date, nullptr);
}

//-------------------------------------------------------------------------------------------------
ReportControls::ReportControls(const TextFile &tf, int *start_line, bool *found_nml,
                               const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  ReportControls(policy_in, wrap)
{
  NamelistEmulator t_nml = reportInput(tf, start_line, found_nml, policy, wrap);
  if (t_nml.getKeywordStatus("syntax") != InputStatus::MISSING) {
    setOutputSyntax(t_nml.getString("syntax"));
  }
  if (t_nml.getKeywordStatus("scope") != InputStatus::MISSING) {
    setOutputScope(t_nml.getString("scope"));
  }
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
    report_syntax = OutputSyntax::MATPLOTLIB;
  }
  else if (strcmpCased(report_layout_in, "matrix_prg", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "matrix_program", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "matrix_algebra", CaseSensitivity::NO)) {
    report_syntax = OutputSyntax::MATRIX_PRG;
  }
  else if (strcmpCased(report_layout_in, "standalone", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "stormm", CaseSensitivity::NO)) {
    report_syntax = OutputSyntax::STANDALONE;
  }
}

} // namespace namelist
} // namespace namelist
