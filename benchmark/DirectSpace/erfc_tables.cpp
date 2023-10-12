#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/log_scale_spline.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/report_table.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::double4;
using stormm::data_types::float4;
#endif
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Place a pair of points and separate them by a set distance.  Return the result as a six-element
// vector orders point 1 {x, y, z}, point 2 {x, y, z}.
//
// Arguments:
//   xrs:       The random number generator used to create the initial locations of the points.
//              Modified as it churns through random numbers.
//   r_target:  The target distance between the two points
//-------------------------------------------------------------------------------------------------
std::vector<double> placePointPair(Xoshiro256ppGenerator *xrs, const double r_target) {
  std::vector<double> result = uniformRand(xrs, 6, 10.0);
  addScalarToVector(&result, -5.0);
  double dx = result[3] - result[0];
  double dy = result[4] - result[1];
  double dz = result[5] - result[2];
  double r2 = (dx * dx) + (dy * dy) + (dz * dz);
  bool refresh;
  do {
    refresh = (r2 < 1.0e-4);
    if (refresh) {
      for (int i = 0; i < 6; i++) {
        result[i] = uniformRand(xrs, 10.0) - 5.0;
        dx = result[3] - result[0];
        dy = result[4] - result[1];
        dz = result[5] - result[2];
        r2 = (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
  } while (refresh);

  // Determine the current separation and expand or contract the distance between the points as
  // needed.
  const double midx = 0.5 * (result[0] + result[3]);
  const double midy = 0.5 * (result[1] + result[4]);
  const double midz = 0.5 * (result[2] + result[5]);
  const double r_init = sqrt(r2);
  const double rscale = 0.5 * r_target / r_init;
  result[0] = midx - (dx * rscale);
  result[1] = midy - (dy * rscale);
  result[2] = midz - (dz * rscale);
  result[3] = midx + (dx * rscale);
  result[4] = midy + (dy * rscale);
  result[5] = midz + (dz * rscale);
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test a logarithmic spline table in terms of accuracy over a particular range.  Print the result
// for display and processing in popular matrix algebra packages.
//
// Arguments:
//-------------------------------------------------------------------------------------------------
template <typename T4>
void analyzeTable(const LogSplineForm target_form_in, const double ew_coeff,
                  const double min_analysis_range, const double max_analysis_range,
                  const int mantissa_bits_in, const TableIndexing indexing_method_in,
                  const BasisFunctions basis_set_in, const float indexing_offset_in,
                  const std::string &test_environment, Xoshiro256ppGenerator *xrs,
                  const std::string &output_file_name, StopWatch *timer) {

  // Prepare to write the summary
  std::ofstream foutp = openOutputFile(output_file_name, PrintSituation::APPEND, "open the output "
                                       "for additional test recording");
  const char comment_guard = commentSymbol(OutputSyntax::MATRIX_PKG);
  printProtectedText("Potential form:  " + getEnumerationName(target_form_in) +
                     "\nIndexing method: " + getEnumerationName(indexing_method_in) +
                     "\nBasis functions: " + getEnumerationName(basis_set_in) + "\n",
                     comment_guard, &foutp, 120);

  // Parse the input for the testing environment
  const bool test_natural = strcmpCased(test_environment, "natural", CaseSensitivity::NO);
  
  // Begin building the report file.  The variable name will be useful for reporting timings of
  // each logarithmic spline table construction.
  std::string var_name, unit_str, title_str;
  switch (target_form_in) {
  case LogSplineForm::ELEC_PME_DIRECT:
    var_name += "u_";
    unit_str = "kcal/mol";
    title_str = "PME Energy";
    break;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    var_name += "ux_";
    unit_str = "kcal/mol";
    title_str = "PME Energy (with Exclusion)";
    break;
  case LogSplineForm::DELEC_PME_DIRECT:
    var_name += "du_";
    unit_str = "kcal/mol-A";
    title_str = "PME Force";
    break;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    var_name += "dux_";
    unit_str = "kcal/mol-A";
    title_str = "PME Force (with Exclusion)";
    break;
  case LogSplineForm::CUSTOM:
    break;
  }
  switch (indexing_method_in) {
  case TableIndexing::ARG:
    var_name += "r_";
    title_str += ", Arg";
    break;
  case TableIndexing::SQUARED_ARG:
    var_name += "r2_";
    title_str += ", Sq. Arg";
    break;
  case TableIndexing::ARG_OFFSET:
    var_name += "ro_";
    title_str += ", Arg Offset";
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    var_name += "r2o_";
    title_str += ", Sq. Arg Offset";
    break;
  }
  switch (basis_set_in) {
  case BasisFunctions::MIXED_FRACTIONS:
    var_name += "frac";
    title_str += ", Fraction Series";
    break;
  case BasisFunctions::POLYNOMIAL:
    var_name += "poly";
    title_str += ", Polynomial";
    break;
  }

  // Compute the spline tables.
  const double kcoul = stormm::symbols::charmm_gromacs_bioq;
  const int table_construction = timer->addCategory("Log Spline Table, " + var_name);
  timer->assignTime(0);
  LogScaleSpline<T4> lgsp(target_form_in, ew_coeff, kcoul, mantissa_bits_in, 4096.0, 0.015625,
                          indexing_method_in, basis_set_in, indexing_offset_in);
  timer->assignTime(table_construction);
  const double dscr = 0.00390625;
  const int npts = (max_analysis_range - min_analysis_range) / dscr;
  std::vector<double> rpts(npts), dbl_eval(npts);
  std::vector<double> flt_eval_means(npts), flt_eval_stdev(npts), flt_eval_mue(npts);
  std::vector<double> spl_eval_means(npts), spl_eval_stdev(npts), spl_eval_mue(npts);
  const int ntrials = (test_natural) ? 32 : 1;
  std::vector<double> flt_eval(ntrials), spl_eval(ntrials);
  const double bfac = 2.0 * ew_coeff / sqrt(stormm::symbols::pi);
  const float kcoulf = kcoul;
  const float bfacf = bfac;
  const float ew_coeff_f = ew_coeff;
  for (int i = 0; i < npts; i++) {
    const double r = min_analysis_range + (static_cast<double>(i) * dscr);
    rpts[i] = r;
    const double ewr = ew_coeff * r;
    switch (target_form_in) {
    case LogSplineForm::ELEC_PME_DIRECT:
      dbl_eval[i] = kcoul * erfc(ew_coeff * r) / r;
      break;
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      dbl_eval[i] = kcoul * (erfc(ew_coeff * r) - 1.0) / r;
      break;
    case LogSplineForm::DELEC_PME_DIRECT:
      dbl_eval[i] = -kcoul * ((bfac * exp(-ewr * ewr)) + (erfc(ew_coeff * r) / r)) / (r * r);
      break;
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      dbl_eval[i] = -kcoul * ((bfac * exp(-ewr * ewr)) + ((erfc(ew_coeff * r) - 1.0) / r)) /
                    (r * r);
      break;
    case LogSplineForm::CUSTOM:
      break;
    }

    // For any given value of r, there are many orientations of two particles that could create
    // the proper distance, and do so with different combinations of displacments along the
    // Cartesian axes.  Find a series of these points.
    for (int j = 0; j < ntrials; j++) {

      // Create a pair of points in a random orientation with a distance r between them.  Compute
      // displacements in float in order evaluate the spline as well as the single-precision
      // analytic result.
      float sp_r2;
      if (test_natural) {
        const std::vector<double> coords = placePointPair(xrs, r);
        const float sp_dx = static_cast<float>(coords[3]) - static_cast<float>(coords[0]);
        const float sp_dy = static_cast<float>(coords[4]) - static_cast<float>(coords[1]);
        const float sp_dz = static_cast<float>(coords[5]) - static_cast<float>(coords[2]);    
        sp_r2 = (sp_dx * sp_dx) + (sp_dy * sp_dy) + (sp_dz * sp_dz);
      }
      else {
        sp_r2 = r * r;
      }
      const float sp_r  = sqrtf(sp_r2);
      const float sp_ewr = ew_coeff_f * sp_r;
      switch (target_form_in) {
      case LogSplineForm::ELEC_PME_DIRECT:
        flt_eval[j] = kcoulf * erfcf(sp_ewr) / sp_r;
        break;
      case LogSplineForm::ELEC_PME_DIRECT_EXCL:
        flt_eval[j] = kcoulf * (erfcf(sp_ewr) - 1.0f) / sp_r;
        break;
      case LogSplineForm::DELEC_PME_DIRECT:
        flt_eval[j] = -kcoulf * ((bfacf * expf(-sp_ewr * sp_ewr)) + (erfcf(sp_ewr) / sp_r)) /
                       sp_r2;
        break;
      case LogSplineForm::DELEC_PME_DIRECT_EXCL:
        flt_eval[j] = -kcoulf * ((bfacf * expf(-sp_ewr * sp_ewr)) +
                                 ((erfcf(sp_ewr) - 1.0f) / sp_r)) / sp_r2;
        break;
      case LogSplineForm::CUSTOM:
        break;
      }
      flt_eval_mue[i] += fabs(flt_eval[j] - dbl_eval[i]);

      // Evaluate the spline based on the way it is indexed.
      switch (indexing_method_in) {
      case TableIndexing::ARG:
        spl_eval[j] = lgsp.evaluate(sp_r);
        break;
      case TableIndexing::SQUARED_ARG:
        spl_eval[j] = lgsp.evaluate(sp_r2);
        break;
      case TableIndexing::ARG_OFFSET:
        spl_eval[j] = lgsp.evaluate(sp_r + indexing_offset_in);
        break;
      case TableIndexing::SQ_ARG_OFFSET:
        spl_eval[j] = lgsp.evaluate(sp_r2 + indexing_offset_in);
        break;
      }
      spl_eval_mue[i] += fabs(spl_eval[j] - dbl_eval[i]);
    }
    if (test_natural) {
      flt_eval_means[i] = mean(flt_eval);
      flt_eval_stdev[i] = variance(flt_eval, VarianceMethod::STANDARD_DEVIATION);
      spl_eval_means[i] = mean(spl_eval);
      spl_eval_stdev[i] = variance(spl_eval, VarianceMethod::STANDARD_DEVIATION);
      flt_eval_mue[i] /= static_cast<double>(ntrials);
      spl_eval_mue[i] /= static_cast<double>(ntrials);
    }
    else {
      flt_eval_means[i] = flt_eval[0];
      flt_eval_stdev[i] = 0.0;
      spl_eval_means[i] = spl_eval[0];
      spl_eval_stdev[i] = 0.0;
    }
  }

  // Smooth the data and contribute to the output.
  const int ave_window = 8;
  const int plot_pts = npts / ave_window;
  const double ave_wt = 1.0 / static_cast<double>(ave_window);
  std::vector<double> test_data(plot_pts * 3);
  for (int i = 0; i < plot_pts; i++) {
    double r_ave = 0.0;
    double spl_mue = 0.0;
    double flt_mue = 0.0;
    for (int j = 0; j < ave_window; j++) {
      r_ave += rpts[(i * ave_window) + j];
      spl_mue += spl_eval_mue[i];
      flt_mue += flt_eval_mue[i];
    }
    r_ave *= ave_wt;
    spl_mue = log10(spl_mue * ave_wt);
    flt_mue = log10(flt_mue * ave_wt);
    test_data[i] = r_ave;
    test_data[i +       plot_pts] = spl_mue;
    test_data[i + (2 * plot_pts)] = flt_mue;
  }

  // One final edit to the variable name, irrelevant to the table of timings
  if (test_natural) {
    var_name += "_natural";
    title_str += ", Natural Process";
  }
  ReportTable test_tab(test_data, { "Distance, A", "Mean Unsigned Error, Spline",
                                    "Mean Unsigned Error, FP32 Analytic" },
                       std::vector<int>(3, 12), var_name, 120);
  test_tab.printTable(&foutp, OutputSyntax::MATRIX_PKG);
  std::string plot_command("figure;\nhold on;\n");
  plot_command += "plot(" + var_name + "(:,1), " + var_name +
                  "(:,3), 'color', [ 0.9 0.3 0.0 ], 'linewidth', 4);";
  plot_command += "plot(" + var_name + "(:,1), " + var_name +
                  "(:,2), 'color', [ 0.1 0.1 0.1 ], 'linewidth', 4);";
  plot_command += "xlabel('Distance, A');\nylabel('log10 Mean Unsigned Error, " + unit_str +
                  "');\n";
  plot_command += "title('" + title_str + "');\n";
  plot_command += "legend('Analytic FP32', 'Spline');\n";
  plot_command += "set(gca, 'fontsize', 36);\n";
  foutp.write(plot_command.data(), plot_command.size());  
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Create a Hybrid object to engage the GPU and absorb any bootup time into "miscellaneous"
  if (oe.getDisplayTimingsOrder()) {
    Hybrid<int> gpu_trigger(1);
    timer.assignTime(0);
  }
  
  // Take in additional inputs
  double cutoff = 10.0;
  double dsum_tol = 1.0e-7;
  int mantissa_bits = 5;
  int igseed = 322029317;
  std::string output_file_name("erfc_table_results.m");
  for (int i = 0; i < argc; i++) {
    if (i < argc - 1 && strcmpCased(argv[i], "-cut", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
          verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC)) {
        cutoff = stod(std::string(argv[i + 1]));
      }
      else {
        rtErr("The -cut optional keyword must be followed by a real number.", "erfc_tables");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-dsum_tol", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
          verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC)) {
        cutoff = stod(std::string(argv[i + 1]));
      }
      else {
        rtErr("The -dsum_tol optional keyword must be followed by a real number.", "erfc_tables");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-mnbits", CaseSensitivity::NO)) {
      bool problem = false;
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        mantissa_bits = stoi(std::string(argv[i + 1]));
      }
      else {
        problem = true;
      }
      if (mantissa_bits < 0) {
        problem = true;
      }
      if (problem) {
        rtErr("The -mnbits optional keyword must be followed by a positive integer.",
              "erfc_tables");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-igseed", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        igseed = stoi(std::string(argv[i + 1]));
      }
      else {
        rtErr("The -igseed optional keyword must be followed by an integer.", "erfc_tables");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-o", CaseSensitivity::YES)) {
      output_file_name = std::string(argv[i + 1]);
    }
  }

  // Initialize the random number generator.
  Xoshiro256ppGenerator xrs(igseed);

  // Initialize the output.
  std::ofstream foutp = openOutputFile(output_file_name, PrintSituation::OVERWRITE, "prime the "
                                       "output file for printing");
  std::string primer("%% Clear the decks\nclear all\nclose all\n\n");
  foutp.write(primer.data(), primer.size());
  foutp.close();
  
  // Compute the Ewald coefficient.
  const double ew_coeff = ewaldCoefficient(cutoff, dsum_tol);
  const std::vector<TableIndexing> tidx_methods = { TableIndexing::SQUARED_ARG,
                                                    TableIndexing::ARG,
                                                    TableIndexing::SQ_ARG_OFFSET,
                                                    TableIndexing::ARG_OFFSET };
  const std::vector<LogSplineForm> lsfrm_methods = { LogSplineForm::DELEC_PME_DIRECT,
                                                     LogSplineForm::DELEC_PME_DIRECT_EXCL };
  for (size_t i = 0; i < tidx_methods.size(); i++) {
    float idx_offset;
    switch (tidx_methods[i]) {
    case TableIndexing::ARG:
    case TableIndexing::SQUARED_ARG:
      idx_offset = 0.0;
      break;
    case TableIndexing::ARG_OFFSET:
    case TableIndexing::SQ_ARG_OFFSET:
      idx_offset = 0.0625;
      break;
    }
    for (size_t j = 0; j < lsfrm_methods.size(); j++) {
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                           tidx_methods[i], BasisFunctions::POLYNOMIAL, idx_offset, "natural",
                           &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                           tidx_methods[i], BasisFunctions::POLYNOMIAL, idx_offset, "lab",
                           &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                           tidx_methods[i], BasisFunctions::MIXED_FRACTIONS, idx_offset, "natural",
                           &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                           tidx_methods[i], BasisFunctions::MIXED_FRACTIONS, idx_offset, "lab",
                           &xrs, output_file_name, &timer);
    }
  }
  
  // Create some cell grids out of periodic systems and compute direct-space interactions.  Test
  // the spline tables in "real-world" applications.
  const double half_cutoff = 0.5 * cutoff;
  const std::vector<std::string> systems = { "trpcage_in_water", "drug_example", "ubiquitin" };
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm(base_top_path, "top", systems, base_crd_path, "inpcrd", systems);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    PhaseSpaceWriter psw = ps.data();
    double cell_a, cell_b, cell_c;
    hessianNormalWidths(psw.invu, &cell_a, &cell_b, &cell_c);
    const int na_cell = cell_a / half_cutoff;
    const int nb_cell = cell_b / half_cutoff;
    const int nc_cell = cell_c / half_cutoff;
    
  }
  
  // Print results
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
}
