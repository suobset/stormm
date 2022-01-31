#include <cstdlib>
#include <stdio.h>
#include <time.h>
#ifdef OMNI_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "Cuda/hybrid.h"
#include "Parsing/parse.h"
#include "Parsing/file_listing.h"
#include "Reporting/code_dox.h"
#include "Topology/topology.h"
#include "Random/random.h"
#include "UnitTesting/unit_test.h"

using omni::topology::Prmtop;
using omni::random::RandomNumberGenerator;

// CHECK
#include <ctime>
#include "DataTypes/common_types.h"
#include "DataTypes/omni_vector_types.h"
#include "Reporting/error_format.h"
using omni::cuda::HybridFormat;
using omni::cuda::Hybrid;
using omni::docs::ObjectIdentifier;
using omni::docs::searchFileForObject;
using omni::errors::rt_err;
using omni::errors::rtWarn;
using omni::parse::listFilesInPath;
using omni::parse::SearchStyle;
using omni::testing::Approx;
using omni::testing::ComparisonType;
using omni::ulint;
#ifndef OMNI_USE_CUDA
using omni::double2;
using omni::int4;
#endif
// END CHECK

int main() {

#if 0
  std::string topname("/home/dcerutti/amber20_src/AmberTools/test/mdgx/TIP4P/wat288.top");
  Prmtop P1(topname);
  topname = "/home/dcerutti/amber20_src/test/tip5p/prmtop";
  Prmtop P2(topname);
  topname = "/home/dcerutti/amber20_src/test/chamber/md_engine/dhfr_cmap/prmtop";
  Prmtop P3(topname);
#endif

#ifdef OMNI_USE_CUDA
  Hybrid<double> H1_unified(10240 * 10240, "H1_unified", HybridFormat::UNIFIED);
  Hybrid<double> H2_decoupled(10240 * 10240, "H1_decoupled", HybridFormat::DECOUPLED);
  Hybrid<double> H3_expedited(10240 * 10240, "H1_expedited", HybridFormat::EXPEDITED);
  struct timespec start, stop;
  double accum;
  cudaDeviceSynchronize();
  clock_gettime(CLOCK_REALTIME, &start);
  for (int i = 0; i < 10; i++) {
    H3_expedited.upload();
  }
  cudaDeviceSynchronize();
  clock_gettime(CLOCK_REALTIME, &stop);
  accum = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Expedited upload time: %12.6lf\n", accum / 10.0);
  clock_gettime(CLOCK_REALTIME, &start);
  for (int i = 0; i < 10; i++) {
    H2_decoupled.upload();
  }
  cudaDeviceSynchronize();
  clock_gettime(CLOCK_REALTIME, &stop);
  accum = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) / 1000000000.0;
  printf("Decoupled upload time: %12.6lf\n", accum / 10.0);
  for (int i = 0; i < H3_expedited.size(); i++) {
    H2_decoupled.putHost(static_cast<double>(i) * 0.1, i);
    H3_expedited.putHost(static_cast<double>(i) * 0.1, i);
  }
  H2_decoupled.upload(10, 10);
  H3_expedited.upload(10, 10);
  for (int i = 10; i < 20; i++) {
    H2_decoupled.putHost(static_cast<double>(i) * 0.001, i);
    H3_expedited.putHost(static_cast<double>(i) * 0.001, i);
  }
  H2_decoupled.download(13, 4);
  H3_expedited.download(13, 4);
  cudaDeviceSynchronize();
  printf("\n");
  for (int i = 0; i < 30; i += 5) {
    printf("  %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf\n", H2_decoupled.readHost(i),
	   H2_decoupled.readHost(i + 1), H2_decoupled.readHost(i + 2),
	   H2_decoupled.readHost(i + 3), H2_decoupled.readHost(i + 4));
  }
  printf("\n");
  for (int i = 0; i < 30; i += 5) {
    printf("  %9.5lf %9.5lf %9.5lf %9.5lf %9.5lf\n", H3_expedited.readHost(i),
	   H3_expedited.readHost(i + 1), H3_expedited.readHost(i + 2),
	   H3_expedited.readHost(i + 3), H3_expedited.readHost(i + 4));
  }

#endif

  // Check random number generator
  RandomNumberGenerator my_rngA, my_rngC;
  RandomNumberGenerator my_rngB(891826);
  printf("Randoms = [\n");
  double sum_of_x = 0.0;
  for (int i = 0; i < 50; i++) {
    double x = my_rngC.gaussianRandomNumber();
    printf("  %12.8lf %12.8lf %12.8lf\n", my_rngA.uniformRandomNumber(),
	   my_rngB.uniformRandomNumber(), my_rngC.uniformRandomNumber());
    sum_of_x += x;
  }
  printf("];\n\n");

  // Check real-valued number comparisons
  if (7.5 == Approx(7.5)) {
    printf("Comparison A works!\n");
  }
  if (7.5 != Approx(7.8)) {
    printf("Comparison B works!\n");
  }
  if (7.5 == Approx(7.6, ComparisonType::RELATIVE, 0.05)) {
    printf("Comparison C works!\n");
  }
  if (7.5 != Approx(7.9, ComparisonType::RELATIVE, 0.05)) {
    printf("Comparison D works!\n");
  }
  if (7.5 == Approx(7.6).margin(0.11)) {
    printf("Comparison E works!\n");
  }
  std::vector<double> v1 = {7.5, 8.0, 9.5};
  std::vector<double> v2 = {7.6, 8.0, 9.6};
  if (v1 == Approx(v2, ComparisonType::MEAN_UNSIGNED_ERROR, 0.1)) {
    printf("Comparison F works!\n");
  }
  float f1 = 6.40;
  float f2 = 6.51;
  if (f1 == Approx(f2).margin(0.12)) {
    printf("Comparison G works!\n");
  }
  std::vector<float> fv1 = {3.5, 5.8, 9.2};
  std::vector<float> fv2 = {3.4, 5.8, 9.3};
  if (fv1 == Approx(fv2, ComparisonType::MEAN_UNSIGNED_ERROR, 0.1)) {
    printf("Comparison H works!\n");
  }
  fv2[0] = 4.4;
  if (fv1 != Approx(fv2, ComparisonType::RELATIVE_RMS_ERROR, 0.02)) {
    printf("Comparison I works!\n");
  }
  std::vector<int> iv1 = {35, 58, 92};
  std::vector<int> iv2 = {44, 58, 93};
  if (iv1 != Approx(iv2, ComparisonType::RELATIVE_RMS_ERROR, 0.02)) {
    printf("Comparison J works!\n");
  }

  // Check file parsing
  std::vector<std::string> T3 = listFilesInPath("/Users/davecerutti/Omni/src/",
                                                SearchStyle::RECURSIVE);

  // Check object searching
  ObjectIdentifier myO = searchFileForObject("Hybrid",
                                             "/Users/davecerutti/Omni/src/Cuda/hybrid.h");
  printf("Found %d instances of Hybrid in %s\n", myO.getInstanceCount(),
         myO.getFileName().c_str());

  // CHECK
#ifdef OMNI_USE_CUDA
  cudaDeviceSynchronize();
#endif
  omni::cuda::mem_balance_sheet.printMemoryProfile();
#ifdef OMNI_USE_CUDA
  cudaDeviceSynchronize();
#endif
  rtWarn("This is a strange and long-winded warning.  I don't know why one needs a warning right "
          "here at all, except to show that there is working formatting for a warning.", "main");
  rt_err("This is a very long error message.  So long, in fact, that it cannot easily fit on a "
         "single line.\n\nThat is why I have a parser to help things fit on a screen in a pretty "
         "arrangement.", "NoClass");
  // END CHECK

  return 0;
}
