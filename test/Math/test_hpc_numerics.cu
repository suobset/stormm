// -*-c++-*-
#include "../../src/UnitTesting/unit_test.h"

using namespace omni::numerics;
using namespace omni::testing;

#include "../../src/Potential/accumulation.i"

// Constants to guide the accumulation
constexpr int ncycles = 10000;
constexpr int fp_bits = 23;

//-------------------------------------------------------------------------------------------------
__global__ __launch_bounds__(large_block_size, 1) kSumFloatStreams(llint* result) {

}

//-------------------------------------------------------------------------------------------------
void sumFloatStreams(Hybrid *acc_array) {
  const int nslots = acc_array->size();
  llint* acc_host = acc_array->data(HybridTargetLevel::HOST);
  for (int i = 0; i < large_block_size; i++) {
    acc_host[i] = 0LL;
  }
  for (int i = 0; i < 10000; i++) {

  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Test split fixed-precision accumulation");
  section(1);
  Hybrid<llint> acc_array(large_block_size, "Accumulators", HybridFormat::EXPEDITED);
  llint* acc_devc = acc_array.data(HybridTargetLevel::DEVICE);
  kSumFloatStreams<<<1, large_block_size>>>(acc_devc);
  
  
  // Print results
  printTestSummary(oe.getVerbosity());

  return 0;

}
