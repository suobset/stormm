#include <cmath>
#include "../../src/DataTypes/common_types.h"
#include "../../src/Reporting/error_format.h"
#include "namelist_element.h"
#include "nml_random.h"

namespace omni {
namespace namelist {

using omni::data_types::uint;
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator randomInput(const TextFile &tf, int *start_line) {
  NamelistEmulator t_nml("random", CaseSensitivity::AUTOMATIC, ExceptionResponse::DIE,
                         "Namelist containing parameters for the random number generator and its "
                         "updates throughout a simulation.");
  t_nml.addKeyword(NamelistElement("igseed", NamelistType::INTEGER,
                                   std::to_string(default_random_seed)));
  t_nml.addKeyword(NamelistElement("igstreams", NamelistType::INTEGER,
                                   std::to_string(default_random_streams)));
  t_nml.addKeyword(NamelistElement("igstride", NamelistType::INTEGER,
                                   std::to_string(default_random_stride)));
  t_nml.addKeyword(NamelistElement("igwarmup", NamelistType::INTEGER,
                                   std::to_string(default_random_warmup)));
  t_nml.addHelp("igseed", "Seed for the first random number stream state vector (additional "
                "streams will be produced from jumps in the sequence based on the first state, "
                "permitting reproducible, parallel random number creation all dependent on a "
                "single stored state.");
  t_nml.addHelp("igstreams", "The number of random number streams with which to produce a cache "
                "of random numbers.  A thousand streams, all based on jumps from some specified "
                "initial state and each producing 32 random numbers, will create a cache of "
                "32,000 deterministic random numbers that is independent of the architecture or "
                "actual thread count in the program.  The default stream count is set high enough "
                "to provide any modern graphics card with plenty of work for all of its threads "
                "to produce random numbers at the rate that the memory bus can store them.");
  t_nml.addHelp("igstride", "Quantity of random numbers to generate from each stream each time "
                "the generator states are hauled out of RAM (or GPU GMEM).  A reasonably large "
                "stride ensures that the random number creation is efficient: the bandwidth of "
                "storing the random number products dwarfs the bandwidth of recalling and then "
                "replacing the random state vectors themselves.  The total random cache size is "
                "the product of igstride x igstreams x sizeof(real number rpresentation).");
  t_nml.addHelp("igwarmup", "Number of cycles to put the first stream's generator through before "
                "spawning separate streams.  If the random seed is especially small, the initial "
                "state may contain a lot of zeros and will therefore not produce high quality "
                "random numbers for the first couple of dozen cycles.  The speed with which "
                "numbers can be generated makes this warmup cost trivial.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.  All calls to this function should
  // proceed in consecutive calls, to make use of the updates to start_line and avoid reading any
  // instance of this namelist twice or skipping instances of it in the search for some other
  // namelist.  An alternative is to keep an independent counter to track progress through the
  // input file in search for &rst namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::NO, tf.getLineCount());

  // Check input
  int warmup_cycles = t_nml.getIntValue("igwarmup");
  validateRandomSeed(t_nml.getIntValue("igseed"), &warmup_cycles);
  
  return t_nml;
}

//-------------------------------------------------------------------------------------------------
void validateRandomSeed(const int igseed, int *warmup_cycles) {
  if (igseed == 0) {
    rtErr("The random seed cannot be zero.", "validateRandomSeed");
  }
  const int nbits = sizeof(int) * 8;
  const uint useed = igseed;
  int nactive = 0;
  for (int i = 0; i < nbits; i++) {
    nactive += ((useed >> i) & 0x1);
  }
  if (*warmup_cycles < ((32 - nactive) * 4)) {
    *warmup_cycles = (32 - nactive) * 4;
    rtWarn("The random seed " + std::to_string(igseed) + " contains " + std::to_string(nactive) +
           " bits set to 1, for which a warmup of at least " + std::to_string((32 - nactive) * 4) +
           " cycles is recommended.  The warmup period has been extended to accommodate the seed.",
           "validateRandomSeed");
  }
}

} // namespace namelist
} // namespace omni
