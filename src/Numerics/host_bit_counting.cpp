#include "copyright.h"
#include "host_bit_counting.h"

namespace stormm {
namespace numerics {

constexpr int byte_bit_counts[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                                    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                                    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                                    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };
                                    
//-------------------------------------------------------------------------------------------------
int hostPopcs(const ushort x) {
  return byte_bit_counts[(x & 0xff)] + byte_bit_counts[((x >> 8) & 0xff)];
}

//-------------------------------------------------------------------------------------------------
int hostPopc(const uint x) {
  return byte_bit_counts[(x & 0xff)        ] + byte_bit_counts[((x >>  8) & 0xff)] +
         byte_bit_counts[((x >> 16) & 0xff)] + byte_bit_counts[((x >> 24) & 0xff)];
}

//-------------------------------------------------------------------------------------------------
int hostPopcll(const ullint x) {
  return byte_bit_counts[(x & 0xff)        ] + byte_bit_counts[((x >>  8) & 0xff)] +
         byte_bit_counts[((x >> 16) & 0xff)] + byte_bit_counts[((x >> 24) & 0xff)] +
         byte_bit_counts[((x >> 32) & 0xff)] + byte_bit_counts[((x >> 40) & 0xff)] +
         byte_bit_counts[((x >> 48) & 0xff)] + byte_bit_counts[((x >> 56) & 0xff)];
}

//-------------------------------------------------------------------------------------------------
int hostFfss(const ushort x) {
  int result = 15;
  ushort test = (0x1 << 15);
  while (result > 0 && (!(x & test))) {
    test >>= 1;
    result--;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int hostFfs(const uint x) {
  int result = 31;
  uint test = (0x1U << 31);
  while (result > 0 && (!(x & test))) {
    test >>= 1;
    result--;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int hostFfsll(const ullint x) {
  int result = 63;
  ullint test = (0x1LLU << 63);
  while (result > 0 && (!(x & test))) {
    test >>= 1;
    result--;
  }
  return result;
}

} // namespace numerics
} // namespace stormm
