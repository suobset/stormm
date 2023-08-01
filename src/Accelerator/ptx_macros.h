#include "copyright.h"

#ifndef STORMM_PTX_MACROS_H
#define STORMM_PTX_MACROS_H

#ifdef STORMM_USE_CUDA
#define SHFL_DOWN(a, b)  __shfl_down_sync(0xffffffff, a, b, 32)
#define SHFL_XOR(a, b)   __shfl_xor_sync(0xffffffff, a, b, 32)
#define SHFL_UP(a, b)    __shfl_up_sync(0xffffffff, a, b, 32)
#define SHFL(a, b)       __shfl_sync(0xffffffff, a, b, 32)
#define SYNCWARP         __syncwarp()

#define WARP_REDUCE_DOWN(var) {                     \
  var += __shfl_down_sync(0xffffffff, var, 16, 32); \
  var += __shfl_down_sync(0xffffffff, var,  8, 32); \
  var += __shfl_down_sync(0xffffffff, var,  4, 32); \
  var += __shfl_down_sync(0xffffffff, var,  2, 32); \
  var += __shfl_down_sync(0xffffffff, var,  1, 32); \
}

#endif

#ifdef STORMM_USE_HIP
#define SHFL_DOWN(a, b)  __shfl_down(0xffffffffffffffff, a, b, 64)
#define SHFL_XOR(a, b)   __shfl_xor(0xffffffffffffffff, a, b, 64)
#define SHFL_UP(a, b)    __shfl_up(0xffffffffffffffff, a, b, 64)
#define SHFL(a, b)       __shfl(0xffffffffffffffff, a, b, 64)
#define SYNCWARP

#define WARP_REDUCE_DOWN(var) {                     \
  var += __shfl_down_sync(0xffffffff, var, 32, 64); \
  var += __shfl_down_sync(0xffffffff, var, 16, 64); \
  var += __shfl_down_sync(0xffffffff, var,  8, 64); \
  var += __shfl_down_sync(0xffffffff, var,  4, 64); \
  var += __shfl_down_sync(0xffffffff, var,  2, 64); \
  var += __shfl_down_sync(0xffffffff, var,  1, 64); \
}

#endif

#endif
