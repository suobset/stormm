#ifndef OMNI_PTX_MACROS_H
#define OMNI_PTX_MACROS_H

#ifdef OMNI_USE_CUDA
#define SHFL_DOWN(a, b)  __shfl_down_sync(0xffffffff, a, b, 32);
#define SHFL_XOR(a, b)   __shfl_xor_sync(0xffffffff, a, b, 32);
#define SHFL_UP(a, b)    __shfl_up_sync(0xffffffff, a, b, 32);
#endif

#ifdef OMNI_USE_HIP
#define SHFL_DOWN(a, b)  __shfl_down_sync(0xffffffffffffffff, a, b, 32);
#define SHFL_XOR(a, b)   __shfl_xor_sync(0xffffffffffffffff, a, b, 32);
#define SHFL_UP(a, b)    __shfl_up_sync(0xffffffffffffffff, a, b, 32);
#endif

#endif
