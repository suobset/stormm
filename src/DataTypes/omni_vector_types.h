// -*-c++-*-
#ifndef OMNI_VECTOR_TYPES_H
#define OMNI_VECTOR_TYPES_H

#include <cstddef>
#include <typeinfo>
#include <typeindex>
#include "Reporting/error_format.h"

#ifdef OMNI_USE_HPC
#include <vector_types.h>
#endif

namespace omni {
namespace data_types {

#ifndef OMNI_USE_HPC
// If CUDA is not defined, then some of the CUDA POD vector types will need to be delineated in
// order to run the program in CPU code.
struct int2 {
  int x;
  int y;
};

struct int3 {
  int x;
  int y;
  int z;
};

struct int4 {
  int x;
  int y;
  int z;
  int w;
};

struct uint2 {
  unsigned int x;
  unsigned int y;
};

struct uint3 {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};

struct uint4 {
  unsigned int x;
  unsigned int y;
  unsigned int z;
  unsigned int w;
};

struct float2 {
  float x;
  float y;
};

struct float3 {
  float x;
  float y;
  float z;
};

struct float4 {
  float x;
  float y;
  float z;
  float w;
};

struct longlong2 {
  long long int x;
  long long int y;
};

struct longlong3 {
  long long int x;
  long long int y;
  long long int z;
};

struct longlong4 {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
};

struct ulonglong2 {
  unsigned long long int x;
  unsigned long long int y;
};

struct ulonglong3 {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
};

struct ulonglong4 {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
};

struct double2 {
  double x;
  double y;
};

struct double3 {
  double x;
  double y;
  double z;
};

struct double4 {
  double x;
  double y;
  double z;
  double w;
};

struct char2 {
  char x;
  char y;
};

struct char3 {
  char x;
  char y;
  char z;
};

struct char4 {
  char x;
  char y;
  char z;
  char w;
};

struct uchar2 {
  unsigned char x;
  unsigned char y;
};

struct uchar3 {
  unsigned char x;
  unsigned char y;
  unsigned char z;
};

struct uchar4 {
  unsigned char x;
  unsigned char y;
  unsigned char z;
  unsigned char w;
};

struct short2 {
  short int x;
  short int y;
};

struct short3 {
  short int x;
  short int y;
  short int z;
};

struct short4 {
  short int x;
  short int y;
  short int z;
  short int w;
};

struct ushort2 {
  unsigned short int x;
  unsigned short int y;
};

struct ushort3 {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
};

struct ushort4 {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
  unsigned short int w;
};

#endif

// Type definitions for consistent nomenclature
typedef longlong2 llint2;
typedef longlong3 llint3;
typedef longlong4 llint4;
typedef ulonglong2 ullint2;
typedef ulonglong3 ullint3;
typedef ulonglong4 ullint4;

// CUDA- and HIP-defined vector types corresponding to each of the above POD types
static const size_t int2_type_index = std::type_index(typeid(int2)).hash_code();
static const size_t int3_type_index = std::type_index(typeid(int3)).hash_code();
static const size_t int4_type_index = std::type_index(typeid(int4)).hash_code();
static const size_t double2_type_index = std::type_index(typeid(double2)).hash_code();
static const size_t double3_type_index = std::type_index(typeid(double3)).hash_code();
static const size_t double4_type_index = std::type_index(typeid(double4)).hash_code();
static const size_t float2_type_index = std::type_index(typeid(float2)).hash_code();
static const size_t float3_type_index = std::type_index(typeid(float3)).hash_code();
static const size_t float4_type_index = std::type_index(typeid(float4)).hash_code();
static const size_t char2_type_index = std::type_index(typeid(char2)).hash_code();
static const size_t char3_type_index = std::type_index(typeid(char3)).hash_code();
static const size_t char4_type_index = std::type_index(typeid(char4)).hash_code();
static const size_t uchar2_type_index = std::type_index(typeid(uchar2)).hash_code();
static const size_t uchar3_type_index = std::type_index(typeid(uchar3)).hash_code();
static const size_t uchar4_type_index = std::type_index(typeid(uchar4)).hash_code();
static const size_t uint2_type_index = std::type_index(typeid(uint2)).hash_code();
static const size_t uint3_type_index = std::type_index(typeid(uint3)).hash_code();
static const size_t uint4_type_index = std::type_index(typeid(uint4)).hash_code();
static const size_t longlong2_type_index = std::type_index(typeid(longlong2)).hash_code();
static const size_t longlong3_type_index = std::type_index(typeid(longlong3)).hash_code();
static const size_t longlong4_type_index = std::type_index(typeid(longlong4)).hash_code();
static const size_t ulonglong2_type_index = std::type_index(typeid(ulonglong2)).hash_code();
static const size_t ulonglong3_type_index = std::type_index(typeid(ulonglong3)).hash_code();
static const size_t ulonglong4_type_index = std::type_index(typeid(ulonglong4)).hash_code();
static const size_t short2_type_index = std::type_index(typeid(short2)).hash_code();
static const size_t short3_type_index = std::type_index(typeid(short3)).hash_code();
static const size_t short4_type_index = std::type_index(typeid(short4)).hash_code();
static const size_t ushort2_type_index = std::type_index(typeid(ushort2)).hash_code();
static const size_t ushort3_type_index = std::type_index(typeid(ushort3)).hash_code();
static const size_t ushort4_type_index = std::type_index(typeid(ushort4)).hash_code();

/// \brief Test whether some data type is a recognized HPC (CUDA or HIP) vector type.
template <typename T> bool isHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int2_type_index || ct == int3_type_index || ct == int4_type_index ||
          ct == double2_type_index || ct == double3_type_index || ct == double4_type_index ||
          ct == float2_type_index || ct == float3_type_index || ct == float4_type_index ||
          ct == char2_type_index || ct == char3_type_index || ct == char4_type_index ||
          ct == uchar2_type_index || ct == uchar3_type_index || ct == uchar4_type_index ||
          ct == uint2_type_index || ct == uint3_type_index || ct == uint4_type_index ||
          ct == longlong2_type_index || ct == longlong3_type_index || ct == longlong4_type_index ||
          ct == ulonglong2_type_index || ct == ulonglong3_type_index ||
          ct == ulonglong4_type_index || ct == short2_type_index || ct == short3_type_index ||
          ct == short4_type_index || ct == ushort2_type_index || ct == ushort3_type_index ||
          ct == ushort4_type_index);
}

/// \brief Test whether some data types is a signed integral HPC vector type
template <typename T> bool isSignedIntegralHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int2_type_index || ct == int3_type_index || ct == int4_type_index ||
          ct == char2_type_index || ct == char3_type_index || ct == char4_type_index ||
          ct == longlong2_type_index || ct == longlong3_type_index || ct == longlong4_type_index ||
          ct == short2_type_index || ct == short3_type_index || ct == short4_type_index);
}

/// \brief Test whether some data types is an unsigned integral HPC vector type
template <typename T> bool isUnsignedIntegralHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == uint2_type_index || ct == uint3_type_index || ct == uint4_type_index ||
          ct == uchar2_type_index || ct == uchar3_type_index || ct == uchar4_type_index ||
          ct == ulonglong2_type_index || ct == ulonglong3_type_index ||
          ct == ulonglong4_type_index || ct == ushort2_type_index || ct == ushort3_type_index ||
          ct == ushort4_type_index);
}

/// \brief Test whether some data types is a real number HPC vector type
template <typename T> bool isFloatingPointHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == double2_type_index || ct == double3_type_index || ct == double4_type_index ||
          ct == float2_type_index || ct == float3_type_index || ct == float4_type_index);
}

/// \brief Get the size of an HPC vector, in terms of the number of member variables
template <typename T> int getHpcVectorTypeSize() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int2_type_index || ct == uint2_type_index || ct == double2_type_index ||
      ct == float2_type_index || ct == char2_type_index || ct == uchar2_type_index ||
      ct == longlong2_type_index || ct == ulonglong2_type_index || ct == short2_type_index ||
      ct == ushort2_type_index) {
    return 2;
  }
  else if (ct == int3_type_index || ct == uint3_type_index || ct == double3_type_index ||
           ct == float3_type_index || ct == char3_type_index || ct == uchar3_type_index ||
           ct == longlong3_type_index || ct == ulonglong3_type_index || ct == short3_type_index ||
           ct == ushort3_type_index) {
    return 3;
  }
  else if (ct == int4_type_index || ct == uint4_type_index || ct == double4_type_index ||
           ct == float4_type_index || ct == char4_type_index || ct == uchar4_type_index ||
           ct == longlong4_type_index || ct == ulonglong4_type_index || ct == short4_type_index ||
           ct == ushort4_type_index) {
    return 4;
  }
  else {
    rtErr("Unknown data type " + std::string(std::type_index(typeid(T)).name()) + " encountered.",
          "getHpcVectorTypeSize");
  }
  __builtin_unreachable();
}

/// \brief Produce a platform-independent name by which to identify one of the scalar data types.
template <typename T> std::string getOmniHpcVectorTypeName() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int2_type_index) return "int2";
  else if (ct == int3_type_index) return "int3";
  else if (ct == int4_type_index) return "int4";
  else if (ct == uint2_type_index) return "unsigned_int2";
  else if (ct == uint3_type_index) return "unsigned_int3";
  else if (ct == uint4_type_index) return "unsigned_int4";
  else if (ct == double2_type_index) return "double2";
  else if (ct == double3_type_index) return "double3";
  else if (ct == double4_type_index) return "double4";
  else if (ct == float2_type_index) return "float2";
  else if (ct == float3_type_index) return "float3";
  else if (ct == float4_type_index) return "float4";
  else if (ct == char2_type_index) return "char2";
  else if (ct == char3_type_index) return "char3";
  else if (ct == char4_type_index) return "char4";
  else if (ct == uchar2_type_index) return "unsigned_char2";
  else if (ct == uchar3_type_index) return "unsigned_char3";
  else if (ct == uchar4_type_index) return "unsigned_char4";
  else if (ct == longlong2_type_index) return "long_long_int2";
  else if (ct == longlong3_type_index) return "long_long_int3";
  else if (ct == longlong4_type_index) return "long_long_int4";
  else if (ct == ulonglong2_type_index) return "unsigned_long_long_int2";
  else if (ct == ulonglong3_type_index) return "unsigned_long_long_int3";
  else if (ct == ulonglong4_type_index) return "unsigned_long_long_int4";
  else if (ct == short2_type_index) return "short_int2";
  else if (ct == short3_type_index) return "short_int3";
  else if (ct == short4_type_index) return "short_int4";
  else if (ct == ushort2_type_index) return "unsigned_short_int2";
  else if (ct == ushort3_type_index) return "unsigned_short_int3";
  else if (ct == ushort4_type_index) return "unsigned_short_int4";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
          "HPC vector type.", "getOmniHpcVectorTypeName");
  }
  __builtin_unreachable();
}

} // namespace data_types
} // namespace omni

namespace omni {

// Make the HPC vectorized tuple types available through Omni namespaces, if they were defined
// without an overarching HPC library.
#ifndef OMNI_USE_HPC
using data_types::int2;
using data_types::int3;
using data_types::int4;
using data_types::double2;
using data_types::double3;
using data_types::double4;
using data_types::float2;
using data_types::float3;
using data_types::float4;
using data_types::char2;
using data_types::char3;
using data_types::char4;
using data_types::uchar2;
using data_types::uchar3;
using data_types::uchar4;
using data_types::uint2;
using data_types::uint3;
using data_types::uint4;
using data_types::longlong2;
using data_types::longlong3;
using data_types::longlong4;
using data_types::ulonglong2;
using data_types::ulonglong3;
using data_types::ulonglong4;
using data_types::short2;
using data_types::short3;
using data_types::short4;
using data_types::ushort2;
using data_types::ushort3;
using data_types::ushort4;
#endif

// Alternate type definitions need to be carried through regardless of HPC compilation
using data_types::llint2;
using data_types::llint3;
using data_types::llint4;
using data_types::ullint2;
using data_types::ullint3;
using data_types::ullint4;

// Make type index identifiers available in all Omni namespaces
using data_types::int2_type_index;
using data_types::int3_type_index;
using data_types::int4_type_index;
using data_types::double2_type_index;
using data_types::double3_type_index;
using data_types::double4_type_index;
using data_types::float2_type_index;
using data_types::float3_type_index;
using data_types::float4_type_index;
using data_types::char2_type_index;
using data_types::char3_type_index;
using data_types::char4_type_index;
using data_types::uchar2_type_index;
using data_types::uchar3_type_index;
using data_types::uchar4_type_index;
using data_types::uint2_type_index;
using data_types::uint3_type_index;
using data_types::uint4_type_index;
using data_types::longlong2_type_index;
using data_types::longlong3_type_index;
using data_types::longlong4_type_index;
using data_types::ulonglong2_type_index;
using data_types::ulonglong3_type_index;
using data_types::ulonglong4_type_index;
using data_types::short2_type_index;
using data_types::short3_type_index;
using data_types::short4_type_index;
using data_types::ushort2_type_index;
using data_types::ushort3_type_index;
using data_types::ushort4_type_index;

} // namespace omni

#endif
