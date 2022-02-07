// -*-c++-*-
#ifndef OMNI_COMMON_TYPES_H
#define OMNI_COMMON_TYPES_H

#include <string>
#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "Reporting/error_format.h"

namespace omni {
namespace data_types {

/// \brief Integral type casts
/// \{
typedef unsigned int uint;
typedef unsigned long int ulint;
typedef unsigned long int ulong;
typedef long long int llint;
typedef unsigned long long int ullint;
typedef unsigned char uchar;
/// \}

/// \brief Type indices for common POD types, from which Hybrid objects might be composed
/// \{
static const size_t int_type_index = std::type_index(typeid(int)).hash_code();
static const size_t double_type_index = std::type_index(typeid(double)).hash_code();
static const size_t float_type_index = std::type_index(typeid(float)).hash_code();
static const size_t char_type_index = std::type_index(typeid(char)).hash_code();
static const size_t uchar_type_index = std::type_index(typeid(uchar)).hash_code();
static const size_t uint_type_index = std::type_index(typeid(uint)).hash_code();
static const size_t ulint_type_index = std::type_index(typeid(ulint)).hash_code();
static const size_t llint_type_index = std::type_index(typeid(llint)).hash_code();
static const size_t ullint_type_index = std::type_index(typeid(ullint)).hash_code();
static const size_t short_type_index = std::type_index(typeid(short int)).hash_code();
static const size_t ushort_type_index = std::type_index(typeid(short unsigned int)).hash_code();
static const size_t bool_type_index = std::type_index(typeid(bool)).hash_code();
/// \}
  
/// \brief Test whether some data type is a recognized scalar.
template <typename T> bool isScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int_type_index || ct == uint_type_index || ct == double_type_index ||
          ct == float_type_index || ct == char_type_index || ct == uchar_type_index ||
          ct == llint_type_index || ct == ullint_type_index || ct == short_type_index ||
          ct == ushort_type_index || ct == ulint_type_index || ct == bool_type_index);
}

/// \brief Test whether some data type is a recognized (signed) integral scalar.
template <typename T> bool isSignedIntegralScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int_type_index || ct == char_type_index || ct == llint_type_index ||
          ct == short_type_index);
}

/// \brief Test whether some data type is a recognized (unsigned) integral scalar.
template <typename T> bool isUnsignedIntegralScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == uint_type_index || ct == ulint_type_index || ct == uchar_type_index ||
          ct == ullint_type_index || ct == ushort_type_index || ct == bool_type_index);
}

/// \brief Test whether some data type is a recongized floating point number.
template <typename T> bool isFloatingPointScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == double_type_index || ct == float_type_index);
}

/// \brief Produce a platform-independent name by which to identify one of the scalara data types.
template <typename T> std::string getOmniScalarTypeName() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int_type_index) return "int";
  else if (ct == uint_type_index) return "unsigned_int";
  else if (ct == ulint_type_index) return "unsigned_long_int";
  else if (ct == double_type_index) return "double";
  else if (ct == float_type_index) return "float";
  else if (ct == char_type_index) return "char";
  else if (ct == uchar_type_index) return "unsigned_char";
  else if (ct == llint_type_index) return "long_long_int";
  else if (ct == ullint_type_index) return "unsigned_long_long_int";
  else if (ct == short_type_index) return "short_int";
  else if (ct == ushort_type_index) return "unsigned_short_int";
  else if (ct == bool_type_index) return "bool";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
	  "scalar type.", "getOmniScalarTypeName");
  }
  __builtin_unreachable();
}
  
} // namespace data_types
} // namespace omni

namespace omni {
using data_types::ulint;
using data_types::ulong;
using data_types::llint;
using data_types::ullint;
using data_types::uchar;
using data_types::int_type_index;
using data_types::double_type_index;
using data_types::float_type_index;
using data_types::char_type_index;
using data_types::uchar_type_index;
using data_types::ulint_type_index;
using data_types::llint_type_index;
using data_types::ullint_type_index;
using data_types::short_type_index;
using data_types::ushort_type_index;
using data_types::bool_type_index;
} // namespace omni

#endif
