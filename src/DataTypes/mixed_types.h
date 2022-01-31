// -*-c++-*-
#ifndef OMNI_MIXED_TYPES_H
#define OMNI_MIXED_TYPES_H

namespace omni {
namespace data_types {

/// \brief A combination of integer and double data.  Like other mixed tuples in this library,
///        this will maintain the x, y, z, and w member variable naming conventions.
struct CombineIDp {
  int x;
  double y;
};

} // namespace data_types
} // namespace omni

namespace omni {
using data_types::CombineIDp;
}

#endif
