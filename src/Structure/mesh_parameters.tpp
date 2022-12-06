// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParamKit<T>::MeshParamKit(const int na_in, const int nb_in, const int nc_in,
                              const int95_t orig_x_in, const int95_t orig_y_in,
                              const int95_t orig_z_in, const T scale_in, const T inv_scale_in,
                              const int scale_bits_in, const T* umat_in, const T* invu_in,
                              const T* widths_in, const int95_t* fp_invu_in) :
    na{na_in}, nb{nb_in}, nc{nc_in}, orig_x{orig_x_in}, orig_y{orig_y_in}, orig_z{orig_z_in},
    scale{scale_in}, inv_scale{inv_scale_in}, scale_bits{scale_bits_in},
    umat{ umat_in[0], umat_in[1], umat_in[2], umat_in[3], umat_in[4], umat_in[5], umat_in[6],
          umat_in[7], umat_in[8] },
    invu{ invu_in[0], invu_in[1], invu_in[2], invu_in[3], invu_in[4], invu_in[5], invu_in[6],
          invu_in[7], invu_in[8] },
    widths{ widths_in[0], widths_in[1], widths_in[2] },
    fp_invu{ fp_invu_in[0], fp_invu_in[1], fp_invu_in[2], fp_invu_in[3], fp_invu_in[4],
             fp_invu_in[5], fp_invu_in[6], fp_invu_in[7], fp_invu_in[8] }
{}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshOrigin() const {
  if (isFloatingPointScalarType<Tcoord>()) {
    std::vector<Tcoord> result(3);
    result[0] = hostInt95ToDouble(origin_x) * inverse_scale_factor;
    result[1] = hostInt95ToDouble(origin_y) * inverse_scale_factor;
    result[2] = hostInt95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else if (isFloatingPointHpcVectorType<Tcoord>()) {
    rtErr("The mesh coordinate origin is available as an HPC vector type through the "
          "getMeshOriginAsTuple() function.", "MeshParameter", "getMeshOrigin");
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.", "MeshParameter", "getMeshOrigin");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T3> T3 MeshParameters::getMeshOriginAsTuple() const {
  if (isFloatingPointScalarType<T3>()) {
    rtErr("The mesh coordinate origin is available as a 3-element array of the desired floating "
          "point type through the getMeshOrigin() function.", "MeshParameter",
          "getMeshOriginAsTuple");
  }
  else if (isFloatingPointHpcVectorType<T3>()) {
    T3 result;
    result.x = hostInt95ToDouble(origin_x) * inverse_scale_factor;
    result.y = hostInt95ToDouble(origin_y) * inverse_scale_factor;
    result.z = hostInt95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.  The fixed-precision representation is "
          "not available as a tuple.", "MeshParameter", "getMeshOriginAsTuple");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
std::vector<Tcoord> MeshParameters::getMeshElementVector(const UnitCellAxis dim) const {
  std::vector<Tcoord> result(3);
  const int icol = static_cast<int>(dim);
  if (isFloatingPointScalarType<Tcoord>()) {
    for (int i = 0; i < 3; i++) {
      result[i] = hostInt95ToDouble(fp_element_invu[i + (3 * icol)]) * inverse_scale_factor;
    }
  }
  else {
    rtErr("Mesh element vectors can only be returned as " + getStormmScalarTypeName<float>() +
          " or " + getStormmScalarTypeName<double>() + ".  For the fixed-precision "
          "representation, use getMeshElementVectorFP().", "MeshParameters",
          "getMeshElementVector");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
std::vector<Tcoord> MeshParameters::getMeshElementVector(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return getMeshElementVector<Tcoord>(UnitCellAxis::A);
  case CartesianDimension::Y:
    return getMeshElementVector<Tcoord>(UnitCellAxis::B);
  case CartesianDimension::Z:
    return getMeshElementVector<Tcoord>(UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T3>
T3 MeshParameters::getMeshElementVectorAsTuple(const UnitCellAxis dim) const {
  T3 result;
  const int icol = static_cast<int>(dim);
  if (isFloatingPointHpcVectorType<T3>()) {
    result.x = hostInt95ToDouble(fp_element_invu[(3 * icol)    ]) * inverse_scale_factor;
    result.y = hostInt95ToDouble(fp_element_invu[(3 * icol) + 1]) * inverse_scale_factor;
    result.z = hostInt95ToDouble(fp_element_invu[(3 * icol) + 2]) * inverse_scale_factor;
  }
  else {
    rtErr("Mesh element vectors can only be returned as " + getStormmScalarTypeName<float3>() +
          " or " + getStormmScalarTypeName<double3>() + " tuples.  To get one of the vectors as "
          "a real-valued vector, use getMeshElementVector().  For the fixed-precision "
          "representation, use getMeshElementVectorFP().", "MeshParameters",
          "getMeshElementVector");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T3>
T3 MeshParameters::getMeshElementVectorAsTuple(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::A);
  case CartesianDimension::Y:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::B);
  case CartesianDimension::Z:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshTransform() const {
  std::vector<Tcoord> result(9);
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (ct == double_type_index) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = element_umat[i];
    }
  }
  else if (ct == float_type_index) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = sp_element_umat[i];
    }
  }
  else {
    rtErr("The transformation matrix into element space is only available in single- or double-"
          "precision floating point numbers.", "MeshParameters", "getMeshTransform");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshInverseTransform() const {
  std::vector<Tcoord> result(9);
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (isFloatingPointScalarType<Tcoord>()) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = hostInt95ToDouble(fp_element_invu[i]) * inverse_scale_factor;
    }
  }
  else {
    rtErr("The inverse transformation matrix (the column matrix of element vectors) is only "
          "available in fixed-precision format or single- or double-precision floating point "
          "numbers.  To get the fixed-precision format, use getMeshInverseTransformAsFP().",
          "MeshParameters", "getMeshTransform");
  }
  return result;
}

} // namespace structure
} // namespace stormm
