#include "copyright.h"
#include "Math/matrix_ops.h"
#include "puremesh.h"

namespace stormm {
namespace structure {

using math::invertSquareMatrix;

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in, const int scale_bits_in) :
    na{na_in}, nb{nb_in}, nc{nc_in}, origin_x{doubleToInt95(origin_x_in)},
    origin_y{doubleToInt95(origin_y_in)}, origin_z{doubleToInt95(origin_z_in)},
    scale_bits{scale_bits_in}, scale_factor{pow(2.0, scale_bits)},
    inverse_scale_factor{1.0 / scale_factor}, element_umat{}, sp_element_umat{}, element_invu{},
    sp_element_invu{}, fp_element_invu{}
{}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in,
                               const std::vector<double> &element_vectors,
                               const int scale_bits_in) :
    MeshParameters(na_in, nb_in, nc_in, origin_x_in, origin_y_in, origin_z_in, scale_bits_in)
{
  if (element_vectors.size() == 9LLU) {
    for (int i = 0; i < 9; i++) {
      element_invu[i] = element_vectors[i];
      sp_element_invu[i] = element_vectors[i];
      fp_element_invu[i] = doubleToInt95(element_vectors[i] * scale_factor);
    }
    invertSquareMatrix(element_invu, element_umat, 3);
    for (int i = 0; i < 9; i++) {
      sp_element_umat[i] = element_umat[i];
    }
  }
  else if (element_vectors.size() == 3LLU) {
    for (int i = 0; i < 9; i++) {
      element_umat[i] = 0.0;
      sp_element_umat[i] = 0.0;
      element_invu[i] = 0.0;
      sp_element_invu[i] = 0.0;
      fp_element_invu[i] = { 0LL, 0 };
    }
    for (int i = 0; i < 3; i++) {
      element_umat[4 * i] = 1.0 / element_vectors[i];
      sp_element_umat[4 * i] = element_umat[4 * i];
      element_invu[4 * i] = element_vectors[i];
      sp_element_invu[4 * i] = element_invu[4 * i];
      fp_element_invu[4 * i] = doubleToInt95(element_vectors[i] * scale_factor);
    }
  }
  else {
    rtErr("The mesh element is defined by a 3x3 matrix.  A total of " +
          std::to_string(element_vectors.size()) + " elements were provided.", "MeshParameters");
  }
}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in, const double element_x,
                               const double element_y, const double element_z,
			       const int scale_bits_in) :
  MeshParameters(na_in, nb_in, nc_in, origin_x_in, origin_y_in, origin_z_in,
                 { element_x, element_y, element_z }, scale_bits_in)
{}

//-------------------------------------------------------------------------------------------------
int MeshParameters::getAxisElementCount(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return na;
  case UnitCellAxis::B:
    return nb;
  case UnitCellAxis::C:
    return nc;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MeshParameters::getAxisElementCount(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return na;
  case CartesianDimension::Y:
    return nb;
  case CartesianDimension::Z:
    return nc;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getMeshOrigin(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return int95ToDouble(origin_x) * inverse_scale_factor;
  case CartesianDimension::Y:
    return int95ToDouble(origin_y) * inverse_scale_factor;
  case CartesianDimension::Z:
    return int95ToDouble(origin_z) * inverse_scale_factor;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshOriginAsFixedPrecision() const {
  std::vector<int95_t> result(3);
  result[0] = origin_x;
  result[1] = origin_y;
  result[2] = origin_z;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t MeshParameters::getMeshOriginAsFixedPrecision(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return origin_x;
  case CartesianDimension::Y:
    return origin_y;
  case CartesianDimension::Z:
    return origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MeshParameters::getScalingBits() const {
  return scale_bits;
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getScalingFactor() const {
  return scale_factor;
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getInverseScalingFactor() const {
  return inverse_scale_factor;
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getAxisCoordinates(const UnitCellAxis mesh_axis,
                                                        const CartesianDimension cart_axis) const {
  std::vector<int95_t> result;
  const int cdim = static_cast<int>(cart_axis);
  switch (mesh_axis) {
  case UnitCellAxis::A:
    result.resize(na + 1);
    result[0] = origin_x;
    for (int i = 0; i < na; i++) {
      result[i + 1] = splitFPSum(result[i], fp_element_invu[cdim]);
    }
    break;
  case UnitCellAxis::B:
    result.resize(nb + 1);
    result[0] = origin_y;
    for (int i = 0; i < na; i++) {
      result[i + 1] = splitFPSum(result[i], fp_element_invu[3 + cdim]);
    }
    break;
  case UnitCellAxis::C:
    result.resize(nc + 1);
    result[0] = origin_z;
    for (int i = 0; i < nc; i++) {
      result[i + 1] = splitFPSum(result[i], fp_element_invu[6 + cdim]);
    }
    break;
  }
  return result;
}

} // namespace structure
} // namespace stormm
