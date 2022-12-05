#include <climits>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/fixed_precision.h"
#include "Math/matrix_ops.h"
#include "mesh_parameters.h"

namespace stormm {
namespace structure {

using numerics::min_localpos_scale_bits;
using math::hessianNormalWidths;
using math::invertSquareMatrix;

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in,
                               const std::vector<double> &element_vectors,
                               const int scale_bits_in) :
    na{na_in}, nb{nb_in}, nc{nc_in}, origin_x{0LL, 0}, origin_y{0LL, 0}, origin_z{0LL, 0},
    scale_bits{scale_bits_in}, scale_factor{pow(2.0, scale_bits)},
    inverse_scale_factor{1.0 / scale_factor}, unit_cell{UnitCellType::NONE}, element_umat{},
    sp_element_umat{}, element_invu{}, sp_element_invu{}, fp_element_invu{}
{
  origin_x = doubleToInt95(origin_x_in * scale_factor);
  origin_y = doubleToInt95(origin_y_in * scale_factor);
  origin_z = doubleToInt95(origin_z_in * scale_factor);
  defineElement(element_vectors);
  validateMeshDimensions();
  validateFixedPrecisionBits();
}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters() :
    MeshParameters(1, 1, 1, 0.0, 0.0, 0.0, { 1.0, 1.0, 1.0 })
{}

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
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in, const double element_width,
                               const int scale_bits_in) :
    MeshParameters(na_in, nb_in, nc_in, origin_x_in, origin_y_in, origin_z_in,
                   { element_width, element_width, element_width }, scale_bits_in)
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
std::vector<int95_t> MeshParameters::getMeshOriginAsFP() const {
  std::vector<int95_t> result(3);
  result[0] = origin_x;
  result[1] = origin_y;
  result[2] = origin_z;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t MeshParameters::getMeshOriginAsFP(const CartesianDimension dim) const {
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
UnitCellType MeshParameters::getMeshCellType() const {
  return unit_cell;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshElementVectorAsFP(const UnitCellAxis dim) const {
  std::vector<int95_t> result(3);
  const int icol = static_cast<int>(dim);
  for (int i = 0; i < 3; i++) {
    result[i] = fp_element_invu[i + (3 * icol)];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshInverseTransformAsFP() const {
  std::vector<int95_t> result(9);
  for (int i = 0; i < 9; i++) {
    result[i] = fp_element_invu[i];
  }
  return result;
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

//-------------------------------------------------------------------------------------------------
MeshParamKit<double> MeshParameters::dpData() const {
  return MeshParamKit<double>(na, nb, nc, origin_x, origin_y, origin_z, scale_factor,
                              inverse_scale_factor, scale_bits, element_umat, element_invu, widths,
                              fp_element_invu);
}

//-------------------------------------------------------------------------------------------------
MeshParamKit<float> MeshParameters::spData() const {
  return MeshParamKit<float>(na, nb, nc, origin_x, origin_y, origin_z, scale_factor,
                             inverse_scale_factor, scale_bits, sp_element_umat, sp_element_invu,
                             sp_widths, fp_element_invu);
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setMeshDimension(const int n_in, const UnitCellAxis mesh_axis) {
  if (n_in < 0) {
    rtErr("A mesh dimension of " + std::to_string(n_in) + " along the " +
          getEnumerationName(mesh_axis) + " is invalid.", "MeshParameters", "setMeshDimension");
  }
  switch (mesh_axis) {
  case UnitCellAxis::A:
    na = n_in;
  case UnitCellAxis::B:
    nb = n_in;
  case UnitCellAxis::C:
    nc = n_in;
  }
  validateMeshDimensions();
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setMeshDimension(const std::vector<int> &n_in) {
  if (n_in.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(n_in.size()) + " provided).",
          "MeshParameters", "setMeshDimension");
  }
  na = n_in[0];
  nb = n_in[1];
  nc = n_in[2];
  validateMeshDimensions();  
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const double v, CartesianDimension cart_axis) {
  setOrigin(doubleToInt95(v), cart_axis);
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const int95_t v, CartesianDimension cart_axis) {
  switch (cart_axis) {
  case CartesianDimension::X:
    origin_x = v;
    break;
  case CartesianDimension::Y:
    origin_y = v;
    break;
  case CartesianDimension::Z:
    origin_z = v;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const std::vector<double> &v) {
  if (v.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(v.size()) + " provided).",
          "MeshParameters", "setOrigin");
  }
  std::vector<int95_t> vfp = { doubleToInt95(v[0]), doubleToInt95(v[1]), doubleToInt95(v[2]) };
  setOrigin(vfp);
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const std::vector<int95_t> &v) {
  if (v.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(v.size()) + " provided).",
          "MeshParameters", "setOrigin");
  }
  origin_x = v[0];
  origin_y = v[1];
  origin_z = v[2];
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setScalingBits(const int scale_bits_in) {
  scale_bits = scale_bits_in;
  validateFixedPrecisionBits();
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::defineElement(const std::vector<double> &element_vectors) {
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

  // The determineUnitCellTypeByShape() function will define what would be very small unit cells,
  // in particular 1 x 1 x 1 Angstrom, as "NONE" type.  This is because there must be a unit cell
  // defined for some purposes, and such a unit cell would be an impossibly small simulation.
  // However, for mesh elements, the default "NONE" type unit cell dimensions are quite common.
  // Use a basic off-diagonal check instead.
  if (fabs(element_invu[1]) > constants::tiny || fabs(element_invu[2]) > constants::tiny ||
      fabs(element_invu[3]) > constants::tiny || fabs(element_invu[5]) > constants::tiny ||
      fabs(element_invu[6]) > constants::tiny || fabs(element_invu[7]) > constants::tiny) {
    unit_cell = UnitCellType::TRICLINIC;
  }
  else {
    unit_cell = UnitCellType::ORTHORHOMBIC;
  }

  // Use the Hessian Normal form to compute the number of mesh elements to search in each direction
  // and color all of the necessary elements around each atom.
  hessianNormalWidths(element_invu, &widths[0], &widths[1], &widths[2]);
  for (int i = 0; i < 3; i++) {
    sp_widths[i] = widths[i];
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::validateMeshDimensions() const {
  const llint total_elements = static_cast<llint>(na) * static_cast<llint>(nb) *
                               static_cast<llint>(nc);
  if (total_elements > INT_MAX) {
    rtErr("The total number of elements on the mesh cannot exceed " + std::to_string(INT_MAX) +
          " (currently " + std::to_string(total_elements) + ").", "MeshParameters",
          "validateMeshDimensions");
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::validateFixedPrecisionBits() const {
  if (scale_bits < min_localpos_scale_bits) {
    rtErr("A minimum of " + std::to_string(min_localpos_scale_bits) + " must be stored after the "
          "decimal in fixed-precision mesh representations (" + std::to_string(scale_bits) +
          " specified).", "MeshParameters", "validateFixedPrecisionBits");
  }
}
 
} // namespace structure
} // namespace stormm
