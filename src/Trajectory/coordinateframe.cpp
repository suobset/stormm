#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "amber_ascii.h"
#include "coordinateframe.h"
#include "trajectory_enumerators.h"

namespace omni {
namespace trajectory {

using math::extractBoxDimensions;
using math::roundUp;
using parse::TextFile;
  
//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const int natom_in, const UnitCellType unit_cell_in,
                                             const double* xcrd_in, const double* ycrd_in,
                                             const double* zcrd_in, const double* umat_in,
                                             const double* invu_in, const double* boxdim_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter::CoordinateFrameWriter(const int natom_in, const UnitCellType unit_cell_in,
                                             double* xcrd_in, double* ycrd_in, double* zcrd_in,
                                             double* umat_in, double* invu_in, double* boxdim_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const int natom_in, const UnitCellType unit_cell_in) :
    file_name{std::string("")},
    frame_number{0},
    atom_count{natom_in},
    unit_cell{unit_cell_in},
    x_coordinates{HybridKind::POINTER, "x_coordinates"},
    y_coordinates{HybridKind::POINTER, "y_coordinates"},
    z_coordinates{HybridKind::POINTER, "z_coordinates"},
    box_space_transform{HybridKind::POINTER, "box_transform"},
    inverse_transform{HybridKind::POINTER, "inv_transform"},
    box_dimensions{HybridKind::POINTER, "box_dimensions"},
    storage{HybridKind::ARRAY, "frame_data"}
{
  allocate();
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const int natom_in, const UnitCellType unit_cell_in,
                                 const double* xcrd_in, const double* ycrd_in,
                                 const double* zcrd_in, const double* umat_in,
                                 const double* invu_in, const double* boxdim_in) :
    CoordinateFrame(natom_in, unit_cell_in)
{
  // No method yet for putting data into a Hybrid object based on a pointer to an array
  double* xtmp = x_coordinates.data();
  double* ytmp = y_coordinates.data();
  double* ztmp = z_coordinates.data();
  for (int i = 0; i < natom_in; i++) {
    xtmp[i] = xcrd_in[i];
    ytmp[i] = ycrd_in[i];
    ztmp[i] = zcrd_in[i];
  }

  // Create transformation matrices if available, or initialize to the identity matrix
  double* utmp = box_space_transform.data();
  double* invutmp = inverse_transform.data();
  double* bdimtmp = box_dimensions.data();
  if (umat_in != nullptr) {
    for (int i = 0; i < 9; i++) {
      utmp[i] = umat_in[i];
      invutmp[i] = invu_in[i];
    }

    // Extract the box dimensions directly from the inverse transformation matrix if necessary
    if (boxdim_in == nullptr) {
      extractBoxDimensions(&bdimtmp[0], &bdimtmp[1], &bdimtmp[2], &bdimtmp[3], &bdimtmp[4],
                           &bdimtmp[5], invu_in);
    }
    else {
      for (int i = 0; i < 6; i++) {
        bdimtmp[i] = boxdim_in[i];
      }
    }
  }
  else {
    for (int i = 0; i < 9; i++) {
      utmp[i]    = static_cast<double>((i & 0x3) == 0);
      invutmp[i] = static_cast<double>((i & 0x3) == 0);
    }
    for (int i = 0; i < 3; i++) {
      bdimtmp[i    ] = 0.0;
      bdimtmp[i + 3] = 0.5 * symbols::pi;
    }
  }
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const std::string &file_name_in,
                                 const CoordinateFileKind file_kind, const int frame_number_in) :
    CoordinateFrame()
{
  buildFromFile(file_name_in, file_kind, frame_number_in);
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(PhaseSpace *ps) :
    CoordinateFrame(ps->getAtomCount(), ps->getUnitCellType())
{
  file_name = ps->getFileName();
  allocate();
  PhaseSpaceWriter psr = ps->data(HybridTargetLevel::HOST);
  double* xptr = x_coordinates.data();
  double* yptr = y_coordinates.data();
  double* zptr = z_coordinates.data();
  for (int i = 0; i < atom_count; i++) {
    xptr[i] = psr.xcrd[i];
    yptr[i] = psr.ycrd[i];
    zptr[i] = psr.zcrd[i];
  }
  double* utmp = box_space_transform.data();
  double* invutmp = inverse_transform.data();
  double* bdimtmp = box_dimensions.data();
  for (int i = 0; i < 9; i++) {
    utmp[i]    = psr.umat[i];
    invutmp[i] = psr.invu[i];
  }
  for (int i = 0; i < 6; i++) {
    bdimtmp[i] = psr.boxdim[i];
  }
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const PhaseSpace &ps) :
    CoordinateFrame(ps.getAtomCount(), ps.getUnitCellType())
{
  file_name = ps.getFileName();
  allocate();
  PhaseSpaceReader psr = ps.data(HybridTargetLevel::HOST);
  double* xptr = x_coordinates.data();
  double* yptr = y_coordinates.data();
  double* zptr = z_coordinates.data();
  for (int i = 0; i < atom_count; i++) {
    xptr[i] = psr.xcrd[i];
    yptr[i] = psr.ycrd[i];
    zptr[i] = psr.zcrd[i];
  }
  double* utmp = box_space_transform.data();
  double* invutmp = inverse_transform.data();
  double* bdimtmp = box_dimensions.data();
  for (int i = 0; i < 9; i++) {
    utmp[i]    = psr.umat[i];
    invutmp[i] = psr.invu[i];
  }
  for (int i = 0; i < 6; i++) {
    bdimtmp[i] = psr.boxdim[i];
  }
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const CoordinateFrame &original) :
    file_name{original.file_name},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    box_space_transform{original.box_space_transform},
    inverse_transform{original.inverse_transform},
    box_dimensions{original.box_dimensions},
    storage{original.storage}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame& CoordinateFrame::operator=(const CoordinateFrame &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy the file name (if applicable) and atom count
  file_name = other.file_name;
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;

  // Copy the Hybrid objects to preserve tags and the proper kinds.  Copying the data storage
  // object carries over all of the other object's contents.
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  box_space_transform = other.box_space_transform;
  inverse_transform = other.inverse_transform;
  box_dimensions = other.box_dimensions;
  storage = other.storage;

  // Allocating again sets the internal POINTER-kind Hybrid objects
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(CoordinateFrame &&original) :
    file_name{std::move(original.file_name)},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    box_space_transform{std::move(original.box_space_transform)},
    inverse_transform{std::move(original.inverse_transform)},
    box_dimensions{std::move(original.box_dimensions)},
    storage{std::move(original.storage)}
{
  // No repair of the pointers is necessary, as they still point to valid data
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame& CoordinateFrame::operator=(CoordinateFrame &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  file_name = std::move(other.file_name);
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  box_space_transform = std::move(other.box_space_transform);
  inverse_transform = std::move(other.inverse_transform);
  box_dimensions = std::move(other.box_dimensions);
  storage = std::move(other.storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::buildFromFile(const std::string &file_name_in,
                                    const CoordinateFileKind file_kind, const int frame_number) {
  file_name = file_name_in;

  // Try to detect the file format if it is not already specified.  If it remains UNKNOWN, that
  // will ultimately lead to an error.
  CoordinateFileKind actual_kind = file_kind;
  if (file_kind == CoordinateFileKind::UNKNOWN) {
    actual_kind = detectCoordinateFileKind(file_name);
  }  
  switch (actual_kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      // The number of atoms must be known a-priori in order to read from a .crd trajectory file.
      if (atom_count == 0) {
        rtErr("A number of atoms matching the trajectory must be known prior to reading a .crd "
              "file.", "CoordinateFrame", "buildFromFile");
      }
      allocate();
      TextFile tf(file_name);
      readAmberCrdFormat(tf, &x_coordinates, &y_coordinates, &z_coordinates, &box_space_transform,
                         &inverse_transform, frame_number);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
    {
      TextFile tf(file_name);
      atom_count = getAmberRestartAtomCount(tf);
      allocate();
      getAmberInputCoordinates(tf, &x_coordinates, &y_coordinates, &z_coordinates,
                               &box_space_transform, &inverse_transform, &box_dimensions);
    }
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      TextFile tf(file_name);
      atom_count = getAmberRestartAtomCount(tf);
      allocate();
      getAmberInputCoordinates(tf, &x_coordinates, &y_coordinates, &z_coordinates,
                               &box_space_transform, &inverse_transform, &box_dimensions);
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " could not be understood.", "CoordinateFrame",
          "buildFromFile");
  }

  // Interpret the box transformation
  unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
}

//-------------------------------------------------------------------------------------------------
int CoordinateFrame::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
UnitCellType CoordinateFrame::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
std::vector<double>
CoordinateFrame::getInterlacedCoordinates(const TrajectoryKind kind,
                                          const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(0, atom_count, kind, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double>
CoordinateFrame::getInterlacedCoordinates(const int low_index, const int high_index,
                                          const TrajectoryKind kind,
                                          const HybridTargetLevel tier) const {

  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "PhaseSpace",
          "getInterlacedCoordinates");
  }
  std::vector<double> result(3 * (high_index - low_index));
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* xptr = x_coordinates.data();
      const double* yptr = y_coordinates.data();
      const double* zptr = z_coordinates.data();
      for (int i = low_index; i < high_index; i++) {
        const int base_idx = 3 * (i - low_index);
        result[base_idx    ] = xptr[i];
        result[base_idx + 1] = yptr[i];
        result[base_idx + 2] = zptr[i];
      }
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> xval = x_coordinates.readDevice(low_index, high_index);
      const std::vector<double> yval = y_coordinates.readDevice(low_index, high_index);
      const std::vector<double> zval = z_coordinates.readDevice(low_index, high_index);
      for (int i = low_index; i < high_index; i++) {
        const int base_idx = 3 * (i - low_index);
        result[base_idx    ] = xval[i];
        result[base_idx + 1] = yval[i];
        result[base_idx + 2] = zval[i];
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getBoxSpaceTransform(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_space_transform.readHost();
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_space_transform.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getInverseTransform(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return inverse_transform.readHost();
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    return inverse_transform.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getBoxDimensions(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_dimensions.readHost();
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_dimensions.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::allocate() {
  const int padded_atom_count  = (atom_count > 0) ? roundUp(atom_count, warp_size_int) : 1;
  const int padded_matrix_size = roundUp(9, warp_size_int);
  storage.resize((3 * padded_atom_count) + (3 * padded_matrix_size));
  x_coordinates.setPointer(&storage,                           0, atom_count);
  y_coordinates.setPointer(&storage,           padded_atom_count, atom_count);
  z_coordinates.setPointer(&storage,       2 * padded_atom_count, atom_count);
  box_space_transform.setPointer(&storage, 3 * padded_atom_count,                             9);
  inverse_transform.setPointer(&storage,  (3 * padded_atom_count) +      padded_matrix_size,  9);
  box_dimensions.setPointer(&storage,     (3 * padded_atom_count) + (2 * padded_matrix_size), 6);
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void CoordinateFrame::upload() {
  storage.upload();
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::download() {
  storage.download();
}
#endif

//-------------------------------------------------------------------------------------------------
const CoordinateFrameReader CoordinateFrame::data(HybridTargetLevel tier) const {
  return CoordinateFrameReader(atom_count, unit_cell, x_coordinates.data(tier),
                               y_coordinates.data(tier), z_coordinates.data(tier),
                               box_space_transform.data(tier), inverse_transform.data(tier),
                               box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter CoordinateFrame::data(HybridTargetLevel tier) {
  return CoordinateFrameWriter(atom_count, unit_cell, x_coordinates.data(tier),
                               y_coordinates.data(tier), z_coordinates.data(tier),
                               box_space_transform.data(tier), inverse_transform.data(tier),
                               box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
const CoordinateFrameReader getCoordinateFrameReader(const PhaseSpace &ps) {
  const PhaseSpaceReader psr = ps.data();
  return CoordinateFrameReader(psr.natom, psr.unit_cell, psr.xcrd, psr.ycrd, psr.zcrd,
                               psr.umat, psr.invu, psr.boxdim);
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter getCoordinateFrameWriter(PhaseSpace *ps) {
  PhaseSpaceWriter psw = ps->data();
  return CoordinateFrameWriter(psw.natom, psw.unit_cell, psw.xcrd, psw.ycrd, psw.zcrd,
                               psw.umat, psw.invu, psw.boxdim);
}

} // namespace trajectory
} // namespace omni
