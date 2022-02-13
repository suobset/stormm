#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "Parsing/textfile.h"
#include "amber_ascii.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace omni {
namespace trajectory {

using card::HybridKind;
using math::roundUp;
using parse::TextFile;

//-------------------------------------------------------------------------------------------------
PhaseSpaceReader::PhaseSpaceReader(const int natom_in, const UnitCellType unit_cell_in,
                                   const double* xcrd_in, const double* ycrd_in,
                                   const double* zcrd_in, const double* umat_in,
                                   const double* invu_in, const double* boxdim_in,
                                   const double* xvel_in, const double* yvel_in,
                                   const double* zvel_in, const double* xfrc_in,
                                   const double* yfrc_in, const double* zfrc_in,
                                   const double* xprv_in, const double* yprv_in,
                                   const double* zprv_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in},
    xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in}, xprv{xprv_in}, yprv{yprv_in}, zprv{zprv_in}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter::PhaseSpaceWriter(const int natom_in, const UnitCellType unit_cell_in,
                                   double* xcrd_in, double* ycrd_in, double* zcrd_in,
                                   double* umat_in, double* invu_in, double* boxdim_in,
                                   double* xvel_in, double* yvel_in, double* zvel_in,
                                   double* xfrc_in, double* yfrc_in, double* zfrc_in,
                                   double* xprv_in, double* yprv_in, double* zprv_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in},
    xfrc{xfrc_in}, yfrc{yfrc_in}, zfrc{zfrc_in}, xprv{xprv_in}, yprv{yprv_in}, zprv{zprv_in}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const int atom_count_in, const UnitCellType unit_cell_in) :
    file_name{std::string("")},
    atom_count{atom_count_in},
    unit_cell{unit_cell_in},
    x_coordinates{HybridKind::POINTER, "x_coordinates"},
    y_coordinates{HybridKind::POINTER, "y_coordinates"},
    z_coordinates{HybridKind::POINTER, "z_coordinates"},
    box_space_transform{HybridKind::POINTER, "box_transform"},
    inverse_transform{HybridKind::POINTER, "inv_transform"},
    box_dimensions{HybridKind::POINTER, "box_dimensions"},
    x_velocities{HybridKind::POINTER, "x_velocities"},
    y_velocities{HybridKind::POINTER, "y_velocities"},
    z_velocities{HybridKind::POINTER, "z_velocities"},
    x_forces{HybridKind::POINTER, "x_forces"},
    y_forces{HybridKind::POINTER, "y_forces"},
    z_forces{HybridKind::POINTER, "z_forces"},
    x_prior_coordinates{HybridKind::POINTER, "x_prior_coords"},
    y_prior_coordinates{HybridKind::POINTER, "y_prior_coords"},
    z_prior_coordinates{HybridKind::POINTER, "z_prior_coords"},
    storage{HybridKind::ARRAY, "phase_space_data"}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const std::string &file_name_in, const CoordinateFileKind file_kind,
                       const int frame_number) :
    PhaseSpace()
{
  buildFromFile(file_name_in, file_kind, frame_number);
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const std::string &file_name_in, const CoordinateFileKind file_kind,
                       const AtomGraph &ag, const int frame_number) :
    PhaseSpace(file_name_in, file_kind, frame_number)
{
  // Check that the coordinates agree with the topology
  if (ag.getAtomCount() != atom_count) {
    rtErr("Coordinates in file " + file_name + " describe " + std::to_string(atom_count) +
          "atoms, whereas topology " + ag.getFileName() + " describes a system of " +
          std::to_string(atom_count) + " atoms.", "PhaseSpace");
  }
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const PhaseSpace &original) :
    file_name{original.file_name},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    box_space_transform{original.box_space_transform},
    inverse_transform{original.inverse_transform},
    box_dimensions{original.box_dimensions},
    x_velocities{original.x_velocities},
    y_velocities{original.y_velocities},
    z_velocities{original.z_velocities},
    x_forces{original.x_forces},
    y_forces{original.y_forces},
    z_forces{original.z_forces},
    x_prior_coordinates{original.x_prior_coordinates},
    y_prior_coordinates{original.y_prior_coordinates},
    z_prior_coordinates{original.z_prior_coordinates},
    storage{original.storage}
{
  // Set the POINTER-kind Hybrids in the new object appropriately.  Even the resize() operation
  // inherent to "allocate" will pass by with little more than a check that the length of the data
  // storage array is already what it should be.
  allocate();
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& PhaseSpace::operator=(const PhaseSpace &other) {
  
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
  x_velocities = other.x_velocities;
  y_velocities = other.y_velocities;
  z_velocities = other.z_velocities;
  x_forces = other.x_forces;
  y_forces = other.y_forces;
  z_forces = other.z_forces;
  x_prior_coordinates = other.x_prior_coordinates;
  y_prior_coordinates = other.y_prior_coordinates;
  z_prior_coordinates = other.z_prior_coordinates;
  storage = other.storage;

  // Allocating again sets the internal POINTER-kind Hybrid objects
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(PhaseSpace &&original) :
    file_name{std::move(original.file_name)},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    box_space_transform{std::move(original.box_space_transform)},
    inverse_transform{std::move(original.inverse_transform)},
    box_dimensions{std::move(original.box_dimensions)},
    x_velocities{std::move(original.x_velocities)},
    y_velocities{std::move(original.y_velocities)},
    z_velocities{std::move(original.z_velocities)},
    x_forces{std::move(original.x_forces)},
    y_forces{std::move(original.y_forces)},
    z_forces{std::move(original.z_forces)},
    x_prior_coordinates{std::move(original.x_prior_coordinates)},
    y_prior_coordinates{std::move(original.y_prior_coordinates)},
    z_prior_coordinates{std::move(original.z_prior_coordinates)},
    storage{std::move(original.storage)}
{
  // No repair of the pointers is necessary, as they still point to valid data
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& PhaseSpace::operator=(PhaseSpace &&other) {

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
  x_velocities = std::move(other.x_velocities);
  y_velocities = std::move(other.y_velocities);
  z_velocities = std::move(other.z_velocities);
  x_forces = std::move(other.x_forces);
  y_forces = std::move(other.y_forces);
  z_forces = std::move(other.z_forces);
  x_prior_coordinates = std::move(other.x_prior_coordinates);
  y_prior_coordinates = std::move(other.y_prior_coordinates);
  z_prior_coordinates = std::move(other.z_prior_coordinates);
  storage = std::move(other.storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::buildFromFile(const std::string &file_name_in, const CoordinateFileKind file_kind,
                               const int frame_number) {
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
              "file.", "PhaseSpace", "buildFromFile");
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
      getAmberRestartVelocities(tf, &x_velocities, &y_velocities, &z_velocities);
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The coordinate file type of " + file_name + " could not be understood.", "PhaseSpace",
          "buildFromFile");
  }

  // Interpret the box transformation
  unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
}

//-------------------------------------------------------------------------------------------------
std::string PhaseSpace::getFileName() const {
  return file_name;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpace::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
UnitCellType PhaseSpace::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const TrajectoryKind kind,
                                                         const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(0, atom_count, kind, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const int low_index,
                                                         const int high_index,
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
      double *xptr, *yptr, *zptr;
      switch (kind) {
      case TrajectoryKind::POSITIONS:
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
      case TrajectoryKind::VELOCITIES:
        {
          const double* xptr = x_velocities.data();
          const double* yptr = y_velocities.data();
          const double* zptr = z_velocities.data();
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            result[base_idx    ] = xptr[i];
            result[base_idx + 1] = yptr[i];
            result[base_idx + 2] = zptr[i];
          }
        }
        break;
      case TrajectoryKind::FORCES:
        {
          const double* xptr = x_forces.data();
          const double* yptr = y_forces.data();
          const double* zptr = z_forces.data();
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            result[base_idx    ] = xptr[i];
            result[base_idx + 1] = yptr[i];
            result[base_idx + 2] = zptr[i];
          }
        }
        break;
      }
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      std::vector<double> xval, yval, zval;
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        xval = x_coordinates.readDevice(low_index, high_index);
        yval = y_coordinates.readDevice(low_index, high_index);
        zval = z_coordinates.readDevice(low_index, high_index);
        break;
      case TrajectoryKind::VELOCITIES:
        xval = x_velocities.readDevice(low_index, high_index);
        yval = y_velocities.readDevice(low_index, high_index);
        zval = z_velocities.readDevice(low_index, high_index);
        break;
      case TrajectoryKind::FORCES:
        xval = x_forces.readDevice(low_index, high_index);
        yval = y_forces.readDevice(low_index, high_index);
        zval = z_forces.readDevice(low_index, high_index);
        break;
      }
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
void PhaseSpace::initializeForces() {
  double* x_ptr = x_forces.data();
  double* y_ptr = y_forces.data();
  double* z_ptr = z_forces.data();
  for (int i = 0; i < atom_count; i++) {
    x_ptr[i] = 0.0;
    y_ptr[i] = 0.0;
    z_ptr[i] = 0.0;
  }
#ifdef OMNI_USE_HPC
  uploadForces();
#endif
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getBoxSpaceTransform(const HybridTargetLevel tier) const {
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
std::vector<double> PhaseSpace::getInverseTransform(const HybridTargetLevel tier) const {
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
const Hybrid<double>* PhaseSpace::getStoragePointer() const {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* PhaseSpace::getStoragePointer() {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::allocate() {
  const int padded_atom_count  = (atom_count > 0) ? roundUp(atom_count, warp_size_int) : 1;
  const int padded_matrix_size = roundUp(9, warp_size_int);
  storage.resize((15 * padded_atom_count) + (3 * padded_matrix_size));
  x_coordinates.setPointer(&storage,                            0, atom_count);
  y_coordinates.setPointer(&storage,            padded_atom_count, atom_count);
  z_coordinates.setPointer(&storage,        2 * padded_atom_count, atom_count);
  box_space_transform.setPointer(&storage,  3 * padded_atom_count, 9);
  inverse_transform.setPointer(&storage,   (3 * padded_atom_count) +      padded_matrix_size, 9);
  box_dimensions.setPointer(&storage,      (3 * padded_atom_count) + (2 * padded_matrix_size), 6);
  const int thus_far = (3 * padded_atom_count) + (3 * padded_matrix_size);
  x_velocities.setPointer(&storage,        thus_far,                           atom_count);
  y_velocities.setPointer(&storage,        thus_far +      padded_atom_count,  atom_count);
  z_velocities.setPointer(&storage,        thus_far + (2 * padded_atom_count), atom_count);
  x_forces.setPointer(&storage,            thus_far + (3 * padded_atom_count), atom_count);
  y_forces.setPointer(&storage,            thus_far + (4 * padded_atom_count), atom_count);
  z_forces.setPointer(&storage,            thus_far + (5 * padded_atom_count), atom_count);
  x_prior_coordinates.setPointer(&storage, thus_far + (6 * padded_atom_count), atom_count);
  y_prior_coordinates.setPointer(&storage, thus_far + (7 * padded_atom_count), atom_count);
  z_prior_coordinates.setPointer(&storage, thus_far + (8 * padded_atom_count), atom_count);
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void PhaseSpace::upload() {
  storage.upload();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadPositions() {
  storage.upload(0, 3 * roundUp(atom_count, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadTransformations() {
  const int stride = roundUp(atom_count, warp_size_int);
  storage.upload(3 * stride, 3 * roundUp(9, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadVelocities() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.upload((3 * stride) + offset, 3 * stride);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadForces() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.upload((6 * stride) + offset, 3 * stride);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadPriorPositions() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.upload((9 * stride) + offset, 3 * stride);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::download() {
  storage.download();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadPositions() {
  storage.download(0, 3 * roundUp(atom_count, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadTransformations() {
  const int stride = roundUp(atom_count, warp_size_int);
  storage.download(3 * stride, 3 * roundUp(9, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadVelocities() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.download((3 * stride) + offset, 3 * stride);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadForces() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.download((6 * stride) + offset, 3 * stride);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadPriorPositions() {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  storage.download((9 * stride) + offset, 3 * stride);
}
#endif

//-------------------------------------------------------------------------------------------------
const PhaseSpaceReader PhaseSpace::data(const HybridTargetLevel tier) const {
  return PhaseSpaceReader(atom_count, unit_cell, x_coordinates.data(tier),
                          y_coordinates.data(tier), z_coordinates.data(tier),
                          box_space_transform.data(tier), inverse_transform.data(tier),
                          box_dimensions.data(tier), x_velocities.data(tier),
                          y_velocities.data(tier), z_velocities.data(tier), x_forces.data(tier),
                          y_forces.data(tier), z_forces.data(tier),
                          x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                          z_prior_coordinates.data(tier));
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter PhaseSpace::data(const HybridTargetLevel tier) {
  return PhaseSpaceWriter(atom_count, unit_cell, x_coordinates.data(tier),
                          y_coordinates.data(tier), z_coordinates.data(tier),
                          box_space_transform.data(tier), inverse_transform.data(tier),
                          box_dimensions.data(tier), x_velocities.data(tier),
                          y_velocities.data(tier), z_velocities.data(tier), x_forces.data(tier),
                          y_forces.data(tier), z_forces.data(tier),
                          x_prior_coordinates.data(tier), y_prior_coordinates.data(tier),
                          z_prior_coordinates.data(tier));
}

//-------------------------------------------------------------------------------------------------
UnitCellType determineUnitCellTypeByShape(const double* inv_ptr) {
  if (inv_ptr[0] < 1.01 && inv_ptr[4] < 1.01 && inv_ptr[8] < 1.01) {

    // No unit cell could be this small, and the identity matrix is what the transforms are set to
    // when there are no boundary conditions.  It is an inocuous matrix that can be put through
    // various operations without leading to catastrophes, although code should check the unit_cell
    // member variable to avoid meaningless calculations.
    return UnitCellType::NONE;
  }
  else if (fabs(inv_ptr[1]) + fabs(inv_ptr[2]) + fabs(inv_ptr[3]) + fabs(inv_ptr[5]) +
           fabs(inv_ptr[6]) + fabs(inv_ptr[7]) < constants::small) {
    return UnitCellType::ORTHORHOMBIC;
  }
  else {
    return UnitCellType::TRICLINIC;
  }
  __builtin_unreachable();
}

} // namespace trajectory
} // namespace omni
