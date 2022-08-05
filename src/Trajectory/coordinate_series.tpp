// -*-c++-*-
namespace stormm {
namespace trajectory {

using numerics::default_trajpos_scale_bits;
using math::roundUp;

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T>::CoordinateSeriesWriter(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in,
                                                  const int gpos_bits_in,
                                                  const double gpos_scale_in,
                                                  const double inv_gpos_scale_in, T* xcrd_in,
                                                  T* ycrd_in, T* zcrd_in, double* umat_in,
                                                  double* invu_in, double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, gpos_bits{gpos_bits_in},
    gpos_scale{gpos_scale_in}, inv_gpos_scale{inv_gpos_scale_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesReader<T>::CoordinateSeriesReader(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in,
                                                  const int gpos_bits_in,
                                                  const double gpos_scale_in,
                                                  const double inv_gpos_scale_in, const T* xcrd_in,
                                                  const T* ycrd_in, const T* zcrd_in,
                                                  const double* umat_in, const double* invu_in,
                                                  const double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, gpos_bits{gpos_bits_in},
    gpos_scale{gpos_scale_in}, inv_gpos_scale{inv_gpos_scale_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesReader<T>::CoordinateSeriesReader(const CoordinateSeriesWriter<T> &csw) :
    natom{csw.natom}, nframe{csw.nframe}, unit_cell{csw.unit_cell}, gpos_bits{csw.gpos_bits},
    gpos_scale{csw.gpos_scale}, inv_gpos_scale{csw.inv_gpos_scale}, xcrd{csw.xcrd}, ycrd{csw.ycrd},
    zcrd{csw.zcrd}, umat{csw.umat}, invu{csw.invu}, boxdim{csw.boxdim}
{}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const int natom_in, const int nframe_in,
                                      const UnitCellType unit_cell_in,
                                      const int globalpos_scale_bits_in) :
    atom_count{natom_in}, frame_count{nframe_in}, frame_capacity{nframe_in},
    globalpos_scale_bits{globalpos_scale_bits_in}, unit_cell{unit_cell_in},
    globalpos_scale{pow(2.0, globalpos_scale_bits)},
    inverse_globalpos_scale{1.0 / globalpos_scale},
    x_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_x_coords"},
    y_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_y_coords"},
    z_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_z_coords"},
    box_space_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                         "cser_umat"},
    inverse_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                       "cser_invu"},
    box_dimensions{static_cast<size_t>(nframe_in) * roundUp<size_t>(6, warp_size_zu),
                   "cser_boxdims"}
{
  // Limits on the valid data types
  if (isFloatingPointScalarType<T>() == false && isSignedIntegralScalarType<T>() == false) {
    rtErr("A CoordinateSeries object is only valid with a signed integer or real number scalar "
          "data type.", "CoordinateSeries");
  }
  if (globalpos_scale_bits < 0) {
    rtErr("A fixed precision representation cannot have a negative bit count (" +
          std::to_string(globalpos_scale_bits) + ") after the decimal.", "CoordinateSeries");
  }
  if (isFloatingPointScalarType<T>() && globalpos_scale_bits != 0) {
    rtErr("A real-numbered representation of coordinates is incompatible with a fixed-precision "
          "representation of " + std::to_string(globalpos_scale_bits) + " bits after the decimal.",
          "CoordinateSeries");
  }
  if (isSignedIntegralScalarType<T>() && globalpos_scale_bits == 0) {
    globalpos_scale_bits = default_trajpos_scale_bits;
    globalpos_scale = pow(2.0, globalpos_scale_bits);
    inverse_globalpos_scale = 1.0 / globalpos_scale;
  }
  allocate(frame_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const std::string &file_name, const int atom_count_in,
                                      const CoordinateFileKind file_kind,
                                      const std::vector<int> &frame_numbers,
                                      const int replica_count, const UnitCellType unit_cell_in,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(atom_count_in, replica_count * static_cast<int>(frame_numbers.size()),
                     unit_cell_in, globalpos_scale_bits_in)
{
  importFromFile(file_name, file_kind, frame_numbers, replica_count, 0);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(PhaseSpace *ps, const int nframe_in,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(ps->getAtomCount(), nframe_in, ps->getUnitCellType(), globalpos_scale_bits_in)
{
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(ps, 0, ps->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const PhaseSpace &ps, const int nframe_in,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(ps.getAtomCount(), nframe_in, ps.getUnitCellType(), globalpos_scale_bits_in)
{
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(ps, 0, ps.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(CoordinateFrame *cf, const int nframe_in,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(cf->getAtomCount(), nframe_in, cf->getUnitCellType(), globalpos_scale_bits_in)
{
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(cf, 0, cf->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const CoordinateFrame &cf, const int nframe_in,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(cf.getAtomCount(), nframe_in, cf.getUnitCellType(), globalpos_scale_bits_in)
{
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(cf, 0, cf.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Toriginal>
CoordinateSeries<T>::CoordinateSeries(const CoordinateSeries<Toriginal> &original,
                                      const int globalpos_scale_bits_in) :
    CoordinateSeries(original.getAtomCount(), original.getFrameCount(), original.getUnitCellType(),
                     globalpos_scale_bits_in)
{
  for (int i = 0; i < frame_count; i++) {
    const CoordinateFrame cf = original.exportFrame(i);
    importCoordinateSet(cf, 0, cf.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int CoordinateSeries<T>::getAtomCount() const {
  return atom_count;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
int CoordinateSeries<T>::getFrameCount() const {
  return frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int CoordinateSeries<T>::getFrameCapacity() const {
  return frame_capacity;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
UnitCellType CoordinateSeries<T>::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int CoordinateSeries<T>::getFixedPrecisionBits() const {
  if (isFloatingPointScalarType<T>()) {
    rtWarn("A CoordinateSeries object with a real-numbered data representation does not have a "
           "meaningful fixed-precision scaling factor.", "CoordinateSeries",
           "getFixedPrecisionBits");
  }
  return globalpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Treport> std::vector<Treport>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index,
                                              const int globalpos_bits_out,
                                              const HybridTargetLevel tier) const {
  return getInterlacedCoordinates<Treport>(frame_index, 0, atom_count, globalpos_bits_out, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Treport> std::vector<Treport>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index, const int low_index,
                                              const int high_index, const int globalpos_bits_out,
                                              const HybridTargetLevel tier) const {

  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "CoordinateSeries",
          "getInterlacedCoordinates");
  }
  const int actual_bits_out = (globalpos_bits_out < 0) ? globalpos_scale_bits : globalpos_bits_out;
  std::vector<Treport> result(3 * (high_index - low_index));
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const T* xptr = x_coordinates.data();
      const T* yptr = y_coordinates.data();
      const T* zptr = z_coordinates.data();
      const size_t frame_offset = natom_zu * fidx_zu;
      if (isSignedIntegralScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        if (globalpos_scale_bits >= actual_bits_out) {
          llint divisor = 1LL;
          for (int i = actual_bits_out; i < globalpos_scale_bits; i++) {
            divisor *= 2LL;
          }
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            const size_t access_idx = frame_offset + static_cast<size_t>(i);
            result[base_idx    ] = static_cast<llint>(xptr[access_idx]) / divisor;
            result[base_idx + 1] = static_cast<llint>(yptr[access_idx]) / divisor;
            result[base_idx + 2] = static_cast<llint>(zptr[access_idx]) / divisor;
          }
        }
        else {
          const int shift_left = actual_bits_out - globalpos_scale_bits;
          llint multiplier = 1LL;
          for (int i = globalpos_scale_bits; i < actual_bits_out; i++) {
            multiplier *= 2LL;
          }
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            const size_t access_idx = frame_offset + static_cast<size_t>(i);
            result[base_idx    ] = static_cast<llint>(xptr[access_idx]) * multiplier;
            result[base_idx + 1] = static_cast<llint>(yptr[access_idx]) * multiplier;
            result[base_idx + 2] = static_cast<llint>(zptr[access_idx]) * multiplier;
          }
        }
      }
      else if (isSignedIntegralScalarType<T>() && isFloatingPointScalarType<Treport>()) {
        const Treport multiplier = 1.0 / pow(2.0, globalpos_scale_bits);
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ]  = xptr[access_idx];
          result[base_idx + 1]  = yptr[access_idx];
          result[base_idx + 2]  = zptr[access_idx];
          result[base_idx    ] *= multiplier;
          result[base_idx + 1] *= multiplier;
          result[base_idx + 2] *= multiplier;
        }
      }
      else if (isFloatingPointScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        const Treport multiplier = pow(2.0, actual_bits_out);
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ]  = llround(xptr[access_idx] * multiplier);
          result[base_idx + 1]  = llround(yptr[access_idx] * multiplier);
          result[base_idx + 2]  = llround(zptr[access_idx] * multiplier);
        }
      }
      else {
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ] = xptr[access_idx];
          result[base_idx + 1] = yptr[access_idx];
          result[base_idx + 2] = zptr[access_idx];
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t llim = (fidx_zu * natom_zu) + static_cast<size_t>(low_index);
      const size_t hlim = (fidx_zu * natom_zu) + static_cast<size_t>(high_index);
      const std::vector<T> xval = x_coordinates.readDevice(llim, hlim);
      const std::vector<T> yval = y_coordinates.readDevice(llim, hlim);
      const std::vector<T> zval = z_coordinates.readDevice(llim, hlim);
      const int irange = high_index - low_index;
      if (isSignedIntegralScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        if (globalpos_scale_bits >= actual_bits_out) {
          const int shift_right = globalpos_scale_bits - actual_bits_out;
          for (int i = 0; i < irange; i++) {
            result[(3 * i)    ] = (static_cast<llint>(xval[i]) >> shift_right);
            result[(3 * i) + 1] = (static_cast<llint>(yval[i]) >> shift_right);
            result[(3 * i) + 2] = (static_cast<llint>(zval[i]) >> shift_right);
          }
        }
        else {
          const int shift_left = actual_bits_out - globalpos_scale_bits;
          for (int i = 0; i < high_index - low_index; i++) {
            result[(3 * i)    ] = (static_cast<llint>(xval[i]) << shift_left);
            result[(3 * i) + 1] = (static_cast<llint>(yval[i]) << shift_left);
            result[(3 * i) + 2] = (static_cast<llint>(zval[i]) << shift_left);
          }
        }
      }
      else if (isSignedIntegralScalarType<T>() && isFloatingPointScalarType<Treport>()) {
        const Treport multiplier = 1.0 / pow(2.0, globalpos_scale_bits);
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ]  = xval[i];
          result[(3 * i) + 1]  = yval[i];
          result[(3 * i) + 2]  = zval[i];
          result[(3 * i)    ] *= multiplier;
          result[(3 * i) + 1] *= multiplier;
          result[(3 * i) + 2] *= multiplier;
        }
      }
      else if (isFloatingPointScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        const Treport multiplier = pow(2.0, globalpos_scale_bits);
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ]  = llround(xval[i] * multiplier);
          result[(3 * i) + 1]  = llround(yval[i] * multiplier);
          result[(3 * i) + 2]  = llround(zval[i] * multiplier);
        }
      }
      else {
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ] = xval[i];
          result[(3 * i) + 1] = yval[i];
          result[(3 * i) + 2] = zval[i];
        }
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getBoxSpaceTransform(const int frame_index,
                                                              const HybridTargetLevel tier) const {
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(9, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_space_transforms.readHost(read_start, 9);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_space_transforms.readDevice(read_start, 9);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getInverseTransform(const int frame_index,
                                                             const HybridTargetLevel tier) const {
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(9, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return inverse_transforms.readHost(read_start, 9);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return inverse_transforms.readDevice(read_start, 9);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getBoxDimensions(const int frame_index,
                                                          const HybridTargetLevel tier) const {
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(6, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_dimensions.readHost(read_start, 6);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_dimensions.readDevice(read_start, 6);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateFrame CoordinateSeries<T>::exportFrame(const int frame_index,
                                                 const HybridTargetLevel tier) const {
  CoordinateFrame result(atom_count, unit_cell);
  CoordinateFrameWriter resultw = result.data();
  const size_t natom_zu = atom_count;
  const size_t fidx_zu  = static_cast<size_t>(frame_index); 
  const size_t frame_offset = roundUp<size_t>(atom_count, warp_size_zu) * fidx_zu;
  const size_t xfrm_offset  = roundUp<size_t>(9, warp_size_zu) * fidx_zu;
  const size_t bdim_offset  = roundUp<size_t>(6, warp_size_zu) * fidx_zu;
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const T* xptr = x_coordinates.data();
      const T* yptr = y_coordinates.data();
      const T* zptr = z_coordinates.data();
      const double* umat_ptr = box_space_transforms.data();
      const double* invu_ptr = inverse_transforms.data();
      const double* bdim_ptr = box_dimensions.data();
      if (globalpos_scale_bits == 0) {
        for (size_t i = 0; i < natom_zu; i++) {
          resultw.xcrd[i] = xptr[frame_offset + i];
          resultw.ycrd[i] = yptr[frame_offset + i];
          resultw.zcrd[i] = zptr[frame_offset + i];
        }
      }
      else {
        const double conv_factor = inverse_globalpos_scale;
        for (size_t i = 0; i < natom_zu; i++) {
          resultw.xcrd[i] = static_cast<double>(xptr[frame_offset + i]) * conv_factor;
          resultw.ycrd[i] = static_cast<double>(yptr[frame_offset + i]) * conv_factor;
          resultw.zcrd[i] = static_cast<double>(zptr[frame_offset + i]) * conv_factor;
        }
      }
      for (size_t i = 0; i < 9LLU; i++) {
        resultw.umat[i] = umat_ptr[xfrm_offset + i];
        resultw.invu[i] = invu_ptr[xfrm_offset + i];
      }
      for (size_t i = 0; i < 6LLU; i++) {
        resultw.boxdim[i] = bdim_ptr[bdim_offset + i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<T> tmp_xcrd = x_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_ycrd = y_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_zcrd = z_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<double> tmp_umat = box_space_transforms.readDevice(xfrm_offset, 9);
      const std::vector<double> tmp_invu = inverse_transforms.readDevice(xfrm_offset, 9);
      const std::vector<double> tmp_bdim = box_dimensions.readDevice(bdim_offset, 6);
      if (globalpos_scale_bits == 0) {
        for (size_t i = 0; i < natom_zu; i++) {
          resultw.xcrd[i] = tmp_xcrd[i];
          resultw.ycrd[i] = tmp_ycrd[i];
          resultw.zcrd[i] = tmp_zcrd[i];
        }
      }
      else {
        const double conv_factor = inverse_globalpos_scale;
        for (size_t i = 0; i < natom_zu; i++) {
          resultw.xcrd[i] = static_cast<double>(tmp_xcrd[i]) * conv_factor;
          resultw.ycrd[i] = static_cast<double>(tmp_ycrd[i]) * conv_factor;
          resultw.zcrd[i] = static_cast<double>(tmp_zcrd[i]) * conv_factor;
        }
      }
      for (size_t i = 0; i < 9LLU; i++) {
        resultw.umat[i] = tmp_umat[i];
        resultw.invu[i] = tmp_invu[i];
      }
      for (size_t i = 0; i < 6LLU; i++) {
        resultw.boxdim[i] = tmp_bdim[i];
      }
    }
    break;
#endif
  }
  
  // Set the frame number based on the position in the series
  result.setFrameNumber(frame_index);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::exportToFile(const std::string &file_name, const CoordinateFileKind kind,
                                       const PrintSituation expectation, const int low_index,
                                       const int high_index, const HybridTargetLevel tier) const {
  if (low_index < 0 || low_index >= frame_count || high_index < low_index ||
      high_index >= frame_count) {
    rtErr("The frame index range " + std::to_string(low_index) + " to " +
          std::to_string(high_index) + " is invalid for a series with " +
          std::to_string(frame_count) + " frames.", "CoordinateSeries", "exportToFile");
  }
  const int actual_high_index = (high_index > low_index) ? high_index : frame_count;
  const PrintSituation aexp = adjustTrajectoryOpeningProtocol(expectation, kind,
                                                              "CoordinateSeries", "exportToFile");
  const DataFormat style = getTrajectoryFormat(kind);
  const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
  std::ofstream foutp;
  foutp = openOutputFile(file_name, aexp, "Open an output file for writing CoordinateSeries "
                         "contents.", style);
  if (fi_exists == false) {
    initializeTrajectory(&foutp, kind, atom_count);
  }
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
    if (kind == CoordinateFileKind::AMBER_INPCRD && actual_high_index - low_index != 1) {
      rtErr("The " + getCoordinateFileKindName(kind) + " file format requires one and only "
            "one frame.  It cannot accept a series of " +
            std::to_string(actual_high_index - low_index) + " frames.", "CoordinateSeries",
            "exportToFile");
    }
    for (int i = low_index; i < actual_high_index; i++) {
      const CoordinateFrame cf = exportFrame(i, tier);
      const CoordinateFrameReader cfr = cf.data();
      writeFrame(&foutp, file_name, kind, atom_count, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                 nullptr, nullptr, nullptr, cfr.unit_cell, cfr.boxdim);
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    rtErr("A restart file cannot be written based on a CoordinateSeries.  The object will not be "
          "able to store both coordinates and velocities needed for checkpointing.",
          "CoordinateFrame", "exportToFile");
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("A trajectory format must be specified.", "CoordinateSeries", "exportToFile");
  }
  foutp.close();
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T> void CoordinateSeries<T>::upload() {
  x_coordinates.upload();
  y_coordinates.upload();
  z_coordinates.upload();
  box_space_transforms.upload();
  inverse_transforms.upload();
  box_dimensions.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void CoordinateSeries<T>::download() {
  x_coordinates.download();
  y_coordinates.download();
  z_coordinates.download();
  box_space_transforms.download();
  inverse_transforms.download();
  box_dimensions.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T> CoordinateSeries<T>::data(const HybridTargetLevel tier) {
  return CoordinateSeriesWriter<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale,
                                   x_coordinates.data(tier), y_coordinates.data(tier),
                                   z_coordinates.data(tier), box_space_transforms.data(tier),
                                   inverse_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> CoordinateSeries<T>::data(const HybridTargetLevel tier) const {
  return CoordinateSeriesReader<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale,
                                   x_coordinates.data(tier), y_coordinates.data(tier),
                                   z_coordinates.data(tier), box_space_transforms.data(tier),
                                   inverse_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrameReader &cfr,
                                              const int atom_start, const int atom_end,
                                              const int frame_index) {
  const int actual_atom_end = (atom_end > atom_start) ? atom_end : cfr.natom;
  if (actual_atom_end - atom_start != atom_count) {
    rtErr("A CoordinateSeries with frames of " + std::to_string(atom_count) + " atoms cannot "
          "accept " + std::to_string(actual_atom_end - atom_start) + " atoms from a system with " +
          std::to_string(cfr.natom) + " atoms.", "CoordinateSeries", "importCoordinateSet");
  }
  
  // Compute the actual frame index that shall be written, and the upper limit of the atoms that
  // will be written.  If the frame index exceeds the current capacity, make more capacity.
  const int actual_frame_index = (frame_index == -1) ? frame_count : frame_index;
  const size_t atom_offset = static_cast<size_t>(actual_frame_index) *
                             roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  const size_t frame_limit = atom_offset + static_cast<size_t>(actual_atom_end - atom_start);

  // Exapnd the capacity to at least 125% of its original size to avoid continuous resizing of the
  // arrays as more frames are added.
  if (actual_frame_index >= frame_capacity) {
    allocate(((actual_frame_index * 5) + 3) / 4);
  }
  T* xptr = x_coordinates.data();
  T* yptr = y_coordinates.data();
  T* zptr = z_coordinates.data();
  double* umat_ptr = box_space_transforms.data();
  double* invu_ptr = inverse_transforms.data();
  double* bdim_ptr = box_dimensions.data();
  size_t jcon = atom_start;
  if (globalpos_scale_bits == 0) {
    for (size_t j = atom_offset; j < frame_limit; j++) {
      xptr[j] = cfr.xcrd[jcon];
      yptr[j] = cfr.ycrd[jcon];
      zptr[j] = cfr.zcrd[jcon];
      jcon++;
    }
  }
  else {
    const double lgpos_scale = globalpos_scale;
    for (size_t j = atom_offset; j < frame_limit; j++) {
      xptr[j] = llround(cfr.xcrd[jcon] * lgpos_scale);
      yptr[j] = llround(cfr.ycrd[jcon] * lgpos_scale);
      zptr[j] = llround(cfr.zcrd[jcon] * lgpos_scale);
      jcon++;
    }
  }
  const size_t xfrm_offset = actual_frame_index * roundUp(static_cast<size_t>(9), warp_size_zu);
  const size_t bdim_offset = actual_frame_index * roundUp(static_cast<size_t>(6), warp_size_zu);
  const size_t xfrm_limit = xfrm_offset + 9LLU;
  const size_t bdim_limit = bdim_offset + 9LLU;
  jcon = 0LLU;
  for (size_t j = xfrm_offset; j < xfrm_limit; j++) {
    umat_ptr[j] = cfr.umat[jcon];
    invu_ptr[j] = cfr.invu[jcon];
    jcon++;
  }
  jcon = 0LLU;
  for (size_t j = bdim_offset; j < bdim_limit; j++) {
    bdim_ptr[j] = cfr.boxdim[jcon];
    jcon++;
  }

  // Update the number of actual frames, if appropriate
  if (actual_frame_index >= frame_count) {
    frame_count = actual_frame_index + 1;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrameReader &cfr,
                                              const int frame_index) {
  importCoordinateSet(cfr, 0, cfr.natom, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrameWriter &cfw,
                                              const int atom_start, const int atom_end,
                                              const int frame_index) {
  importCoordinateSet(CoordinateFrameReader(cfw), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrameWriter &cfw,
                                              const int frame_index) {
  importCoordinateSet(CoordinateFrameReader(cfw), 0, cfw.natom, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrame &cf, const int atom_start,
                                              const int atom_end, const int frame_index) {
  importCoordinateSet(cf.data(), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrame *cf, const int frame_index) {
  importCoordinateSet(cf->data(), 0, cf->getAtomCount(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const PhaseSpace &ps, const int atom_start,
                                              const int atom_end, const int frame_index) {
  importCoordinateSet(CoordinateFrameReader(ps), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const PhaseSpace *ps, const int frame_index) {
  importCoordinateSet(CoordinateFrameReader(ps), 0, ps->getAtomCount(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importFromFile(const std::string &file_name,
                                         const CoordinateFileKind file_kind,
                                         const std::vector<int> &frame_numbers,
                                         const int replica_count, const int frame_index_start) {
  
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
              "file.", "CoordinateSeries", "buildFromFile");
      }
      TextFile tf(file_name);
      std::vector<int> actual_frame_numbers;
      if (frame_numbers.size() == 0LLU) {
        const size_t nframe_detected = countAmberCrdFormatFrames(tf, atom_count, unit_cell);
        actual_frame_numbers.resize(nframe_detected);
        for (int i = 0; i < nframe_detected; i++) {
          actual_frame_numbers[i] = i;
        }
      }
      else {
        actual_frame_numbers.insert(actual_frame_numbers.end(), frame_numbers.begin(),
                                    frame_numbers.end());
      }
      const int nimports = actual_frame_numbers.size();
      const int orig_frame_count = frame_count;
      if (frame_index_start < 0) {
        resize((nimports * replica_count) + frame_count);
      }
      else {
        resize((nimports * replica_count) + frame_index_start);
      }

      // The file must be read in as double precision, despite the multiple modes that a
      // CoordinateSeries can operate in.  Buffer the data in a CoordinateFrame object, feeding
      // its pointers to a low-level overload of the reading function.  Transfer the results to
      // the CoordinateSeries.
      CoordinateFrame tmp_cf(atom_count, unit_cell);
      CoordinateFrameWriter tmp_cfw = tmp_cf.data();
      const CoordinateFrameReader tmp_cfr(tmp_cfw);
      T* xptr = x_coordinates.data();
      T* yptr = y_coordinates.data();
      T* zptr = z_coordinates.data();
      double* umat_ptr = box_space_transforms.data();
      double* invu_ptr = inverse_transforms.data();
      double* bdim_ptr = box_dimensions.data();
      for (int i = 0; i < nimports; i++) {
        readAmberCrdFormat(tf, tmp_cfw.xcrd, tmp_cfw.ycrd, tmp_cfw.zcrd, tmp_cfw.natom, unit_cell,
                           tmp_cfw.umat, tmp_cfw.invu, tmp_cfw.boxdim, actual_frame_numbers[i]);
        for (int j = 0; j < replica_count; j++) {
          if (frame_index_start < 0) {
            importCoordinateSet(tmp_cfr, orig_frame_count + i + (j * nimports));
          }
          else {
            importCoordinateSet(tmp_cfr, frame_index_start + i + (j * nimports));            
          }
        }
      }
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      // The number of atoms need not be known when first reading an Amber ASCII input coordinate
      // or restart file--it is present in those files, and any NetCDF coordinate file.  However,
      // if the number of atoms currently in the CoordinateSeries is not zero and there are one or
      // more frames, assume that it is a legitimate number and enforce that any subsequent files
      // match that number of atoms.
      TextFile tf(file_name);
      const int test_atom_count = getAmberRestartAtomCount(tf);
      if (atom_count > 0 && frame_count > 0 && test_atom_count != atom_count) {
        rtErr("The atom count detected in " + file_name + " (" + std::to_string(test_atom_count) +
              " does not agree with the atom count of " + std::to_string(frame_count) +
              " frames of the existing series (" + std::to_string(atom_count) + " atoms).",
              "CoordinateSeries", "importFromFile");
      }
      else {
        atom_count = test_atom_count;
      }
      if (frame_index_start < 0) {
        resize(replica_count + frame_count);
      }
      else {
        resize(replica_count + frame_index_start);
      }
      CoordinateFrame tmp_cf(atom_count, unit_cell);
      CoordinateFrameWriter tmp_cfw = tmp_cf.data();
      const CoordinateFrameReader tmp_cfr(tmp_cfw);
      T* xptr = x_coordinates.data();
      T* yptr = y_coordinates.data();
      T* zptr = z_coordinates.data();
      double* umat_ptr = box_space_transforms.data();
      double* invu_ptr = inverse_transforms.data();
      double* bdim_ptr = box_dimensions.data();
      getAmberInputCoordinates(tf, tmp_cfw.xcrd, tmp_cfw.ycrd, tmp_cfw.zcrd, tmp_cfw.natom,
                               tmp_cfw.umat, tmp_cfw.invu, tmp_cfw.boxdim);
      const int orig_frame_count = frame_count;
      for (int i = 0; i < replica_count; i++) {
        if (frame_index_start < 0) {
          importCoordinateSet(tmp_cfr, orig_frame_count + i);
        }
        else {
          importCoordinateSet(tmp_cfr, frame_index_start + i);
        }
      }
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " could not be understood.", "CoordinateSeries",
          "importFromFile");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::reserve(const int new_frame_capacity) {
  allocate(new_frame_capacity);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::shrinkToFit() {
  const size_t natom_zu = roundUp<size_t>(atom_count, warp_size_zu);
  const size_t fc_zu    = static_cast<size_t>(frame_count);
  const size_t xfrm_zu  = roundUp<size_t>(9, warp_size_zu);
  const size_t bdim_zu  = roundUp<size_t>(6, warp_size_zu);
  x_coordinates.resize(natom_zu * fc_zu);
  y_coordinates.resize(natom_zu * fc_zu);
  z_coordinates.resize(natom_zu * fc_zu);
  box_space_transforms.resize(xfrm_zu * fc_zu);
  inverse_transforms.resize(xfrm_zu * fc_zu);
  box_dimensions.resize(xfrm_zu * fc_zu);
  x_coordinates.shrinkToFit();
  y_coordinates.shrinkToFit();
  z_coordinates.shrinkToFit();
  box_space_transforms.shrinkToFit();
  inverse_transforms.shrinkToFit();
  box_dimensions.shrinkToFit();
  frame_capacity = frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count) {
  allocate(new_frame_count);
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrameReader &cfr,
                                 const int atom_start, const int atom_end) {
  allocate(new_frame_count);
  const int orig_frame_count = frame_count;
  for (int i = orig_frame_count; i < new_frame_count; i++) {
    importCoordinateSet(cfr, atom_start, atom_end, i);
  }
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrameWriter &cfw,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrame &cf,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, CoordinateFrame *cf,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const PhaseSpace &ps,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, PhaseSpace *ps,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrameReader &cfr, const int atom_start,
                                   const int atom_end) {
  const int orig_frame_count = frame_count;
  if (frame_count >= frame_capacity) {
    allocate(((frame_count * 5) + 3) / 4);
  }
  importCoordinateSet(cfr, atom_start, atom_end, orig_frame_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrameWriter &cfw, const int atom_start,
                                   const int atom_end) {
  pushBack(CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrame &cf, const int atom_start,
                                   const int atom_end) {
  pushBack(cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(CoordinateFrame *cf, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const PhaseSpace &ps, const int atom_start,
                                   const int atom_end) {
  pushBack(CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(PhaseSpace *ps, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::allocate(const int new_frame_capacity) {
  if (new_frame_capacity > frame_capacity) {
    frame_capacity = new_frame_capacity;
    const size_t fc_zu = static_cast<size_t>(frame_capacity);
    const size_t total_atoms = fc_zu * roundUp(static_cast<size_t>(atom_count), warp_size_zu);
    const size_t total_xfrm = fc_zu * roundUp<size_t>(9, warp_size_zu);
    const size_t total_bdim = fc_zu * roundUp<size_t>(6, warp_size_zu);
    
    // Allocate space for the new capacity.  This will reallocate each array, but one by one in
    // order to avoid keeping what could be nearly double the amount of data at any one time.
    x_coordinates.resize(total_atoms);
    y_coordinates.resize(total_atoms);
    z_coordinates.resize(total_atoms);
    box_space_transforms.resize(total_xfrm);
    inverse_transforms.resize(total_xfrm);
    box_dimensions.resize(total_bdim);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig, typename Tnew>
CoordinateSeries<Tnew> changeCoordinateSeriesType(const CoordinateSeries<Torig> &cs,
                                                  const int globalpos_scale_bits_in) {
  const CoordinateSeriesReader<Torig> csr = cs.data();
  CoordinateSeries<Tnew> result(csr.natom, csr.nframe, csr.unit_cell, globalpos_scale_bits_in);
  CoordinateSeriesWriter<Tnew> resultw = result.data();
  const size_t natom_zu        = csr.natom;
  const size_t padded_natom_zu = roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t padded_xfrm_zu  = roundUp<size_t>(9, warp_size_zu);
  const size_t padded_bdim_zu  = roundUp<size_t>(6, warp_size_zu);
  const size_t fc_zu           = static_cast<size_t>(csr.nframe);
  if (isFloatingPointScalarType<Torig>() && isFloatingPointScalarType<Tnew>()) {
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = csr.xcrd[pos];
        resultw.ycrd[pos] = csr.ycrd[pos];
        resultw.zcrd[pos] = csr.zcrd[pos];
      }
    }
  }
  else if (isFloatingPointScalarType<Torig>() && isSignedIntegralScalarType<Tnew>()) {
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = llround(csr.xcrd[pos] * resultw.gpos_scale);
        resultw.ycrd[pos] = llround(csr.ycrd[pos] * resultw.gpos_scale);
        resultw.zcrd[pos] = llround(csr.zcrd[pos] * resultw.gpos_scale);
      }
    }
  }
  else if (isFloatingPointScalarType<Tnew>() && isSignedIntegralScalarType<Torig>()) {
    const double conv_factor = pow(2.0, -globalpos_scale_bits_in);
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = static_cast<double>(csr.xcrd[pos]) * csr.inv_gpos_scale;
        resultw.ycrd[pos] = static_cast<double>(csr.ycrd[pos]) * csr.inv_gpos_scale;
        resultw.zcrd[pos] = static_cast<double>(csr.zcrd[pos]) * csr.inv_gpos_scale;
      }
    }
  }
  else if (isSignedIntegralScalarType<Tnew>() && isSignedIntegralScalarType<Torig>()) {
    if (globalpos_scale_bits_in == csr.gpos_bits) {
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          resultw.xcrd[pos] = static_cast<Tnew>(csr.xcrd[pos]);
          resultw.ycrd[pos] = static_cast<Tnew>(csr.ycrd[pos]);
          resultw.zcrd[pos] = static_cast<Tnew>(csr.zcrd[pos]);
        }
      }
    }
    else if (globalpos_scale_bits_in > csr.gpos_bits) {
      const int forward_shift = globalpos_scale_bits_in - csr.gpos_bits;
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          const llint ixval = csr.xcrd[pos];
          const llint iyval = csr.ycrd[pos];
          const llint izval = csr.zcrd[pos];
          resultw.xcrd[pos] = static_cast<Tnew>(ixval << forward_shift);
          resultw.ycrd[pos] = static_cast<Tnew>(iyval << forward_shift);
          resultw.zcrd[pos] = static_cast<Tnew>(izval << forward_shift);
        }
      }
    }
    else {
      const int backward_shift = globalpos_scale_bits_in - csr.gpos_bits;
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          const llint ixval = csr.xcrd[pos];
          const llint iyval = csr.ycrd[pos];
          const llint izval = csr.zcrd[pos];
          resultw.xcrd[pos] = static_cast<Tnew>(ixval >> backward_shift);
          resultw.ycrd[pos] = static_cast<Tnew>(iyval >> backward_shift);
          resultw.zcrd[pos] = static_cast<Tnew>(izval >> backward_shift);
        }
      }
    }
  }
  for (size_t i = 0; i < fc_zu; i++) {
    for (size_t j = 0; j < 9LLU; j++) {
      const size_t pos = (i * padded_xfrm_zu) + j;
      resultw.umat[pos] = csr.umat[pos];
      resultw.invu[pos] = csr.invu[pos];
    }
    for (size_t j = 0; j < 6LLU; j++) {
      const size_t pos = (i * padded_bdim_zu) + j;
      resultw.boxdim[pos] = csr.boxdim[pos];
    }
  }
  
  return result;
}

} // namespace trajectory
} // namespace stormm
