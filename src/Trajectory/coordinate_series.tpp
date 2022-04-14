// -*-c++-*-
namespace omni {
namespace trajectory {

using math::roundUp;

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T>::CoordinateSeriesWriter(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in, T* xcrd_in,
                                                  T* ycrd_in, T* zcrd_in, double* umat_in,
                                                  double* invu_in, double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesReader<T>::CoordinateSeriesReader(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in,
                                                  const T* xcrd_in, const T* ycrd_in,
                                                  const T* zcrd_in, const double* umat_in,
                                                  const double* invu_in, const double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const int natom_in, const int nframe_in,
                                      const UnitCellType unit_cell_in) :
    atom_count{natom_in}, frame_count{nframe_in}, frame_capacity{nframe_in},
    unit_cell{unit_cell_in},
    x_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_x_coords"},
    y_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_y_coords"},
    z_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_z_coords"},
    box_space_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                         "cser_umat"},
    inverse_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                       "cser_invu"},
    box_dimensions{static_cast<size_t>(nframe_in) * roundUp<size_t>(6, warp_size_zu),
                   "cser_boxdims"}
{
  allocate(frame_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const std::string &file_name, const int atom_count_in,
                                      const CoordinateFileKind file_kind,
                                      const std::vector<int> &frame_numbers,
                                      const int replica_count, const UnitCellType unit_cell_in) :
    CoordinateSeries(atom_count_in, replica_count * static_cast<int>(frame_numbers.size()),
                     unit_cell_in)
{
  importFromFile(file_name, file_kind, frame_numbers, replica_count, 0);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(PhaseSpace *ps, const int nframe_in) :
    CoordinateSeries(ps->getAtomCount(), nframe_in, ps->getUnitCellType())
{
  allocate(frame_count);
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(ps, 0, ps->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const PhaseSpace &ps, const int nframe_in) :
    CoordinateSeries(ps.getAtomCount(), nframe_in, ps.getUnitCellType())
{
  allocate(frame_count);
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(ps, 0, ps.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(CoordinateFrame *cf, const int nframe_in) :
    CoordinateSeries(cf->getAtomCount(), nframe_in, cf->getUnitCellType())
{
  allocate(frame_count);
  for (int i = 0; i < frame_count; i++) {
    importCoordinateSet(cf, 0, cf->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const CoordinateFrame &cf, const int nframe_in) :
    CoordinateSeries(cf.getAtomCount(), nframe_in, cf.getUnitCellType())
{
  allocate(frame_count);
  for (int i = 0; i < frame_count; i++) {
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
template <typename T> std::vector<T>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index,
                                              const HybridTargetLevel tier) const {

  return getInterlacedCoordinates(frame_index, 0, atom_count, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index, const int low_index,
                                              const int high_index,
                                              const HybridTargetLevel tier) const {

  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "CoordinateSeries",
          "getInterlacedCoordinates");
  }
  std::vector<double> result(3 * (high_index - low_index));
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* xptr = x_coordinates.data();
      const double* yptr = y_coordinates.data();
      const double* zptr = z_coordinates.data();
      const size_t frame_offset = natom_zu * fidx_zu;
      for (int i = low_index; i < high_index; i++) {
        const int base_idx = 3 * (i - low_index);
        const size_t access_idx = frame_offset + static_cast<size_t>(i);
        result[base_idx    ] = xptr[access_idx];
        result[base_idx + 1] = yptr[access_idx];
        result[base_idx + 2] = zptr[access_idx];
      }
    }
    break;
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t llim = (fidx_zu * natom_zu) + static_cast<size_t>(low_index);
      const size_t hlim = (fidx_zu * natom_zu) + static_cast<size_t>(high_index);
      const std::vector<double> xval = x_coordinates.readDevice(llim, hlim);
      const std::vector<double> yval = y_coordinates.readDevice(llim, hlim);
      const std::vector<double> zval = z_coordinates.readDevice(llim, hlim);
      const size_t frame_offset = natom_zu * fidx_zu;
      for (int i = 0; i < high_index - low_index; i++) {
        result[(3 * i)    ] = xval[i];
        result[(3 * i) + 1] = yval[i];
        result[(3 * i) + 2] = zval[i];
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
#ifdef OMNI_USE_HPC
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
#ifdef OMNI_USE_HPC
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
#ifdef OMNI_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_dimensions.readDevice(read_start, 6);
    break;
#endif
  }
  __builtin_unreachable();
}

#ifdef OMNI_USE_HPC
//-------------------------------------------------------------------------------------------------
void CoordinateSeries::upload() {
  x_coordinates.upload();
  y_coordinates.upload();
  z_coordinates.upload();
  box_space_transforms.upload();
  inverse_transforms.upload();
  box_dimensions.upload();
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::download() {
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
  return CoordinateSeriesWriter<T>(atom_count, frame_count, unit_cell, x_coordinates.data(tier),
                                   y_coordinates.data(tier), z_coordinates.data(tier),
                                   box_space_transforms.data(tier), inverse_transforms.data(tier),
                                   box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> CoordinateSeries<T>::data(const HybridTargetLevel tier) const {
  return CoordinateSeriesReader<T>(atom_count, frame_count, unit_cell, x_coordinates.data(tier),
                                   y_coordinates.data(tier), z_coordinates.data(tier),
                                   box_space_transforms.data(tier), inverse_transforms.data(tier),
                                   box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importCoordinateSet(const CoordinateFrameReader &cfr,
                                              const int atom_start, const int atom_end,
                                              const int frame_index) {
  const int actual_atom_end = (atom_end > atom_start) ? atom_end : cfr.natom;

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
  for (size_t j = atom_offset; j < frame_limit; j++) {
    xptr[j] = cfr.xcrd[jcon];
    yptr[j] = cfr.ycrd[jcon];
    zptr[j] = cfr.zcrd[jcon];
    jcon++;
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
  importCoordinateSet(CoordinateFrameReader(CoordinateFrameWriter(ps)), 0, ps->getAtomCount(),
                      frame_index);
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

} // namespace trajectory
} // namespace omni
