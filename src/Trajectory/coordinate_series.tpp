// -*-c++-*-
namespace omni {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const int natom_in, const int nframe_in,
                                   const UnitCellType unit_cell_in) :
    atom_count{natom_in}, frame_count{nframe_in}, frame_capacity{nframe_in},
    frame_numbers{nframe_in, "cser_frame_nums"},
    x_coordinates{HybridKind::POINTER, "cser_x_coords"},
    y_coordinates{HybridKind::POINTER, "cser_y_coords"},
    z_coordinates{HybridKind::POINTER, "cser_z_coords"},
    box_space_transforms{HybridKind::POINTER, "cser_umat"},
    inverse_transforms{HybridKind::POINTER, "cser_invu"},
    box_dimensions{HybridKind::POINTER, "cser_boxdims"},
    particle_data{static_cast<size_t>(natom_in) * static_cast<size_t>(nframe_in),
                  "cser_particle_data"},
    box_data{nframe_in * ((2 * roundUp(9, warp_size_int)) + roundUp(6, warp_size_int)),
             "cser_box_data"}
{}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const std::string &file_name,
                                   const CoordinateFileKind file_kind,
                                   const std::vector<int> &frame_numbers_in,
                                   const int replica_count_in, const int atom_count_in) :
    CoordinateSeries(atom_count_in, replica_count_in * static_cast<int>(frame_numbers_in.size()))
{

}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const TextFile &tf, const CoordinateFileKind file_kind,
                                   const std::vector<int> &frame_numbers_in) {
}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(PhaseSpace *ps, const int nframe_in) {
}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const PhaseSpace &ps, const int nframe_in) {

}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(CoordinateFrame *cf, const int nframe_in) {

}

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const CoordinateFrame &cf, const int nframe_in) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::importCoordinateSet(const CoordinateFrameReader &cfr,
                                           const int atom_start, const int atom_end,
                                           const int frame_index) {

  // Compute the actual frame index that shall be written, and the upper limit of the atoms that
  // will be written.  If the frame index exceeds the current capacity, make more capacity.
  const int actual_frame_index = (frame_index == -1) ? frame_count : frame_index;
  const size_t atom_offset = static_cast<size_t>(actual_frame_index) *
                             roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  const size_t frame_limit = atom_offset + static_cast<size_t>(atom_end - atom_start);

  // Exapnd the capacity to at least 125% of its original size to avoid continuous resizing of the
  // arrays as more frames are added.
  if (actual_frame_index >= frame_capacity) {
    frame_capacity = ((actual_frame_index * 5) + 3) / 4;
    reallocate();
  }

  T* xptr = x_coordinates->data();
  T* yptr = y_coordinates->data();
  T* zptr = z_coordinates->data();
  double* umat_ptr = box_space_transforms->data();
  double* invu_ptr = inverse_transforms->data();
  double* bdim_ptr = box_dimensions->data();
  size_t jcon = atom_start;
  for (size_t j = atom_offset; j < frame_limit; j++) {
    xptr[j] = cfr.xcrd[jcon];
    yptr[j] = cfr.ycrd[jcon];
    zptr[j] = cfr.zcrd[jcon];
    jcon++;
  }
  const size_t xfrm_offset = i * roundUp(static_cast<size_t>(9), warp_size_zu);
  const size_t bdim_offset = i * roundUp(static_cast<size_t>(6), warp_size_zu);
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
void CoordinateSeries::importCoordinateSet(const CoordinateFrameReader &cfr,
                                           const int frame_index) {
  importCoordinateSet(cfr, 0, cfr.atom_count, frame_index);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::importFromFile(const std::string &file_name,
                                      const CoordinateFileKind file_kind,
                                      const std::vector<int> &frame_numbers,
                                      const int replica_count_in, const int frame_index_start) {

  // Check that there is sufficient space in the object.  If there is not, allocate more.
  if (frame_count < 
  
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
      allocate();
      TextFile tf(file_name);
      if (frame_numbers.size() == 0LLU) {
        const size_t nframe_detected = countAmberCrdFormatFrames(tf, atom_count, unit_cell);

        // The file must be read in as double precision, despite the multiple modes that a
        // CoordinateSeries can operate in.  Buffer the data in a CoordinateFrame object, feeding
        // its pointers to a low-level overload of the reading function.  Transfer the results to
        // the CoordinateSeries.
        CoordinateFrame tmp_cf(atom_count, unit_cell);
        CoordinateFrameWriter tmp_cfw = tmp_cf.data();
        T* xptr = x_coordinates->data();
        T* yptr = y_coordinates->data();
        T* zptr = z_coordinates->data();
        double* umat_ptr = box_space_transforms->data();
        double* invu_ptr = inverse_transforms->data();
        double* bdim_ptr = box_dimensions->data();
        for (size_t i = 0; i < nframe_detected; i++) {
          readAmberCrdFormat(tf, tmp_cfw.xcrd, tmp_cfw.ycrd, tmp_cfw.zcrd, unit_cell,
                             tmp_cfw.umat, tmp_cfw.invu, tmp_cfw.boxdim, i);
        }
      }
    }
    break;
  case CoordinateFileKind::AMBER_CRD:
  }
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::allocate() {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::reallocate() {
  
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const CoordinateFrameReader &cfr,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const CoordinateFrameWriter &cfw,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const CoordinateFrame &cf,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const CoordinateFrame *cf,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const PhaseSpace &ps,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count_in, const PhaseSpace *ps,
                              const int atom_start, const int atom_end) {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::reserve() {

}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack() {

}

} // namespace trajectory
} // namespace omni
