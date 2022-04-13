// -*-c++-*-
namespace omni {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const int natom_in, const int nframe_in,
                                   const UnitCellType unit_cell_in) :
    atom_count{natom_in}, frame_count{nframe_in}, frame_capacity{nframe_in},
    frame_numbers{nframe_in, "cser_frame_nums"},
    x_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_x_coords"},
    y_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_y_coords"},
    z_coordinates{static_cast<size_t>(nframe_in) * static_cast<size_t>(natom_in), "cser_z_coords"},
    box_space_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                         "cser_umat"},
    inverse_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                       "cser_invu"},
    box_dimensions{static_cast<size_t>(nframe_in) * roundUp<size_t>(6, warp_size_zu),
                   "cser_boxdims"},
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
                                   const std::vector<int> &frame_numbers_in,
                                   const int replica_count_in, const int atom_count_in) :
{}

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
    allocate(((actual_frame_index * 5) + 3) / 4);
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
void CoordinateSeries::reserve(const int new_frame_capacity) {
  allocate(new_frame_capacity);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count) {
  allocate(new_frame_count);
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, const CoordinateFrameReader &cfr,
                              const int atom_start, const int atom_end) {
  allocate(new_frame_count);
  const int orig_frame_count = frame_count;
  for (int i = orig_frame_count; i < new_frame_count; i++) {
    importCoordinateSet(cfr, atom_start, atom_end, i);    
  }
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, const CoordinateFrameWriter &cfw,
                              const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, const CoordinateFrame &cf,
                              const int atom_start, const int atom_end) {
  resize(new_frame_count, cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, CoordinateFrame *cf, const int atom_start,
                              const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, const PhaseSpace &ps,
                              const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::resize(const int new_frame_count, PhaseSpace *ps,
                              const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(const CoordinateFrameReader &cfr, const int atom_start,
                                const int atom_end) {
  if (frame_count >= frame_capacity) {
    allocate(((frame_count * 5) + 3) / 4);
  }
  importCoordinateSet(cfr, atom_start, atom_end, frame_count);
  frame_count++;
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(const CoordinateFrameWriter &cfw, const int atom_start,
                                const int atom_end) {
  pushBack(CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(const CoordinateFrame &cf, const int atom_start,
                                const int atom_end) {
  pushBack(cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(CoordinateFrame *cf, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(const PhaseSpace &ps, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::pushBack(PhaseSpace *ps, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
void CoordinateSeries::allocate(const int new_frame_capacity) {
  if (new_frame_capacity > frame_capacity) {
    frame_capacity = new_frame_capacity;
    const size_t fc_zu = static_cast<size_t>(frame_capacity);
    const size_t total_atoms = fc_zu * roundUp(static_cast<size_t>(atom_count), warp_size_zu);
    const size_t total_xfrm = fc_zu * roundUp<size_t>(9, warp_size_zu);
    const size_t total_bdim = fc_zu * roundUp<size_t>(6, warp_size_zu);

    // Allocate space for the new capacity.  This will reallocate each array, but one by one in
    // order to avoid keeping what could be nearly double the amount of data at any one time.
    frame_numbers.resize(frame_capacity);
    x_coordinates.resize(total_atoms);
    y_coordinates.resize(total_atoms);
    z_coordinates.resize(total_atoms);
    box_space_transforms.resize(total_xfrm);
    inverse_transforms.resize(total_xfrm);
    box_dimensions.resize(total_bdim);

    // Reset the array sizes to the current storage 
    const size_t fn_zu = static_cast<size_t>(frame_count);
    const size_t active_atoms = fn_zu * roundUp(static_cast<size_t>(atom_count), warp_size_zu);
    const size_t active_xfrm = fn_zu * roundUp<size_t>(9, warp_size_zu);
    const size_t active_bdim = fn_zu * roundUp<size_t>(6, warp_size_zu);
    frame_numbers.resize(frame_count);
    x_coordinates.resize(active_atoms);
    y_coordinates.resize(active_atoms);
    z_coordinates.resize(active_atoms);
    box_space_transforms.resize(active_xfrm);
    inverse_transforms.resize(active_xfrm);
    box_dimensions.resize(active_bdim);
  }
}

} // namespace trajectory
} // namespace omni
