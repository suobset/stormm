// -*-c++-*-
namespace omni {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
CoordinateSeries::CoordinateSeries(const int natom_in, const int nframe_in,
                                   const UnitCellType unit_cell_in) :
    file_name{std::string("")}, atom_count{natom_in}, frame_count{nframe_in},
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
CoordinateSeries::CoordinateSeries(const std::string &file_name_in,
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
void CoordinateSeries::buildFromFile(const std::string &file_name_in,
                                     const CoordinateFileKind file_kind,
                                     const std::vector<int> &frame_numbers,
                                     const int replica_count_in, const int atom_count_in) {
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
          const size_t atom_offset = i * roundUp(static_cast<size_t>(atom_count), warp_size_zu);
          const size_t frame_limit = atom_offset + static_cast<size_t>(atom_count);
          size_t jcon = 0LLU;
          for (size_t j = atom_offset; j < frame_limit; j++) {
            xptr[j] = tmp_cfw.xcrd[jcon];
            yptr[j] = tmp_cfw.ycrd[jcon];
            zptr[j] = tmp_cfw.zcrd[jcon];
            jcon++;
          }
          const size_t xfrm_offset = i * roundUp(static_cast<size_t>(9), warp_size_zu);
          const size_t bdim_offset = i * roundUp(static_cast<size_t>(6), warp_size_zu);
          const size_t xfrm_limit = xfrm_offset + 9LLU;
          const size_t bdim_limit = bdim_offset + 9LLU;
          jcon = 0LLU;
          for (size_t j = xfrm_offset; j < xfrm_limit; j++) {
            umat_ptr[j] = tmp_cfw.umat[jcon];
            invu_ptr[j] = tmp_cfw.invu[jcon];
            jcon++;
          }
          jcon = 0LLU;
          for (size_t j = bdim_offset; j < bdim_limit; j++) {
            bdim_ptr[j] = tmp_cfw.boxdim[jcon];
            jcon++;
          }
        }
      }
    }
  case CoordinateFileKind::AMBER_CRD:
}
  
} // namespace trajectory
} // namespace omni
