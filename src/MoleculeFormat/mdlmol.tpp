// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const T* xcrd, const T* ycrd, const T* zcrd,
               const double inv_scale, const int molecule_index) :
  MdlMol(ExceptionResponse::SILENT)
{
  const AtomGraph *ag = chemfe->getTopologyPointer();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  allocate(cdk, nbk, molecule_index);

  // Fill in the coordinates
  impartCoordinates(xcrd, ycrd, zcrd, inv_scale, cdk, molecule_index);

  // Fill in additional arrays
  transferTopologicalDetails(chemfe, molecule_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const T* xcrd, const T* ycrd, const T* zcrd,
               const int molecule_index) :
  MdlMol(chemfe, xcrd, ycrd, zcrd, 1.0, molecule_index)
{
  // Check that the data type was a floating-point type
  if (isFloatingPointScalarType<T>() == false) {
    rtErr("Data from integral types is for fixed-precision representations and must be "
          "accompanied by a scaling factor to take the values back into real, internal units.",
          "MdlMol");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMol::impartCoordinates(const T* xcrd, const T* ycrd, const T* zcrd,
                               const double scale_factor) {
  if (scale_factor != 1.0) {
    for (int i = 0; i < atom_count; i++) {
      coordinates[i].x = xcrd[i] * scale_factor;
      coordinates[i].y = ycrd[i] * scale_factor;
      coordinates[i].z = zcrd[i] * scale_factor;
    }
  }
  else {
    for (int i = 0; i < atom_count; i++) {
      coordinates[i].x = xcrd[i];
      coordinates[i].y = ycrd[i];
      coordinates[i].z = zcrd[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMol::impartCoordinates(const T* xcrd, const T* ycrd, const T* zcrd,
                               const double scale_factor, const ChemicalDetailsKit &cdk,
                               const int molecule_index) {
  if (molecule_index < 0 || molecule_index >= cdk.nmol) {
    rtErr("Molecule index " + std::to_string(molecule_index) + " is invalid for a topology with " +
          std::to_string(cdk.nmol) + " molecules.", "MdlMol", "impartCoordinates");
  }
  const int llim = cdk.mol_limits[molecule_index];
  const int hlim = cdk.mol_limits[molecule_index + 1];
  if (hlim - llim != atom_count) {
    rtErr("Molecule " + std::to_string(molecule_index) + " with " + std::to_string(hlim - llim) +
          " atoms is incompatible with an MDL MOL object of " +	std::to_string(atom_count) +
          " atoms.", "MdlMol", "impartCoordinates");
  }
  if (scale_factor != 1.0) {
    for (int i = llim; i < hlim; i++) {
      const int iatom = cdk.mol_contents[i];
      coordinates[i - llim].x = xcrd[iatom] * scale_factor;
      coordinates[i - llim].y = ycrd[iatom] * scale_factor;
      coordinates[i - llim].z = zcrd[iatom] * scale_factor;
    }
  }
  else {
    for (int i = llim; i < hlim; i++) {
      const int iatom = cdk.mol_contents[i];
      coordinates[i - llim].x = xcrd[iatom];
      coordinates[i - llim].y = ycrd[iatom];
      coordinates[i - llim].z = zcrd[iatom];
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMol::impartCoordinates(const CoordinateSeriesReader<T> &csr, const int frame_index,
                               const HybridTargetLevel tier) {
  if (frame_index < 0 || frame_index >= csr.nframe) {
    rtErr("Invalid frame index request " + std::to_string(frame_index) + " for a series of " +
          std::to_string(csr.nframe) + " frames.", "MdlMol", "impartCoordinates");
  }
  checkAtomCount(csr.natom);
  const size_t padded_natom = roundUp(csr.natom, warp_size_int);
  const size_t offset = static_cast<size_t>(padded_natom) * static_cast<size_t>(frame_index);
  impartCoordinates(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset], csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMol::impartCoordinates(const CoordinateSeriesWriter<T> &csw, const int frame_index,
                               const HybridTargetLevel tier) {
  if (frame_index < 0 || frame_index >= csw.nframe) {
    rtErr("Invalid frame index request " + std::to_string(frame_index) + " for a series of " +
          std::to_string(csw.nframe) + " frames.", "MdlMol", "impartCoordinates");
  }
  checkAtomCount(csw.natom);
  const size_t padded_natom = roundUp(csw.natom, warp_size_int);
  const size_t offset = static_cast<size_t>(padded_natom) * static_cast<size_t>(frame_index);
  impartCoordinates(&csw.xcrd[offset], &csw.ycrd[offset], &csw.zcrd[offset], csw.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MdlMol::impartCoordinates(const CoordinateSeries<T> &cs, const int frame_index,
                               const HybridTargetLevel tier) {
  impartCoordinates(cs.data(tier), frame_index);
}
  
} // namespace structure
} // namespace stormm
