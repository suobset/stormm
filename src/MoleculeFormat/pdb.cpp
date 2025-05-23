#include "copyright.h"
#include "Chemistry/periodic_table.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "FileManagement/file_util.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Parsing/parse.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinate_util.h"
#include "Trajectory/trajectory_enumerators.h"
#include "molecule_file_io.h"
#include "pdb.h"

namespace stormm {
namespace structure {

using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::elemental_symbols;
using constants::CartesianDimension;
using constants::CaseSensitivity;
using constants::PrecisionModel;
using diskutil::openOutputFile;
using parse::strcmpCased;
using parse::TextFile;
using parse::TextFileReader;
using parse::uppercase;
using stmath::computeBoxTransform;
using stmath::roundUp;
using synthesis::CondensateReader;
using synthesis::PsSynthesisReader;
using topology::ChemicalDetailsKit;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::determineUnitCellType;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
using trajectory::TrajectoryKind;
  
//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const int atom_count_in, const int model_count_in) :
    file_name{std::string("")}, atom_count{atom_count_in}, model_count{model_count_in},
    padded_atom_count{roundUp(atom_count_in, warp_size_int)}, ter_card_count{0},
    unit_cell{UnitCellType::NONE},
    atom_classes{}, atom_numbers{}, atom_names{}, residue_names{}, chain_names{},
    residue_numbers{}, x_coordinates{}, y_coordinates{}, z_coordinates{}, box_dimensions{},
    occupancies{}, b_factors{}, atom_symbols{}, atom_formal_charges{}, ter_card_locations{},
    anisotropy_xx{}, anisotropy_xy{}, anisotropy_xz{}, anisotropy_yy{}, anisotropy_yz{},
    anisotropy_zz{},
    print_formal_charges{false}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const std::string &file_name_in, const int model_index,
         const std::vector<char> &alternate_prefs) :
    Pdb(0, 1)
{
  file_name = file_name_in;
  TextFile tf(file_name);
  TextFileReader tfr = tf.data();

  // Scan for the number of models in the file.  Check that the requested entry is reasonable.
  const std::vector<int2> model_line_lims = findEntryLimits(tf, std::string("MODEL"));
  model_count = 1;
  if (model_index >= static_cast<int>(model_line_lims.size())) {
    rtErr("Model entry " + std::to_string(model_index) + " is invalid for a PDB file with " +
          std::to_string(model_line_lims.size()) + " detected models.", "Pdb");
  }
  else if (model_index >= 0) {

    // Read the specified model from the PDB file.
    const int line_start = model_line_lims[model_index].x;
    const int line_end   = model_line_lims[model_index].y;
    atom_count = countAtomRecords(tf, line_start, line_end);
    padded_atom_count = roundUp(atom_count, warp_size_int);
    allocate();
    readAtomRecords(tf, line_start, line_end, 0);
    seekUnitCellDimensions(tf, line_start, line_end, 0);
    if (box_dimensions[3] > 1.0e-6) {
      unit_cell = determineUnitCellType<double>(box_dimensions);
    }
  }
  else {

    // Read all models in the PDB file.
    model_count = model_line_lims.size();
    const int xfrm_stride = roundUp(6, warp_size_int);
    for (int i = 0; i < model_count; i++) {
      const int line_start = model_line_lims[i].x;
      const int line_end   = model_line_lims[i].y;
      if (i == 0) {
        atom_count = countAtomRecords(tf, line_start, line_end);
        padded_atom_count = roundUp(atom_count, warp_size_int);
        allocate();
      }
      else {
        const int i_natom = countAtomRecords(tf, line_start, line_end);
        if (atom_count != i_natom) {
          rtErr("Model entry " + std::to_string(i) + " contains " + std::to_string(i_natom) +
                " atoms, whereas " + std::to_string(atom_count) + " are expected based on prior "
                "models read from PDB file " + file_name_in + ".", "Pdb");
        }
      }
      readAtomRecords(tf, line_start, line_end, i * padded_atom_count);
      seekUnitCellDimensions(tf, line_start, line_end, i * xfrm_stride);
      if (box_dimensions[(i * xfrm_stride) + 3] > 1.0e-6) {
        unit_cell = determineUnitCellType<double>(box_dimensions);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const std::string &file_name_in, const std::vector<int> &model_indices,
         const std::vector<char> &alternate_prefs) :
    Pdb(0, 1)
{
  file_name = file_name_in;
  TextFile tf(file_name);
  TextFileReader tfr = tf.data();

  // Scan for the number of models in the file
  const std::vector<int2> model_line_lims = findEntryLimits(tf, std::string("MODEL"));
  model_count = model_indices.size();
  allocate();
  for (int i = 0; i < model_count; i++) {
    if (model_indices[i] < 0 || model_indices[i] >= model_line_lims.size()) {
      rtErr("Model entry " + std::to_string(model_indices[i]) + " is invalid for a PDB file "
            "with " + std::to_string(model_line_lims.size()) + " detected models.", "Pdb");
    }
    const int line_start = model_line_lims[model_indices[i]].x;
    const int line_end   = model_line_lims[model_indices[i]].y;
    if (i == 0) {
      atom_count = countAtomRecords(tf, line_start, line_end);
      padded_atom_count = roundUp(atom_count, warp_size_int);
    }
    else {
      const int i_natom = countAtomRecords(tf, line_start, line_end);
      if (atom_count != i_natom) {
        rtErr("Model entry " + std::to_string(model_indices[i]) + " contains " +
              std::to_string(i_natom) + " atoms, whereas " + std::to_string(atom_count) +
              " are expected based on prior models read from PDB file " + file_name_in + ".",
              "Pdb");
      }
    }
    readAtomRecords(tf, line_start, line_end, i * padded_atom_count, alternate_prefs);
  }
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const AtomGraph *ag, const CoordinateFrame *cf) :
    Pdb(ag->getAtomCount(), 1)
{
  loadTopologicalData(*ag);
  loadCoordinates(cf);
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const AtomGraph &ag, const CoordinateFrame &cf) :
    Pdb(ag.getAtomCount(), 1)
{
  loadTopologicalData(ag.getSelfPointer());
  loadCoordinates(cf);
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const AtomGraph *ag, const PhaseSpace *ps) :
    Pdb(ag->getAtomCount(), 1)
{
  loadTopologicalData(*ag);
  loadCoordinates(ps);
}

//-------------------------------------------------------------------------------------------------
Pdb::Pdb(const AtomGraph &ag, const PhaseSpace &ps) :
    Pdb(ag.getAtomCount(), 1)
{
  loadTopologicalData(ag.getSelfPointer());
  loadCoordinates(ps);
}

//-------------------------------------------------------------------------------------------------
int Pdb::getAtomCount() const {
  return atom_count;
}
  
//-------------------------------------------------------------------------------------------------
int Pdb::getModelCount() const {
  return model_count;
}

//-------------------------------------------------------------------------------------------------
bool Pdb::getFormalChargePrinting() const {
  return print_formal_charges;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame Pdb::exportCoordinateFrame(const int model_index,
                                           const HybridFormat layout) const {
  UnitCellType unit_cell;
  const int xfrm_offset = model_index * roundUp(6, warp_size_int);
  const int mdl_offset = padded_atom_count * model_index;
  CoordinateFrame result(atom_count,
                         determineUnitCellType<double>(&box_dimensions.data()[xfrm_offset]),
                         layout);
  std::vector<double> mdl_umat, mdl_invu;
  computeTransforms(model_index, &mdl_umat, &mdl_invu);
  switch (layout) {
  case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
#endif
    {
      CoordinateFrameWriter rw = result.data();
      for (int i = 0; i < atom_count; i++) {
        rw.xcrd[i] = x_coordinates[mdl_offset + i];
        rw.ycrd[i] = y_coordinates[mdl_offset + i];
        rw.zcrd[i] = z_coordinates[mdl_offset + i];
      }
      for (int i = 0; i < 9; i++) {
        rw.umat[i] = mdl_umat[i];
        rw.invu[i] = mdl_invu[i];
      }
      for (int i = 0; i < 6; i++) {
        rw.boxdim[i] = box_dimensions[xfrm_offset + i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::DEVICE_ONLY:
    {
      Hybrid<double>* dvc_x = result.getCoordinateHandle(CartesianDimension::X);
      Hybrid<double>* dvc_y = result.getCoordinateHandle(CartesianDimension::Y);
      Hybrid<double>* dvc_z = result.getCoordinateHandle(CartesianDimension::Z);
      Hybrid<double>* dvc_dims = result.getBoxDimensionsHandle();
      Hybrid<double>* dvc_umat = result.getBoxTransformHandle();
      Hybrid<double>* dvc_invu = result.getInverseTransformHandle();
      dvc_x->putDevice(&x_coordinates[mdl_offset], 0, atom_count);
      dvc_y->putDevice(&y_coordinates[mdl_offset], 0, atom_count);
      dvc_z->putDevice(&z_coordinates[mdl_offset], 0, atom_count);
      dvc_dims->putDevice(&box_dimensions[xfrm_offset], 0, 6);
      dvc_umat->putDevice(mdl_umat);
      dvc_invu->putDevice(mdl_invu);
    }
    break;
#endif
  }
#ifdef STORMM_USE_HPC
  // The upload() member function of the underlying Hybrid class will handle the format cases
  result.upload();
#endif
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace Pdb::exportPhaseSpace(const int model_index, const HybridFormat layout) const {
  UnitCellType unit_cell;
  const int xfrm_offset = model_index * roundUp(6, warp_size_int);
  const int mdl_offset = padded_atom_count * model_index;
  PhaseSpace result(atom_count, determineUnitCellType<double>(&box_dimensions.data()[xfrm_offset]),
                    layout);
  std::vector<double> mdl_umat, mdl_invu;
  computeTransforms(model_index, &mdl_umat, &mdl_invu);
  switch (layout) {
  case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
#endif
    {
      PhaseSpaceWriter rw = result.data();
      for (int i = 0; i < atom_count; i++) {
        rw.xcrd[i] = x_coordinates[mdl_offset + i];
        rw.ycrd[i] = y_coordinates[mdl_offset + i];
        rw.zcrd[i] = z_coordinates[mdl_offset + i];
        rw.xalt[i] = x_coordinates[mdl_offset + i];
        rw.yalt[i] = y_coordinates[mdl_offset + i];
        rw.zalt[i] = z_coordinates[mdl_offset + i];
      }
      for (int i = 0; i < 9; i++) {
        rw.umat[i] = mdl_umat[i];
        rw.invu[i] = mdl_invu[i];
        rw.umat_alt[i] = mdl_umat[i];
        rw.invu_alt[i] = mdl_invu[i];
      }
      for (int i = 0; i < 6; i++) {
        rw.boxdim[i] = box_dimensions[xfrm_offset + i];
        rw.boxdim_alt[i] = box_dimensions[xfrm_offset + i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::DEVICE_ONLY:
    {
      Hybrid<double>* dvc_x = result.getCoordinateHandle(CartesianDimension::X,
                                                         TrajectoryKind::POSITIONS,
                                                         CoordinateCycle::WHITE);
      Hybrid<double>* dvc_y = result.getCoordinateHandle(CartesianDimension::Y,
                                                         TrajectoryKind::POSITIONS,
                                                         CoordinateCycle::WHITE);
      Hybrid<double>* dvc_z = result.getCoordinateHandle(CartesianDimension::Z,
                                                         TrajectoryKind::POSITIONS,
                                                         CoordinateCycle::WHITE);
      Hybrid<double>* dvc_xalt = result.getCoordinateHandle(CartesianDimension::X,
                                                            TrajectoryKind::POSITIONS,
                                                            CoordinateCycle::BLACK);
      Hybrid<double>* dvc_yalt = result.getCoordinateHandle(CartesianDimension::Y,
                                                            TrajectoryKind::POSITIONS,
                                                            CoordinateCycle::BLACK);
      Hybrid<double>* dvc_zalt = result.getCoordinateHandle(CartesianDimension::Z,
                                                            TrajectoryKind::POSITIONS,
                                                            CoordinateCycle::BLACK);
      Hybrid<double>* dvc_dims = result.getBoxDimensionsHandle();
      Hybrid<double>* dvc_umat = result.getBoxTransformHandle();
      Hybrid<double>* dvc_invu = result.getInverseTransformHandle();
      Hybrid<double>* dvc_dims_alt = result.getBoxDimensionsHandle();
      Hybrid<double>* dvc_umat_alt = result.getBoxTransformHandle();
      Hybrid<double>* dvc_invu_alt = result.getInverseTransformHandle();
      dvc_x->putDevice(&x_coordinates[mdl_offset], 0, atom_count);
      dvc_y->putDevice(&y_coordinates[mdl_offset], 0, atom_count);
      dvc_z->putDevice(&z_coordinates[mdl_offset], 0, atom_count);
      dvc_dims->putDevice(&box_dimensions[xfrm_offset], 0, 6);
      dvc_umat->putDevice(mdl_umat);
      dvc_invu->putDevice(mdl_invu);
    }
    break;
#endif
  }
#ifdef STORMM_USE_HPC
  // The upload() member function of the underlying Hybrid class will handle the format cases
  result.upload();
#endif
  return result;
}

//-------------------------------------------------------------------------------------------------
void Pdb::writeToFile(const std::string &file_name, const PrintSituation expectation,
                      const std::vector<int> &model_indices) const {
  std::ofstream foutp = openOutputFile(file_name, expectation, "write the contents of a Pdb class "
                                       "object to a file");
  if (model_indices.size() > 0) {
    for (size_t i = 0; i < model_indices.size(); i++) {
      writeModel(&foutp, model_indices[i]);
    }
  }
  else {
    for (size_t i = 0; i < model_count; i++) {
      writeModel(&foutp, i);
    }
  }

  // Print the closing records
  char buffer[128];
  snprintf(buffer, 128, "END\n");
  foutp.write(buffer, strlen(buffer));
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadTopologicalData(const AtomGraph *ag) {
  if (ag->getAtomCount() != atom_count) {
    rtErr("The specified topology contains " + std::to_string(ag->getAtomCount()) + " atoms, but "
          "the object expects " + std::to_string(atom_count) + ".", "Pdb", "loadTopologicalData");
  }
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const std::vector<int> res_idx = ag->getResidueIndex();
  const int last_solute_atom = ag->getLastSoluteAtom();
  for (int i = 0; i < atom_count; i++) {
    atom_numbers[i] = cdk.atom_numbers[i];
    atom_names[i] = cdk.atom_names[i];
    residue_numbers[i] = cdk.res_numbers[i];
    residue_names[i] = cdk.res_names[res_idx[i]];
    atom_symbols[i] = uppercase(elemental_symbols[cdk.z_numbers[i]]);
    chain_names[i] = (i <= last_solute_atom) ? 'A' : ' ';
    atom_formal_charges[i] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadTopologicalData(const AtomGraph &ag) {
  loadTopologicalData(ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const llint* xcrd_in, const llint* ycrd_in, const llint* zcrd_in,
                          const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                          const double gpos_scale, const int model_index) {
  validateModelIndex(model_index, "loadCoordinates");
  const double inv_gpos_scale = 1.0 / gpos_scale;
  for (int i = 0; i < atom_count; i++) {
    const int ip = (model_index * padded_atom_count) + i;
    x_coordinates[ip] = hostInt95ToDouble(xcrd_in[i], xcrd_ovrf[i]) * inv_gpos_scale;
    y_coordinates[ip] = hostInt95ToDouble(ycrd_in[i], ycrd_ovrf[i]) * inv_gpos_scale;
    z_coordinates[ip] = hostInt95ToDouble(zcrd_in[i], zcrd_ovrf[i]) * inv_gpos_scale;
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const CoordinateFrame *cf, const int model_index) {
  const CoordinateFrameReader cfr = cf->data();
  loadCoordinates<double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const CoordinateFrame &cf, const int model_index) {
  const CoordinateFrameReader cfr = cf.data();
  loadCoordinates<double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const PhaseSpace *ps, const int model_index) {
  const PhaseSpaceReader psr = ps->data();
  loadCoordinates<double>(psr.xcrd, psr.ycrd, psr.zcrd, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const PhaseSpace &ps, const int model_index) {
  const PhaseSpaceReader psr = ps.data();
  loadCoordinates<double>(psr.xcrd, psr.ycrd, psr.zcrd, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                          const int model_index) {
  validateModelIndex(model_index, "loadCoordinates");
  const PsSynthesisReader poly_psr = poly_ps->data();
  if (poly_psr.atom_counts[system_index] != atom_count) {
    rtErr("Unable to load a system with " + std::to_string(poly_psr.atom_counts[system_index]) +
          " atoms into an object allocated for " + std::to_string(atom_count) + ".", "Pdb",
          "loadCoordinates");
  }
  const size_t offset = poly_psr.atom_starts[system_index];
  loadCoordinates(&poly_psr.xcrd[offset], &poly_psr.ycrd[offset], &poly_psr.zcrd[offset],
                  &poly_psr.xcrd_ovrf[offset], &poly_psr.ycrd_ovrf[offset],
                  &poly_psr.zcrd_ovrf[offset], poly_psr.gpos_scale, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                          const int model_index) {
  loadCoordinates(poly_ps.getSelfPointer(), system_index, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const Condensate *cdns, const int system_index, const int model_index) {
  validateModelIndex(model_index, "loadCoordinates");
  const CondensateReader cdnsr = cdns->data();
  if (system_index >= cdnsr.system_count) {
    rtErr("System " + std::to_string(system_index) + " is invalid for a Condensate with " +
          std::to_string(cdnsr.system_count) + " systems.", "Pdb", "loadCoordinates");
  }
  if (cdnsr.atom_counts[system_index] != atom_count) {
    rtErr("Unable to load a system with " + std::to_string(cdnsr.atom_counts[system_index]) +
          " atoms into an object allocated for " + std::to_string(atom_count) + ".", "Pdb",
          "loadCoordinates");
  }
  const size_t offset = cdnsr.atom_starts[system_index];
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    loadCoordinates<double>(&cdnsr.xcrd[offset], &cdnsr.ycrd[offset], &cdnsr.zcrd[offset],
                            model_index);
    break;
  case PrecisionModel::SINGLE:
    loadCoordinates<float>(&cdnsr.xcrd_sp[offset], &cdnsr.ycrd_sp[offset], &cdnsr.zcrd_sp[offset],
                           model_index);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::loadCoordinates(const Condensate &cdns, const int system_index, const int model_index) {
  loadCoordinates(cdns.getSelfPointer(), system_index, model_index);
}

//-------------------------------------------------------------------------------------------------
void Pdb::setFormalChargePrinting(bool print_formal_charges_in) {
  print_formal_charges = print_formal_charges_in;
}
  
//-------------------------------------------------------------------------------------------------
void Pdb::allocate() {
  atom_classes.resize(atom_count);
  atom_numbers.resize(atom_count);
  atom_names.resize(atom_count);
  residue_names.resize(atom_count);
  chain_names.resize(atom_count);
  residue_numbers.resize(atom_count);
  x_coordinates.resize(padded_atom_count * model_count);
  y_coordinates.resize(padded_atom_count * model_count);
  z_coordinates.resize(padded_atom_count * model_count);
  const int xfrm_length = roundUp(6, warp_size_int) * model_count;
  box_dimensions.resize(xfrm_length);
  for (int i = 0; i < xfrm_length; i++) {
    box_dimensions[i] = 0.0;
  }
  occupancies.resize(padded_atom_count * model_count);
  b_factors.resize(padded_atom_count * model_count);
  atom_symbols.resize(atom_count);
  atom_formal_charges.resize(atom_count);
  anisotropy_xx.resize(padded_atom_count * model_count);
  anisotropy_xy.resize(padded_atom_count * model_count);
  anisotropy_xz.resize(padded_atom_count * model_count);
  anisotropy_yy.resize(padded_atom_count * model_count);
  anisotropy_yz.resize(padded_atom_count * model_count);
  anisotropy_zz.resize(padded_atom_count * model_count);
}
  
//-------------------------------------------------------------------------------------------------
void Pdb::computeTransforms(const int model_index, std::vector<double> *umat,
                            std::vector<double> *invu) const {
  umat->resize(9);
  invu->resize(9);
  double* umat_ptr = umat->data();
  double* invu_ptr = invu->data();
  const int xfrm_offset = model_index * roundUp(6, warp_size_int);
  const int mdl_offset = padded_atom_count * model_index;
  if (fabs(box_dimensions[mdl_offset + 3]) < 1.0e-6 ||
      fabs(box_dimensions[mdl_offset + 4]) < 1.0e-6 || 
      fabs(box_dimensions[mdl_offset + 5]) < 1.0e-6) {
    for (int i = 0; i < 9; i++) {
      umat_ptr[i] = ((i & 0x3) == 0);
      invu_ptr[i] = umat_ptr[i];
    }
  }
  else {
    computeBoxTransform(&box_dimensions[xfrm_offset], umat->data(), invu->data());
  }  
}

//-------------------------------------------------------------------------------------------------
int Pdb::countAtomRecords(const TextFile &tf, const int line_start, const int line_end) const {
  const TextFileReader tfr = tf.data();
  const int actual_line_end = (line_end < 0 ||
                               line_end >= tfr.line_count) ? tfr.line_count : line_end;
  int result = 0;
  for (int i = line_start; i < actual_line_end; i++) {
    const size_t l_zero = tfr.line_limits[i];
    result += (tfr.line_lengths[i] >= 54 &&
               (tfr.text[l_zero] == 'A' && tfr.text[l_zero + 1] == 'T' &&
                tfr.text[l_zero + 2] == 'O' && tfr.text[l_zero + 3] == 'M') ||
               (tfr.text[l_zero] == 'H' && tfr.text[l_zero + 1] == 'E' &&
                tfr.text[l_zero + 2] == 'T' && tfr.text[l_zero + 3] == 'A' &&
                tfr.text[l_zero + 4] == 'T' && tfr.text[l_zero + 5] == 'M'));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void Pdb::readAtomRecords(const TextFile &tf, const int line_start, const int line_end,
                          const int atom_start, const std::vector<char> &alternate_prefs) {
  const TextFileReader tfr = tf.data();
  int atmcon = atom_start;
  for (int i = line_start; i < line_end; i++) {
    const size_t lz = tfr.line_limits[i];

    // Check for TER cards if this is the first model
    if (atom_start == 0 && tfr.line_lengths[i] >= 26 &&
        tfr.text[lz] == 'T' && tfr.text[lz + 1] == 'E' && tfr.text[lz + 2] == 'R') {
      ter_card_locations.push_back(atmcon);
    }

    // 
    if (tfr.line_lengths[i] >= 54 &&
        (tfr.text[lz    ] == 'A' && tfr.text[lz + 1] == 'T' && tfr.text[lz + 2] == 'O' &&
         tfr.text[lz + 3] == 'M') ||
        (tfr.text[lz    ] == 'H' && tfr.text[lz + 1] == 'E' && tfr.text[lz + 2] == 'T' &&
         tfr.text[lz + 3] == 'A' && tfr.text[lz + 4] == 'T' && tfr.text[lz + 5] == 'M')) {

      // Check that the atom's alternate location indicator is consistent with specifications.
      bool alt_loc_match = false;
      const size_t nalt = alternate_prefs.size();
      const char alt_loc_id = tfr.text[lz + 16];
      for (size_t j = 0; j < nalt; j++) {
        alt_loc_match = (alt_loc_match || alt_loc_id == alternate_prefs[j]);
      }
      if (alt_loc_match == false) {
        continue;
      }

      // Check that the coordinates are properly formatted.
      if (tfr.text[lz + 34] != '.' || tfr.text[lz + 42] != '.' || tfr.text[lz + 50] != '.') {
        rtErr("The file " + tf.getFileName() + " was submitted for coordinate extraction in PDB "
              "format, but does not appear to conform to the format on line " +
              std::to_string(i + 1) + ".", "extractPdbCoordinates");
      }

      // Determine the atom classification.
      if (atom_start == 0) {
        if (tfr.text[lz] == 'A') {
          atom_classes[atmcon] = PdbAtomKind::ATOM;
        }
        else {
          atom_classes[atmcon] = PdbAtomKind::HETATM;
        }
      }

      // Transcribe the atom and residue numbers
      char4 tmp_atom_name, tmp_resi_name;
      char *endptr;
      char buffer[9];
      buffer[5] = '\0';
      for (int j = 0; j < 5; j++) {
        buffer[j] = tfr.text[lz + 6 + j];
      }
      if (atom_start == 0) {
        atom_numbers[atmcon] = strtol(buffer, &endptr, 10);
      }
      buffer[4] = '\0';
      for (int j = 0; j < 4; j++) {
        buffer[j] = tfr.text[lz + 22 + j];
      }
      if (atom_start == 0) {
        residue_numbers[atmcon] = strtol(buffer, &endptr, 10);
      }
      
      // Transcribe the atom, residue, and chain names
      tmp_atom_name.x = tfr.text[lz + 12];
      tmp_atom_name.y = tfr.text[lz + 13];
      tmp_atom_name.z = tfr.text[lz + 14];
      tmp_atom_name.w = tfr.text[lz + 15];
      int iter = 0;
      while (tmp_atom_name.x == ' ' && iter < 4) {
        tmp_atom_name.x = tmp_atom_name.y;
        tmp_atom_name.y = tmp_atom_name.z;
        tmp_atom_name.z = tmp_atom_name.w;
        tmp_atom_name.w = ' ';
        iter++;
      }
      tmp_resi_name.x = tfr.text[lz + 17];
      tmp_resi_name.y = tfr.text[lz + 18];
      tmp_resi_name.z = tfr.text[lz + 19];
      tmp_resi_name.w = ' ';
      if (atom_start == 0) {
        atom_names[atmcon] = tmp_atom_name;
        residue_names[atmcon] = tmp_resi_name;
        chain_names[atmcon] = tfr.text[lz + 21];
      }

      // Transcribe numeric data for the atom: coordinates, occupancy, and the (isotropic)
      // temperature factor
      double xyz[3];
      double ob[2] = { 1.0, 0.0 };
      buffer[8]	= '\0';
      for (int j = 0; j	< 3; j++) {
        const size_t kmin = 30 + (8 * j);
        const size_t kmax = kmin + 8;
        for (size_t k = kmin; k < kmax; k++) {
          buffer[k - kmin] = tfr.text[lz + k];
        }
        xyz[j] = strtod(buffer, &endptr);
      }
      if (tfr.line_lengths[i] >= 66) {
        buffer[6] = '\0';
        for (int j = 0; j < 2; j++) {
          const size_t kmin = 54 + (6 * j);
          const size_t kmax = kmin + 6;
          for (size_t k = kmin; k < kmax; k++) {
            buffer[k - kmin] = tfr.text[lz + k];
          }
          ob[j] = strtod(buffer, &endptr);
        }
      }
      x_coordinates[atmcon] = xyz[0];
      y_coordinates[atmcon] = xyz[1];
      z_coordinates[atmcon] = xyz[2];
      occupancies[atmcon] = ob[0];
      b_factors[atmcon] = ob[1];

      // Transcribe the atom symbol and charge
      if (atom_start == 0) {
        if (tfr.line_lengths[i] >= 78) {
          char2 tmp_atom_symbol = { ' ', ' ' };
          tmp_atom_symbol.x = tfr.text[lz + 76];
          tmp_atom_symbol.y = tfr.text[lz + 77];
          atom_symbols[atmcon] = tmp_atom_symbol;
        }
        else {
          atom_symbols[atmcon] = { ' ', ' ' };
        }
        if (tfr.line_lengths[i] >= 80) {
          buffer[0] = tfr.text[lz + 78];
          buffer[1] = tfr.text[lz + 79];
          buffer[2] = '\0';
          atom_formal_charges[atmcon] = strtol(buffer, &endptr, 10);
        }
        else {
          atom_formal_charges[atmcon] = 0;
        }
      }
      
      // Look into anisotropic temperature factor records
      if (i + 1 < tfr.line_count && tfr.line_lengths[i + 1] >= 70) {
        const size_t lz_aniso = tfr.line_limits[i + 1];
        if (tfr.text[lz_aniso    ] == 'A' && tfr.text[lz_aniso + 1] == 'N' &&
            tfr.text[lz_aniso + 2] == 'I' && tfr.text[lz_aniso + 3] == 'S' &&
            tfr.text[lz_aniso + 4] == 'O' && tfr.text[lz_aniso + 5] == 'U') {

          // Check that the atom data matches.
          bool match = true;
          for (size_t j = 6; j < 27; j++) {
            match = (match && tfr.text[lz + j] == tfr.text[lz_aniso + j]);
          }
          if (match) {
            int anisou[6];
            buffer[7] = '\0';
            for (int j = 0; j < 6; j++) {
              const size_t kmin = 28 + (7 * j);
              const size_t kmax = kmin + 7;
              for (size_t k = kmin; k < kmax; k++) {
                buffer[k - kmin] = tfr.text[lz + k];
              }
              anisou[j] = strtol(buffer, &endptr, 10);
            }
            anisotropy_xx[atmcon] = anisou[0];
            anisotropy_yy[atmcon] = anisou[1];
            anisotropy_zz[atmcon] = anisou[2];
            anisotropy_xy[atmcon] = anisou[3];
            anisotropy_xz[atmcon] = anisou[4];
            anisotropy_yz[atmcon] = anisou[5];
          }
          i++;
        }
      }
      atmcon++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::seekUnitCellDimensions(const TextFile &tf, const int line_start, const int line_end,
                                 const int box_start) {
  const TextFileReader tfr = tf.data();
  for (int i = line_start; i < line_end; i++) {
    const size_t lz = tfr.line_limits[i];
    if (tfr.line_limits[i + 1] - lz >= 54 &&
        strcmpCased(&tfr.text[lz], "CRYST1", CaseSensitivity::YES)) {
      char len_buffer[10], ang_buffer[8];
      len_buffer[9] = '\0';
      ang_buffer[7] = '\0';
      char *endptr;
      for (int j = 0; j < 3; j++) {
        const int lzp_len = lz + 6 + (j * 9);
        for (int k = 0; k < 9; k++) {
          len_buffer[k] = tfr.text[lzp_len + k];
        }
        box_dimensions[box_start + j] = strtod(len_buffer, &endptr) * symbols::pi / 180.0;
        const int lzp_ang = lz + 33 + (j * 7);
        for (int k = 0; k < 7; k++) {
          ang_buffer[k] = tfr.text[lzp_ang + k];
        }
        box_dimensions[box_start + 3 + j] = strtod(ang_buffer, &endptr);
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::validateModelIndex(const int model_index, const char* caller) const {
  if (model_index >= model_count) {
    rtErr("Model index " + std::to_string(model_index) + " is invalid for a Pdb object with " +
          std::to_string(model_count) + " models.", "Pdb", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void Pdb::writeModel(std::ofstream *foutp, const int model_index) const {
  validateModelIndex(model_index);
  char buffer[128];
  const int xfrm_offset = model_index * roundUp(6, warp_size_int);
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    snprintf(buffer, 128, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
             box_dimensions[xfrm_offset], box_dimensions[xfrm_offset + 1],
             box_dimensions[xfrm_offset + 2], box_dimensions[xfrm_offset + 3],
             box_dimensions[xfrm_offset + 4], box_dimensions[xfrm_offset + 5]);
    foutp->write(buffer, strlen(buffer));
    break;
  }
  const int atom_llim = padded_atom_count * model_index;
  const int atom_hlim = atom_llim + atom_count;
  std::vector<bool> has_ter(atom_count, false);
  for (int i = 0; i < ter_card_count; i++) {
    has_ter[ter_card_locations[i]] = true;
  }
  for (int i = atom_llim; i < atom_hlim; i++) {
    const size_t base_i = i - atom_llim;
    switch (atom_classes[base_i]) {
    case PdbAtomKind::ATOM:
      snprintf(buffer, 128, "ATOM  ");
      break;
    case PdbAtomKind::HETATM:
      snprintf(buffer, 128, "HETATM");
      break;
    }

    // Format the atom name for PDB conventions
    char4 tmp_atom_name = atom_names[base_i];
    if (tmp_atom_name.w == ' ') {
      tmp_atom_name.w = tmp_atom_name.z;
      tmp_atom_name.z = tmp_atom_name.y;
      tmp_atom_name.y = tmp_atom_name.x;
      tmp_atom_name.x = ' ';
    }
    snprintf(&buffer[6], 122, "%5d %c%c%c%c %c%c%c%c%c%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf      "
             "     %c%c", atom_numbers[base_i], tmp_atom_name.x, tmp_atom_name.y, tmp_atom_name.z,
             tmp_atom_name.w, residue_names[base_i].x, residue_names[base_i].y,
             residue_names[base_i].z, residue_names[base_i].w, chain_names[base_i],
             residue_numbers[base_i], x_coordinates[i], y_coordinates[i], z_coordinates[i],
             occupancies[i], b_factors[i], atom_symbols[base_i].x, atom_symbols[base_i].y);
    const int slen = strlen(buffer);
    if (print_formal_charges) {
      snprintf(&buffer[slen], 128 - slen, "%2d\n", atom_formal_charges[base_i]);
    }
    else {
      snprintf(&buffer[slen], 128 - slen, "\n");
    }
    foutp->write(buffer, strlen(buffer));
    if (has_ter[base_i]) {
      snprintf(buffer, 128, "TER   %5d      %c%c%c%c%c%4d \n", atom_numbers[base_i],
               residue_names[base_i].x, residue_names[base_i].y, residue_names[base_i].z,
               residue_names[base_i].w, chain_names[base_i], residue_numbers[base_i]);
      foutp->write(buffer, strlen(buffer));
    }
  }
  snprintf(buffer, 128, "ENDMDL\n");
  foutp->write(buffer, strlen(buffer));
}

} // namespace structure
} // namespace stormm
