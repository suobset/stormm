#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::tiny;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::math::computeBoxTransform;
using omni::random::Xoshiro256ppGenerator;
using omni::math::minValue;
using omni::math::maxAbsValue;
using omni::math::maxValue;
using omni::math::matrixVectorMultiply;
using omni::math::sum;
using omni::symbols::pi;
using omni::topology::AtomGraph;
using omni::topology::UnitCellType;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFrameWriter;
using namespace omni::structure;
using namespace omni::testing;

int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Basic checks on re-imaging calculations");

  // Section 2
  section("Check distance, angle, and dihedral calculations");
  
  // Get a realistic system
  const char osc = osSeparator();
  const std::string base_top_path = oe.getOmniSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getOmniSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string drug_top_path = base_top_path + osc + "drug_example.top";
  const std::string drug_crd_path = base_crd_path + osc + "drug_example.inpcrd";
  const bool files_exist = (getDrivePathType(drug_top_path) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_path) == DrivePathType::FILE);
  AtomGraph drug_ag = (files_exist) ? AtomGraph(drug_top_path) : AtomGraph();
  CoordinateFrame drug_cf = (files_exist) ? CoordinateFrame(drug_crd_path,
                                                            CoordinateFileKind::AMBER_INPCRD) :
                                            CoordinateFrame();
  if (files_exist == false) {
    rtWarn("Files for a drug molecule, in water and inside a periodic box, were not found.  Check "
           "the $OMNI_SOURCE environment variable to ensure that " + drug_top_path + " and " +
           drug_crd_path + " become valid paths.  Some tests will be skipped",
           "test_local_arrangement");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;

  // Image a single number (needed for NMR restraints with periodicity, for example)
  section(1);
  const std::vector<double> point_samples = {  0.75,  0.25,  0.35, -0.65, -1.87,  2.33 };
  const double trial_range = 0.5;
  const std::vector<double> primu_answer  = {  0.25,  0.25,  0.35,  0.35,  0.13,  0.33 };
  const std::vector<double> minim_answer  = { -0.25, -0.25, -0.15, -0.15,  0.13, -0.17 };
  std::vector<double> primu_result(point_samples.size()), minim_result(point_samples.size());
  for (size_t i = 0; i < point_samples.size(); i++) {
    primu_result[i] = imageValue(point_samples[i], trial_range, ImagingMethod::PRIMARY_UNIT_CELL);
    minim_result[i] = imageValue(point_samples[i], trial_range, ImagingMethod::MINIMUM_IMAGE);
  }
  check(primu_result, RelationalOperator::EQUAL, Approx(primu_answer).margin(tiny), "Value "
        "imaging into a [0, range) interval does not work as expected.");
  check(minim_result, RelationalOperator::EQUAL, Approx(minim_answer).margin(tiny), "Value "
        "imaging into a [-0.5 * range, 0.5 * range) interval does not work as expected.");
  
  // Create a fake rectilinear system
  const std::vector<double> rectilinear_box = { 15.0, 24.0, 18.0, 0.5 * pi, 0.5  * pi, 0.5 * pi };
  const std::vector<double> rectilinear_x_crd = {  5.1,  0.7,  7.9, 27.9 };
  const std::vector<double> rectilinear_y_crd = {  3.8, -3.4, 13.5, -9.4 };
  const std::vector<double> rectilinear_z_crd = { -4.9,  1.6, -3.1, -38.2 };
  std::vector<double> rectilinear_umat(9), rectilinear_invu(9);
  computeBoxTransform(rectilinear_box, rectilinear_umat.data(), rectilinear_invu.data());
  double x1 = rectilinear_x_crd[0];
  double y1 = rectilinear_y_crd[0];
  double z1 = rectilinear_z_crd[0];
  imageCoordinates(&x1, &y1, &z1, rectilinear_umat.data(), rectilinear_invu.data(),
                   UnitCellType::ORTHORHOMBIC, ImagingMethod::MINIMUM_IMAGE);
  check(std::vector<double>({ x1, y1, z1}), RelationalOperator::EQUAL,
        std::vector<double>({ 5.1, 3.8, -4.9 }), "Re-imaging a single point in a rectilinear box "
        "by the minimum image convention fails.");
  imageCoordinates(&x1, &y1, &z1, rectilinear_umat.data(), rectilinear_invu.data(),
                   UnitCellType::ORTHORHOMBIC, ImagingMethod::PRIMARY_UNIT_CELL);
  check(std::vector<double>({ x1, y1, z1}), RelationalOperator::EQUAL,
        std::vector<double>({ 5.1, 3.8, 13.1 }), "Re-imaging a single point in a rectilinear box "
        "into the primary unit cell (first octant) fails.");
  std::vector<double> rect_x_copy(rectilinear_x_crd);
  std::vector<double> rect_y_copy(rectilinear_y_crd);
  std::vector<double> rect_z_copy(rectilinear_z_crd);
  imageCoordinates(&rect_x_copy, &rect_y_copy, &rect_z_copy, rectilinear_umat.data(),
                   rectilinear_invu.data(), UnitCellType::ORTHORHOMBIC,
                   ImagingMethod::MINIMUM_IMAGE);
  check(rect_x_copy, RelationalOperator::EQUAL, std::vector<double>({ 5.1, 0.7, -7.1, -2.1 }),
        "Re-imaging rectilinear Cartesian X coordinates by the minimum image convention fails.");  
  check(rect_y_copy, RelationalOperator::EQUAL, std::vector<double>({ 3.8, -3.4, -10.5, -9.4 }),
        "Re-imaging rectilinear Cartesian Y coordinates by the minimum image convention fails.");  
  check(rect_z_copy, RelationalOperator::EQUAL, std::vector<double>({ -4.9,  1.6, -3.1, -2.2 }),
        "Re-imaging rectilinear Cartesian Z coordinates by the minimum image convention fails.");

  // Create a dense spread of points that will need lots of re-imaging in the rectilinear box
  const int npts = 500;
  std::vector<double> dense_x_crd(npts), dense_y_crd(npts), dense_z_crd(npts);
  Xoshiro256ppGenerator xsr(918733245);
  for (int i = 0; i < npts; i++) {
    dense_x_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
    dense_y_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
    dense_z_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
  }
  std::vector<double> dense_x_copy(dense_x_crd);
  std::vector<double> dense_y_copy(dense_y_crd);
  std::vector<double> dense_z_copy(dense_z_crd);
  imageCoordinates(&dense_x_copy, &dense_y_copy, &dense_z_copy, rectilinear_umat.data(),
                   rectilinear_invu.data(), UnitCellType::ORTHORHOMBIC,
                   ImagingMethod::MINIMUM_IMAGE);
  std::vector<double> box_x_disp(npts), box_y_disp(npts), box_z_disp(npts);
  for (int i = 0; i < npts; i++) {
    box_x_disp[i] = (dense_x_copy[i] - dense_x_crd[i]) * rectilinear_umat[0];
    box_y_disp[i] = (dense_y_copy[i] - dense_y_crd[i]) * rectilinear_umat[4];
    box_z_disp[i] = (dense_z_copy[i] - dense_z_crd[i]) * rectilinear_umat[8];
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  check(box_x_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the X dimension.");
  check(box_y_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the Y dimension.");
  check(box_z_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the Z dimension.");
  std::vector<double> min_rect_crd = { minValue(dense_x_copy), minValue(dense_y_copy),
                                       minValue(dense_z_copy) };
  std::vector<double> max_rect_crd = { maxValue(dense_x_copy), maxValue(dense_y_copy),
                                       maxValue(dense_z_copy) };
  std::vector<double> min_rect_ans = { -0.5 * rectilinear_box[0], -0.5 * rectilinear_box[1],
                                       -0.5 * rectilinear_box[2] };
  std::vector<double> max_rect_ans = { 0.5 * rectilinear_box[0], 0.5 * rectilinear_box[1],
                                       0.5 * rectilinear_box[2] };
  check(min_rect_crd, RelationalOperator::GE, Approx(min_rect_ans).margin(tiny), "Re-imaging by "
        "the minimum image convention puts some particles outside the minimum expected range.");
  check(max_rect_crd, RelationalOperator::LT, Approx(max_rect_ans).margin(tiny), "Re-imaging by "
        "the minimum image convention puts some particles outside the maximum expected range.");
  imageCoordinates(&dense_x_copy, &dense_y_copy, &dense_z_copy, rectilinear_umat.data(),
                   rectilinear_invu.data(), UnitCellType::ORTHORHOMBIC,
                   ImagingMethod::PRIMARY_UNIT_CELL);
  for (int i = 0; i < npts; i++) {
    box_x_disp[i] = (dense_x_copy[i] - dense_x_crd[i]) * rectilinear_umat[0];
    box_y_disp[i] = (dense_y_copy[i] - dense_y_crd[i]) * rectilinear_umat[4];
    box_z_disp[i] = (dense_z_copy[i] - dense_z_crd[i]) * rectilinear_umat[8];
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  check(box_x_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the X dimension.");
  check(box_y_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the Y dimension.");
  check(box_z_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the Z dimension.");
  min_rect_crd = { minValue(dense_x_copy), minValue(dense_y_copy), minValue(dense_z_copy) };
  max_rect_crd = { maxValue(dense_x_copy), maxValue(dense_y_copy), maxValue(dense_z_copy) };
  min_rect_ans = { 0.0, 0.0, 0.0 };
  max_rect_ans = { rectilinear_box[0], rectilinear_box[1], rectilinear_box[2] };
  check(min_rect_crd, RelationalOperator::GE, Approx(min_rect_ans).margin(tiny), "Re-imaging into "
        "the primary unit cell puts some particles outside the minimum expected range.");
  check(max_rect_crd, RelationalOperator::LT, Approx(max_rect_ans).margin(tiny), "Re-imaging into "
        "the primary unit cell puts some particles outside the maximum expected range.");
  
  // Create and test a fake triclinic system
  const std::vector<double> triclinic_box = { 16.0, 17.0, 15.0, 0.6 * pi, 0.53 * pi, 0.55 * pi };
  std::vector<double> triclinic_umat(9), triclinic_invu(9);
  computeBoxTransform(triclinic_box, triclinic_umat.data(), triclinic_invu.data());
  for (int i = 0; i < npts; i++) {
    dense_x_copy[i] = dense_x_crd[i];
    dense_y_copy[i] = dense_y_crd[i];
    dense_z_copy[i] = dense_z_crd[i];
  }
  imageCoordinates(&dense_x_copy, &dense_y_copy, &dense_z_copy, triclinic_umat.data(),
                   triclinic_invu.data(), UnitCellType::TRICLINIC, ImagingMethod::MINIMUM_IMAGE);
  for (int i = 0; i < npts; i++) {
    const double dx = dense_x_copy[i] - dense_x_crd[i];
    const double dy = dense_y_copy[i] - dense_y_crd[i];
    const double dz = dense_z_copy[i] - dense_z_crd[i];
    box_x_disp[i] = (triclinic_umat[0] * dx) + (triclinic_umat[3] * dy) + (triclinic_umat[6] * dz);
    box_y_disp[i] = (triclinic_umat[1] * dx) + (triclinic_umat[4] * dy) + (triclinic_umat[7] * dz);
    box_z_disp[i] = (triclinic_umat[2] * dx) + (triclinic_umat[5] * dy) + (triclinic_umat[8] * dz);
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  std::vector<double> non_integral_disp = { maxAbsValue(box_x_disp), maxAbsValue(box_y_disp),
                                            maxAbsValue(box_z_disp) };
  check(non_integral_disp, RelationalOperator::EQUAL,
        Approx(std::vector<double>(3, 0.0)).margin(tiny), "Re-imaging in a triclinic unit cell by "
        "the minimum image convention moves particles by non-integral box lengths.");
  std::vector<double> frac_tric_x(npts), frac_tric_y(npts), frac_tric_z(npts);
  std::vector<double> tmp_crd(3);
  std::vector<double> tmp_frac(3, 0.0);
  for (int i = 0; i < npts; i++) {
    tmp_crd[0] = dense_x_copy[i];
    tmp_crd[1] = dense_y_copy[i];
    tmp_crd[2] = dense_z_copy[i];
    matrixVectorMultiply(triclinic_umat.data(), tmp_crd.data(), tmp_frac.data(), 3, 3, 1.0, 1.0,
                         0.0);
    frac_tric_x[i] = tmp_frac[0];
    frac_tric_y[i] = tmp_frac[1];
    frac_tric_z[i] = tmp_frac[2];
  }
  std::vector<double> min_tric_frac = { minValue(frac_tric_x), minValue(frac_tric_y),
                                        minValue(frac_tric_z) };
  std::vector<double> max_tric_frac = { maxValue(frac_tric_x), maxValue(frac_tric_y),
                                        maxValue(frac_tric_z) };
  check(min_tric_frac, RelationalOperator::GE, std::vector<double>({ -0.5, -0.5, -0.5 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  check(max_tric_frac, RelationalOperator::LT, std::vector<double>({ 0.5, 0.5, 0.5 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  for (int i = 0; i < npts; i++) {
    dense_x_copy[i] = dense_x_crd[i];
    dense_y_copy[i] = dense_y_crd[i];
    dense_z_copy[i] = dense_z_crd[i];
  }
  imageCoordinates(&dense_x_copy, &dense_y_copy, &dense_z_copy, triclinic_umat.data(),
                   triclinic_invu.data(), UnitCellType::TRICLINIC,
                   ImagingMethod::PRIMARY_UNIT_CELL);
  for (int i = 0; i < npts; i++) {
    const double dx = dense_x_copy[i] - dense_x_crd[i];
    const double dy = dense_y_copy[i] - dense_y_crd[i];
    const double dz = dense_z_copy[i] - dense_z_crd[i];
    box_x_disp[i] = (triclinic_umat[0] * dx) + (triclinic_umat[3] * dy) + (triclinic_umat[6] * dz);
    box_y_disp[i] = (triclinic_umat[1] * dx) + (triclinic_umat[4] * dy) + (triclinic_umat[7] * dz);
    box_z_disp[i] = (triclinic_umat[2] * dx) + (triclinic_umat[5] * dy) + (triclinic_umat[8] * dz);
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  non_integral_disp = { maxAbsValue(box_x_disp), maxAbsValue(box_y_disp),
                        maxAbsValue(box_z_disp) };
  check(non_integral_disp, RelationalOperator::EQUAL,
        Approx(std::vector<double>(3, 0.0)).margin(tiny), "Re-imaging to the primary unit cell in "
        "triclinic system moves particles by non-integral box lengths.");
  for (int i = 0; i < npts; i++) {
    tmp_crd[0] = dense_x_copy[i];
    tmp_crd[1] = dense_y_copy[i];
    tmp_crd[2] = dense_z_copy[i];
    matrixVectorMultiply(triclinic_umat.data(), tmp_crd.data(), tmp_frac.data(), 3, 3, 1.0, 1.0,
                         0.0);
    frac_tric_x[i] = tmp_frac[0];
    frac_tric_y[i] = tmp_frac[1];
    frac_tric_z[i] = tmp_frac[2];
  }
  min_tric_frac[0] = minValue(frac_tric_x);
  min_tric_frac[1] = minValue(frac_tric_y);
  min_tric_frac[2] = minValue(frac_tric_z);
  max_tric_frac[0] = maxValue(frac_tric_x);
  max_tric_frac[1] = maxValue(frac_tric_y);
  max_tric_frac[2] = maxValue(frac_tric_z);
  check(min_tric_frac, RelationalOperator::GE, std::vector<double>({ 0.0, 0.0, 0.0 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  check(max_tric_frac, RelationalOperator::LT, std::vector<double>({ 1.0, 1.0, 1.0 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");  

  // Check internal coordinate computations
  section(2);
  check(distance(0, 1, drug_cf), RelationalOperator::EQUAL, Approx(1.3656697990).margin(tiny),
        "The distance between two atoms is not computed correctly.");
  CoordinateFrameWriter cfw = drug_cf.data();
  cfw.xcrd[9] +=       cfw.boxdim[0];
  cfw.ycrd[9] += 2.0 * cfw.boxdim[1];
  cfw.zcrd[9] -= 3.0 * cfw.boxdim[2];
  check(distance(9, 10, drug_cf), RelationalOperator::EQUAL, Approx(1.5113186295).margin(tiny),
        "The distance between two atoms is not computed correctly after one of them has been "
        "pushed into another unit cell image.");
  check(angle(13, 46, 47, drug_cf), RelationalOperator::EQUAL, Approx(1.9124059543).margin(tiny),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  cfw.xcrd[22] -= 4.0 * cfw.boxdim[0];
  cfw.ycrd[22] +=       cfw.boxdim[1];
  cfw.zcrd[22] -= 2.0 * cfw.boxdim[2];
  check(angle(40, 22, 50, drug_cf), RelationalOperator::EQUAL, Approx(1.8791980388).margin(tiny),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  check(dihedral_angle(9, 10, 11, 12, drug_cf), RelationalOperator::EQUAL,
        Approx(-3.1069347104).margin(tiny), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");
  check(dihedral_angle(20, 9, 10, 11, drug_cf), RelationalOperator::EQUAL,
        Approx(-1.6233704738).margin(tiny), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
