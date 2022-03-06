#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/omni_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Structure/virtual_site_handling.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/unit_test.h"

using omni::constants::small;
using omni::constants::tiny;
using omni::data_types::double3;
using omni::diskutil::DrivePathType;
using omni::diskutil::getDrivePathType;
using omni::diskutil::osSeparator;
using omni::errors::rtWarn;
using omni::math::computeBoxTransform;
using omni::random::Xoshiro256ppGenerator;
using omni::math::angleBetweenVectors;
using omni::math::crossProduct;
using omni::math::dot;
using omni::math::magnitude;
using omni::math::maxAbsValue;
using omni::math::maxValue;
using omni::math::matrixVectorMultiply;
using omni::math::minValue;
using omni::math::normalize;
using omni::math::perpendicularComponent;
using omni::math::pointPlaneDistance;
using omni::math::sum;
using omni::symbols::pi;
using omni::topology::AtomGraph;
using omni::topology::UnitCellType;
using omni::topology::VirtualSiteKind;
using omni::trajectory::CoordinateFileKind;
using omni::trajectory::CoordinateFrame;
using omni::trajectory::CoordinateFrameWriter;
using namespace omni::structure;
using namespace omni::testing;

//-------------------------------------------------------------------------------------------------
// Scramble the positions of virtual sites in a system.  Scatter frame atoms between box images
// but do not alter their re-imaged relative positions.
//
// Arguments:
//   ps:   Coordinates and box dimensions of the system
//   ag:   System topology (identifies virtual sites)
//   xsr:  Random number generator
//-------------------------------------------------------------------------------------------------
void scrambleSystemCoordinates(PhaseSpace *ps, const AtomGraph &ag, Xoshiro256ppGenerator *xsr) {
  PhaseSpaceWriter psw = ps->data();

  // Scramble and recover the virtual site locations
  for (int i = 0; i < ag.getAtomCount(); i++) {
    if (ag.getAtomicNumber(i) == 0) {
      psw.xcrd[i] += xsr->gaussianRandomNumber();
      psw.ycrd[i] += xsr->gaussianRandomNumber();
      psw.zcrd[i] += xsr->gaussianRandomNumber();
    }
    else {
      const double indx = floor(8.0 * (xsr->uniformRandomNumber() - 0.5));
      const double indy = floor(8.0 * (xsr->uniformRandomNumber() - 0.5));
      const double indz = floor(8.0 * (xsr->uniformRandomNumber() - 0.5));
      psw.xcrd[i] += (psw.invu[0] * indx) + (psw.invu[3] * indy) + (psw.invu[6] * indz);
      psw.ycrd[i] += (psw.invu[1] * indx) + (psw.invu[4] * indy) + (psw.invu[7] * indz);
      psw.zcrd[i] += (psw.invu[2] * indx) + (psw.invu[5] * indy) + (psw.invu[8] * indz);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Center a system on its first atom and re-image it according to the minimum image convention.
//
// Arguments:
//   ps:   Coordinates and box dimensions of the system
//-------------------------------------------------------------------------------------------------
void centerAndReimageSystem(PhaseSpace *ps) {
  PhaseSpaceWriter psw = ps->data();
  const double dx = psw.xcrd[0];
  const double dy = psw.ycrd[0];
  const double dz = psw.zcrd[0];
  for (int i = 0; i < psw.natom; i++) {
    psw.xcrd[i] -= dx;
    psw.ycrd[i] -= dy;
    psw.zcrd[i] -= dz;
  }
  imageCoordinates(ps, ImagingMethod::MINIMUM_IMAGE);
}

//-------------------------------------------------------------------------------------------------
// Check the virtual site placement of a system using code independent of the virtual site
// positioning functions.  Each virtual site gets up to three double-precision results, testing
// its characteristics against expected results.  Tests for the structure will be run internally.
//
// Arguments:
//   ps:        Coordinates of the system
//   ag:        System topology
//   answers:   Packed array of expected results
//   do_tests:  Indicator of whether it is worth pursuing the tests
//-------------------------------------------------------------------------------------------------
void checkVirtualSiteMetrics(const PhaseSpace &ps, const AtomGraph &ag,
                             const std::vector<double3> &answers, const TestPriority do_tests) {
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double3> result(vsk.nsite, { 0.0, 0.0, 0.0 });
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom = vsk.frame1_idx[i];
    const int frame2_atom = vsk.frame2_idx[i];
    const std::vector<double> p_vs = { psr.xcrd[vsite_atom] - psr.xcrd[parent_atom],
                                       psr.ycrd[vsite_atom] - psr.ycrd[parent_atom],
                                       psr.zcrd[vsite_atom] - psr.zcrd[parent_atom] };
    const std::vector<double> p_f2 = { psr.xcrd[frame2_atom] - psr.xcrd[parent_atom],
                                       psr.ycrd[frame2_atom] - psr.ycrd[parent_atom],
                                       psr.zcrd[frame2_atom] - psr.zcrd[parent_atom] };
    const double p_vs_distance = sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                      (p_vs[2] * p_vs[2]));
    const double p_f2_distance = sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                      (p_f2[2] * p_f2[2]));
    int frame3_atom, frame4_atom;
    double p_f3_distance, p_f4_distance;
    std::vector<double> p_f3(3), p_f4(3);
    const VirtualSiteKind kind = static_cast<VirtualSiteKind>(vsk.vs_types[i]);

    // Get details of the the third frame atom
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
    case VirtualSiteKind::FIXED_4:
      frame3_atom = vsk.frame3_idx[i];
      p_f3 = { psr.xcrd[frame3_atom] - psr.xcrd[parent_atom],
               psr.ycrd[frame3_atom] - psr.ycrd[parent_atom],
               psr.zcrd[frame3_atom] - psr.zcrd[parent_atom] };
      break;
    }

    // Get details of the the fourth frame atom
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      break;
    case VirtualSiteKind::FIXED_4:
      frame4_atom = vsk.frame4_idx[i];
      p_f4 = { psr.xcrd[frame4_atom] - psr.xcrd[parent_atom],
               psr.ycrd[frame4_atom] - psr.ycrd[parent_atom],
               psr.zcrd[frame4_atom] - psr.zcrd[parent_atom] };
      break;
    }

    // Evaluate each frame
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
      {
        // Check the ratio of distances between the virtual site and its parent atom
        result[i].x = p_vs_distance / p_f2_distance;

        // Check co-linearity and fill in the second slot
        result[i].y = dot(p_f2, p_vs) / (magnitude(p_f2) * magnitude(p_vs));
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        // Check the distance between the virtual site and its parent atom
        result[i].x = p_vs_distance;

        // Check co-linearity and fill in the second slot
        result[i].y = dot(p_f2, p_vs) / (magnitude(p_f2) * magnitude(p_vs));
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        // Check the distance between the virtual site and its parent atom, comparing it to the
        // expected distance.
        std::vector<double> overall_displacement(3);
        for (int j = 0; j < 3; j++) {
          overall_displacement[j] = (p_f2[j] * vsk.dim1[i]) + (p_f3[j] * vsk.dim2[i]);
        }
        result[i].x = magnitude(p_vs) / magnitude(overall_displacement);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);

        // Check the angle between the virtual site and the frame atoms.
        const std::vector<double> f2_f3 = { psr.xcrd[frame3_atom] - psr.xcrd[frame2_atom],
                                            psr.ycrd[frame3_atom] - psr.ycrd[frame2_atom],
                                            psr.zcrd[frame3_atom] - psr.zcrd[frame2_atom] };
        const std::vector<double> f2_p = { psr.xcrd[parent_atom] - psr.xcrd[frame2_atom],
                                           psr.ycrd[parent_atom] - psr.ycrd[frame2_atom],
                                           psr.zcrd[parent_atom] - psr.zcrd[frame2_atom] };
        const std::vector<double> f23_t_f2p = perpendicularComponent(f2_f3, f2_p);
        const std::vector<double> pvs_t_f2p = perpendicularComponent(p_vs, f2_p);
        result[i].z = angleBetweenVectors(p_f2, p_vs) + angleBetweenVectors(f23_t_f2p, pvs_t_f2p);
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the distance between this point and the plane, comparing it to the distance
        // prescribed by the ratio of the cross product between the vectors defining the plane.
        std::vector<double> f2_x_f3(3);
        crossProduct(p_f2, p_f3, &f2_x_f3);
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs) - (magnitude(f2_x_f3) * vsk.dim1[i] *
                                                              vsk.dim2[i] * vsk.dim3[i]);

        // Compute the distance of the virtual site from the axis perpendicular to the plane
        // running through the parent atom.
        std::vector<double> vs_t = perpendicularComponent(p_vs, f2_x_f3);
        result[i].z = magnitude(vs_t);
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Compile a list of all frame types found in this system.
  int n_unique = 0;
  std::vector<int> unique_frm;
  std::vector<bool> coverage(vsk.nsite);
  for (int i = 0; i < vsk.nsite; i++) {
    if (coverage[i]) {
      continue;
    }
    for (int j = i; j < vsk.nsite; j++) {
      coverage[j] = (coverage[j] || vsk.vs_types[j] == vsk.vs_types[i]);
    }
    unique_frm.push_back(vsk.vs_types[i]);
    n_unique++;
  }
  std::string frame_type_list;
  for (int i = 0; i < n_unique; i++) {
    frame_type_list += getVirtualSiteFrameName(static_cast<VirtualSiteKind>(unique_frm[i]));
    if (i < n_unique - 1) {
      frame_type_list += ", ";
    }
  }
  
  // Compare the first result
  std::vector<double> rvec(vsk.nsite), avec(vsk.nsite);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].x;
    avec[i] = answers[i].x;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(small), "Metric 1, measuring "
        "relative and absolute distances between virtual particles and their parent atoms, fails "
        "in a system described by topology " + ag.getFileName() + ".  Virtual site types present "
        "in this system: " + frame_type_list + ".", do_tests);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].y;
    avec[i] = answers[i].y;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(small), "Metric 2, measuring dot "
        "products and expected distances to test colinearity or point-plane distances, fails in a "
        "system described by topology " + ag.getFileName() + ".  Virtual site types present in "
        "this system: " + frame_type_list + ".", do_tests);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].z;
    avec[i] = answers[i].z;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(small), "Metric 3, measuring "
        "miscellaneous quantities like the distance between the virtual site and a normal vector "
        "(Out-3) or the angle made by a virtual site with its nearby frame atoms (FAD-3), fails "
        "in a system described by topology " + ag.getFileName() + ".  Virtual site types present "
        "in this system: " + frame_type_list + ".", do_tests);        
}

//-------------------------------------------------------------------------------------------------
int main(int argc, char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);

  // Section 1
  section("Basic checks on re-imaging calculations");

  // Section 2
  section("Check distance, angle, and dihedral calculations");

  // Section 3
  section("Virtual site placement");
  
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
  check(distance(0, 1, drug_cf), RelationalOperator::EQUAL, Approx(1.3656697990).margin(small),
        "The distance between two atoms is not computed correctly.");
  CoordinateFrameWriter cfw = drug_cf.data();
  cfw.xcrd[9] +=       cfw.boxdim[0];
  cfw.ycrd[9] += 2.0 * cfw.boxdim[1];
  cfw.zcrd[9] -= 3.0 * cfw.boxdim[2];
  check(distance(9, 10, drug_cf), RelationalOperator::EQUAL, Approx(1.5113186295).margin(small),
        "The distance between two atoms is not computed correctly after one of them has been "
        "pushed into another unit cell image.");
  check(angle(13, 46, 47, drug_cf), RelationalOperator::EQUAL, Approx(1.9124059543).margin(small),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  cfw.xcrd[22] -= 4.0 * cfw.boxdim[0];
  cfw.ycrd[22] +=       cfw.boxdim[1];
  cfw.zcrd[22] -= 2.0 * cfw.boxdim[2];
  check(angle(40, 22, 50, drug_cf), RelationalOperator::EQUAL, Approx(1.8791980388).margin(small),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  check(dihedral_angle(9, 10, 11, 12, drug_cf), RelationalOperator::EQUAL,
        Approx(-3.1069347104).margin(small), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");
  check(dihedral_angle(20, 9, 10, 11, drug_cf), RelationalOperator::EQUAL,
        Approx(-1.6233704738).margin(small), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");

  // Check the placement of virtual sites, frame type by frame type.  Scramble and recover the
  // virtual site locations, then place the entire system back in the primary unit cell and
  // check the geometries.
  section(3);
  const std::string brbz_top_path = base_top_path + osc + "bromobenzene_vs.top";
  const std::string brbz_crd_path = base_crd_path + osc + "bromobenzene_vs.inpcrd";
  const std::string stro_top_path = base_top_path + osc + "stereo_L1_vs.top";
  const std::string stro_crd_path = base_crd_path + osc + "stereo_L1_vs.inpcrd";
  const std::string symm_top_path = base_top_path + osc + "symmetry_L1_vs.top";
  const std::string symm_crd_path = base_crd_path + osc + "symmetry_L1_vs.inpcrd";
  const bool vsfi_exist = (getDrivePathType(brbz_top_path) == DrivePathType::FILE &&
                           getDrivePathType(brbz_crd_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_top_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_crd_path) == DrivePathType::FILE &&
                           getDrivePathType(stro_top_path) == DrivePathType::FILE &&
                           getDrivePathType(stro_crd_path) == DrivePathType::FILE);
  const TestPriority do_vs_tests = (vsfi_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph brbz_ag = (vsfi_exist) ? AtomGraph(brbz_top_path) : AtomGraph();
  PhaseSpace brbz_ps = (vsfi_exist) ? PhaseSpace(brbz_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  const int n_brbz_vs = brbz_ag.getVirtualSiteCount();
  std::vector<int> brbz_frame_type_answer(n_brbz_vs);
  brbz_frame_type_answer[0] = static_cast<int>(VirtualSiteKind::FIXED_2);
  for (int i = 1; i < n_brbz_vs; i++) {
    brbz_frame_type_answer[i] = static_cast<int>(VirtualSiteKind::OUT_3);
  }
  std::vector<int> brbz_detected_frame_types(n_brbz_vs);
  for (int i = 0; i < n_brbz_vs; i++) {
    brbz_detected_frame_types[i] = static_cast<int>(brbz_ag.getVirtualSiteFrameType(i));
  }
  check(brbz_frame_type_answer, RelationalOperator::EQUAL, brbz_detected_frame_types,
        "The bromobenzene system, containing a mixture of FIXED_2 and OUT_3 virtual sites, did "
        "not correctly report its virtual site content.", do_vs_tests);  
  PhaseSpaceWriter brbz_psw = brbz_ps.data();
  scrambleSystemCoordinates(&brbz_ps, brbz_ag, &xsr);
  placeVirtualSites(&brbz_ps, brbz_ag);
  centerAndReimageSystem(&brbz_ps);
  const std::vector<double3> brbz_answers = { {  0.800000000, -1.000000000,  0.000000000 },
                                              {  0.434257916,  1.000000000,  0.000000000 },
                                              {  0.434257916,  1.000000000,  0.000000000 },
                                              {  0.429810021,  1.000000000,  0.000000000 },
                                              {  0.429810021,  1.000000000,  0.000000000 },
                                              {  0.441946790,  1.000000000,  0.000000000 },
                                              {  0.441946790,  1.000000000,  0.000000000 },
                                              {  0.440347769,  1.000000000,  0.000000000 },
                                              {  0.440347769,  1.000000000,  0.000000000 },
                                              {  0.438562255,  1.000000000,  0.000000000 },
                                              {  0.438562255,  1.000000000,  0.000000000 },
                                              {  0.421362406,  1.000000000,  0.000000000 },
                                              {  0.421362406,  1.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(brbz_ps, brbz_ag, brbz_answers, do_vs_tests);
  AtomGraph stro_ag = (vsfi_exist) ? AtomGraph(stro_top_path) : AtomGraph();
  PhaseSpace stro_ps = (vsfi_exist) ? PhaseSpace(stro_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  PhaseSpaceWriter stro_psw = stro_ps.data();
  scrambleSystemCoordinates(&stro_ps, stro_ag, &xsr);
  placeVirtualSites(&stro_ps, stro_ag);
  centerAndReimageSystem(&stro_ps);
  const std::vector<double3> stro_answers = { {  1.000000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  1.000000000,  0.000000000 },
                                              {  0.300000000, -1.000000000,  0.000000000 },
                                              {  0.250000000,  0.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(stro_ps, stro_ag, stro_answers, do_vs_tests);
  AtomGraph symm_ag = (vsfi_exist) ? AtomGraph(symm_top_path) : AtomGraph();
  PhaseSpace symm_ps = (vsfi_exist) ? PhaseSpace(symm_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  PhaseSpaceWriter symm_psw = symm_ps.data();
  scrambleSystemCoordinates(&symm_ps, symm_ag, &xsr);
  placeVirtualSites(&symm_ps, symm_ag);
  centerAndReimageSystem(&symm_ps);
  const std::vector<double3> symm_answers = { {  0.500000000,  0.000000000,  2.094395121 },
                                              {  0.500000000,  0.000000000,  2.094395115 },
                                              {  0.500000000,  0.000000000,  2.094395115 },
                                              {  0.500000000,  0.000000000,  2.094395115 },
                                              {  0.500000000,  0.000000000,  2.094395100 },
                                              {  0.500000000,  0.000000000,  2.094395100 },
                                              {  0.450000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  0.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(symm_ps, symm_ag, symm_answers, do_vs_tests);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
}
