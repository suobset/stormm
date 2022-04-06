#include <cmath>
#include "Math/rounding.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_abstracts.h"
#include "isomerization.h"
#include "local_arrangement.h"
#include "structure_enumerators.h"

namespace omni {
namespace structure {

using math::roundUp;
using parse::char4ToString;
using topology::NonbondedKit;

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(double* xcrd, double* ycrd, double* zcrd, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  const int natom = moving_atoms.size();
  double center_x, center_y, center_z;
  std::vector<double> x_moves(natom), y_moves(natom), z_moves(natom);

  // Center and image the coordinates
  center_x = xcrd[atom_j];
  center_y = ycrd[atom_j];
  center_z = zcrd[atom_j];
  xcrd[atom_i] -= center_x;
  ycrd[atom_i] -= center_y;
  zcrd[atom_i] -= center_z;
  xcrd[atom_j] = 0.0;
  ycrd[atom_j] = 0.0;
  zcrd[atom_j] = 0.0;
  for (int i = 0; i < natom; i++) {
    const int mk = moving_atoms[i];
    double mvx = xcrd[mk];
    double mvy = ycrd[mk];
    double mvz = zcrd[mk];
    xcrd[mk] -= center_x;
    ycrd[mk] -= center_y;
    zcrd[mk] -= center_z;
    x_moves[mk] = mvx - xcrd[mk];
    y_moves[mk] = mvy - ycrd[mk];
    z_moves[mk] = mvz - zcrd[mk];
  }

  // Define the vector of rotation, then the matrix
  const double dx = xcrd[atom_j] - xcrd[atom_i];
  const double dy = ycrd[atom_j] - ycrd[atom_i];
  const double dz = zcrd[atom_j] - zcrd[atom_i];
  const double cos_ra = cos(rotation_angle);
  const double sin_ra = sin(rotation_angle);
  const std::vector<double> rmat = { cos_ra + (dx * dx * (1.0 - cos_ra)),
                                     (dy * dx * (1.0 - cos_ra)) + (dz * sin_ra),
                                     (dz * dx * (1.0 - cos_ra)) - (dy * sin_ra),
                                     (dx * dy * (1.0 - cos_ra)) - (dz * sin_ra),
                                     cos_ra + (dy * dy * (1.0 - cos_ra)),
                                     (dz * dy * (1.0 - cos_ra)) + (dx * sin_ra),
                                     (dx * dz * (1.0 - cos_ra)) + (dy * sin_ra),
                                     (dy * dz * (1.0 - cos_ra)) - (dx * sin_ra),
                                     cos_ra + (dz * dz * (1.0 - cos_ra)) };

  // Loop over all moving particles and rotate about the vector
  for (int i = 0; i < natom; i++) {
    const int mk = moving_atoms[i];
    const double nx = (rmat[0] * xcrd[mk]) + (rmat[3] * ycrd[mk]) + (rmat[6] * zcrd[mk]);
    const double ny = (rmat[1] * xcrd[mk]) + (rmat[4] * ycrd[mk]) + (rmat[7] * zcrd[mk]);
    const double nz = (rmat[2] * xcrd[mk]) + (rmat[5] * ycrd[mk]) + (rmat[8] * zcrd[mk]);
    xcrd[mk] = nx;
    ycrd[mk] = ny;
    zcrd[mk] = nz;
  }

  // Restore the original imaging and location of the bond atoms and moving atoms
  xcrd[atom_j] = center_x;
  ycrd[atom_j] = center_y;
  zcrd[atom_j] = center_z;
  xcrd[atom_i] += center_x;
  ycrd[atom_i] += center_y;
  zcrd[atom_i] += center_z;
  for (int i = 0; i < natom; i++) {

    // The centering for each movable atom has already been folded into the imaging move
    const int mk = moving_atoms[i];
    xcrd[mk] += x_moves[mk];
    ycrd[mk] += y_moves[mk];
    zcrd[mk] += z_moves[mk];
  }
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrame *cf, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cf->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrameWriter cfw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cfw.xcrd, cfw.ycrd, cfw.zcrd, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpace *ps, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(ps->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceWriter psw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(psw.xcrd, psw.ycrd, psw.zcrd, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceSynthesis *psynth, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  rotateAboutBond(psynth->data(), system_index, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PsSynthesisWriter psynthw, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  const int system_offset = psynthw.atom_starts[system_index];
  const int nmove = moving_atoms.size();
  std::vector<int> local_atoms(nmove);
  for (int i = 0; i < nmove; i++) {
    local_atoms[i] = i + 2;
  }
  std::vector<double> xcrd(nmove + 2), ycrd(nmove + 2), zcrd(nmove + 2);
  const double inv_scl = psynthw.inv_gpos_scale;
  const longlong4 tmpcrdi = psynthw.xyz_qlj[atom_i + system_offset];
  xcrd[0] = static_cast<double>(tmpcrdi.x) * inv_scl;
  ycrd[0] = static_cast<double>(tmpcrdi.y) * inv_scl;
  zcrd[0] = static_cast<double>(tmpcrdi.z) * inv_scl;
  const longlong4 tmpcrdj = psynthw.xyz_qlj[atom_j + system_offset];
  xcrd[1] = static_cast<double>(tmpcrdj.x) * inv_scl;
  ycrd[1] = static_cast<double>(tmpcrdj.y) * inv_scl;
  zcrd[1] = static_cast<double>(tmpcrdj.z) * inv_scl;
  for (int i = 0; i < nmove; i++) {
    const longlong4 tmpcrd = psynthw.xyz_qlj[moving_atoms[i] + system_offset];
    xcrd[i + 2] = static_cast<double>(tmpcrd.x) * inv_scl;
    ycrd[i + 2] = static_cast<double>(tmpcrd.y) * inv_scl;
    zcrd[i + 2] = static_cast<double>(tmpcrd.z) * inv_scl;
  }
  const int xfrm_offset = roundUp(9, warp_size_int) * system_index;
  rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, rotation_angle);
  const double pos_scl = psynthw.gpos_scale;
  for (int i = 0; i < nmove; i++) {
    const int ixyz_dest = moving_atoms[i] + system_offset;
    psynthw.xyz_qlj[ixyz_dest].x = round(xcrd[i + 2] * pos_scl);
    psynthw.xyz_qlj[ixyz_dest].y = round(ycrd[i + 2] * pos_scl);
    psynthw.xyz_qlj[ixyz_dest].z = round(zcrd[i + 2] * pos_scl);
  }
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(double* xcrd, double* ycrd, double* zcrd, const double* umat,
                      const double* invu, const UnitCellType unit_cell, const int chiral_center,
                      const std::vector<int> &moving_atoms, const int root_a, const int root_b) {

  // Find the bisector of the root_a : chiral_center : root_b angle, then rotate the moving atoms
  // 180 degrees about this bond.
  
  
}
  
} // namespace structure
} // namespace omni
