#include <cmath>
#include "Constants/symbol_values.h"
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
  }
  
  // Define the vector of rotation, then the matrix
  double dx = xcrd[atom_j] - xcrd[atom_i];
  double dy = ycrd[atom_j] - ycrd[atom_i];
  double dz = zcrd[atom_j] - zcrd[atom_i];
  const double invdr = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
  dx *= invdr;
  dy *= invdr;
  dz *= invdr;
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
    xcrd[mk] += center_x;
    ycrd[mk] += center_y;
    zcrd[mk] += center_z;
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
  rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, rotation_angle);
  const double pos_scl = psynthw.gpos_scale;
  for (int i = 0; i < nmove; i++) {
    const int ixyz_dest = moving_atoms[i] + system_offset;
    psynthw.xyz_qlj[ixyz_dest].x = llround(xcrd[i + 2] * pos_scl);
    psynthw.xyz_qlj[ixyz_dest].y = llround(ycrd[i + 2] * pos_scl);
    psynthw.xyz_qlj[ixyz_dest].z = llround(zcrd[i + 2] * pos_scl);
  }
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceSynthesis *psynth, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  rotateAboutBond(psynth->data(), system_index, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateSeriesWriter<double> csw, const int frame_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms, 
                     const double rotation_angle) {
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  rotateAboutBond(&csw.xcrd[fidx_zu * natom_zu], &csw.ycrd[fidx_zu * natom_zu],
                  &csw.zcrd[fidx_zu * natom_zu], atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(double* xcrd, double* ycrd, double* zcrd, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int natom, const int root_a, const int root_b) {

  switch (protocol) {
  case ChiralInversionProtocol::ROTATE:
    {
      // Find the bisector of the root_a : chiral_center : root_b angle.  Shift the root_b atom to
      // lie along the line of the bisector, rotate the moving atoms 180 degrees about this "bond,"
      // then replace the root_b atom.
      const double orig_bx = xcrd[root_b];
      const double orig_by = ycrd[root_b];
      const double orig_bz = zcrd[root_b];
      const double midpoint_x = orig_bx - xcrd[root_a];
      const double midpoint_y = orig_by - ycrd[root_a];
      const double midpoint_z = orig_bz - zcrd[root_a];
      xcrd[root_b] = midpoint_x;
      ycrd[root_b] = midpoint_y;
      zcrd[root_b] = midpoint_z;
      rotateAboutBond(xcrd, ycrd, zcrd, root_b, chiral_center, moving_atoms, symbols::pi);
      xcrd[root_b] = orig_bx;
      ycrd[root_b] = orig_by;
      zcrd[root_b] = orig_bz;
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    for (int i = 0; i < natom; i++) {
      xcrd[i] = -xcrd[i];
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrame *cf, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {
  flipChiralCenter(cf->data(), chiral_center, protocol, moving_atoms, root_a, root_b);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrameWriter cfw, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {
  flipChiralCenter(cfw.xcrd, cfw.ycrd, cfw.zcrd, chiral_center, protocol, moving_atoms, cfw.natom,
                   root_a, root_b);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpace *ps, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {
  flipChiralCenter(ps->data(), chiral_center, protocol, moving_atoms, root_a, root_b);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceWriter psw, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {
  flipChiralCenter(psw.xcrd, psw.ycrd, psw.zcrd, chiral_center, protocol, moving_atoms, psw.natom,
                   root_a, root_b);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PsSynthesisWriter psynthw, const int system_index, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {

  // Prepare an array of local atom indices to indicate moving atoms.
  switch (protocol) {
  case ChiralInversionProtocol::ROTATE:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int nmove = moving_atoms.size();
      std::vector<int> local_atoms(nmove);
      for (int i = 0; i < nmove; i++) {
        local_atoms[i] = i + 2;
      }
      std::vector<double> xcrd(nmove + 2), ycrd(nmove + 2), zcrd(nmove + 2);
      const double inv_scl = psynthw.inv_gpos_scale;

      // Prepare double-precision representations of the relevant coordinates, including the
      // midpoint of the two root atoms for the heaviest chiral branches and the chiral atom
      // itself.
      const longlong4 tmp_roota = psynthw.xyz_qlj[root_a + system_offset];
      const longlong4 tmp_rootb = psynthw.xyz_qlj[root_b + system_offset];
      const longlong4 tmp_ccen  = psynthw.xyz_qlj[chiral_center + system_offset];
      xcrd[0] = static_cast<double>(tmp_roota.x + ((tmp_rootb.x - tmp_roota.x) / 2LL)) * inv_scl;
      ycrd[0] = static_cast<double>(tmp_roota.y + ((tmp_rootb.y - tmp_roota.y) / 2LL)) * inv_scl;
      zcrd[0] = static_cast<double>(tmp_roota.z + ((tmp_rootb.z - tmp_roota.z) / 2LL)) * inv_scl;
      xcrd[1] = static_cast<double>(tmp_ccen.x) * inv_scl;
      ycrd[1] = static_cast<double>(tmp_ccen.y) * inv_scl;
      zcrd[1] = static_cast<double>(tmp_ccen.z) * inv_scl;
      for (int i = 0; i < nmove; i++) {
        const longlong4 tmpcrd = psynthw.xyz_qlj[moving_atoms[i] + system_offset];
        xcrd[i + 2] = static_cast<double>(tmpcrd.x) * inv_scl;
        ycrd[i + 2] = static_cast<double>(tmpcrd.y) * inv_scl;
        zcrd[i + 2] = static_cast<double>(tmpcrd.z) * inv_scl;
      }
      rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, symbols::pi);
      const double pos_scl = psynthw.gpos_scale;
      for (int i = 0; i < nmove; i++) {
        const int ixyz_dest = moving_atoms[i] + system_offset;
        psynthw.xyz_qlj[ixyz_dest].x = llround(xcrd[i + 2] * pos_scl);
        psynthw.xyz_qlj[ixyz_dest].y = llround(ycrd[i + 2] * pos_scl);
        psynthw.xyz_qlj[ixyz_dest].z = llround(zcrd[i + 2] * pos_scl);
      }
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int natom = psynthw.atom_counts[system_index];
      for (int i = 0; i < natom; i++) {
        const longlong4 tmp_crd = psynthw.xyz_qlj[i + system_offset];
        const longlong4 tmp_rfl = { -tmp_crd.x, tmp_crd.y, tmp_crd.z, tmp_crd.w };
        psynthw.xyz_qlj[i + system_offset] = tmp_rfl;
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceSynthesis *psynth, const int system_index, const int chiral_center,
                      const ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      const int root_a, const int root_b) {
  flipChiralCenter(psynth->data(), system_index, chiral_center, protocol, moving_atoms, root_a,
                   root_b);
}

} // namespace structure
} // namespace omni
