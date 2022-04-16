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
void flipChiralCenter(double* xcrd, double* ycrd, double* zcrd, const int center_idx,
                      const int natom, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      // Unpack the appropriate inversion group
      const int root_a = inversion_groups[center_idx].root_atom;
      const int root_b = inversion_groups[center_idx].pivot_atom;
      const int ccen   = chiral_centers[center_idx];
      
      // Find the bisector of the root_a : chiral_center : root_b angle.  Shift the root_b atom to
      // lie along the line of the bisector, rotate the moving atoms 180 degrees about this "bond,"
      // then replace the root_b atom.
      const double orig_bx = xcrd[root_b];
      const double orig_by = ycrd[root_b];
      const double orig_bz = zcrd[root_b];
      double dbx = orig_bx - xcrd[ccen];
      double dby = orig_by - ycrd[ccen];
      double dbz = orig_bz - zcrd[ccen];
      const double invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
      dbx = xcrd[ccen] + (dbx * invrb);
      dby = ycrd[ccen] + (dby * invrb);
      dbz = zcrd[ccen] + (dbz * invrb);
      double dax = xcrd[root_a] - xcrd[ccen];
      double day = ycrd[root_a] - ycrd[ccen];
      double daz = zcrd[root_a] - zcrd[ccen];
      const double invra = 1.0 / sqrt((dax * dax) + (day * day) + (daz * daz));
      dax = xcrd[ccen] + (dax * invra);
      day = ycrd[ccen] + (day * invra);
      daz = zcrd[ccen] + (daz * invra);
      const double midpoint_x = 0.5 * (dbx + dax);
      const double midpoint_y = 0.5 * (dby + day);
      const double midpoint_z = 0.5 * (dbz + daz);
      xcrd[root_b] = midpoint_x;
      ycrd[root_b] = midpoint_y;
      zcrd[root_b] = midpoint_z;
      rotateAboutBond(xcrd, ycrd, zcrd, root_b, chiral_centers[center_idx],
                      inversion_groups[center_idx].rotatable_atoms, symbols::pi);
      xcrd[root_b] = orig_bx;
      ycrd[root_b] = orig_by;
      zcrd[root_b] = orig_bz;
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      // Find the molecule home of the present center and flip those atoms only.
      for (int i = 0; i < natom; i++) {
        xcrd[i] = -xcrd[i];
      }
      
      // All centers have been flipped by the reflection.  Loop over all other chiral centers and
      // perform chiral inversions in order to flip the others back, where possible.  Infinite
      // recursion is limited by the fact that at most one chiral center in a molecule can be given
      // the protocol "REFLECT," and reflection only triggers subsequent rotations.
      const int nchirals = chiral_protocols.size();
      for (int i = 0; i < nchirals; i++) {
        if (chiral_protocols[i] == ChiralInversionProtocol::ROTATE) {
          flipChiralCenter(xcrd, ycrd, zcrd, i, natom, chiral_centers, chiral_protocols,
                           inversion_groups);
        }
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrame *cf, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter(cf->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrameWriter cfw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter(cfw.xcrd, cfw.ycrd, cfw.zcrd, center_idx, cfw.natom, chiral_centers,
                   chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpace *ps, const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter(ps->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceWriter psw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter(psw.xcrd, psw.ycrd, psw.zcrd, center_idx, psw.natom, chiral_centers,
                   chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PsSynthesisWriter psynthw, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {

  // Prepare an array of local atom indices to indicate moving atoms.
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int nmove = inversion_groups[center_idx].rotatable_atoms.size();
      std::vector<int> local_atoms(nmove);
      for (int i = 0; i < nmove; i++) {
        local_atoms[i] = i + 2;
      }
      std::vector<double> xcrd(nmove + 2), ycrd(nmove + 2), zcrd(nmove + 2);
      const double inv_scl = psynthw.inv_gpos_scale;

      // Prepare double-precision representations of the relevant coordinates, including the
      // midpoint of the two root atoms for the heaviest chiral branches and the chiral atom
      // itself.
      const int root_a = inversion_groups[center_idx].root_atom;
      const int root_b = inversion_groups[center_idx].pivot_atom;
      const longlong4 tmp_roota = psynthw.xyz_qlj[root_a + system_offset];
      const longlong4 tmp_rootb = psynthw.xyz_qlj[root_b + system_offset];
      const longlong4 tmp_ccen  = psynthw.xyz_qlj[chiral_centers[center_idx] + system_offset];
      const double ccenx = static_cast<double>(tmp_ccen.x) * inv_scl;
      const double cceny = static_cast<double>(tmp_ccen.y) * inv_scl;
      const double ccenz = static_cast<double>(tmp_ccen.z) * inv_scl;
      double dbx = static_cast<double>(tmp_rootb.x - tmp_ccen.x) * inv_scl;
      double dby = static_cast<double>(tmp_rootb.y - tmp_ccen.y) * inv_scl;
      double dbz = static_cast<double>(tmp_rootb.z - tmp_ccen.z) * inv_scl;
      const double invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      double dax = static_cast<double>(tmp_roota.x - tmp_ccen.x) * inv_scl;
      double day = static_cast<double>(tmp_roota.y - tmp_ccen.y) * inv_scl;
      double daz = static_cast<double>(tmp_roota.z - tmp_ccen.z) * inv_scl;
      const double invra = 1.0 / sqrt((dax * dax) + (day * day) + (daz * daz));
      dax = ccenx + (dax * invra);
      day = cceny + (day * invra);
      daz = ccenz + (daz * invra);
      const double midpoint_x = 0.5 * (dbx + dax);
      const double midpoint_y = 0.5 * (dby + day);
      const double midpoint_z = 0.5 * (dbz + daz);
      xcrd[0] = midpoint_x;
      ycrd[0] = midpoint_y;
      zcrd[0] = midpoint_z;
      xcrd[1] = ccenx;
      ycrd[1] = cceny;
      zcrd[1] = ccenz;
      const int* moving_atoms = inversion_groups[center_idx].rotatable_atoms.data();
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

      // All centers have been flipped by the reflection.  Loop over all other chiral centers and
      // perform chiral inversions in order to flip the others back, where possible.
      const int nchiral = chiral_protocols.size();
      for (int i = 0; i < nchiral; i++) {
        if (chiral_protocols[i] == ChiralInversionProtocol::ROTATE) {
          flipChiralCenter(psynthw, system_index, i, chiral_centers, chiral_protocols,
                           inversion_groups);
        }
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceSynthesis *psynth, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<RotatorGroup> &inversion_groups) {
  flipChiralCenter(psynth->data(), system_index, center_idx, chiral_centers, chiral_protocols,
                   inversion_groups);
}

} // namespace structure
} // namespace omni
