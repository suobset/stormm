#include <cmath>
#include "Constants/symbol_values.h"
#include "isomerization.h"
#include "local_arrangement.h"
#include "structure_enumerators.h"

namespace omni {
namespace structure {

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
  xcrd[0] = static_cast<double>(psynthw.xcrd[atom_i + system_offset]) * inv_scl;
  ycrd[0] = static_cast<double>(psynthw.ycrd[atom_i + system_offset]) * inv_scl;
  zcrd[0] = static_cast<double>(psynthw.zcrd[atom_i + system_offset]) * inv_scl;
  xcrd[1] = static_cast<double>(psynthw.xcrd[atom_j + system_offset]) * inv_scl;
  ycrd[1] = static_cast<double>(psynthw.ycrd[atom_j + system_offset]) * inv_scl;
  zcrd[1] = static_cast<double>(psynthw.zcrd[atom_j + system_offset]) * inv_scl;
  for (int i = 0; i < nmove; i++) {
    const int ixyz_orig = moving_atoms[i] + system_offset;
    xcrd[i + 2] = static_cast<double>(psynthw.xcrd[ixyz_orig]) * inv_scl;
    ycrd[i + 2] = static_cast<double>(psynthw.ycrd[ixyz_orig]) * inv_scl;
    zcrd[i + 2] = static_cast<double>(psynthw.zcrd[ixyz_orig]) * inv_scl;
  }
  rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, rotation_angle);
  const double pos_scl = psynthw.gpos_scale;
  for (int i = 0; i < nmove; i++) {
    const int ixyz_dest = moving_atoms[i] + system_offset;
    psynthw.xcrd[ixyz_dest] = llround(xcrd[i + 2] * pos_scl);
    psynthw.ycrd[ixyz_dest] = llround(ycrd[i + 2] * pos_scl);
    psynthw.zcrd[ixyz_dest] = llround(zcrd[i + 2] * pos_scl);
  }
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceSynthesis *psynth, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  rotateAboutBond(psynth->data(), system_index, atom_i, atom_j, moving_atoms, rotation_angle);
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
  flipChiralCenter<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, center_idx, chiral_centers,
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
  flipChiralCenter<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, center_idx, chiral_centers,
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
      const llint tmp_roota_x = psynthw.xcrd[root_a + system_offset];
      const llint tmp_roota_y = psynthw.ycrd[root_a + system_offset];
      const llint tmp_roota_z = psynthw.zcrd[root_a + system_offset];
      const llint tmp_rootb_x = psynthw.xcrd[root_b + system_offset];
      const llint tmp_rootb_y = psynthw.ycrd[root_b + system_offset];
      const llint tmp_rootb_z = psynthw.zcrd[root_b + system_offset];
      const llint tmp_ccen_x  = psynthw.xcrd[chiral_centers[center_idx] + system_offset];
      const llint tmp_ccen_y  = psynthw.ycrd[chiral_centers[center_idx] + system_offset];
      const llint tmp_ccen_z  = psynthw.zcrd[chiral_centers[center_idx] + system_offset];
      const double ccenx = static_cast<double>(tmp_ccen_x) * inv_scl;
      const double cceny = static_cast<double>(tmp_ccen_y) * inv_scl;
      const double ccenz = static_cast<double>(tmp_ccen_z) * inv_scl;
      double dbx = static_cast<double>(tmp_rootb_x - tmp_ccen_x) * inv_scl;
      double dby = static_cast<double>(tmp_rootb_y - tmp_ccen_y) * inv_scl;
      double dbz = static_cast<double>(tmp_rootb_z - tmp_ccen_z) * inv_scl;
      const double invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      double dax = static_cast<double>(tmp_roota_x - tmp_ccen_x) * inv_scl;
      double day = static_cast<double>(tmp_roota_y - tmp_ccen_y) * inv_scl;
      double daz = static_cast<double>(tmp_roota_z - tmp_ccen_z) * inv_scl;
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
        const llint tmpcrd_x = psynthw.xcrd[moving_atoms[i] + system_offset];
        const llint tmpcrd_y = psynthw.ycrd[moving_atoms[i] + system_offset];
        const llint tmpcrd_z = psynthw.zcrd[moving_atoms[i] + system_offset];
        xcrd[i + 2] = static_cast<double>(tmpcrd_x) * inv_scl;
        ycrd[i + 2] = static_cast<double>(tmpcrd_y) * inv_scl;
        zcrd[i + 2] = static_cast<double>(tmpcrd_z) * inv_scl;
      }
      rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, symbols::pi);
      const double pos_scl = psynthw.gpos_scale;
      for (int i = 0; i < nmove; i++) {
        const int ixyz_dest = moving_atoms[i] + system_offset;
        psynthw.xcrd[ixyz_dest] = llround(xcrd[i + 2] * pos_scl);
        psynthw.ycrd[ixyz_dest] = llround(ycrd[i + 2] * pos_scl);
        psynthw.zcrd[ixyz_dest] = llround(zcrd[i + 2] * pos_scl);
      }
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int natom = inversion_groups[center_idx].rotatable_atoms.size();
      for (int i = 0; i < natom; i++) {
        const int atom_idx = inversion_groups[center_idx].rotatable_atoms[i] + system_offset;
        psynthw.xcrd[atom_idx] = -psynthw.xcrd[atom_idx];
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
