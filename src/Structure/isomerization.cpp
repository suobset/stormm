#include <cmath>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "isomerization.h"
#include "local_arrangement.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using data_types::int95_t;
  
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
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter(cf->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrameWriter cfw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, center_idx, chiral_centers,
                                   chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpace *ps, const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter(ps->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceWriter psw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, center_idx, chiral_centers,
                                   chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PsSynthesisWriter psynthw, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {

  // Prepare an array of local atom indices to indicate moving atoms.
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int nmove = inversion_groups[center_idx].getMovingAtomCount();
      std::vector<int> local_atoms(nmove);
      for (int i = 0; i < nmove; i++) {
        local_atoms[i] = i + 2;
      }
      std::vector<double> xcrd(nmove + 2), ycrd(nmove + 2), zcrd(nmove + 2);
      const double inv_scl = psynthw.inv_gpos_scale;

      // Prepare double-precision representations of the relevant coordinates, including the
      // midpoint of the two root atoms for the heaviest chiral branches and the chiral atom
      // itself.
      const int root_idx = inversion_groups[center_idx].getRootAtom() + system_offset;
      const int pivt_idx = inversion_groups[center_idx].getPivotAtom() + system_offset;
      const int chir_idx = chiral_centers[center_idx] + system_offset;
      const int95_t tmp_roota_x = { psynthw.xcrd[root_idx], psynthw.xcrd_ovrf[root_idx] };
      const int95_t tmp_roota_y = { psynthw.ycrd[root_idx], psynthw.ycrd_ovrf[root_idx] };
      const int95_t tmp_roota_z = { psynthw.zcrd[root_idx], psynthw.zcrd_ovrf[root_idx] };
      const int95_t tmp_rootb_x = { psynthw.xcrd[pivt_idx], psynthw.xcrd_ovrf[pivt_idx] };
      const int95_t tmp_rootb_y = { psynthw.ycrd[pivt_idx], psynthw.ycrd_ovrf[pivt_idx] };
      const int95_t tmp_rootb_z = { psynthw.zcrd[pivt_idx], psynthw.zcrd_ovrf[pivt_idx] };
      const int95_t tmp_ccen_x  = { psynthw.xcrd[chir_idx], psynthw.xcrd_ovrf[chir_idx] };
      const int95_t tmp_ccen_y  = { psynthw.ycrd[chir_idx], psynthw.ycrd_ovrf[chir_idx] };
      const int95_t tmp_ccen_z  = { psynthw.zcrd[chir_idx], psynthw.zcrd_ovrf[chir_idx] };
      const double ccenx = int95ToDouble(tmp_ccen_x) * inv_scl;
      const double cceny = int95ToDouble(tmp_ccen_y) * inv_scl;
      const double ccenz = int95ToDouble(tmp_ccen_z) * inv_scl;
      double dbx = int95ToDouble(splitFPSum(tmp_rootb_x, -tmp_ccen_x.x, -tmp_ccen_x.y)) * inv_scl;
      double dby = int95ToDouble(splitFPSum(tmp_rootb_y, -tmp_ccen_y.x, -tmp_ccen_y.y)) * inv_scl;
      double dbz = int95ToDouble(splitFPSum(tmp_rootb_z, -tmp_ccen_z.x, -tmp_ccen_z.y)) * inv_scl;
      const double invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      double dax = int95ToDouble(splitFPSum(tmp_roota_x, -tmp_ccen_x.x, -tmp_ccen_x.y)) * inv_scl;
      double day = int95ToDouble(splitFPSum(tmp_roota_y, -tmp_ccen_y.x, -tmp_ccen_y.y)) * inv_scl;
      double daz = int95ToDouble(splitFPSum(tmp_roota_z, -tmp_ccen_z.x, -tmp_ccen_z.y)) * inv_scl;
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
      const int* moving_atoms = inversion_groups[center_idx].getMovingAtoms().data();
      for (int i = 0; i < nmove; i++) {
        const int atom_idx = moving_atoms[i] + system_offset;
        xcrd[i + 2] = int95ToDouble(psynthw.xcrd[atom_idx], psynthw.xcrd_ovrf[atom_idx]) * inv_scl;
        ycrd[i + 2] = int95ToDouble(psynthw.ycrd[atom_idx], psynthw.ycrd_ovrf[atom_idx]) * inv_scl;
        zcrd[i + 2] = int95ToDouble(psynthw.zcrd[atom_idx], psynthw.zcrd_ovrf[atom_idx]) * inv_scl;
      }
      rotateAboutBond(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms, symbols::pi);
      const double pos_scl = psynthw.gpos_scale;
      for (int i = 0; i < nmove; i++) {
        const int ixyz_dest = moving_atoms[i] + system_offset;
        const int95_t ixcrd = doubleToInt95(xcrd[i + 2] * pos_scl);
        const int95_t iycrd = doubleToInt95(ycrd[i + 2] * pos_scl);
        const int95_t izcrd = doubleToInt95(zcrd[i + 2] * pos_scl);
        psynthw.xcrd[ixyz_dest] = ixcrd.x;
        psynthw.ycrd[ixyz_dest] = iycrd.x;
        psynthw.zcrd[ixyz_dest] = izcrd.x;
        psynthw.xcrd_ovrf[ixyz_dest] = ixcrd.y;
        psynthw.ycrd_ovrf[ixyz_dest] = iycrd.y;
        psynthw.zcrd_ovrf[ixyz_dest] = izcrd.y;
      }
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int natom = inversion_groups[center_idx].getMovingAtomCount();
      const int* moving_atoms = inversion_groups[center_idx].getMovingAtoms().data();
      for (int i = 0; i < natom; i++) {
        const int atom_idx = moving_atoms[i] + system_offset;
        psynthw.xcrd[atom_idx]      = -psynthw.xcrd[atom_idx];
        psynthw.xcrd_ovrf[atom_idx] = -psynthw.xcrd_ovrf[atom_idx];
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
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter(psynth->data(), system_index, center_idx, chiral_centers, chiral_protocols,
                   inversion_groups);
}

} // namespace structure
} // namespace stormm
