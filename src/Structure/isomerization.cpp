#include <cmath>
#include "isomerization.h"

namespace omni {
namespace structure {
  
//-------------------------------------------------------------------------------------------------
void rotateAboutBond(double* xcrd, double* ycrd, double* zcrd, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {

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
                                     cos_ra + (dy * dy * (1.0 -cos_ra)),
                                     (dz * dy * (1.0 - cos_ra)) + (dx * sin_ra),
                                     (dx * dz * (1.0 - cos_ra)) + (dy * sin_ra),
                                     (dy * dz * (1.0 - cos_ra)) - (dx * sin_ra),
                                     cos_ra + (dz * dz * (1.0 -cos_ra)) };

  // Loop over all moving particles and rotate about the vector
  const int natom = moving_atoms.size();
  for (int i = 0; i < natom; i++) {
    const int mk = moving_atoms[i];
    const double nx = (rmat[0] * xcrd[mk]) + (rmat[3] * ycrd[mk]) + (rmat[6] * zcrd[mk]);
    const double ny = (rmat[1] * xcrd[mk]) + (rmat[4] * ycrd[mk]) + (rmat[7] * zcrd[mk]);
    const double nz = (rmat[2] * xcrd[mk]) + (rmat[5] * ycrd[mk]) + (rmat[8] * zcrd[mk]);
    xcrd[mk] = nx;
    ycrd[mk] = ny;
    zcrd[mk] = nz;
  }
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrame *cf, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cf->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrameWriter *cfw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cfw->xcrd, cfw->ycrd, cfw->zcrd, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpace *ps, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(ps->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceWriter *cfw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(psw->xcrd, psw->ycrd, psw->zcrd, atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceSynthesis *psynth, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  rotateAboutBond(psynth->data(), atom_i + system_offset, atom_j + system_offset,
                  moving_atoms_offset, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PsSynthesisWriter *psynthw, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const double rotation_angle) {
  const int system_offset = psynthw->atom_starts[system_index];
  const int system_natom = psynthw->atom_counts[system_index];
  const int nmove = moving_atoms.size();
  std::vector<int> moving_atoms_offset(nmove);
  for (int i = 0; i < nmove; i++) {
    if (moving_atoms[i] >= system_natom) {
      rtErr("System " + std::to_string(system_index + 1) + " has " + std::to_string(system_natom) +
            " atoms, but index " + moving_atoms[i] + " was requested.", "rotateAboutBond");
    }
    moving_atoms_offset[i] = moving_atoms[i] + system_offset;
  }
  if (atom_i < 0 || atom_i >= system_natom) {
    rtErr("System " + std::to_string(system_index + 1) + " has " + std::to_string(system_natom) +
          " atoms, but index " + atom_i + " was requested.", "rotateAboutBond");
  }
  rotateAboutBond(psynthw->xcrd, psynthw->ycrd, psynthw->zcrd, atom_i + offset, atom_j + offset,
                  moving_atoms, rotation_angle);
}

} // namespace structure
} // namespace omni
