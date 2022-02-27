#include "isomerization.h"

namespace omni {
namespace structure {

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(double* xcrd, double* ycrd, double* zcrd, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {

  // Define the vector of rotation

  // Loop over all moving particles and rotate about the vector
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
  rotateAboutBond(psynthw->xcrd, psynthw->ycrd, psynthw->zcrd, atom_i + offset, atom_j + offset,
                  moving_atoms, rotation_angle);
}

} // namespace structure
} // namespace omni
