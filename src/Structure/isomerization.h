// -*-c++-*-
#ifndef OMNI_ISOMERIZATION_H
#define OMNI_ISOMERIZATION_H

#include "Chemistry/chemistry_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Synthesis/phasespace_synthesis.h"

namespace omni {
namespace structure {

using chemistry::ChiralInversionProtocol;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
  
/// \brief Rotate a molecule such that selected internal coordinates achieve particular values.
///        The internal coordinates must define a rotatable bond (bonds in rings are prohibited
///        from rotation under this procedure).  This operation takes place only on the CPU.  The
///        results must be uploaded to the GPU but can serve to seed coordinates for more detailed
///        manipulations.
///
/// Overloaded:
///   - Take the X, Y, and Z particle coordinates directly
///   - Take a CoordinateFrame or writeable abstract for the coordinates
///   - Take a PhaseSpace object or writeable abstract for the coordinates
///
/// \param xcrd            Cartesian X coordinates of particles
/// \param ycrd            Cartesian Y coordinates of particles
/// \param zcrd            Cartesian Z coordinates of particles
/// \param cf              Coordinates to manipulate (will be modified upon return)
/// \param cfw             CoordianteFrame writer abstract
/// \param ps              PhaseSpace object with coordinates, velocities, and forces (the
///                        coordinates will be modified upon return)
/// \param psw             PhaseSpace writer abstract
/// \param psynth          Collection of systems, with coordinates, velocities, and forces
///                        represented in fixed precision
/// \param psynthw         PhaseSpaceSynthesis writer abstract
/// \param atom_i          Root of the rotatable bond
/// \param atom_j          Second atom of the rotatable bond.
/// \param moving_atoms    Atoms branching from atom_j and distal to atom_i.  These will rotate.
/// \param rotation_angle  The angle about which to rotate atoms indexed in moving_atoms.  The
///                        angle is oriented by the right hand rule with one's hand on atom_i and
///                        thumb pointed towards atom_j.
/// \{
void rotateAboutBond(double* xcrd, double* ycrd, double* zcrd, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(CoordinateFrame *cf, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(CoordinateFrameWriter cfw, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(PhaseSpace *ps, int atom_i, int atom_j, const std::vector<int> &moving_atoms,
                     double rotation_angle);

void rotateAboutBond(PhaseSpaceWriter psw, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(PsSynthesisWriter psynthw, int system_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(PhaseSpaceSynthesis *psynth, int system_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(CoordinateSeriesWriter<double> csw, int frame_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);
/// \}

/// \brief Rotate two branches of a chiral center 180 degrees, so as to invert the center.  
void flipChiralCenter(double* xcrd, double* ycrd, double* zcrd, int chiral_center,
                      ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      int natom, int root_a, int root_b);

void flipChiralCenter(CoordinateFrame *cf, int chiral_center, ChiralInversionProtocol protocol,
                      const std::vector<int> &moving_atoms, int root_a, int root_b);

void flipChiralCenter(CoordinateFrameWriter cfw, int chiral_center,
                      ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      int root_a, int root_b);

void flipChiralCenter(PhaseSpace *ps, int chiral_center, ChiralInversionProtocol protocol,
                      const std::vector<int> &moving_atoms, int root_a, int root_b);

void flipChiralCenter(PhaseSpaceWriter psw, int chiral_center, ChiralInversionProtocol protocol,
                      const std::vector<int> &moving_atoms, int root_a, int root_b);

void flipChiralCenter(PsSynthesisWriter psynthw, int system_index, int chiral_center,
                      ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      int root_a, int root_b);

void flipChiralCenter(PhaseSpaceSynthesis *psynth, int system_index, int chiral_center,
                      ChiralInversionProtocol protocol, const std::vector<int> &moving_atoms,
                      int root_a, int root_b);

} // namespace structure
} // namespace omni

#endif
