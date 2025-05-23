// -*-c++-*-
#ifndef STORMM_NBEMULATOR_H
#define STORMM_NBEMULATOR_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "Namelists/nml_emulate.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Topology/atomgraph.h"
#include "energy_enumerators.h"
#include "nbemulator_util.h"

namespace stormm {
namespace energy {
  
using namelist::EmulatorControls;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SynthesisCacheMap;
using synthesis::SystemCache;

/// \brief In order to make computed binding energies agree with experiment, it may be necessary to
///        include fitted potentials between a choice of atom types.  The means of choosing the
///        atom types and the constraints on fitting are important, and incorporated as methods of
///        this class.
class NBEmulator {
public:

  /// \brief Instances of the class may be read from a formatted file (a set of predetermined
  ///        potentials) or computed based on a synthesis of coordinates, topologies, chemical
  ///        features and chosen atom types (a set of fitted potentials).  In the latter case,
  ///        methods of the class can export the fitted potentials to a formatted file.
  ///
  /// \param file_name      Name of a file containing fitted parameters
  /// \param poly_ps        A collection of systems with known coordinates, corresponding to
  ///                       entries of target_nrg
  /// \param poly_ag        Toplogies for the systems in poly_ps, used to evaluate the molecular
  ///                       mechanics energies of each system
  /// \param chemfe_v       Array of chemical features for each unique topology in the coordinate
  ///                       synthesis, presented in the same order as the unique topologies of the
  ///                       coordinate synthesis
  /// \param mm_context     A series of molecular mechanics terms that the fitted model is expected
  ///                       to operate alongside
  /// \param atom_sources   List of definitions for the atom types giving rise to pair potentials
  /// \param basis_width    The width of the support in every basis function.  In all cases the
  ///                       basis function is a cubic spline with zero derivative at the end
  ///                       points, a value of 1.0 nearer to the origin (zero separation between
  ///                       two particles) and a value of 0.0 further from the origin.  The basis
  ///                       function therefore resembles a smoothed step function with demarcations
  ///                       set by this width and values in the basis_origins array.
  /// \param basis_origins  The lower limit of the support for each basis function.  The length of
  ///                       this vector indicates the total number of basis functions in every pair
  ///                       potential between any two atom types.
  /// \param clusters       Some target energies may reflect not the overall energy of a particular
  ///                       system but the difference in energy for a combination of systems.  This
  ///                       vector of vectors lists the various combinations.  If the vector is
  ///                       left empty, each system in the synthesis will be expected to correspond
  ///                       to one value of the target_nrg_in vector.  Otherwise, one element of
  ///                       clusters must correspond to each element of target_nrg_in.
  /// \param force_zeros    A list of structures and atoms within those structures wherein the
  ///                       forces on atoms due to all contributions (fitted potentials plus
  ///                       molecular mechanics terms) are expected to be zero.
  /// \{
  NBEmulator(const std::string &file_name);

  NBEmulator(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
             const SystemCache &sysc, const SynthesisCacheMap &scmap,
             const EmulatorControls &emulcon);
  /// \}

private:

  /// The stiffness matrix will not be allocated if potentials are read from a file and the object
  /// is to be applied in subsequent energy calculations.  Otherwise, the matrix is an intricate
  /// construction, with each column assigned a particular variable.  Descriptions of each column
  /// are stored in the column_descriptors array, descriptions of each row in the row_descriptors
  /// array.
  Hybrid<double> stiffness_matrix;

  /// Descriptions of the contents of each row of the stiffness matrix.  Like the matrix itself,
  /// this array is only allocated and filled if potentials are being fitted.
  std::vector<std::string> row_descriptors;
  
  /// Descriptions of the contents of each column of the stiffness matrix.  Like the matrix itself,
  /// this array is only allocatged and filled if potentials are being fitted.
  std::vector<std::string> column_descriptors;

  /// Target energies (or forces) for solving the equation Ax = b in the least squares sense,
  /// where A is the stiffness matrix.  The potentials are fitted to minimize the difference
  /// between calculated molecular mechanics plus potential energies and these values.
  Hybrid<double> target_nrg;
};

} // namespace energy
} // namespace stormm

#endif
