// -*-c++-*-
#ifndef CONFORMER_ANALYSIS_H
#define CONFORMER_ANALYSIS_H

#include <vector>
#include "copyright.h"
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/Namelists/nml_conformer.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Structure/local_arrangement.h"

namespace conf_app {
namespace analysis {

using stormm::energy::ScoreCard;
using stormm::energy::StaticExclusionMask;
using stormm::namelist::ConformerControls;
using stormm::namelist::FilesControls;
using stormm::namelist::default_minimize_clash_r0;
using stormm::namelist::default_minimize_clash_ratio;
using stormm::structure::ClashReport;
using stormm::structure::MdlMol;
using stormm::synthesis::Condensate;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::PsSynthesisReader;
using stormm::synthesis::SystemCache;
    
/// \brief Compute the angles made by each rotatable bond.  This will compute the angle based on
///        the highest Z-number atoms attached to either end of the bond (not counting the atoms at
///        the other end of the bond).  If there is a tie, the atom with the lower topological
///        index will be preferred.

/// \brief Filter structures to obtain the ones with the best energies separated by some minimum
///        positional RMSD threshold.
///
/// \param poly_ps  The synthesis of coordinates, holding all of the energy-minimized structures
/// \param emin     The minimum energy values and history
/// \param confcon  Conformer namelist user input
std::vector<int> filterMinimizedStructures(const PhaseSpaceSynthesis &poly_ps,
                                           const std::vector<StaticExclusionMask> &poly_ps_masks,
                                           const ScoreCard &emin,
                                           const ConformerControls &confcon);

/// \brief Print the best structures for each system, grouped by topology, in accord with user
///        input.
///
/// \param poly_ps       Synthesis of topologies from each
/// \param best_confs    List of the best configurations
/// \param emin          Energy tracking with histories from each minimization
/// \param sc            The cache of all systems read from disk
/// \param sdf_recovery  Vector of recovered MdlMol objects for each system reach from disk
/// \param fcon          Control data from a &files namelist in the user input
void printResults(const PhaseSpaceSynthesis &poly_ps, const std::vector<int> &best_confs,
                  const ScoreCard &emin, const SystemCache &sc,
                  const std::vector<MdlMol> &sdf_recovery, const FilesControls &fcon);

} // namespace analysis
} // namespace conf_app

#endif
