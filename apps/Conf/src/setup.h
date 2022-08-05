// -*-c++-*-
#ifndef CONFORMER_SETUP_H
#define CONFORMER_SETUP_H

#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/Namelists/user_settings.h"

namespace conf_app {
namespace setup {

using stormm::namelist::UserSettings;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SystemCache;
using stormm::testing::StopWatch;

/// \brief Enumerate coarse-grained conformer sampling approaches: full if the combinatorial
///        permutations are few enough, randomized without replacement if the combinatorial
///        permutations are up to 32 times larger, and randomized with replacement if the
///        combinatorial permutations are huge.
enum class SamplingStrategy {
  FULL,     ///< Sample all permutations of conformers
  LIMITED,  ///< Randomly sample conformers without replacement
  SPARSE    ///< Randomly sample conformers with no regard to whether replacement could occur
};

/// \brief Expand the initial list of systems into a complete list of initial states for the
///        population of conformers which conformer.omni will minimize in search of the
///        lowest-energy states.  This is done by parsing each topology into its chemical
///        details and then enumerating the rotatable bonds, cis-trans isomers, and chiral centers
///        which it could sample.
///
/// \param ui  User input settings, obtained from the input deck
/// \param sc  Cache of topologies and initial structures.  The list of topologies will not be
///            expanded by this procedure, but the list of structures will undergo a radical
///            expansion.
/// \param tm  Timer to record the wall time spent on various setup procedures
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     StopWatch *tm);

} // namespace setup
} // namespace conf_app

#endif
