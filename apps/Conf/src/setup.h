// -*-c++-*-
#ifndef CONFORMER_SETUP_H
#define CONFORMER_SETUP_H

#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Random/random.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/Namelists/user_settings.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::ChiralInversionProtocol;
using stormm::chemistry::IsomerPlan;
using stormm::namelist::UserSettings;
using stormm::random::Xoshiro256ppGenerator;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SystemCache;
using stormm::testing::StopWatch;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;

/// \brief Enumerate coarse-grained conformer sampling approaches: full if the combinatorial
///        permutations are few enough, randomized without replacement if the combinatorial
///        permutations are up to 32 times larger, and randomized with replacement if the
///        combinatorial permutations are huge.
enum class SamplingStrategy {
  FULL,     ///< Sample all permutations of conformers
  LIMITED,  ///< Randomly sample conformers without replacement
  SPARSE    ///< Randomly sample conformers with no regard to whether replacement could occur
};

/// \brief Determine whether two ways of isomerizing a molecule are linked, whether by sharing an
///        atom in the chiral center or rotatable bond, or by having any of those atoms be, in
///        turn, bonded to one another.
///
/// \param isomerizers  List of atom groups involved in each isomer creation for the molecule
/// \param permi        The first origin of isomerization and thus permutations
/// \param permj        The second origin of isomerization and thus permutations
/// \param nbk          Contains the non-bonded exclusions, including bonded exclusions
bool permutationsAreLinked(const std::vector<IsomerPlan> &isomerizers, int permi, int permj,
                           const NonbondedKit<double> &nbk);
  
/// \brief Return an estimate of the number of all permutations covering isomerizations involving
///        the same atoms or atoms that are bonded to one another.
///
/// \param limits       The number of possible settings for each permutation
/// \param isomerizers  Groups of atoms affecting different conformers and isomers, including bond
///                     rotations
/// \param ag           Topology for the molecule of interest
double computeLocalPermutations(const std::vector<int> &limits,
                                const std::vector<IsomerPlan> &isomerizers, const AtomGraph *ag);
  
/// \brief Expand the initial list of systems into a complete list of initial states for the
///        population of conformers which conformer.omni will minimize in search of the
///        lowest-energy states.  This is done by parsing each topology into its chemical
///        details and then enumerating the rotatable bonds, cis-trans isomers, and chiral centers
///        which it could sample.
///
/// \param ui   User input settings, obtained from the input deck
/// \param sc   Cache of topologies and initial structures.  The list of topologies will not be
///             expanded by this procedure, but the list of structures will undergo a radical
///             expansion.
/// \param xrs  Random number generator (the master generator will guide the CPU-based,
///             coarse-grained conformer selection)
/// \param tm   Timer to record the wall time spent on various setup procedures
PhaseSpaceSynthesis expandConformers(const UserSettings &ui, const SystemCache &sc,
                                     Xoshiro256ppGenerator *xrs, StopWatch *tm);

} // namespace setup
} // namespace conf_app

#endif
