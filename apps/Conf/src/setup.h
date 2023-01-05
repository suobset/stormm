// -*-c++-*-
#ifndef CONFORMER_SETUP_H
#define CONFORMER_SETUP_H

#include "../../../src/Chemistry/atommask.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Math/tickcounter.h"
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/Namelists/nml_conformer.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Random/random.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/UnitTesting/stopwatch.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::AtomMask;
using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::ChiralInversionProtocol;
using stormm::chemistry::IsomerPlan;
using stormm::constants::ExceptionResponse;
using stormm::energy::StaticExclusionMask;
using stormm::math::TickCounter;
using stormm::namelist::ConformerControls;
using stormm::namelist::UserSettings;
using stormm::random::Xoshiro256ppGenerator;
using stormm::structure::MdlMol;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SystemCache;
using stormm::testing::StopWatch;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::CoordinateFrame;

/// \brief Timings strings to track various setup procedures
/// \{
constexpr char tm_input_parsing[] = "Input parsing";
constexpr char tm_feature_detection[] = "Chemical feature detection";
constexpr char tm_coordinate_expansion[] = "Coordinate expansion";
/// \}
  
/// \brief Enumerate coarse-grained conformer sampling approaches: full if the combinatorial
///        permutations are few enough, randomized without replacement if the combinatorial
///        permutations are up to 32 times larger, and randomized with replacement if the
///        combinatorial permutations are huge.
enum class SamplingStrategy {
  FULL,     ///< Sample all permutations of conformers
  LIMITED,  ///< Randomly sample conformers without replacement
  SPARSE    ///< Randomly sample conformers with no regard to whether replacement could occur
};
  
/// \brief Get the core mask for one of the topologies in the cache of systems.  This will compute
///        the mask for one example of the systems using the topology, and as such should not
///        be used with core masks defined by coordinate-related details.
///
/// \param conf_input   The &conformer namelist user input
/// \param sdf_example  An MDL MOL structure containing coordinates and possibly a data item
///                     defining topology-specific core mask information
/// \param ps           An example of the coordinates, used if the MDL MOL entry is not available
/// \param ag           The topology of interest
/// \param chemfe       Chemical features computed for the topology of interest
/// \param policy       Indicate what to do if errors are encountered when making the mask (passed
///                     down from command-line input)
AtomMask getCoreMask(const ConformerControls &conf_input, const MdlMol &sdf_example,
                     const PhaseSpace &ps, const AtomGraph *ag, const ChemicalFeatures &chemfe,
                     const ExceptionResponse policy);

/// \brief Add implicit solvent models to topologies.  Set the position-restrained cores of each
///        proto-conformer in the system cache and return the list of such masks.
///
/// \param ui            User input settings, obtained from the input deck
/// \param sc            Cache of topologies and initial structures.  The topologies are modified
///                      temporarily to immobilize certain atoms when determining rotatable bond
///                      groups, then returned to their original mobility states.  Implicit solvent
///                      conditions are added to all topologies.
/// \param sdf_recovery  Recovered SD file entries corresponding to each system (or placeholders
///                      if the coordinates did not originate from an SD file)
/// \param tm            Timer to record the wall time spent on various setup procedures
void setGenerativeConditions(const UserSettings &ui, SystemCache *sc,
                             const std::vector<MdlMol> &sdf_recovery, StopWatch *tm);
  
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

/// \brief Create a conformation based on some original coordinate set and a prescription for how
///        to permute it.
/// 
/// \param cf                       The original coordinates, which will be cloned and manipulated
/// \param ag                       System topology
/// \param exclusion_mask           Mask of excluded non-bonded interactions
/// \param isomerizers              List of atom groups that move with each change of an
///                                 independent state enumerated by ptrack
/// \param ptrack                   Control object for enumerating each permuation (may be modified
///                                 iternaly, but will be returned in its original input state) 
/// \param chiral_centers           A list of the molecule's chiral centers
/// \param chiral_center_plans      Methods for inverting (or instructions not to invert) each
///                                 chiral center
/// \param chiral_variable_indices  Indices of each chiral center in the larger list of all
///                                 rotatable bonds, cis-trans invertible bonds, and invertible
///                                 chiral centers
/// \param invertible_groups        Groups of atoms that move upon inverting a chiral center
/// \param max_seeding_attempts     Number of allowed attempts for correcting a clashing initial
///                                 configuration
/// \param self_clash_ratio         Minimum ratio of interparticle distance to their respective
///                                 pair van-der Waals sigma values, before the interaction is
///                                 considered a clash
CoordinateFrame forgeConformation(const CoordinateFrame &cf, const AtomGraph *ag,
                                  const StaticExclusionMask &exclusion_mask,
                                  const std::vector<IsomerPlan> &isomerizers,
                                  TickCounter<double> *ptrack,
                                  const std::vector<int> &chiral_centers,
                                  const std::vector<ChiralInversionProtocol> &chiral_center_plans,
                                  const std::vector<int> &chiral_variable_indices,
                                  const std::vector<IsomerPlan> &invertible_groups,
                                  const int max_seeding_attempts, const double self_clash_ratio);

/// \brief Count the numbr of conformers that each system might conceivably produce, based on the
///        number of rotatable bonds, cis-trans isomers, and invertible chiral centers.  These
///        counts will serve to lay out space for the coordinate synthesis that can later be
///        populated and re-populated.
///
/// \param conf_input  Conformer input controls
/// \param sc          The cache of systems, supply ChemicalFeatures for each of them
/// \param tm          Timer to record the wall time spent on various setup procedures
std::vector<int> calculateReplicaCounts(const ConformerControls &conf_input,
                                        const SystemCache &sc, StopWatch *tm);

/// \brief Create the sandbox for working with conformers.  All coordinates for each system in this
///        workspace will be initialized from the systems cache, to then be modified.
///
/// \param replica_counts  Numbers of each system to include in the synthesis, based on rotamer
///                        counts determined in the calculateReplicaCounts() above.
/// \param sc              The cache of systems, supply ChemicalFeatures for each of them
/// \param tm              Timer to record the wall time spent on various setup procedures
PhaseSpaceSynthesis buildReplicaWorkspace(const std::vector<int> &replica_counts,
                                          const SystemCache &sc, StopWatch *tm);
  
/// \brief Expand the initial list of systems into a complete list of initial states for the
///        population of conformers which conformer.omni will minimize in search of the
///        lowest-energy states.  This is done by parsing each topology into its chemical
///        details and then enumerating the rotatable bonds, cis-trans isomers, and chiral centers
///        which it could sample.
///
/// \param poly_ps          Pre-allocated synthesis of coordinates (modified and returned)
/// \param replica_counts   Numbers of replicas anticipated for each system
/// \param ui               User input settings, obtained from the input deck
/// \param sc               Cache of topologies and initial structures.  The list of topologies will
///                         not be expanded by this procedure, but the list of structures will
///                         undergo a radical expansion.
/// \param exclusion_masks  List of exclusion masks for each unique topology in poly_ps
/// \param xrs              Random number generator (the master generator will guide the CPU-based,
///                         coarse-grained conformer selection)
/// \param tm               Timer to record the wall time spent on various setup procedures
void expandConformers(PhaseSpaceSynthesis *poly_ps, const std::vector<int> &replica_counts,
                      const UserSettings &ui, const SystemCache &sc,
                      const std::vector<StaticExclusionMask> &exclusion_masks,
                      Xoshiro256ppGenerator *xrs, StopWatch *tm);

} // namespace setup
} // namespace conf_app

#endif
