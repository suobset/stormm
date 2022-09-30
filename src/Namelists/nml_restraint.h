// -*-c++-*-
#ifndef STORMM_NML_RST_H
#define STORMM_NML_RST_H

#include "copyright.h"
#include "input.h"
#include "namelist_emulator.h"
#include "namelist_enumerations.h"
#include "Chemistry/chemical_features.h"
#include "Parsing/textfile.h"
#include "Restraints/bounded_restraint.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace namelist {

using chemistry::ChemicalFeatures;
using restraints::BoundedRestraint;
using topology::AtomGraph;
using trajectory::CoordinateFrameReader;
  
/// \brief Object to encapsulate and dispense restraint information collected from a single
///        &restraint namelist.
class RestraintControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param policy_in   Requested error handling behavior
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \{
  RestraintControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  RestraintControls(const TextFile &tf, int *start_line,
                    ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

  /// \brief Get the restraint specified by this namelist.  There are no other getter functions
  ///        for individual details, nor setters to manually edit the result.  Restraints can be
  ///        built by many means, namelists being only one, and editing should be done on the
  ///        more tractable BoundedRestraint form.
  ///
  /// Overloaded:
  ///   - Produce a restraint based on atom index specifications
  ///   - Produce a restraint based on atom masks (probably from user input)
  ///
  /// \param ag      The system topology (needed to check the atom indexing--this member function
  ///                can only be called once the topology to which the restraint shall apply has
  ///                been constructed, which may be well after user input has been transcribed)
  /// \param chemfe  Chemical features of the system (needed to evaluate atom masks)
  /// \param cfr     Coordinates of the system (needed to evaluate atom masks)
  /// \{
  BoundedRestraint getRestraint(const AtomGraph *ag, const ChemicalFeatures *chemfe,
                                const CoordinateFrameReader &cfr) const;

  BoundedRestraint getRestraint(const AtomGraph *ag) const;
  /// \}  

private:
  ExceptionResponse policy;  ///< Protocol to follow in the event of bad input data
  bool restraint_is_valid;   ///< Indicator that the restraints passes various quality checks
  std::string system;        ///< A system label that can be matched to various entries based on
                             ///<   the -sys keyword in &files namelist input
  std::string ensemble;      ///< Indicates a collection of restraints to apply throughout any
                             ///<   specified systems.  This will supercede any atom indices or
                             ///<   masks specified in the same &restraint namelist.
  std::string mask_i;        ///< Atom I (first atom) mask
  std::string mask_j;        ///< Atom J (second atom) mask
  std::string mask_k;        ///< Atom K (third atom) mask
  std::string mask_l;        ///< Atom L (fourth atom) mask
  int atom_i;                ///< Atom I topological index
  int atom_j;                ///< Atom J topological index
  int atom_k;                ///< Atom K topological index
  int atom_l;                ///< Atom L topological index
  int order;                 ///< Order of the restraint: 0 = no restraint, 1 = positional
                             ///<   restraint, 2 = distance restraint, 3 = angle restraint,
                             ///<   4 = dihedral restraint
  int initiation_step;       ///< Step at which restraint application begins
  int maturation_step;       ///< Step at which restraint application becomes mature
  double initial_k2;         ///< Initial (or static) left-hand stiffness constant
  double initial_k3;         ///< Initial (or static) right-hand stiffness constant
  double initial_r1;         ///< Initial (or static) leftmost harmonic boundary of k2's support
  double initial_r2;         ///< Initial (or static) rightmost harmonic boundary of k2's support
  double initial_r3;         ///< Initial (or static) leftmost harmonic boundary of k3's support
  double initial_r4;         ///< Initial (or static) rightmost harmonic boundary of k3's support
  double mature_k2;          ///< Long-term value of k2, but only if k2 is not static
  double mature_k3;          ///< Long-term value of k3, but only if k3 is not static
  double mature_r1;          ///< Long term value of r1, but only if r1 is not static
  double mature_r2;          ///< Long-term value of r2, but only if r2 is not static
  double mature_r3;          ///< Long-term value of r3, but only if r3 is not static
  double mature_r4;          ///< Long-term value of r4, but only if r4 is not static
  double3 initial_crd;       ///< Initial (or static) Cartesian tether for a positional restraint
  double3 mature_crd;        ///< Long-term Cartesian tether for a positional restraint

  /// \brief Get a verdict on whether this object specifies atoms by masks or by indices.
  RestraintAnchoring getAtomSpecification() const;

  /// \brief Enforce that all atoms of a restraint be specified either with topological indices or
  ///        atom masks.  If there are mixed specifications, the indices can be converted to atom
  ///        masks without any need to see the actual topology, so do that.
  void enforceSpecification();

  /// \brief Check critical stiffness and equilibrium constants for time-dependent restraints.
  ///        Exit with error if requested, otherwise try to fix the values based on initial
  ///        settings.
  ///
  /// \param nstep1_found   Indicator of whether there is time dependency.  Any differences between
  ///                       the fact that nstep1 and nstep2 (the initiation and maturation time
  ///                       steps) have been set will have already been addressed.
  /// \param r2a_found      Indicator of whether the final r2 parameter has been specified
  /// \param r3a_found      Indicator of whether the final r3 parameter has been specified
  /// \param k2a_found      Indicator of whether the final k2 parameter has been specified
  /// \param k3a_found      Indicator of whether the final k3 parameter has been specified
  /// \param starting_line  Line of the input file at which the &restraint namelist begins
  /// \param filename       Name of the input file where a problematic &restraint namelist may have
  ///                       been found
  void checkFinalRestraintSettings(bool nstep1_found, bool r2a_found, bool r3a_found,
                                   bool k2a_found, bool k3a_found, int starting_line,
                                   const std::string &filename);
};
 
/// \brief Produce a namelist for specifying an NMR restraint, equivalent to the eponymous namelist
///        in sander of pmemd.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is often intended to be repeatable)
/// \param policy      Reaction to exceptions encountered during namelist reading
NamelistEmulator restraintInput(const TextFile &tf, int *start_line,
                                ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace namelist
} // namespace stormm

#endif
