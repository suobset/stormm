// -*-c++-*-
#ifndef OMNI_NML_RST_H
#define OMNI_NML_RST_H

#include "input.h"
#include "namelist_emulator.h"
#include "namelist_enumerations.h"
#include "Parsing/textfile.h"
#include "Restraints/bounded_restraint.h"

namespace omni {
namespace namelist {

using restraints::BoundedRestraint;
  
/// \brief Object to encapsulate and dispense restraint information collected from a single
///        &restraint namelist.
struct RestraintControls {

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

  /// \brief Get a verdict on whether this object specifies atoms by masks or by indices.
  RestraintAnchoring getAtomSpecification() const;

  /// \brief Get the restraint specified by this namelist.  There are no other getter functions
  ///        for individual details, nor setters to manually edit the result.  Restraints can be
  ///        built by many means, namelists being only one, and editing should be done on the
  ///        more tractable BoundedRestraint form.
  BoundedRestraint getRestraint() const;

private:
  ExceptionResponse policy;  ///< Protocol to follow in the event of bad input data
  std::string input_file;    ///< Name of the input file containing the &restraint namelist
  int starting_line;         ///< Starting line for the &restraint namelist in the input file
                             ///<   (numbering starts at line 1)
  std::string mask_i;        ///< Atom I (first atom) mask
  std::string mask_j;        ///< Atom J (second atom) mask
  std::string mask_k;        ///< Atom K (third atom) mask
  std::string mask_l;        ///< Atom L (fourth atom) maskk
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

  /// \brief Validate the information in this namelist--can it make a valid restraint?
  void validateRestraint() const;  
};
 
/// \brief Produce a namelist for specifying an NMR restraint, equivalent to the eponymous namelist
///        in sander of pmemd.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
/// \param policy      Reaction to exceptions encountered during namelist reading
NamelistEmulator restraintInput(const TextFile &tf, int *start_line,
                                ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace namelist
} // namespace omni

#endif
