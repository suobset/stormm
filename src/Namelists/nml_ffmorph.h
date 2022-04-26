// -*-c++-*-
#ifndef OMNI_NML_FFMORPH_H
#define OMNI_NML_FFMORPH_H

#include "Constants/behavior.h"
#include "ForceField/forcefield_element.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace omni {
namespace namelist {

using modeling::ForceFieldElement;
  
/// \brief Object to encapsulate force field morphing operations.  Take information from an input
///        file or a series of setters and validate each piece of data as it appears with private
///        member functions.
class FFMorphControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &ffmorph namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \{
  FFMorphControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  FFMorphControls(const TextFile &tf, int *start_line, bool *found_nml,
                  ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

private:
  ExceptionResponse policy;       ///< Set the behavior when bad inputs are encountered.  DIE =
                                  ///<   abort program, WARN = warn the user, and likely reset to
                                  ///<   the default value if one is available, SILENT = do not
                                  ///<   warn the user, but also likely reset to the default value
                                  ///<   if one is available.

  /// Harmonic bond parameters
  std::vector<ForceFieldElement> harmonic_bonds;

  /// Harmonic angle parameters, three-point angle bending with harmonic penalties based on
  /// deviations of the I-J-K angle from an ideal value
  std::vector<ForceFieldElement> harmonic_angles;

  /// Cosine-based dihedral parameters (proper or improper)
  std::vector<ForceFieldElement> cosine_dihedrals;

  /// Urey-Bradley angle parameters, two-point stretching terms with harmonic penalties based on
  /// deviations of the I-K angle from an ideal separation
  std::vector<ForceFieldElement> urey_bradley_angles;  

  /// CHARMM improper dihedral parameters
  std::vector<ForceFieldElement> charmm_impropers;

  /// CMAP bicubic spline surface terms
  std::vector<ForceFieldElement> cmap_surfaces;

  /// Lennard-Jones parameters, or other van-der Waals approximations
  std::vector<ForceFieldElement> van_der_waals_properties;

  /// Charge parameters
  std::vector<ForceFieldElement> charge_properties;

  /// Virtual site frames--the non-bonded properties of the virtual site are controlled by
  /// entries in van_der_waals_properties and charge_properties, respectively.
  std::vector<ForceFieldElement> virtual_sites;
};
  
} // namespace namelist
} // namespave omni

#endif
