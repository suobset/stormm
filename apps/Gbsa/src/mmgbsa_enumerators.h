// -*-c++-*-
#ifndef MMGBSA_UTIL_H
#define MMGBSA_UTIL_H

#include <string>
#include <vector>
#include "../../../src/copyright.h"

namespace mmgbsa {

//-------------------------------------------------------------------------------------------------
// Enumerate the different types of problems in the problem set below.
//-------------------------------------------------------------------------------------------------
enum class MMGBSAProblemKind {
  LIGAND,    // The problem is an instance of one of the ligands
  RECEPTOR,  // The problem is an instance of the _apo_ receptor
  COMPLEX    // The problem is a complex of one of the ligand poses and the receptor
};

//-------------------------------------------------------------------------------------------------
// Enumerate the different weighting schemes for averaging results within each unique ligand.  The
// number of poses averaged in each case may be subject to further restrictions.
//-------------------------------------------------------------------------------------------------
enum class WeightingScheme {
  BOLTZMANN,  // Apply a Boltzmann weight to each refined pose, based on its total energy relative
              //   to others
  FLAT        // Give every refined pose the same weight in the final everage
};

//-------------------------------------------------------------------------------------------------
// Enumerate the choices for reporting energies.
//-------------------------------------------------------------------------------------------------
enum class EnergyReport {
  COMPLETE,           // Report all energy results for every instance of each ligand pose
  COMPOUND_AVERAGES,  // Report the averages, subject to the weighting scheme above, for each
                      //   uniqueligand in the test set
  BEST_RESULTS        // Report only the best results out of all ligands
};

//-------------------------------------------------------------------------------------------------
// Enumerate the choices for reporting docked pose structures.  The binding site carveout and
// ligand will be included.
//-------------------------------------------------------------------------------------------------
enum class PrintedPoseFormat {
  PDB,   // Report structures in Protein Data Bank format
  AMBER  // Report structures as Amber coordinates with unique Amber topology files for the
         //   reduced ligand / binding site system
};
  
//-------------------------------------------------------------------------------------------------
// Return a humand-readable string indicating one enumeration of the problem type.
//
// Arguments:
//   input:  The problem type enumeration to translate
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(MMGBSAProblemKind input);
std::string getEnumerationName(WeightingScheme input);
std::string getEnumerationName(EnergyReport input);
std::string getEnumerationName(PrintedPoseFormat input);

//-------------------------------------------------------------------------------------------------
// Translate a human-readable (and likely user-supplied) character string into one of the
// enumerations for the WeightingScheme class.
//
// Arguments:
//   input:  The string to translate
//-------------------------------------------------------------------------------------------------
WeightingScheme translateWeightingScheme(const std::string &input);

//-------------------------------------------------------------------------------------------------
// Translate a human-readable (and likely user-supplied) character string into an enumeration to
// guide the number of poses per ligand on which to report energies and structures.
//
// Arguments:
//   input:  The string to translate
//-------------------------------------------------------------------------------------------------
EnergyReport translateEnergyReport(const std::string &input);

//-------------------------------------------------------------------------------------------------
// Translate a human-readable (and likely user-supplied) character string into an enumeration of
// the output structure type.
//
// Arguments:
//   input:  The string to translate
//-------------------------------------------------------------------------------------------------
PrintedPoseFormat translatePrintedPoseFormat(const std::string &input);
  
} // namespace mmgbsa

#endif

