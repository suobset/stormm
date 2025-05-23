// -*-c++-*-
#ifndef MMGBSA_NAMELIST_H
#define MMGBSA_NAMELIST_H

#include <string>
#include <vector>
#include "../../../src/copyright.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Namelists/input.h"
#include "../../../src/Namelists/namelist_emulator.h"
#include "../../../src/Namelists/namelist_element.h"
#include "../../../src/Parsing/parsing_enumerators.h"
#include "../../../src/Parsing/textfile.h"
#include "mmgbsa_enumerators.h"

namespace mmgbsa {

using stormm::constants::ExceptionResponse;
using stormm::namelist::readNamelist;
using stormm::namelist::NamelistEmulator;
using stormm::parse::TextFile;
using stormm::parse::WrapTextSearch;
  
// The default weighting scheme is to have all poses of any given ligand averaged with equal
// weight.  This is appropraite when taking snapshots of a trajectory, where the populations in
// different configurations may already be Boltzmann-weighted.
const char default_weighting_scheme[] = "flat";
const char default_energy_report_depth[] = "compound_averages";
  
// The maximum number of poses to use when taking an average of many poses of the same ligand.
// This default of zero implies that all poses should be processed.
const int default_ligand_max_poses = 0;

// The number of best poses to report likewise takes a value of zero, directing all poses to be
// reported.  The default format for printing those poses is a Protein Data Bank file.  Separate
// files will be printed for each pose of a distinct ligand.
const int default_best_poses = 0;
const char default_structure_format[] = "pdb";
const char default_structure_base_name[] = "pose";

// The temperature at which to compute Boltzmann-weighted averages for the energies of multiple
// receptor, ligand, and complex states.
const double default_mmgbsa_temperature = 298.15;

// Default values for defining the receptor binding site carveout.  This is used to focus the
// calculation to a selection of atoms near the ligand or in critical loops of the target protein.
// The default mask is to select no (additional) atoms.  The rotating group size limit is based on
// the number of (real) atoms in a Tryptophan residue.
const double default_carveout_cutoff = 5.0;
const double default_carveout_halo = 2.0;
const char default_carveout_mask[] = "!:*";
const int default_rotating_group_limit = 18;
  
//-------------------------------------------------------------------------------------------------
// Class to collect user inputs taken from an &mmgbsa namelist
//-------------------------------------------------------------------------------------------------
class MMGBSAControls {
public:

  // The constructor can prepare an object with default settings or read the corresponding
  // namelist to accept user input.                                                         
  //
  // Arguments:
  //   tf:          Input file translated into RAM                                             
  //   start_line:  Line of the input file to begin searching for the &solvent namelist
  //   found_nml:   Indication of whether the namelist was found in the input file
  //   policy_in:   Requested error handling behavior
  //   wrap:        Indicate that the search for a &conformer namelist should carry on from
  //                the beginning of an input file if no such namelist is found starting
  //                from the original starting point
  MMGBSAControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  MMGBSAControls(const TextFile &tf, int *start_line, bool *found_nml,
                 ExceptionResponse policy_in = ExceptionResponse::DIE,
                 WrapTextSearch wrap = WrapTextSearch::NO);

  // As with other control objects, copy and move constructors, plus copy and move
  // assignment operators, can all take their default forms.
  //
  // Arguments:
  //   original:  The original object to copy or move
  //   other:     Another object placed on the right hand side of the assignment statement
  MMGBSAControls(const MMGBSAControls &original) = default;
  MMGBSAControls(MMGBSAControls &&original) = default;
  MMGBSAControls& operator=(const MMGBSAControls &original) = default;
  MMGBSAControls& operator=(MMGBSAControls &&original) = default;

  // Get the number of refined pose energies to report (0 indicates that all poses should be
  // reported).
  int getPoseReportCount() const;
  
  // Get the weighting scheme to be used in averaging the results for producing a final energy
  // for each ligand.
  WeightingScheme getWeightingScheme() const;

  // Get the temperature to be used in a Boltzmann-weighted average of the provided poses.  This
  // also applies if computing dissocation constants.
  double getTemperature() const;

  // Get the energy reporting depth, or comprehensiveness, for the data set.
  EnergyReport getEnergyReportDepth() const;

  // Get the carveout distance, in Angstroms, for selecting receptor atoms which shall be mobile.
  double getCarveoutCutoff() const;
  
  // Get the distance, in Angstroms, for selecting receptor atoms which shall be included as
  // explicit parts of a calculation but held immobile.
  double getCarveoutHalo() const;

  // Indicate whether the user has requested a proximity-based carveout.  The carveout can also be
  // expanded by an atom mask.
  bool proximityCarveout() const;
  
  // Get the user-provided atom mask describing additional atoms to include in the receptor's
  // mobile carveout.
  const std::string& getCarveoutExtraMask() const;

  // Get the size limit (in terms of the number of of real atoms) on rotating groups that can be
  // included in the carveout region if a subset of their atoms lies within the designated
  // distance of the reference ligand or ligands.
  int getRotatingGroupLimit() const;

  // Indicate whether the carveout should include complete residues.
  bool completeResidueCarveout() const;

  // Indicate whether the receptor structure should be taken from the minimized complex structure
  // or from a fresh minimization of the apo receptor.
  bool takeReceptorFromComplex() const;
  
  // Get the user-provided atom mask describing additional atoms to include in the receptor's
  // mobile carveout.
  const std::string& getCarveoutLigandReference() const;

  // Get the indication as to whether an energy decomposition should be printed in the final
  // analysis.
  bool printEnergyDecomposition() const;

  // Get an indication as to whether the final, energy-minimized structures should be printed.
  bool printFinalStructures() const;

  // Get the format for printed structures.
  PrintedPoseFormat getPrintedStructureFormat() const;

  // Get the base name of output structures.
  const std::string& getOutputStructureBase() const;
  
  // Set the number of poses (per ligand) to report.
  //
  // Arguments:
  //   ligand_max_poses_in:  The number of poses to set
  void setPoseReportCount(int ligand_max_poses_in);
  
  // Set the weighting scheme to use when averaging results among poses of an individual ligand.
  //
  // Overloaded:
  //   - Provide the choice as a Standard Template Library string, which will be translated into
  //     the proper enumeration
  //   - Provide the enumeration itself
  //
  // Arguments:
  //   averaging_in:  The desired weighting scheme
  void setWeightingScheme(const std::string &averaging_in);
  void setWeightingScheme(WeightingScheme averaging_in);

  // Set the temperature at which to perform Boltzmann-weighted averages and, if applicable, to
  // compute dissociation constants.
  //
  // Arguments:
  //   temperature_in:  The temperature to set
  void setTemperature(double temperature_in);
  
  // Set the depth of energy reporting.  If any ligands have more than one pose, all poses may be
  // reported while the final energies are reported using some averaging technique.  Overloading
  // follows from setWeightingScheme(), above.
  //
  // Arguments:
  //   report_depth_in:  The desired completeness of energy value reporting
  void setEnergyReportDepth(const std::string &report_depth_in);
  void setEnergyReportDepth(EnergyReport report_depth_in);

  // Set the flag for reporting a decomposition of the molecular mechanics energy components.
  //
  // Arguments:
  //   decompose_energy_in:  Set to true to report an energy decomposition
  void setEnergyDecompositionReporting(bool decompose_energy_in);

  // Set the flag for printing final structures of each ligand / receptor complex.
  //
  // Arguments:
  //   print_structures_in:  Set to true to have final structures printed in a format to be
  //                         defined by setPrintedStructureFormat()
  void setStructurePrinting(bool print_structures_in);

  // Set the format of printed structures.
  //
  // Overloaded:
  //   - Provide the choice as a Standard Template Library string
  //   - Provide the choice as an enumeration
  //
  // Arguments:
  //   structure_format_in:  The desired structure format
  void setFinalStructureFormat(const std::string &structure_format_in);
  void setFinalStructureFormat(PrintedPoseFormat structure_format_in);
  
  // Set the carveout distance-based cutoff.
  //
  // Arguments:
  //   carveout_cutoff_in:  The cutoff to set
  void setCarveoutCutoff(int carveout_cutoff_in);

  // Set the carveout distance-based halo.  Atoms in the halo will be included in non-bonded
  // calculations for accuracy and context, but immobile.
  //
  // Arguments:
  //   carveout_halo_in:  The cutoff to set
  void setCarveoutHalo(int carveout_halo_in);

  // Set the extra atom mask to be used in designing the receptor carveout region.
  //
  // Arguments:
  //   carveout_extra_mask_in:  The mask to use
  void setCarveoutExtraMask(const std::string &carveout_extra_mask_in);

  // Set the size limit (in terms of real atoms) on rotating groups that can be used to extend
  // the receptor's binding site carveout.
  //
  // Arguments:
  //   limit_in:  The rotating group size limit to set
  void setRotatingGroup_limit(int limit_in);
  
private:

  ExceptionResponse policy;    // The course of action to take in response to bad input 
  int ligand_max_poses;        // The number of poses to report, minimum one, or to use in each
                               //   ensemble calculation of a binding energy
  int best_poses;              // The number of best poses to report, in a separate section on
                               //   energies as well as in output structures
  WeightingScheme averaging;   // The method for averaging the energies of multiple structures, if
                               //   ligands have more than one structure under a single topology
  double temperature;          // The temperature at which to compute probabilities based on
                               //   computed energies, and at which to consider Boltzmann-weighted
                               //   averages, depending on other user-supplied settings
  double carveout_cutoff;      // The distance, in Angstroms, from the ligand's or ligands' initial
                               //   pose which determines the mobile atoms of the receptor.  If
                               //   no value is specified and a mask of receptor atoms is provided
                               //   instead, the mask will indicate the only receptor atoms to be
                               //   included.  If no mask is specified, the default value of this
                               //   cutoff will determine the carveout region.
  double carveout_halo;        // Atoms within this range of the carveout region will be immobile
                               //   but explicitly included in calculations.  Units of Angstroms.
  bool proximity_carveout;     // Flag to indicate that the carveout is (at least in part) defined
                               //   by proximity of receptor atoms to ligand atoms
  bool receptor_from_complex;  // Flag to have the "energy minimized" structure of the receptor
                               //   taken from the complex, not from an independent relaxation of
                               //   some original, apo structure of the receptor
  EnergyReport report_depth;   // The extent to which energies shall be reported in the final
                               //   analysis.  At a minimum, averages (one value per ligand) will
                               //   be reported, but if one or more ligands has multiple poses and
                               //   the report depth calls for a raw delivery of all data, two
                               //   sets of numbers may be reported.
  bool decompose_energy;       // Flag to indicate that an energy decomposition should be
                               //   presented.  If this is left unspecified (default FALSE), naming
                               //   energy components in the accompanying &report namelist will, in
                               //   effect, trigger energy printing.
  bool print_structures;       // Flag to indicate that the final, energy-minimized states of each
                               //   ligand and receptor complex should be printed.  The structures
                               //   are limited to the ligand and a subset of the receptor (the
                               //   binding site carveout) defined by various member variables,
                               //   above.
  
  // The format to use when printing sturcture results, e.g. PDB.
  PrintedPoseFormat structure_format;

  // The base name of the output structures is only relevant if structure printing is active.
  std::string structure_base_name;

  // A mask defining atoms in the receptor which shall be included in the carveout region,
  // irrespective of their proximity to the reference ligand or ligands.  If this mask is specified
  // and no carveout cutoff is given, the mask will be the sole definition of the carveout.
  std::string carveout_extra_mask;

  // The carveout can be extended if some rotatable groups have a portion of their atoms (but not
  // all) within the distance-based cutoff.  This is the maximum number of real atoms that a
  // rotatable group can have in order for its entirety to be included in the carveout region.  If
  // this is set to zero (the default), the carveout will still consist of selections from the
  // distance-based criterion and any user-specified atom mask.
  int carveout_rotating_group_limit;

  // Indicate that the carveout region should ultimately comprise all atoms of any residue that
  // the other criteria affect any part of.  This can help include entire loops in the definition
  // of the binding site.
  bool carveout_whole_residues;
  
  // Give the option to specify a particular ligand to take as the reference for generating the
  // carveout's distance-based mask.  If this string reads "all", then every ligand in the data
  // set will be used.  The string may also specify a diminution of the ligand's name, in which
  // case all ligands matching the substring will be taken.
  std::string carveout_ligand_reference;
  
  // Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;

  // Validate a temperature setting.  The temperature cannot be negative.
  void validateTemperature();

  // Validate a cutoff setting, either for the carveout cutoff or halo, by verifying that it is
  // greater than zero.  This function will reference the object's policy to direct its course of
  // action in the event of bad input.
  //
  // Arguments:
  //   cutoff:     The cutoff to set, modified and returned if the check passes
  //   cutoff_in:  The cutoff to check, ensuring that it
  //   caller:     Name of the calling function (for error tracing purposes)
  void validateDistanceCutoff(double* cutoff, double cutoff_in, const char* caller = nullptr);
};

//-------------------------------------------------------------------------------------------------
// A free function to prepare the &mmgbsa namelist, wrapping user commands that modulate the
// MM-GBSA workflow.  While minimization and details of the potential function are handled in
// namelists from the main STORMM libraries, this namelist will handle information like the number
// of best results to report and a preferred weighting scheme.
//
// Arguments:
//   tf          Text of file containing the input deck, read into RAM
//   start_line  Line of the input file at which to begin the scan
//   found       Indicator that the namelist was found in the input file
//   policy      Response to bad inputs
//   wrap        Indicate that the search for a &conformer namelist should carry on from the
//               beginning of an input file if no such namelist is found starting from the
//               original starting point
//-------------------------------------------------------------------------------------------------
NamelistEmulator mmgbsaInput(const TextFile &tf, int *start_line, bool *found,
                             ExceptionResponse policy = ExceptionResponse::DIE,
                             WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace mmgbsa

#endif
