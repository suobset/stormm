#include <fstream>
#include <iostream>
#include <time.h>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parsing_enumerators.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

using diskutil::DataFormat;
using diskutil::openOutputFile;
using parse::NumberFormat;
using parse::stringToChar4Vector;
using parse::polyNumericVector;
using parse::printNumberSeries;
using stmath::maxAbsValue;
  
//-------------------------------------------------------------------------------------------------
void AtomGraph::printToFile(const std::string &output_file, const TopologyKind output_style,
                            const PrintSituation expectation,
                            const ExceptionResponse pr_policy) const {

  // It is not permitted to append any of theavailable topology formats with additional data.
  PrintSituation actual_expectation;
  switch (expectation) {
  case PrintSituation::OVERWRITE:
    actual_expectation = PrintSituation::OVERWRITE;
    break;
  case PrintSituation::OPEN_NEW:
  case PrintSituation::UNKNOWN:
    actual_expectation = PrintSituation::OPEN_NEW;
    break;
  case PrintSituation::APPEND:
    switch (pr_policy) {
    case ExceptionResponse::DIE:
      rtErr("It is not permissible to append a topology file (" + output_file + ") with "
            "additional data.", "AtomGraph", "printToFile");
    case ExceptionResponse::WARN:
      rtWarn("It is not permissible to append a topology file (" + output_file + ") with "
             "additional data.  The file will be opened as new, so long as it does not already "
             "exist.", "AtomGraph", "printToFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    actual_expectation = PrintSituation::OPEN_NEW;
  }

  // Branch for various sub-functions to print each style of topology.    
  switch (output_style) {
  case TopologyKind::AMBER:
    printAmberFormat(output_file, expectation, pr_policy);
    break;
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    rtErr("Printing topologies in " + getEnumerationName(output_style) + " format is not yet "
          "supported.", "AtomGraph", "printToFile");
    break;
  case TopologyKind::UNKNOWN:
    switch (pr_policy) {
    case ExceptionResponse::DIE:
      rtErr("A valid topology format must be specified for " + output_file + ".", "AtomGraph",
            "printToFile");
    case ExceptionResponse::WARN:
      rtWarn("A valid topology format must be specified for " + output_file + ".  The topology "
             "will be printed in " + getEnumerationName(TopologyKind::AMBER) + " format.",
             "AtomGraph", "printToFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    printAmberFormat(output_file, expectation, pr_policy);
  }
}
  
//-------------------------------------------------------------------------------------------------
void AtomGraph::printAmberFormat(const std::string &output_file, const PrintSituation expectation,
                                 const ExceptionResponse pr_policy) const {
  std::ofstream foutp = openOutputFile(output_file, expectation, "write an AtomGraph object in "
                                       "AMBER format", DataFormat::ASCII);
  char date_string[80];
  strftime(date_string, sizeof(date_string), "%Y/%m/%d %H:%M:%S", &date);
  const int nchar_date_string = strlen(date_string);
  for (int i = 0; i < nchar_date_string - 2; i++) {
    date_string[i] = date_string[i + 2];
  }
  for (int i = nchar_date_string - 2; i <= nchar_date_string; i++) { 
    date_string[i] = '\0';
  }
  const std::string vs_line(std::string("%VERSION  VERSION_STAMP = ") + version_stamp +
                            "  DATE = " + std::string(date_string) + "\n");
  foutp.write(vs_line.data(), vs_line.size());
  printAmberSection(&foutp, "TITLE", polyNumericVector(stringToChar4Vector(title, ' ')), 20, 4, 0,
                    NumberFormat::CHAR4, "print the title to an AMBER format topology");
  printAmberSection(&foutp, "POINTERS", polyNumericVector(descriptors.readHost()), 10, 8, 0,
                    NumberFormat::INTEGER, "print general descriptors to an AMBER format "
                    "topology");
  printAmberSection(&foutp, "ATOM_NAME", polyNumericVector(atom_names.readHost()), 20, 4, 0,
                    NumberFormat::CHAR4, "print atom names to an AMBER format topology");
  std::vector<double> tmp_atomic_charges = atomic_charges.readHost();
  for (int i = 0; i < atom_count; i++) {
    tmp_atomic_charges[i] *= 18.2223;
  }
  printAmberSection(&foutp, "CHARGE", polyNumericVector(tmp_atomic_charges), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print atomic partial charges to an AMBER format "
                    "topology");
  printAmberSection(&foutp, "ATOMIC_NUMBER", polyNumericVector(atomic_numbers.readHost()), 10, 8,
                    0, NumberFormat::INTEGER, "print atomic numbers to an AMBER format topology.");
  printAmberSection(&foutp, "MASS", polyNumericVector(atomic_masses.readHost()), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print atomic masses to an AMBER format topology");
  std::vector<int> tmp_atom_type_indices = lennard_jones_indices.readHost(0, atom_count);
  for (int i = 0; i < atom_count; i++) {
    tmp_atom_type_indices[i] += 1;
  }
  printAmberSection(&foutp, "ATOM_TYPE_INDEX", polyNumericVector(tmp_atom_type_indices), 10, 8, 0,
                    NumberFormat::INTEGER, "print atom type indices to an AMBER format topology");

  // Re-derive the excluded atom counts and list
  std::vector<int> tmp_exclusion_counts(atom_count);
  std::vector<int> tmp_exclusion_list(total_exclusions);
  const int* atom_excl_bnds_ptr = atom_exclusion_bounds.data();
  const int* atom_excl_list_ptr = atom_exclusion_list.data();
  int n_zero_excl = 0;
  for (int i = 0; i < atom_count; i++) {
    n_zero_excl += (atom_excl_bnds_ptr[i + 1] == atom_excl_bnds_ptr[i]);
  }
  if (total_exclusions != atom_exclusion_bounds.readHost(atom_count) + n_zero_excl) {
    rtErr("The number of total exclusions (" + std::to_string(total_exclusions) + ") is not "
          "reflected in the upper bound of exclusions for the final atom in the list (" +
          std::to_string(atom_exclusion_bounds.readHost(atom_count)) + ") when printing a "
          "topology to AMBER format file " + output_file + ".", "AtomGraph");
  }
  int econ = 0;
  for (int i = 0; i < atom_count; i++) {
    tmp_exclusion_counts[i] = std::max(atom_excl_bnds_ptr[i + 1] - atom_excl_bnds_ptr[i], 1);
    if (atom_excl_bnds_ptr[i + 1] > atom_excl_bnds_ptr[i]) {
      for (int j = atom_excl_bnds_ptr[i]; j < atom_excl_bnds_ptr[i + 1]; j++) {
        tmp_exclusion_list[econ] = atom_excl_list_ptr[j] + 1;
        econ++;
      }
    }
    else {
      tmp_exclusion_list[econ] = 0;
      econ++;
    }
  }
  printAmberSection(&foutp, "NUMBER_EXCLUDED_ATOMS", polyNumericVector(tmp_exclusion_counts), 10,
                    8, 0, NumberFormat::INTEGER, "print atom type indices to an AMBER format "
                    "topology");
  std::vector<int> tmp_nbparm_map(lj_type_count * lj_type_count);
  int ljidx_con = 1;
  for (int i = 0; i < lj_type_count; i++) {
    for (int j = 0; j <= i; j++) {
      tmp_nbparm_map[(i * lj_type_count) + j] = ljidx_con;
      tmp_nbparm_map[(j * lj_type_count) + i] = ljidx_con;
      ljidx_con++;
    }
  }
  printAmberSection(&foutp, "NONBONDED_PARM_INDEX", polyNumericVector(tmp_nbparm_map), 10,
                    8, 0, NumberFormat::INTEGER, "print the matrix of non-bonded parameter "
                    "indices (easily calculable at runtime, but printed out) to an AMBER format "
                    "topology");
  printAmberSection(&foutp, "RESIDUE_LABEL", polyNumericVector(residue_names.readHost()), 20, 4, 0,
                    NumberFormat::CHAR4, "print a topology's residue names to an AMBER format "
                    "file");
  std::vector<int> tmp_residue_starts = residue_limits.readHost(0, residue_count);
  for (int i = 0; i < residue_count; i++) {
    tmp_residue_starts[i] += 1;
  }
  printAmberSection(&foutp, "RESIDUE_POINTER", polyNumericVector(tmp_residue_starts), 10, 8, 0,
                    NumberFormat::INTEGER, "print the first atoms of each residue to an AMBER "
                    "format topology");
  printAmberSection(&foutp, "BOND_FORCE_CONSTANT", polyNumericVector(bond_stiffnesses.readHost()),
                    5, 16, 8, NumberFormat::SCIENTIFIC, "print the bond stiffnesses to an AMBER "
                    "format topology");
  printAmberSection(&foutp, "BOND_EQUIL_VALUE", polyNumericVector(bond_equilibria.readHost()),
                    5, 16, 8, NumberFormat::SCIENTIFIC, "print the bond equilibrium lengths to "
                    "an AMBER format topology");
  printAmberSection(&foutp, "ANGLE_FORCE_CONSTANT",
                    polyNumericVector(angl_stiffnesses.readHost()), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the bond stiffnesses to an AMBER format "
                    "topology");
  printAmberSection(&foutp, "ANGLE_EQUIL_VALUE", polyNumericVector(angl_equilibria.readHost()),
                    5, 16, 8, NumberFormat::SCIENTIFIC, "print the angle equilibrium lengths to "
                    "an AMBER format topology");
  printAmberSection(&foutp, "DIHEDRAL_FORCE_CONSTANT",
                    polyNumericVector(dihe_amplitudes.readHost()), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the dihedral amplitude parameters to "
                    "an AMBER format topology");
  printAmberSection(&foutp, "DIHEDRAL_PERIODICITY",
                    polyNumericVector(dihe_periodicities.readHost()), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the dihedral periodicity parameters to "
                    "an AMBER format topology");
  printAmberSection(&foutp, "DIHEDRAL_PHASE", polyNumericVector(dihe_phase_angles.readHost()), 5,
                    16, 8, NumberFormat::SCIENTIFIC, "print the dihedral phase angle parameters "
                    "to an AMBER format topology");

  // Recover the 1:4 scaling factors for each dihedral parameter set
  std::vector<double> tmp_scee_factors(dihe_parameter_count);
  std::vector<double> tmp_scnb_factors(dihe_parameter_count);
  for (int i = 0; i < dihe_term_count; i++) {
    const int dihe_parm_idx = dihe_parameter_indices.readHost(i);
    const int attn_parm_idx = dihe14_parameter_indices.readHost(i);
    tmp_scee_factors[dihe_parm_idx] = attn14_elec_factors.readHost(attn_parm_idx);
    tmp_scnb_factors[dihe_parm_idx] = attn14_vdw_factors.readHost(attn_parm_idx);
  }
  printAmberSection(&foutp, "SCEE_SCALE_FACTOR", polyNumericVector(tmp_scee_factors), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the electrostatic 1:4 scaling factors for "
                    "each dihedral parameter set to an AMBER format topology");
  printAmberSection(&foutp, "SCNB_SCALE_FACTOR", polyNumericVector(tmp_scnb_factors), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the van-der Waals 1:4 scaling factors for "
                    "each dihedral parameter set to an AMBER format topology");
  printAmberSection(&foutp, "SOLTY", polyNumericVector(solty_info.readHost()), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the unused atom type information (SOLTY) to "
                    "an AMBER format topology");
  const int ljmat_tri_count = (lj_type_count * (lj_type_count + 1)) / 2;
  std::vector<double> tmp_lj_a_values(ljmat_tri_count);
  std::vector<double> tmp_lj_b_values(ljmat_tri_count);
  std::vector<double> tmp_lj_c_values(ljmat_tri_count);
  std::vector<double> tmp_lj_14_a_values(ljmat_tri_count);
  std::vector<double> tmp_lj_14_b_values(ljmat_tri_count);
  std::vector<double> tmp_lj_14_c_values(ljmat_tri_count);
  int nljv = 0;
  for (int i = 0; i < lj_type_count; i++) {
    for (int j = 0; j <= i; j++) {
      tmp_lj_a_values[nljv] = lj_a_values.readHost((i * lj_type_count) + j);
      tmp_lj_b_values[nljv] = lj_b_values.readHost((i * lj_type_count) + j);
      tmp_lj_c_values[nljv] = lj_c_values.readHost((i * lj_type_count) + j);
      tmp_lj_14_a_values[nljv] = lj_14_a_values.readHost((i * lj_type_count) + j);
      tmp_lj_14_b_values[nljv] = lj_14_b_values.readHost((i * lj_type_count) + j);
      tmp_lj_14_c_values[nljv] = lj_14_c_values.readHost((i * lj_type_count) + j);
      nljv++;
    }
  }
  printAmberSection(&foutp, "LENNARD_JONES_ACOEF", polyNumericVector(tmp_lj_a_values), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the Lennard-Jones A coefficients to an AMBER "
                    "format topology");
  printAmberSection(&foutp, "LENNARD_JONES_BCOEF", polyNumericVector(tmp_lj_b_values), 5, 16, 8,
                    NumberFormat::SCIENTIFIC, "print the Lennard-Jones B coefficients to an AMBER "
                    "format topology");
  bool has_lj_ccoef =  (maxAbsValue(tmp_lj_c_values) > 1.0e-6);
  if (has_lj_ccoef) {
    printAmberSection(&foutp, "LENNARD_JONES_CCOEF", polyNumericVector(tmp_lj_c_values), 5, 16, 8,
                      NumberFormat::SCIENTIFIC, "print the Lennard-Jones C coefficients to an "
                      "AMBER format topology");
  }
  if (has_14_lennard_jones_data) {
    printAmberSection(&foutp, "LENNARD_JONES_14_ACOEF", polyNumericVector(tmp_lj_14_a_values), 5,
                      16, 8, NumberFormat::SCIENTIFIC, "print the 1:4 Lennard-Jones A "
                      "coefficients to an AMBER format topology");
    printAmberSection(&foutp, "LENNARD_JONES_14_BCOEF", polyNumericVector(tmp_lj_14_b_values), 5,
                      16, 8, NumberFormat::SCIENTIFIC, "print the 1:4 Lennard-Jones B "
                      "coefficients to an AMBER format topology");
    if (has_lj_ccoef) {
      printAmberSection(&foutp, "LENNARD_JONES_14_CCOEF", polyNumericVector(tmp_lj_14_c_values), 5,
                        16, 8, NumberFormat::SCIENTIFIC, "print the 1:4 Lennard-Jones C "
                        "coefficients to an AMBER format topology");
    }
  }

  // Loop back over the valence terms to assemble those with and without hydrogen.
  const ChemicalDetailsKit cdk = this->getChemicalDetailsKit();
  const ValenceKit vk = this->getDoublePrecisionValenceKit();
  std::vector<int> tmp_bond_with_h(3 * bond_term_with_hydrogen);
  std::vector<int> tmp_bond_without_h(3 * bond_term_without_hydrogen);
  int n_wh = 0;
  int n_woh = 0;
  for (int pos = 0; pos < bond_term_count; pos++) {
    const int atom_i = vk.bond_i_atoms[pos];
    const int atom_j = vk.bond_j_atoms[pos];
    const int b_parm = vk.bond_param_idx[pos];

    // One possible edge case is that there is a bond between a hydrogen and a virtual site.  That
    // will be counted as a "bond including hydrogen" but it is not a bond in any chemical sense.
    if (cdk.z_numbers[atom_i] == 1 || cdk.z_numbers[atom_j] == 1) {
      if (n_wh < bond_term_with_hydrogen) {
        tmp_bond_with_h[(3 * n_wh)    ] = 3 * atom_i;
        tmp_bond_with_h[(3 * n_wh) + 1] = 3 * atom_j;
        tmp_bond_with_h[(3 * n_wh) + 2] = b_parm + 1;
      }
      n_wh++;
    }
    else {
      if (n_woh < bond_term_without_hydrogen) {
        tmp_bond_without_h[(3 * n_woh)    ] = 3 * atom_i;
        tmp_bond_without_h[(3 * n_woh) + 1] = 3 * atom_j;
        tmp_bond_without_h[(3 * n_woh) + 2] = b_parm + 1;
      }
      n_woh++;
    }
  }
  if (n_wh != bond_term_with_hydrogen || n_woh != bond_term_without_hydrogen) {
    rtErr("The numbers of bonds counted with and without hydrogen (" + std::to_string(n_wh) +
          " and " + std::to_string(n_woh) + ") do not agree with the topology's internal records "
          "(" + std::to_string(bond_term_with_hydrogen) + " and " +
          std::to_string(bond_term_without_hydrogen) + ").  File to print: " + output_file + ".",
          "AtomGraph", "printAmberFormat");
  }
  printAmberSection(&foutp, "BONDS_INC_HYDROGEN", polyNumericVector(tmp_bond_with_h), 10, 8, 0,
                    NumberFormat::INTEGER, "print the atom and parameter indexing for bonds with "
                    "hydrogen to an AMBER format topology");
  printAmberSection(&foutp, "BONDS_WITHOUT_HYDROGEN", polyNumericVector(tmp_bond_without_h), 10, 8,
                    0, NumberFormat::INTEGER, "print the atom and parameter indexing for bonds "
                    "with hydrogen to an AMBER format topology");
  std::vector<int> tmp_angl_with_h(4 * angl_term_with_hydrogen);
  std::vector<int> tmp_angl_without_h(4 * angl_term_without_hydrogen);
  n_wh = 0;
  n_woh = 0;
  for (int pos = 0; pos < angl_term_count; pos++) {
    const int atom_i = vk.angl_i_atoms[pos];
    const int atom_j = vk.angl_j_atoms[pos];
    const int atom_k = vk.angl_k_atoms[pos];
    const int a_parm = vk.angl_param_idx[pos];

    // One possible edge case is that there is a bond between a hydrogen and a virtual site.  That
    // will be counted as a "bond including hydrogen" but it is not a bond in any chemical sense.
    if (cdk.z_numbers[atom_i] == 1 || cdk.z_numbers[atom_j] == 1 || cdk.z_numbers[atom_k] == 1) {
      if (n_wh < angl_term_with_hydrogen) {
        tmp_angl_with_h[(4 * n_wh)    ] = 3 * atom_i;
        tmp_angl_with_h[(4 * n_wh) + 1] = 3 * atom_j;
        tmp_angl_with_h[(4 * n_wh) + 2] = 3 * atom_k;
        tmp_angl_with_h[(4 * n_wh) + 3] = a_parm + 1;
      }
      n_wh++;
    }
    else {
      if (n_woh < angl_term_without_hydrogen) {
        tmp_angl_without_h[(4 * n_woh)    ] = 3 * atom_i;
        tmp_angl_without_h[(4 * n_woh) + 1] = 3 * atom_j;
        tmp_angl_without_h[(4 * n_woh) + 2] = 3 * atom_k;
        tmp_angl_without_h[(4 * n_woh) + 3] = a_parm + 1;
      }
      n_woh++;
    }
  }
  if (n_wh != angl_term_with_hydrogen || n_woh != angl_term_without_hydrogen) {
    rtErr("The numbers of bond angles counted with and without hydrogen (" + std::to_string(n_wh) +
          " and " + std::to_string(n_woh) + ") do not agree with the topology's internal records "
          "(" + std::to_string(angl_term_with_hydrogen) + " and " +
          std::to_string(angl_term_without_hydrogen) + ").  File to print: " + output_file + ".",
          "AtomGraph", "printAmberFormat");
  }
  printAmberSection(&foutp, "ANGLES_INC_HYDROGEN", polyNumericVector(tmp_angl_with_h), 10, 8, 0,
                    NumberFormat::INTEGER, "print the atom and parameter indexing for bond angles "
                    "with hydrogen to an AMBER format topology");
  printAmberSection(&foutp, "ANGLES_WITHOUT_HYDROGEN", polyNumericVector(tmp_angl_without_h), 10,
                    8, 0, NumberFormat::INTEGER, "print the atom and parameter indexing for bond "
                    "angles with hydrogen to an AMBER format topology");
  std::vector<int> tmp_dihe_with_h(5 * dihe_term_with_hydrogen);
  std::vector<int> tmp_dihe_without_h(5 * dihe_term_without_hydrogen);
  n_wh = 0;
  n_woh = 0;
  for (int pos = 0; pos < dihe_term_count; pos++) {
    const int atom_i = vk.dihe_i_atoms[pos];
    const int atom_j = vk.dihe_j_atoms[pos];
    const int atom_k = vk.dihe_k_atoms[pos];
    const int atom_l = vk.dihe_l_atoms[pos];
    const int h_parm = vk.dihe_param_idx[pos];

    // One possible edge case is that there is a bond between a hydrogen and a virtual site.  That
    // will be counted as a "bond including hydrogen" but it is not a bond in any chemical sense.
    if (cdk.z_numbers[atom_i] == 1 || cdk.z_numbers[atom_j] == 1 || cdk.z_numbers[atom_k] == 1 ||
        cdk.z_numbers[atom_l] == 1) {
      if (n_wh < dihe_term_with_hydrogen) {
        tmp_dihe_with_h[(5 * n_wh)    ] = 3 * atom_i;
        tmp_dihe_with_h[(5 * n_wh) + 1] = 3 * atom_j;
        tmp_dihe_with_h[(5 * n_wh) + 2] = 3 * atom_k;
        tmp_dihe_with_h[(5 * n_wh) + 3] = 3 * atom_l;
        tmp_dihe_with_h[(5 * n_wh) + 4] = h_parm + 1;
        switch (static_cast<TorsionKind>(vk.dihe_modifiers[pos].w)) {
        case TorsionKind::PROPER:
          break;
        case TorsionKind::IMPROPER:
          tmp_dihe_with_h[(5 * n_wh) + 3] *= -1;
          break;
        case TorsionKind::PROPER_NO_14:
          tmp_dihe_with_h[(5 * n_wh) + 2] *= -1;
          break;
        case TorsionKind::IMPROPER_NO_14:
          tmp_dihe_with_h[(5 * n_wh) + 2] *= -1;
          tmp_dihe_with_h[(5 * n_wh) + 3] *= -1;
          break;
        }
      }
      n_wh++;
    }
    else {
      if (n_woh < dihe_term_without_hydrogen) {
        tmp_dihe_without_h[(5 * n_woh)    ] = 3 * atom_i;
        tmp_dihe_without_h[(5 * n_woh) + 1] = 3 * atom_j;
        tmp_dihe_without_h[(5 * n_woh) + 2] = 3 * atom_k;
        tmp_dihe_without_h[(5 * n_woh) + 3] = 3 * atom_l;
        tmp_dihe_without_h[(5 * n_woh) + 4] = h_parm + 1;
        switch (static_cast<TorsionKind>(vk.dihe_modifiers[pos].w)) {
        case TorsionKind::PROPER:
          break;
        case TorsionKind::IMPROPER:
          tmp_dihe_without_h[(5 * n_woh) + 3] *= -1;
          break;
        case TorsionKind::PROPER_NO_14:
          tmp_dihe_without_h[(5 * n_woh) + 2] *= -1;
          break;
        case TorsionKind::IMPROPER_NO_14:
          tmp_dihe_without_h[(5 * n_woh) + 2] *= -1;
          tmp_dihe_without_h[(5 * n_woh) + 3] *= -1;
          break;
        }
      }
      n_woh++;
    }
  }
  if (n_wh != dihe_term_with_hydrogen || n_woh != dihe_term_without_hydrogen) {
    rtErr("The numbers of dihedral terms counted with and without hydrogen (" +
          std::to_string(n_wh) + " and " + std::to_string(n_woh) + ") do not agree with the "
          "topology's internal records (" + std::to_string(dihe_term_with_hydrogen) + " and " +
          std::to_string(dihe_term_without_hydrogen) + ").  File to print: " + output_file + ".",
          "AtomGraph", "printAmberFormat");
  }
  printAmberSection(&foutp, "DIHEDRALS_INC_HYDROGEN", polyNumericVector(tmp_dihe_with_h), 10, 8, 0,
                    NumberFormat::INTEGER, "print the atom and parameter indexing for bond angles "
                    "with hydrogen to an AMBER format topology");
  printAmberSection(&foutp, "DIHEDRALS_WITHOUT_HYDROGEN", polyNumericVector(tmp_dihe_without_h),
                    10, 8, 0, NumberFormat::INTEGER, "print the atom and parameter indexing for "
                    "dihedral terms with hydrogen to an AMBER format topology");
  printAmberSection(&foutp, "EXCLUDED_ATOMS_LIST", polyNumericVector(tmp_exclusion_list),
                    10, 8, 0, NumberFormat::INTEGER, "print the excluded atoms list to an AMBER "
                    "format topology");
  if (hbond_10_12_parameter_count > 0) {
    printAmberSection(&foutp, "HBOND_ACOEF", polyNumericVector(hbond_a_values.readHost()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding A "
                      "coefficients to an AMBER format topology");
    printAmberSection(&foutp, "HBOND_BCOEF", polyNumericVector(hbond_b_values.readHost()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding B "
                      "coefficients to an AMBER format topology");
    printAmberSection(&foutp, "HBCUT", polyNumericVector(hbond_cutoffs.readHost()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding cutoff "
                      "parameters to an AMBER format topology");
  }
  else {
    printAmberSection(&foutp, "HBOND_ACOEF", polyNumericVector(std::vector<double>()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding A "
                      "coefficients to an AMBER format topology");
    printAmberSection(&foutp, "HBOND_BCOEF", polyNumericVector(std::vector<double>()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding B "
                      "coefficients to an AMBER format topology");
    printAmberSection(&foutp, "HBCUT", polyNumericVector(std::vector<double>()),
                      5, 16, 8, NumberFormat::SCIENTIFIC, "print the hydrogen bonding cutoff "
                      "parameters to an AMBER format topology");
  }
  printAmberSection(&foutp, "AMBER_ATOM_TYPE", polyNumericVector(atom_types.readHost()),
                    20, 4, 0, NumberFormat::CHAR4, "print the atom type names to an AMBER format "
                    "topology");
  printAmberSection(&foutp, "TREE_CHAIN_CLASSIFICATION",
                    polyNumericVector(tree_symbols.readHost()), 20, 4, 0, NumberFormat::CHAR4,
                    "print the tree chain classifications to an AMBER format topology");
  printAmberSection(&foutp, "JOIN_ARRAY", polyNumericVector(tree_joining_info.readHost()),
                    10, 8, 0, NumberFormat::INTEGER, "print the tree joining information to an "
                    "AMBER format topology");
  printAmberSection(&foutp, "IROTAT", polyNumericVector(last_rotator_info.readHost()),
                    10, 8, 0, NumberFormat::INTEGER, "print the rotator information to an AMBER "
                    "format topology");
  printAmberSection(&foutp, "RADIUS_SET", polyNumericVector(atomic_pb_radii.readHost()),
                    5, 16, 8, NumberFormat::SCIENTIFIC, "print the atomic implicit solvent radii "
                    "to an AMBER format topology");
  printAmberSection(&foutp, "SCREEN", polyNumericVector(gb_screening_factors.readHost()),
                    5, 16, 8, NumberFormat::SCIENTIFIC, "print the implicit solvent screening "
                    "factors to an AMBER format topology");
  printAmberSection(&foutp, "IPOL",
                    polyNumericVector(std::vector<int>(1, static_cast<int>(use_polarization))),
                    1, 8, 0, NumberFormat::INTEGER, "print the polarization setting to an "
                    "AMBER format topology");
  switch (use_polarization) {
  case PolarizationSetting::OFF:
    break;
  case PolarizationSetting::ON:
    printAmberSection(&foutp, "POLARIZABILITY",
                      polyNumericVector(atomic_polarizabilities.readHost()), 5, 16, 8,
                      NumberFormat::SCIENTIFIC, "print the atomic polarizability coefficients to "
                      "an AMBER format topology");
    break;
  }
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::printAmberSection(std::ofstream *foutp, const std::string &section_name,
                                  const std::vector<PolyNumeric> &values,
                                  const int values_per_line, const int value_width,
                                  const int decimal_places, const NumberFormat fmt,
                                  const std::string &task) const {
  std::string header("%FLAG " + section_name + "\n%FORMAT(" + std::to_string(values_per_line));
  switch (fmt) {
  case NumberFormat::SCIENTIFIC:
    header += "E" + std::to_string(value_width) + "." + std::to_string(decimal_places) + ")\n";
    break;
  case NumberFormat::STANDARD_REAL:
    header += "F" + std::to_string(value_width) + "." + std::to_string(decimal_places) + ")\n";
    break;
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    header += "I" + std::to_string(value_width) + ")\n";
    break;
  case NumberFormat::CHAR4:
    header += "a" + std::to_string(value_width) + ")\n";
    break;
  }
  foutp->write(header.data(), header.size());
  if (values.size() > 0) {
    printNumberSeries(foutp, values, values_per_line, value_width, decimal_places, fmt,
                      "printAmberSection", task);
  }
  else {
    foutp->write("\n", 1);
  }
}

} // namespace topology
} // namespace stormm
