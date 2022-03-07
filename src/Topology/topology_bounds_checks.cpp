#include <string>
#include "Reporting/error_format.h"
#include "topology_bounds_checks.h"

namespace omni {
namespace topology {

//-------------------------------------------------------------------------------------------------
void atomValidityCheck(const int index, const int max_atoms, const char* class_caller,
                       const char* method_caller) {
  if (index >= max_atoms) {
    rtErr("Atom index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_atoms) + " atoms.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void atomValidityCheck(const int low_index, const int high_index, const int max_atoms,
                       const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_atoms || low_index >= high_index) {
    rtErr("Atom range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_atoms) + " atoms.",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void residueValidityCheck(const int index, const int max_residues, const char* class_caller,
                          const char* method_caller) {
  if (index >= max_residues) {
    rtErr("Residue index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_residues) + " residues.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void residueValidityCheck(const int low_index, const int high_index, const int max_residues,
                          const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_residues || low_index >= high_index) {
    rtErr("Residue range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_residues) + " residues.",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void moleculeValidityCheck(const int index, const int max_molecules, const char* class_caller,
                           const char* method_caller) {
  if (index >= max_molecules) {
    rtErr("Molecule index " + std::to_string(index) + " is invalid for a topology with " +
          std::to_string(max_molecules) + " molecules.", class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void moleculeValidityCheck(const int low_index, const int high_index, const int max_molecules,
                           const char* class_caller, const char* method_caller) {
  if (low_index < 0 || high_index > max_molecules || low_index >= high_index) {
    rtErr("Molecule range " + std::to_string(low_index) + " - " + std::to_string(high_index) +
          " is invalid for a topology with " + std::to_string(max_molecules) + " molecules.",
          class_caller, method_caller);
  }
}

}
}
