// -*-c++-*-
#ifndef STORMM_IMPLICIT_SOLVENT_WORKSPACE_H
#define STORMM_IMPLICIT_SOLVENT_WORKSPACE_H

#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"

namespace stormm {
namespace synthesis {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::PrecisionModel;

/// \brief A simple abstract for the implicit solvent workspace.  There are not readers and writers
///        as the only useful application involves this workspace being writeable.
struct ISWorkspaceKit {

  /// \brief The constructor takes a straight list of arguments.
  ISWorkspaceKit(int fp_bits_in, llint* psi_in, int* psi_overflow_in, llint* sum_deijda_in,
                 int* sum_deijda_overflow_in);

  /// \brief The copy and move constructors are defaulted, assignment operators implicitly deleted
  ///        due to const members in this object.
  ///
  /// \param original  The object to copy or move
  /// \{
  ISWorkspaceKit(const ISWorkspaceKit &original) = default;
  ISWorkspaceKit(ISWorkspaceKit &&original) = default;
  /// \}

  const int fp_bits;           ///< Fixed-precision bits after the decimal
  llint* psi;                  ///< 
  int*   psi_overflow;
  llint* sum_deijda;
  int*   sum_deijda_overflow;
};
  
/// \brief A small collection of arrays to manage temporary accumulators for computing Born radii
///        and their derivatives.
class ImplicitSolventWorkspace {
public:
  
  /// \brief The object contains arrays to store Born radii and derivatives to very high precision.
  ///
  /// Overloaded:
  ///   - Accept an explicit bit count for storing Born radii and derivative information in
  ///     fixed-precision format.
  ///   - Accept a precision model and apply an automatic scaling factor for the Born radii and
  ///     derivatives.
  ///
  /// \param atom_starts  Starting positions of the atoms for each system in the coordinate arrays
  ///                     of the related AtomGraphSynthesis and PhaseSpaceSynthesis
  /// \param atom_counts  Atom counts in each system of the related syntheses
  /// \param bit_count    The number of bits to keep after the decimal (units of Angstroms for the
  ///                     radii, or kcal/mol-Angstroms per change in radius)
  /// \param prec         The precision model to prepare for (automates selection of bit_count)
  /// \{
  ImplicitSolventWorkspace(const Hybrid<int> &atom_starts, const Hybrid<int> &atom_counts,
                           int bit_count);

  ImplicitSolventWorkspace(const Hybrid<int> &atom_starts, const Hybrid<int> &atom_counts,
                           PrecisionModel prec);
  /// \}

  /// \brief Copy and move constructors, as well as the respective assignment operators, take the
  ///        usual form for an object with no POINTER-kind Hybrid objects or other pointers to
  ///        repair.
  ///
  /// \param original  The object to copy or move
  /// \param other     Right-hand side object for assignment operators
  /// \{
  ImplicitSolventWorkspace(const ImplicitSolventWorkspace &original) = default;
  ImplicitSolventWorkspace(ImplicitSolventWorkspace &&original) = default;
  ImplicitSolventWorkspace& operator=(const ImplicitSolventWorkspace &original) = default;
  ImplicitSolventWorkspace& operator=(ImplicitSolventWorkspace &&original) = default;
  /// \}
  
  /// \brief Get the precision in which Born radii and sums of the per-atom energy derivatives
  ///        are stored.
  int getFixedPrecisionBits() const;

  /// \brief Get the abstract, containing pointers to data on the host or device.
  ///
  /// \param tier  Level at which to retrieve pointers (the CPU host, or GPU device)
  ISWorkspaceKit data(HybridTargetLevel tier = HybridTargetLevel::HOST);

#ifdef STORMM_USE_HPC
  /// \brief Upload the object's contents from the host to the GPU device.  This may be useful in
  ///        debugging scenarios.
  void upload();

  /// \brief Download the object's contents from the GPU device to the host.  This is also useful
  ///        in debugging scenarios.
  void download();
#endif

private:
  int fp_bits;                      ///< Fixed-precision bits after the decimal
  Hybrid<llint> psi;                ///< The accumulates psi values for ecah atom, which are then
                                    ///<   transformed, in place, to become the effective GB radii
  Hybrid<int> psi_overflow;         ///< Overflow accumulators for psi
  Hybrid<llint> sum_deijda;         ///< Sums of the energy derivative for each atom i making pair
                                    ///<   contributions with all other atoms j with respect to an
                                    ///<   order parameter
  Hybrid<int> sum_deijda_overflow;  ///< Overflow accumulators for aum_deijda
  Hybrid<llint> llint_data;         ///< ARRAY-kind Hybrid targeted by the preceding POITNER-kind
                                    ///<   objects to store all long long integer data
  Hybrid<llint> int_data;           ///< ARRAY-kind Hybrid targeted by the preceding POITNER-kind
                                    ///<   objects to store all overflow data
};

} // namespace synthesis
} // namespace stormm

#endif
