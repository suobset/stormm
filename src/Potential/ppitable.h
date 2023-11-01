// -*-c++-*-
#ifndef STORMM_PPI_TABLE_H
#define STORMM_PPI_TABLE_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Math/log_scale_spline.h"
#include "Math/math_enumerators.h"
#include "energy_enumerators.h"
#include "pme_util.h"

namespace stormm {
namespace energy {
  
using constants::ExceptionResponse;
using constants::PrecisionModel;
using symbols::amber_ancient_bioq;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using stmath::BasisFunctions;
using stmath::LogScaleSpline;
using stmath::LogSplineForm;
using stmath::LogSplineTable;
using stmath::TableIndexing;
  
/// \brief A tabulated non-bonded potential, with or without exclusions, to be used in the context
///        of particle-particle, particle-mesh calculations.  The key is to create two tables, one
///        for the non-excluded form of the interaction and the other for the excluded form.  The
///        entries for each table will then be concatenated, such that all non-excluded
///        interactions are contiguous and then all exclude interactions are contiguous.  The
///        offset for accessing an excluded interaction based on an index calculated from a
///        particle-particle distance is stored alongside the tabulated splines in the abstract.
class PPITable {
public:

  /// \brief The constuctor can accept all of the arguments that might be useful for making a
  ///        LosScaleSpline, or a LogScaleSpline itself.  Tables for both the otential and the
  ///        derivative will be computed.
  /// \{
  PPITable(NonbondedTheme theme_in = NonbondedTheme::ELECTROSTATIC,
           PrecisionModel precision_in = PrecisionModel::SINGLE,
           BasisFunctions basis_set_in = BasisFunctions::POLYNOMIAL,
           TableIndexing indexing_method_in = TableIndexing::SQUARED_ARG,
           double cutoff_in = default_pme_cutoff, double dsum_tol_in = default_dsum_tol,
           int mantissa_bits_in = 5, double coulomb_in = amber_ancient_bioq,
           double min_spl_compute_range_in = 0.015625, double min_offset_in = 0.125);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
           double cutoff_in = default_pme_cutoff);
  
  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
           const LogScaleSpline<T4> &spl_c);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b);

  template <typename T4>
  PPITable(const LogScaleSpline<T4> &spl_a);
  /// \}
  
  /// \brief The presence of POINTER-kind Hybrid objects implies pointer repairs that require the
  ///        copy and move constructors as well as assignment operators to be spelled out.
  ///
  /// \param original  THe original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  PPITable(const PPITable &original);
  PPITable(PPITable &&original);
  PPITable& operator=(const PPITable &original);
  PPITable& operator=(PPITable &&original);
  /// \}

  /// \brief Get the non-bonded potential described by an object of this class.
  NonbondedTheme getTheme() const;

  /// \brief Get the precision model for this table.
  PrecisionModel getPrecision() const;

  /// \brief Get the cutoff on particle-particle interactions.
  double getCutoff() const;

  /// \brief Get the maximum range of the spline.
  double getMaximumRange() const;

  /// \brief Get the direct sum tolerance.
  double getDirectSumTolerance() const;

  /// \brief Get the number of bits of the mantissa used for table indexing.
  int getBitStride() const;

  /// \brief Get the Ewald coefficient used to perform the switching between short- and
  ///        long-ranged potentials.
  double getEwaldCoefficient() const;

  /// \brief Get the Gaussian RMS sigma parameter (one half the inverse of the Ewald coefficient)
  ///        used to perform the switching between short- and long-ranged potentials.
  double getGaussianWidth() const;
  
private:

  NonbondedTheme theme;           ///< The type of non-bonded potential encoded in the tabulated
                                  ///<   coefficients.
  PrecisionModel precision;       ///< Precision in which this object expresses its table
                                  ///<   coefficients. This indicates whether coeffs for sp_coeffs
                                  ///<   holds any data.
  BasisFunctions basis_set;       ///< The type of basis functions to use in spline approximations
  TableIndexing indexing_method;  ///< The indexing method that converts the argument of the
                                  ///<   underlying function into a spline table element
  double cutoff;                  ///< Maximum inter-particle distance for which interactions are
                                  ///<   valid
  double max_range;               ///< The maximum interparticle distance for which the table can
                                  ///<   produce a result
  double dsum_tol;                ///< The direct sum tolerance, a proportion of each interaction
                                  ///<   that is discarded beginning at the particle-particle pair
                                  ///<   cutoff
  double ew_coeff;                ///< The Ewald coefficient that governs the potential splitting
                                  ///<   between short- and long-ranged components
  int mantissa_bits;              ///< Number of mantissa bits used to index into various tables
  double coulomb;                 ///< The object's take on Coulomb's constant
  int exclusion_offset;           ///< The first index of the table of excluded interactions.  The
                                  ///<   series of non-excluded interactions appears first in each
                                  ///<   table, whether of potentials or derivatives.

  // The main arrays, in full (double) precision
  Hybrid<double> energy;     ///< Double-precision coefficients for the energy, without exclusions
  Hybrid<double> force;      ///< Double-precision coefficients for the force, without exclusions

  /// \brief Double-precision coefficients for the energy, with excluded 1:4 interactions
  Hybrid<double> energy_with_exclusions;

  /// \brief Double-precision coefficients for the force, with excluded 1:4 interactions
  Hybrid<double> force_with_exclusions;

  // The main arrays, in single-precision representations
  Hybrid<float> sp_energy;   ///< Single-precision coefficients for the energy
  Hybrid<float> sp_force;    ///< Single-precision coefficients for the force between particles

  /// \brief Single-precision coefficients for the energy, with excluded 1:4 interactions
  Hybrid<double> sp_energy_with_exclusions;

  /// \brief Single-precision coefficients for the force, with excluded 1:4 interactions
  Hybrid<double> sp_force_with_exclusions;

  // Storage for all coefficients (these are ARRAY-kind Hybrids, unlike the POINTER-kind Hybrids
  // above)
  Hybrid<double> coeffs;     ///< Coefficients of the tables computed in double precision.  This
                             ///<   is an ARRAY-kind Hybrid targeted by the  double-precision
                             ///<   POINTER-kind Hybrid arrays above.
  Hybrid<float> sp_coeffs;   ///< Coefficients of the tables computed in single precision

  /// \brief Find a potential function, without exclusions, among the spline tables provided.
  ///
  /// \param spl_a  The first spline table provided
  /// \param spl_b  The second spline table provided
  /// \param spl_c  The third spline table provided
  /// \param spl_d  The fourth spline table provided
  template <typename T4>
  const LogScaleSpline<T4>& findNonExclPotential(const LogScaleSpline<T4> &spl_a,
                                                 const LogScaleSpline<T4> &spl_b,
                                                 const LogScaleSpline<T4> &spl_c,
                                                 const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a potential function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findExclPotential(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b,
                                              const LogScaleSpline<T4> &spl_c,
                                              const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, without exclusions, among the spline tables provided.
  ///        Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findNonExclForce(const LogScaleSpline<T4> &spl_a,
                                             const LogScaleSpline<T4> &spl_b,
                                             const LogScaleSpline<T4> &spl_c,
                                             const LogScaleSpline<T4> &spl_d) const;

  /// \brief Find a force function, with non-bonded 1:4 exclusions, among the spline tables
  ///        provided.  Descriptions of input parameters follow from findNonExclPotential() above.
  template <typename T4>
  const LogScaleSpline<T4>& findExclForce(const LogScaleSpline<T4> &spl_a,
                                          const LogScaleSpline<T4> &spl_b,
                                          const LogScaleSpline<T4> &spl_c,
                                          const LogScaleSpline<T4> &spl_d) const;

  /// \brief Check a bit in a mask to signify the presence of a particular, necessary function.
  ///
  /// \param spl_x  A spline containing one function relevant to the PPITable
  template <typename T4> uint checkPriority(const LogScaleSpline<T4> &spl_x) const;
  
  /// \brief Determine the highest-priority missing functional form.
  ///
  /// \param holdings  A bitmask checked 1 for the availability of relevant functions: 0x1 for
  ///                  the potential, 0x2 for the excluded potential, 0x4 for the force, 0x8 for
  ///                  the excluded force
  template <typename T4> LogSplineForm findMissingForm(const uint holdings) const;
  
  /// \brief Construct the necessary splines for this aggregate table, based on one to three
  ///        other splines which specify the type.  The priority system is first to find the
  ///        non-excluded energy, then the excluded energy, then the non-excluded force, and
  ///        finally the excluded force (between particles).  Various overloads work given up to
  ///        three available logarithmic spline tables.
  ///
  /// \param spl_a  The first of up to three known spline tables
  /// \param spl_b  The second of up to three known spline tables
  /// \param spl_c  The third and last known spline table
  /// \{
  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                      const LogScaleSpline<T4> &spl_b,
                                      const LogScaleSpline<T4> &spl_c) const;

  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a,
                                       const LogScaleSpline<T4> &spl_b) const;

  template <typename T4>
  LogScaleSpline<T4> getTablePriority(const LogScaleSpline<T4> &spl_a) const;
  /// \}
};

} // namespace energy
} // namespace stormm

#endif
