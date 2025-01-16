// -*-c++-*-
#ifndef STORMM_AG_REFINEMENT_H
#define STORMM_AG_REFINEMENT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/series_ops.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Structure/structure_enumerators.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

using constants::ExceptionResponse;
using parse::WildCardKind;
using stmath::minValue;
using stmath::incrementingSeries;
using structure::ApplyConstraints;

/// \brief A class for tracking the correspondence of two parameter sets.  The parameter sets
///        must be defined as a series of arrays with a 1:1 index mapping among them (i.e. contents
///        of the kth indices of the first, second, and third arrays define the kth parameter set).
template <typename T> class ParameterUnion {
public:
  
  /// \brief The constructor takes the parameters for each of two parameter sets one array at a
  ///        time.
  ///
  /// Overloaded:
  ///   - Provide the inputs as C-style arrays
  ///   - Provide the inputs as Standard Template Library vectors
  ///
  /// \param comp_xa       The first component of parameters in the first set
  /// \param comp_ya       The first (required) component of parameters in the first set
  /// \param comp_za       The third (optional) component of parameters in the first set
  /// \param aparam_count  The number of parameters in the first set
  /// \param comp_xb       The first component of parameters in the second set
  /// \param comp_yb       The second (required) component of parameters in the second set
  /// \param comp_zb       The third (optional) component of parameters in the second set
  /// \param bparam_count  The number of parameters in the second set
  /// \param match_tol     Tolernace by which to declare two parameters match
  /// \{
  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za, const T* comp_wa,
                          const T* comp_va, int aparam_count, const T* comp_xb, const T* comp_yb,
                          const T* comp_zb, const T* comp_wb, const T* comp_vb, int bparam_count,
                          double match_tol = constants::small);

  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                          int aparam_count_in, const T* comp_xb, const T* comp_yb,
                          const T* comp_zb, int bparam_count_in,
                          double match_tol = constants::small);

  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, int aparam_count_in,
                          const T* comp_xb, const T* comp_yb, int bparam_count_in,
                          double match_tol = constants::small);

  explicit ParameterUnion(const T* comp_xa, int aparam_count_in, const T* comp_xb,
                          int bparam_count_in, double match_tol = constants::small);

  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za, const T* comp_wa,
                          const T* comp_va, int aparam_count_in);

  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, const T* comp_za,
                          int aparam_count_in);

  explicit ParameterUnion(const T* comp_xa, const T* comp_ya, int aparam_count_in);

  explicit ParameterUnion(const T* comp_xa, int aparam_count_in);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          const std::vector<T> &comp_za, const std::vector<T> &comp_wa,
                          const std::vector<T> &comp_va, const std::vector<T> &comp_xb,
                          const std::vector<T> &comp_yb, const std::vector<T> &comp_zb,
                          const std::vector<T> &comp_wb, const std::vector<T> &comp_vb,
                          double match_tol = constants::small);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          const std::vector<T> &comp_za, const std::vector<T> &comp_xb,
                          const std::vector<T> &comp_yb, const std::vector<T> &comp_zb,
                          double match_tol = constants::small);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          const std::vector<T> &comp_xb, const std::vector<T> &comp_yb,
                          double match_tol = constants::small);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          double match_tol);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          const std::vector<T> &comp_za, const std::vector<T> &comp_wa,
                          const std::vector<T> &comp_va);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya,
                          const std::vector<T> &comp_za);

  explicit ParameterUnion(const std::vector<T> &comp_xa, const std::vector<T> &comp_ya);

  explicit ParameterUnion(const std::vector<T> &comp_xa);
  /// \}

  /// \brief The default copy and move constructors as well as assignment operators are all valid,
  ///        as there are only Standard Template Library components in the class.
  ///
  /// \param original  The original object to copy or move
  /// \param other     An existing object placed on the right-hand side of the assignment statement
  /// \{
  ParameterUnion(const ParameterUnion &original) = default;
  ParameterUnion(ParameterUnion &&original) = default;
  ParameterUnion& operator=(const ParameterUnion &original) = default;
  ParameterUnion& operator=(ParameterUnion &&original) = default;
  /// \}

  /// \brief Get the number of parameters in one of the input sets.
  ///
  /// \param set_index  The index of the set of interest.  There is no other description, such as
  ///                   which topology the set came from, stored within the object.  This object
  ///                   serves some constructors of the AtomGraph class, and therefore cannot
  ///                   include things for which it is a dependency.
  int getInputSetParameterCount(int set_index) const;
  
  /// \brief Get the number of unique parameter sets.
  int getUniqueParameterCount() const;

  /// \brief Get the union of one of the parameter components.
  ///
  /// \param comp_idx  Index of the component to retrieve.  Acceptable values include 0, 1, or 2.
  std::vector<T> getUnion(int comp_idx) const;

  /// \brief Get the correspondence between the input parameter sets and the consensus tables
  ///        within the object.
  ///
  /// Overloaded:
  ///   - Provide the input set index (alone) to obtain a list of correspondences for every
  ///     parameter in the input set
  ///   - Provide the input set index as well as the index of a parameter within that set to
  ///     obtain a specific parameter's place in the consensus tables
  ///
  /// \param set_index    Index of the input parameter set
  /// \param param_index  The index of a particular parameter within the set
  /// \{
  const std::vector<int>& getCorrespondence(int set_index) const;
  int getCorrespondence(int set_index, int param_index) const;
  /// \}

  /// \brief Add a set of parameters to the existing union.
  ///
  /// Overloaded:
  ///   - Provide C-style arrays with a trusted length
  ///   - Provide Standard Template Library vectors for the input
  ///   
  ///
  /// \param comp_x  The first component of the parameter set to include
  /// \param comp_y  The (required) second component of the parameter set to include
  /// \param comp_z  The (optional) third component of the parameter set to include
  /// \param nparm   The trusted length of comp_x, comp_y, and comp_z
  /// \param tol     The tolerance for declaring a new parameter to match an existing one
  /// \{
  void addSet(const T* comp_x, const T* comp_y, const T* comp_z, const T* comp_w, const T* comp_v,
              int nparm, double tol = constants::small);

  void addSet(const T* comp_x, const T* comp_y, const T* comp_z, int nparm,
              double tol = constants::small);

  void addSet(const T* comp_x, const T* comp_y, int nparm, double tol = constants::small);

  void addSet(const T* comp_x, int nparm, double tol = constants::small);

  void addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
              const std::vector<T> &comp_z, const std::vector<T> &comp_w,
              const std::vector<T> &comp_v, double tol = constants::small);

  void addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
              const std::vector<T> &comp_z, double tol = constants::small);

  void addSet(const std::vector<T> &comp_x, const std::vector<T> &comp_y,
              double tol = constants::small);
  /// \}
  
private:
  int total_parameters;              ///< The total number of unique parameters in the union
  int component_count;               ///< The number of components in each parameter
  int set_count;                     ///< The number of input sets processed by the objectcd
  std::vector<T> union_parameter_x;  ///< The union of the first components of each parameter set
  std::vector<T> union_parameter_y;  ///< The union of the second components of each parameter set
  std::vector<T> union_parameter_z;  ///< The union of the optional third components of each
                                     ///<   parameter set
  std::vector<T> union_parameter_w;  ///< The union of the optional fourth components of each
                                     ///<   parameter set
  std::vector<T> union_parameter_v;  ///< The union of the optional fifth components of each
                                     ///<   parameter set
  
  /// Enumerate which parameters of the unified set the parameters of each input set correspond to.
  /// The first parameter set has a 1:1 correspondence with the first part of the unified set, i.e.
  /// parameter index k in the first parameter set is also parameter index k in the unified set.
  std::vector<std::vector<int>> set_to_consensus_map;

  /// \brief Validate the index of a requested input set.
  ///
  /// \param set_index  The index of the set of interest
  /// \param caller     Name of the calling function, for error tracing purposes
  void validateSetIndex(int set_index, const char* caller = nullptr) const;
};

/// \brief A class to merge two sets of CMAP surfaces.  While the AtomGraphSynthesis has a means
///        for doing this, it is more advanced and fits in the context of a different data
///        structure.  Similar ideas will be applied in the limited case of finding the union of
///        two parameter sets.
class CmapSurfaceUnion {
public:
  
  /// \brief Rather than take ValenceKit abstracts and introduce a dependency on
  ///        atomgraph_abstracts.h, the constructor cherry-picks the relevant details from each
  ///        valence interaction abstract of the two topologies.
  ///
  /// Overloaded:
  ///   - Initialize an empty object
  ///   - Provide the data as raw, C-style arrays
  ///   - Provide the data as Standard Template Library vectors
  ///
  /// \param surf_a     Concatenated array of CMAP surface values for the first topology
  /// \param dim_a      Dimensions of the CMAP surfaces for the first topology (all surfaces are
  ///                   assumed to be square grids).  Summing the squares of the elements of dim_a,
  ///                   as is done internally, provides a prefix sum indicating the bounds of each
  ///                   unique CMAP surface.
  /// \param nmap_a     The number of distinct CMAP surfaces in the first topology, the trusted
  ///                   length of dim_a
  /// \param surf_b     Concatenated array of CMAP surface values for the second topology
  /// \param dim_b      Dimensions of the CMAP surfaces for the second topology
  /// \param nmap_b     The number of distinct CMAP surfaces in the second topology
  /// \param match_tol  The tolerance for declaring that the values at two points of different
  ///                   CMAPs are identical.  All points on both surfaces must align to within this
  ///                   tolerance in order for the surfaces to be considered a match.
  /// \{
  CmapSurfaceUnion();

  CmapSurfaceUnion(const double* surf_a, const int* dim_a, int nmap_a, const double* surf_b,
                   const int* dim_b, int nmap_b, double match_tol = constants::small);

  CmapSurfaceUnion(const double* surf_a, const int* dim_a, int nmap_a);

  CmapSurfaceUnion(const std::vector<double> &surf_a, const std::vector<int> &dim_a,
                   const std::vector<double> &surf_b, const std::vector<int> &dim_b,
                   double match_tol = constants::small);

  CmapSurfaceUnion(const std::vector<double> &surf_a, const std::vector<int> &dim_a);
  /// \}

  /// \brief With only Standard Template Library elements for arrays and scalar member variables
  ///        otherwise, the default copy and move constructors as well as the copy and move
  ///        assignment operators will be adequate.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of an assignment statement
  /// \{
  CmapSurfaceUnion(const CmapSurfaceUnion &orig) = default;
  CmapSurfaceUnion(CmapSurfaceUnion &&orig) = default;
  CmapSurfaceUnion& operator=(const CmapSurfaceUnion &orig) = default;
  CmapSurfaceUnion& operator=(CmapSurfaceUnion &&orig) = default;
  /// \}

  /// \brief Get the number of unique CMAP surfaces.
  int getUniqueSurfaceCount() const;

  /// \brief Get the number of contributing topologies.
  int getContributingTopologyCount() const;
  
  /// \brief Get a CMAP surface based on its index in the consensus parameter set.
  ///
  /// \param surf_index  The index of the surface of interest, with a count beginning at zero
  std::vector<double> getCmapSurface(int surf_index) const;

  /// \brief Get a vector of all surfaces in the union.
  const std::vector<double>& getAllSurfaces() const;

  /// \brief Get the dimension of one or all surfaces in the consensus parameter set.  Each CMAP
  ///        surface is assumed to be a square grid.
  ///
  /// Overloaded:
  ///   - Get all surface dimensions by providing no argument
  ///   - Get a specific CMAP's dimension by specifying a parameter index
  ///
  /// \param surf_index  The index of the CMAP surface of interest
  /// \{
  const std::vector<int>& getSurfaceDimensions() const;
  int getSurfaceDimensions(int surf_index) const;
  /// \}

  /// \brief Access the CMAP surface parameter index correspondence for a particular input set
  ///        (input sets will, for most purposes, be specific topologies).
  ///
  /// Overloaded:
  ///   - Provide the CMAP index within the topology of interest to get its index in the consensus
  ///   - Provide only an indication of the topology to get the mapping for all of its CMAP
  ///     parameters onto the consensus table
  ///
  /// \param set_index   Index of the topology of interest
  /// \param surf_index  Index of the surface of interest within its native topology
  /// \{
  const std::vector<int>& getCorrespondence(int set_index) const;
  int getCorrespondence(int set_index, int surf_index) const;
  /// \}
  
  /// \brief Add another set of CMAP terms to the union.
  ///
  /// \param surf_v  The CMAP surface values with which to populate the object, concatenated
  /// \param dims_v  The dimensions of each surface
  /// \param nmap    The number of distinct CMAP surfaces
  /// \param tol     The tolerance for declaring that the values at two points of different CMAPs
  ///                are identical.  All points on both surfaces must align to within this
  ///                tolerance in order for the surfaces to be considered a match.
  void addSet(const double* surf_v, const int* dim_v, int nmap, double tol = constants::small);
  
private:

  int surface_count;                   ///< The number of unique surfaces in the union of two sets
  int set_count;                       ///< The number of sets covered by the consensus
  std::vector<double> surface_values;  ///< Values of the parameters for each CMAP surface,
                                       ///<   concatenated without padding
  std::vector<int> surface_dims;       ///< Dimensions of each surface (all CMAPs are assumed to
                                       ///<   be square)
  std::vector<int> surface_bounds;     ///< Bounds array on surface_values, above, setting the
                                       ///<   start and stop points for each column-major surface
                                       ///<   matrix.  This is a prefix sum over the squares of
                                       ///<   values in surface_dims.

  /// Enumerate which surfaces of the unified list the CMAP surfaces of each set correspond to.
  /// The first parameter set has a 1:1 correspondence with the first part of the unified set, as
  /// as the case in the unions of other parameters.
  std::vector<std::vector<int>> set_to_consensus_map;

  /// \brief Populate the object based exclusively on one topology's set of CMAP terms.  This can
  ///        be executed in preparation for expanding the set with the other topology's set of
  ///        CMAP terms.
  ///
  /// \param surf_v  The CMAP surface values with which to populate the object, concatenated
  /// \param dims_v  The dimensions of each surface
  /// \param nmap   The number of dinstinct CMAP surfaces
  void populateByOneSet(const double* surf_v, const int* dim_v, int nmap);

  /// \brief Validate the index of a CMAP parameter surface in the consensus set.
  ///
  /// \param surf_index  Index of the surface of interest
  /// \param caller      Name of the calling function (for error tracing purposes)
  void validateParameterIndex(int surf_index, const char* caller) const;

  /// \brief Validate the index of a parameter set comprised by the consensus tables.
  ///
  /// \param set_index  Index of the input parameter set (topology) of interest
  /// \param caller     Name of the calling function (for error tracing purposes)
  void validateSetIndex(int set_index, const char* caller) const;
};
  
/// \brief Unguarded struct to assemble the basic bond, angle, and dihedral indexing from an Amber
///        or other topology file.  All Class I force fields will have, or could have, terms like
///        these.  All descriptions follow from the eponymous member variables in an AtomGraph.
class BasicValenceTable {
public:

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each type of basic valence term
  ///
  /// \param natom_in          The number of atoms in the system (for bounds arrays)
  /// \param nbond_in          The number of bonds to prepare for (atom indices will be filled if
  ///                          more information is provided)
  /// \param nangl_in          The number of angles to prepare for (atom indices will be filled if
  ///                          more information is provided)
  /// \param ndihe_in          The number of dihedrals to prepare for (atom indices will be filled
  ///                          if more information is provided)
  /// \{
  BasicValenceTable();
  BasicValenceTable(int natom_in, int nbond_in, int nangl_in, int ndihe_in,
                    const std::vector<int> &bond_i_atoms_in = {},
                    const std::vector<int> &bond_j_atoms_in = {},
                    const std::vector<int> &bond_param_idx_in = {},
                    const std::vector<int> &angl_i_atoms_in = {},
                    const std::vector<int> &angl_j_atoms_in = {},
                    const std::vector<int> &angl_k_atoms_in = {},
                    const std::vector<int> &angl_param_idx_in = {},
                    const std::vector<int> &dihe_i_atoms_in = {},
                    const std::vector<int> &dihe_j_atoms_in = {},
                    const std::vector<int> &dihe_k_atoms_in = {},
                    const std::vector<int> &dihe_l_atoms_in = {},
                    const std::vector<int> &dihe_param_idx_in = {},
                    const std::vector<char4> &bond_mods_in = {},
                    const std::vector<char4> &angl_mods_in = {},
                    const std::vector<char4> &dihe_mods_in = {});
  /// \}

  /// \brief Detail the relevance of each bond or angle in the system.  Bonds involving virtual
  ///        sites are irrelevant to typical molecular mechanics calculations, as are angles (the
  ///        virtual sites are set at particular geometries based on their frame atoms), although
  ///        dihedrals involving virtual sites may have unique potentials that depend on more than
  ///        just the frame atoms.  When geometric constraints are in effect, bonds involving
  ///        hydrogen atoms are irrelevant to typical molecular mechanics calculations.  Angles
  ///        wherein all three atoms are bonded in a rigid triangle (not, for example, an -NH3
  ///        group where the hydrogens can move around and the H-N-H angles can change, but rigid
  ///        waters where the H-O-H angle is fixed due to an H-H "bond") are also irrelevant.
  ///
  /// \param use_shake   Indicate whether rigid geometry constraints will be put into effect for
  ///                    bonds to individual hydrogen atoms, which can alter how many bond and
  ///                    angle terms the topology evaluates
  /// \param use_settle  Indicate whether the rigid geometry of trigonal water molecules (or other
  ///                    three-atom groups which make an isosceles triangle) will be enforced by
  ///                    the SETTLE algorithm.  This is a stronger constraint than merely
  ///                    constraining individual bonds between hydrogens and the one heavy atom.
  void checkBondAngleRelevance(ApplyConstraints use_shake, ApplyConstraints use_settle,
                               const std::vector<int> &atomic_numbers);
  
  /// \brief Populate the atom assignments for this object.  This can only be called once the
  ///        Atom index arrays have been assigned.
  void makeAtomAssignments();

  // The member variables of this class are public.
  int total_bonds;                        ///< Total number of harmonic bond terms in the topology
  int total_angls;                        ///< Total number of harmonic angle terms in the topology
  int total_dihes;                        ///< Total number of cosine-based proper and improper
                                          ///<   torsions in the topology
  std::vector<int> bond_i_atoms;          ///< First atom of each bond term
  std::vector<int> bond_j_atoms;          ///< Second atom of each bond term
  std::vector<int> bond_param_idx;        ///< Parameter index for each bond term
  std::vector<int> angl_i_atoms;          ///< First atom of each harmonic angle term
  std::vector<int> angl_j_atoms;          ///< Second atom of each harmonic angle term
  std::vector<int> angl_k_atoms;          ///< Third atom of each harmonic angle term
  std::vector<int> angl_param_idx;        ///< Parameter index for each harmonic angle term
  std::vector<int> dihe_i_atoms;          ///< First atom in each cosine-based torsion term
  std::vector<int> dihe_j_atoms;          ///< Second atom in each cosine-based torsion term
  std::vector<int> dihe_k_atoms;          ///< Third atom in each cosine-based torsion term
  std::vector<int> dihe_l_atoms;          ///< Fourth atom in each cosine-based torsion term
  std::vector<int> dihe_param_idx;        ///< Parameter index for each coside-based torsion term
  std::vector<char4> bond_mods;           ///< Notes on each harmonic bond term, i.e. "contains
                                          ///<   hydrogen"
  std::vector<char4> angl_mods;           ///< Notes on each harmonic angle term
  std::vector<char4> dihe_mods;           ///< Notes on each cosine-based torsion term, i.e. "is
                                          ///<   improper, has no associated 1-4 interaction"
  std::vector<int> bond_assigned_atoms;   ///< Each bond term is assigned to an atom, and details
                                          ///<   of the bond assignments for the ith atom are held
                                          ///<   between indices bond_assigned_bounds[i] and
                                          ///<   bond_assigned_bounds[i + 1] of the atoms, index,
                                          ///<   and terms arrays.  The atoms array holds the
                                          ///<   topological indices of atoms in the bond, other
                                          ///<   than the obivous atom i to which the bond is
                                          ///<   assigned.
  std::vector<int> bond_assigned_index;   ///< Parameter indices relevant to bond assignments
  std::vector<int> bond_assigned_terms;   ///< Bond term indices relevant to bond assignments
  std::vector<int> bond_assigned_bounds;  ///< Bounds array stating the limits of bond terms
                                          ///<   assigned to each atom (see the longer explanation
                                          ///<   above)
  std::vector<int> angl_assigned_atoms;   ///< See the description above for the corresponding
                                          ///<   member variable about bonds.  This array has two
                                          ///<   entries for every assignment, so multiply the
                                          ///<   indexing and realize that the relevant atoms for
                                          ///<   bonds assigned to atom i are contained in entries
                                          ///<   2 * angl_assigned_bounds[i] to
                                          ///<   2 * angl_assigned_bounds[i + 1].
  std::vector<int> angl_assigned_index;   ///< Parameter indices relevant to angle assignments
  std::vector<int> angl_assigned_terms;   ///< Angle term indices relevant to angle assignments
  std::vector<int> angl_assigned_bounds;  ///< Bounds array stating the limits of angle terms
                                          ///<   assigned to each atom
  std::vector<int> dihe_assigned_atoms;   ///< See the descriptions about for other atoms relevant
                                          ///<   to bond and angle terms.  There are three atoms
                                          ///<   per assignment in this array.
  std::vector<int> dihe_assigned_index;   ///< Parameter indices relevant to torsion assignments
  std::vector<int> dihe_assigned_terms;   ///< Torsion term indices relevant to torsion assignments
  std::vector<int> dihe_assigned_bounds;  ///< Bounds array for the limits of torsion assignments
                                          ///<   to each atom
  std::vector<char4> bond_assigned_mods;  ///< Notes on bonds assigned to each atom (follows the
                                          ///<   ordering set forth in bond_assigned_bounds)
  std::vector<char4> angl_assigned_mods;  ///< Notes on angles assigned to each atom
  std::vector<char4> dihe_assigned_mods;  ///< Notes on torsions assigned to each atom
  std::vector<uint> bond_relevance;       ///< Notes on whether to evaluate each bond term.  Each
                                          ///<   bond gets a bit in one of the unsigned integer
                                          ///<   elements, 1 for "yes, evaluate the bond" and 0 for
                                          ///<   "no, do not evaluate the bond by the standard MM
                                          ///<   evaluation functions."
  std::vector<uint> angl_relevance;       ///< Notes on whether to evaluate each angle term.  Each
                                          ///<   angle gets a bit, similar to bond_relevance,
                                          ///<   above.
};

/// \brief Unguarded struct to assemble special "CHARMM" force field terms read from an Amber or
///        other topology file.  These include Urey-Bradley harmonic angle terms, CHARMM improper
///        harmonic dihedrals, and CMAP terms.  All descriptions follow from the eponymous member
///        variables in an AtomGraph.
class CharmmValenceTable {
public:
  
  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each type of CHARMM valence term
  ///
  /// \param natom_in  The number of atoms in the system (for bounds arrays)
  /// \param nubrd_in  The number of Urey-Bradley harmonic angles terms to prepare for (atom
  ///                  indices will be filled if enough information is provided)
  /// \param nimpr_in  The number of CHARMM improper dihedrals to prepare for (atom indices will
  ///                  be filled if enough information is provided)
  /// \param ncmap_in  The number of CMAP terms to prepare for (atom indices will be filled if
  ///                  enough information is provided)
  /// \{
  CharmmValenceTable();
  CharmmValenceTable(int natom_in, int nubrd_in, int nimpr_in, int ncmap_in,
                     const std::vector<int> &ubrd_i_atoms_in = {},
                     const std::vector<int> &ubrd_k_atoms_in = {},
                     const std::vector<int> &ubrd_param_idx_in = {},
                     const std::vector<int> &impr_i_atoms_in = {},
                     const std::vector<int> &impr_j_atoms_in = {},
                     const std::vector<int> &impr_k_atoms_in = {},
                     const std::vector<int> &impr_l_atoms_in = {},
                     const std::vector<int> &impr_param_idx_in = {},
                     const std::vector<int> &cmap_i_atoms_in = {},
                     const std::vector<int> &cmap_j_atoms_in = {},
                     const std::vector<int> &cmap_k_atoms_in = {},
                     const std::vector<int> &cmap_l_atoms_in = {},
                     const std::vector<int> &cmap_m_atoms_in = {},
                     const std::vector<int> &cmap_param_idx_in = {});
  /// \}
  
  /// \brief Populate the atom assignments for this object.  This can only be called once the
  ///        Atom index arrays have been assigned.
  void makeAtomAssignments();

  // The member variables of this class are public.
  int total_ub_angles;
  int total_impropers;
  int total_cmaps;
  std::vector<int> ubrd_i_atoms;
  std::vector<int> ubrd_k_atoms;
  std::vector<int> ubrd_param_idx;
  std::vector<int> impr_i_atoms;
  std::vector<int> impr_j_atoms;
  std::vector<int> impr_k_atoms;
  std::vector<int> impr_l_atoms;
  std::vector<int> impr_param_idx;
  std::vector<int> cmap_i_atoms;
  std::vector<int> cmap_j_atoms;
  std::vector<int> cmap_k_atoms;
  std::vector<int> cmap_l_atoms;
  std::vector<int> cmap_m_atoms;
  std::vector<int> cmap_param_idx;
  std::vector<int> ubrd_assigned_atoms;
  std::vector<int> ubrd_assigned_index;
  std::vector<int> ubrd_assigned_terms;
  std::vector<int> ubrd_assigned_bounds;
  std::vector<int> impr_assigned_atoms;
  std::vector<int> impr_assigned_index;
  std::vector<int> impr_assigned_terms;
  std::vector<int> impr_assigned_bounds;
  std::vector<int> cmap_assigned_atoms;
  std::vector<int> cmap_assigned_index;
  std::vector<int> cmap_assigned_terms;
  std::vector<int> cmap_assigned_bounds;
};

/// \brief Unguarded struct to stage 1:4 exclusion screening terms.  This will translate a series
///        of unique dihedral parameter sets, containing 1:4 scaling parameters for every dihedral,
///        into collections of unique 1:4 scaling parameters (screening factors) and a list of
///        indices that each diehdral can follow in order to access their factors from within these
///        parameter collections.  Scaling parameters of zero are an obligatory, unique member of
///        these collections.  This arrangement sets up for force computation in which each
///        dihedral references four cached atoms, then a set of cached parameters for its cosine
///        form, then a pair of cached screening terms for the 1:4 interactions of its I and L
///        atoms.
struct AttenuationParameterSet {

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// \param set_count_in
  AttenuationParameterSet(int set_count_in = 0);
  
  int total_14_sets;                          ///< Total number of unique 1:4 scaling factor pairs
  std::vector<int> dihe14_parameter_indices;  ///< Array with length equal to the number of unique
                                              ///<   dihedral parameter sets in the original
                                              ///<   topology.  Each dihedral parameter set finds
                                              ///<   its original screening factors at these
                                              ///<   indices in the following (likely much shorter)
                                              ///<   arrays.
  std::vector<double> elec_screening_factors; ///< Condensed 1:4 electrostatic attenuation factors
  std::vector<double> vdw_screening_factors;  ///< Condensed 1:4 van-der Waals attenuation factors
};

/// \brief Unguarded struct to condense the atom exclusion lists read from an Amber topology.  This
///        will transform the list of exclusion counts into an capped prefix sum (an exclusive
///        prefix sum in indices 0...N-1, or an inclusive prefix sum in indices 1...N).  Any
///        atoms with no actual exclusions will be listed as having zero exclusions, and blank
///        atoms will be removed from the excluded atoms list.
struct CondensedExclusions {

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each virtual site
  ///
  /// \param natom_in        
  /// \param total_exclusions_in
  /// \{
  CondensedExclusions();
  CondensedExclusions(int natom_in, int total_exclusions_in);
  /// \}

  int total_exclusions;
  std::vector<int> atom_excl_bounds;
  std::vector<int> atom_excl_list;
};

/// \brief Unguarded struct to collect information about virtual sites in the topology.
struct VirtualSiteTable {

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each virtual site 
  ///
  /// \param natom_in     The number of atoms in the system (for an array indicating which atom is
  ///                     which virtual site)
  /// \param vs_count_in  The number of virtual sites in the system (for bounds arrays)
  /// \{
  VirtualSiteTable();

  VirtualSiteTable(int natom_in, int vs_count_in);

  VirtualSiteTable(int natom_in, const std::vector<int> &vs_atoms_in,
                   const std::vector<int> &frame_types_in, const std::vector<int> &frame1_atoms_in,
                   const std::vector<int> &frame2_atoms_in,
                   const std::vector<int> &frame3_atoms_in,
                   const std::vector<int> &frame4_atoms_in,const std::vector<int> &param_idx_in,
                   const std::vector<double> &frame_dim1_in,
                   const std::vector<double> &frame_dim2_in,
                   const std::vector<double> &frame_dim3_in);
/// \}

  int vs_count;                     ///< Number of virtual sites found, a check on the topology
  std::vector<int> vs_numbers;      ///< An array covering all atoms in the topology and stating
                                    ///<   which virtual site they correspond to, containing -1
                                    ///<   entries if the atom is not a virtual site.  This array
                                    ///<   is critical for rapidly tracing the virtual site of
                                    ///<   interest given a topological index, such as finding the
                                    ///<   parent atom of the ith atom in a topology, knowing that
                                    ///<   the ith atom is a virtual site.
  std::vector<int> vs_atoms;        ///< Atoms identified as virtual sites (by having zero mass).
                                    ///<   This compact list includes the topological indices only
                                    ///<   of atoms that hav ebeen identified as virtual sites,
                                    ///<   i.e. 3, 7, 11, ..., 4N + 3 in a box of TIP4P-Ew water.
  std::vector<int> frame_types;     ///< Frame types for each virtual site
  std::vector<int> frame1_atoms;    ///< Parent atoms
  std::vector<int> frame2_atoms;    ///< Second frame atoms
  std::vector<int> frame3_atoms;    ///< Third frame atoms (for some sites)
  std::vector<int> frame4_atoms;    ///< Fourth frame atoms (for the most complex frame type)
  std::vector<int> param_idx;       ///< Parameter indices indicating the set of frame dimensions
                                    ///<   controlling this particular site (this is part of a
                                    ///<   broader effort to cut down on memory bandwidth
                                    ///<   requirements)
  std::vector<double> frame_dim1;   ///< First frame dimension, i.e. distance from VS to parent
  std::vector<double> frame_dim2;   ///< Second frame dimension
  std::vector<double> frame_dim3;   ///< Third frame dimension
};

/// \brief Unguarded struct to help in the construction of a AtomGraph's categorized non-bonded
///        exclusion lists.
struct Map1234 {

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each virtual site
  ///
  /// \param natom_in       The number of atoms in the system (for bounds arrays)
  /// \param nb11_count_in  The (known) number of 1:1 non-bonded exclusions
  /// \param nb12_count_in  The (known) number of 1:2 non-bonded exclusions
  /// \param nb13_count_in  The (known) number of 1:3 non-bonded exclusions
  /// \param nb14_count_in  The (known) number of 1:4 non-bonded exclusions
  /// \{
  Map1234();
  Map1234(int natom_in, int nb11_count_in, int nb12_count_in, int nb13_count_in,
          int nb14_count_in);
  /// \}
  
  std::vector<int> nb11_excl_bounds;  ///< 1:1 exclusion bounds list for each atom (needed for
                                      ///<   virtual sites and atoms with virtual sites)
  std::vector<int> nb11_excl_list;    ///< 1:1 exclusions list (double-counts all exclusions)
  std::vector<int> nb12_excl_bounds;  ///< Bond exclusion bounds list for each atom
  std::vector<int> nb12_excl_list;    ///< Bonded exclusions (double-counts all exclusions)
  std::vector<int> nb13_excl_bounds;  ///< Bond angle exclusion bounds list for each atom
  std::vector<int> nb13_excl_list;    ///< Bond angle exclusions (double-counts all exclusions)
  std::vector<int> nb14_excl_bounds;  ///< 1:4 exclusions bounds list for each atom
  std::vector<int> nb14_excl_list;    ///< 1:4 exclusions (double-counts all exclusions)
};

/// \brief Unguarded struct to hold refined quantities associated with CMAP objects
struct CmapAccessories {

  /// \brief The constructor simply allocates memory, if dimensions are available.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Create an object with pre-allocated memory for each CMAP surface
  ///
  /// \param cmap_dimensions_in  Dimensions of each CMAP surface
  /// \{
  CmapAccessories();
  CmapAccessories(const std::vector<int> &cmap_dimensions_in);
  /// \}
  
  std::vector<double> phi_derivatives;      ///< First derivatives along the CMAP's first axis
  std::vector<double> psi_derivatives;      ///< First derivatives along the CMAP's second axis
  std::vector<double> phi_psi_derivatives;  ///< Cross derivatives at all grid points
  std::vector<int> patch_matrix_bounds;     ///< Bounds of each map's 16-element patch series in
                                            ///<   the patch_matrix_form array
  std::vector<double> patch_matrix_form;    ///< Interlaced, inflated array of CMAP surface values
                                            ///<   and derivatives
                                            ///< The rank 4, square coefficients matrix
                                            ///<   [ a00, a01, ... a33 ] may be computed as a
                                            ///<   matrix multiplication to then obtain the value
                                            ///<   and derivatives of a bicubic spline interpolant
                                            ///<   at an arbitrary point.  Between any four points
                                            ///<   of the grid, however, the coefficients are all
                                            ///<   the same and take up the same amount of space as
                                            ///<   the 16 values used to compute them.  This array
                                            ///<   stores the a00 ... a33 coefficients for every
                                            ///<   grid segment in contiguous stretches of 16
                                            ///<   numbers.
};

/// \brief Encapsulate the six derived parameters for SETTLE computations.  A vector of these
///        structs will be accumulated into a ConstraintTable for any given topology, enumerating
///        each distinct SETTLE-suitable molecule.  In GPU implementations the information will be
///        held by two- and four-tuples.
struct SettleParm {
  double mo;     ///< Mass of the heavy "oxygen" atoms in a SETTLE group
  double mh;     ///< Mass of the light "hydrogen" atoms in a SETTLE group
  double moh;    ///< Combined mass of the heavy and one of the light atoms in a SETTLE group
  double mormt;  ///< Proportional mass of the "oxygen" atoms in a SETTLE group
  double mhrmt;  ///< Proportional mass of "hydrogen" atoms in a SETTLE group
  double ra;     ///< Internal distance measurement (heavy atom - center of mass) of a SETTLE group
  double rb;     ///< Internal distance measurement of a SETTLE group (light atoms - center of mass
                 ///<   along the heavy atom / CoM axis)
  double rc;     ///< Internal distance measurement of a SETTLE group (half the distance between
                 ///<   light atoms)
};

/// \brief Unguarded struct to hold the results of a constraint analysis.  Constraint groups as
///        well as fast, rigid waters will be selected.  All atoms that might participate in
///        constraints will be selected and drawn into groups; whether those groups are applied
///        is a matter for input parameters to decide further on in the program.
struct ConstraintTable {

  /// \brief The constructor handles all allocations based on its input arguments.
  ///
  /// \param atomic_numbers   Atomic numbers of all atoms in the system (0 for virtual sites)
  /// \param atomic_mass      Masses of all atoms, in atomic mass units (Daltons)
  /// \param mol_limits       Starting an ending bounds for reading the contents of each molecule
  ///                         in mol_contents
  /// \param mol_contents     Topological atom indices indicating the contents of each molecule
  /// \param mol_home         Molecules to which each atom belongs, the kth element of the array
  ///                         being the molecule home of atom k in the topology
  /// \param all_nb_excl      Connectivity information can be inferred by examining the 1:2
  ///                         exclusion information, which is reflexive within this object (if
  ///                         atom j excludes atom i, then atom i also excludes atom j)
  /// \param bond_equilibria  Equilibrium lengths of all harmonic bond parameters
  /// \param angl_equilibria  Equilibrium lengths of all harmonic angle parameters
  ConstraintTable(const std::vector<int> &atomic_numbers, const std::vector<double> &atomic_masses,
                  const std::vector<int> &mol_limits, const std::vector<int> &mol_contents,
                  const std::vector<int> &mol_home, const BasicValenceTable &bvt,
                  const Map1234 &all_nb_excl, const std::vector<double> &bond_equilibria,
                  const std::vector<double> &angl_equilibria);

  int cnst_group_count;                   ///< Number of constraint groups found
  int settle_group_count;                 ///< Number of SETTLE-suitable fast rigid waters found
  int cnst_parameter_count;               ///< Number of unique constraint group parameters sets.
                                          ///<   A constraint group parameter set consists of the
                                          ///<   masses of every atom and the bond lengths between
                                          ///<   the central atom and each of its subsidiary
                                          ///<   hydrogens, in order.  Every constraint group of
                                          ///<   unique size, i.e. one central atom with three
                                          ///<   subsidiaries, will have at least one unique
                                          ///<   parameter set.
  int settle_parameter_count;             ///< Number of unique SETTLE parameter sets (a SETTLE
                                          ///<   parameter set consists of the masses of "oxygen"
                                          ///<   and "hydrogen" atoms plus the bond lengths between
                                          ///<   any pair of the three atoms--five numbers in all)
  std::vector<int> cnst_group_list;       ///< Concatenated list of atoms in all constraint groups
  std::vector<int> cnst_group_bounds;     ///< Bounds array for cnst_group_list
  std::vector<int> cnst_group_param_idx;  ///< Indices into the list of constraint group parameter
                                          ///<   sets for each constraint group
  std::vector<int> settle_ox_atoms;       ///< Topological indices of "oxygen" atoms in SETTLE
                                          ///<   groups (this can also include other residues that
                                          ///<   match the "one heavy atom, two symmetric light
                                          ///<   atoms" shape requirements of the SETTLE algorithm)
  std::vector<int> settle_h1_atoms;       ///< Topological indices of the first "hydrogen" atoms in
                                          ///<   SETTLE groups
  std::vector<int> settle_h2_atoms;       ///< Topological indices of the second "hydrogen" atoms
                                          ///<   in SETTLE groups
  std::vector<int> settle_param_idx;      ///< Indices into the lists of SETTLE constraint group
                                          ///<   parameter sets

  /// Parameters for constraint groups--each stretch of the array is one parameter set, with
  /// demarcations provided by the following bounds array.  Each tuple encodes an inverse atom
  /// mass in the x member and a target bond length in the y member.  The order of these tuples
  /// follows the order of atoms in each constraint group: central heavy atom, followed by a
  /// series of bound hydrogen atoms.  The bond length assigned to the central heavy atom is
  /// meaningless.
  std::vector<double2> cnst_parameter_list;

  /// Bounds array for constraint_group_parameters
  std::vector<int> cnst_parameter_bounds;

  // Settle parameter sets
  std::vector<SettleParm> settle_measurements;
};

/// \brief Unguarded struct to handle 1:4 attenuated interactions that are not explicitly
///        handled by the short-ranged non-bond scaling parameters of some dihedral linking the
///        atoms
struct AttenuatedPair {
  int atom_i;           ///< The first atom in the pair (could have been atom I in some dihedral)
  int atom_j;           ///< The first atom in the pair (could have been atom L in some dihedral)
  double lj_scaling;    ///< Scaling factor for van-der Waals interactions
  double elec_scaling;  ///< Scaling factor for electrostatic interactions
};
  
/// \brief Smooth a set of charges subject to the constraint that, when expressed in a particular
///        precision, they will add to the correct integer value.  Ignore if the initial sum is
///        too far from any integral value.  Charges are assumed to be in internal units (atomic
///        units, multiples of the proton charge).
///
/// \param q                      Vector of charges to smooth (modified and returned)
/// \param tmp_charge_parameters  Charge parameters that the system will use (filled and returned)
/// \param tmp_charge_indices     Atom indices into the charge parameters array (filled and
///                               returned)
/// \param q_param_count          The number of unique charge parameters (determined and returned)
/// \param rounding_tol           The maximum tolerance at which to initiate charge rounding.  If
///                               a system has a net charge of 0.0002e, that looks like a roundoff
///                               error in some interpretation of a force field term, so smooth
///                               that over the existing atoms.  If the system has a charge of
///                               -0.9997e, then the system seems to be built with a net charge of
///                               -1, so sprinkle charge increments among the charges to make that
///                               come out to the value that was intended.  If the system has a
///                               value of 0.5837e, then something fundamental may be wrong, but
///                               either way it's not something that STORMM is prepared to redress.
/// \param charge_discretization  The increment with which to discretize charges
/// \param file_name              Optional name of the topology file whose charges are being
///                               smoothed.  Leave blank if this is not about a particular
///                               topology.
void smoothCharges(std::vector<double> *q, std::vector<double> *tmp_charge_parameters,
                   std::vector<int> *tmp_charge_indices, int *charge_parameter_count,
                   double rounding_tol, double precision,
                   const std::string &filename = std::string(""));
  
/// \brief Expand the Lennard-Jones tables based on the non-bonded parameter indices, switching to
///        larger tables that store nearly twice the original data but can be accessed without
///        referencing the non-bonded parameter indices (only the Lennard-Jones indices of each
///        atom are then needed).
///
/// \param lj_a_values     Lennard-Jones A coefficients (expanded and returned)
/// \param lj_b_values     Lennard-Jones B coefficients (expanded and returned)
/// \param lj_c_values     Lennard-Jones C coefficients (expanded and returned)
/// \param lj_14_a_values  Lennard-Jones 1:4 A coefficients (expanded and returned)
/// \param lj_14_b_values  Lennard-Jones 1:4 B coefficients (expanded and returned)
/// \param lj_14_c_values  Lennard-Jones 1:4 C coefficients (expanded and returned)
/// \param hb_a_values     Hydrogen bonding 10-12 A coefficients (expanded and returned)
/// \param hb_b_values     Hydrogen bonding 10-12 B coefficients (expanded and returned)
/// \param hb_c_values     Hydrogen bonding cutoff limits (expanded and returned)
/// \param n_lj_types      The number of Lennard-Jones types
/// \param nb_param_index  The non-bonded parameter indices (this table can be accessed by rows
///                        and columns indexed by a pair of atoms' Lennard-Jones parameter indices)
void expandLennardJonesTables(std::vector<double> *lj_a_values, std::vector<double> *lj_b_values,
                              std::vector<double> *lj_c_values,
                              std::vector<double> *lj_14_a_values,
                              std::vector<double> *lj_14_b_values,
                              std::vector<double> *lj_14_c_values,
                              std::vector<double> *hb_a_values, std::vector<double> *hb_b_values,
                              std::vector<double> *hb_c_values, int n_lj_types,
                              const std::vector<int> &nb_param_index);
  
/// \brief Process the exclusions recorded from an Amber topology.  See the CondensedExclusion
///        helper struct for a description of how the file format can contain confusing information
///        that is best removed.  The results from this function can be stashed in the appropriate
///        segments of the topology's Hybrid data arrays.
///
/// \param raw_counts      The exclusion counts recorded in the Amber prmtop file
/// \param raw_exclusions  The exclusions recorded in the Amber prmtop file
/// \param file_name       The name of the Amber prmtop file (for error reporting)
CondensedExclusions processExclusions(const std::vector<int> &raw_counts,
                                      const std::vector<int> &raw_exclusions,
                                      const std::string &file_name);

/// \brief Back-calculate the exclusions for each atom in a non-redundant format that would be more
///        suitable for printing into a topology file.  The indexing of atoms continues to begin at
///        zero, and exclusions j of atom i will be such that j > i.  All of 1:1, 1:2, 1:3, and 1:4
///        exclusions will be included.  In addition, if an atom has no forward exclusions, it will
///        be recorded as having one exclusion to the zeroth atom (index -1).
///
/// \param nb_excl  Tables of 1:1, 1:2, 1:3, and 1:4 exclusions
CondensedExclusions calculatePrmtopExclusions(const Map1234 nb_excl);

/// \brief The total number of exclusions reported in an Amber prmtop will have explicit zeros
///        listed for atoms that exclude no atoms of higher topological index.  This function will
///        identify all such atoms based on the complete list of exclusions.
///
/// \param nb_excl  Tables of 1:1, 1:2, 1:3, and 1:4 exclusions
int countPrmtopZeroExclusions(const Map1234 nb_excl);

/// \brief Process indexing for "basic" family valence parameters from an Amber topology.  Note the
///        meaning of "basic" family parameters: harmonic bonds based on the distance between atoms
///        I and J; harmonic angles for atoms I, J, and K; and cosine-based dihedrals for atoms I,
///        J, K, and L.  An Amber topology encodes this indexing as two lists of N+1 integer tuples
///        for Nth-order terms: one for terms containing hydrogen and one for terms not containing
///        hydrogen.  The numbers in these tuples reflect Fortran indexing, with some degree of
///        pre-computation that is hereby removed in favor of C / C++ format indices for atoms and
///        parameters.
///
/// \param atom_count          Total number of atoms in the system
/// \param tmp_bond_atoms_h    Lists of bond atoms (and parameters) containing hydrogen
/// \param tmp_bond_atoms_noh  Lists of bond atoms (and parameters) containing only heavy atoms
/// \param tmp_angl_atoms_h    Lists of angle atoms (and parameters) containing hydrogen
/// \param tmp_angl_atoms_noh  Lists of angle atoms (and parameters) containing only heavy atoms
/// \param tmp_dihe_atoms_h    Lists of dihedral atoms (and parameters) containing hydrogen
/// \param tmp_dihe_atoms_noh  Lists of dihedral atoms (and parameters) containing only heavy atoms
BasicValenceTable basicValenceIndexing(int atom_count, const std::vector<int> &tmp_bond_atoms_h,
                                       const std::vector<int> &tmp_bond_atoms_noh,
                                       const std::vector<int> &tmp_angl_atoms_h,
                                       const std::vector<int> &tmp_angl_atoms_noh,
                                       const std::vector<int> &tmp_dihe_atoms_h,
                                       const std::vector<int> &tmp_dihe_atoms_noh);

/// \brief Process indexing for "CHARMM" family force field terms.  These terms originate in
///        the CHARMM force field, but have been imported into Amber force fields to some extent.
///
/// \param atom_count             Total number of atoms in the system
/// \param tmp_ub_atoms           Lists of Urey-Bradley atoms (and parameter indices)
/// \param tmp_charmm_impr_atoms  Lists of CHARMM improper atoms (and parameter indices)
/// \param tmp_cmap_atoms         Lists of CHARMM CMAP atoms (and surface parameter indices)
CharmmValenceTable charmmValenceIndexing(int atom_count, const std::vector<int> &tmp_ub_atoms,
                                         const std::vector<int> &tmp_charmm_impr_atoms,
                                         const std::vector<int> &tmp_cmap_atoms);

/// \brief Condense the 1:4 attenuations from the available basic dihedral parameters into the
///        unguarded struct above.
///
/// \param bvt                       Information on dihedral interactions, among other things
/// \param dihe_elec_screenings      Electrostatic 1:4 screening factors
/// \param dihe_vdw_screenings       van-der Waals 1:4 screening factors
/// \param default_elec14_screening  Not all topologies will offer arrays of 1:4 scaling factors
///                                  to screen the short-range non-bonded interactions.  This
///                                  default, which is not available from the Amber prmtop format
///                                  but rather expected to some in through an MD control namelist,
///                                  will serve when more detailed information is missing.
/// \param default_vdw14_screening   Substitute when 1:4 van-der Waals (Lennard-Jones) scaling
///                                  factors are missing.
AttenuationParameterSet condenseScreeningFactors(const BasicValenceTable &bvt,
                                                 const std::vector<double> &dihe_elec_screenings,
                                                 const std::vector<double> &dihe_vdw_screenings,
                                                 const double default_elec14_screening,
                                                 const double default_vdw14_screening);

/// \brief Process the information relevant to virtual sites.  This function contains a check on
///        the number of virtual sites by independently checking the number of particles with zero
///        mass.  This will unify virtual sites with implicit frame types in the topology's bonded
///        interactions with custom-defined virtual sites.
///
/// \param expected_vsite_count  The expected number of virtual sites (as stated in the topology)
/// \param file_name             Name of the topology file guiding this (for error reporting)
/// \param tmp_masses            Array of atom masses, in a raw form read in from the topology
/// \param bvt                   Table of basic valence parameters, containing arrays to denote
///                              which atoms are bonded to one another
/// \param tmp_bond_equilibria   Parameter array of bond equilibrium lengths
/// \param tmp_angl_equilibria   Parameter array of angle equilibrium spreads
/// \param tmp_atom_names        Name of all atoms in the system, for tracing purposes if a problem
///                              is found during this procedure
/// \param tmp_residue_names     Names of all residues in the system, for tracing purposes if a
///                              problem is found during this procedure
/// \param tmp_residue_limits    Limits (starting atom indices) of all residues if a problem is
///                              found during this procedure
/// \param vsite_custom_frames   Custom virtual site frame indexing data
/// \param vsite_custom_details  Custom virtual site frame dimensions
VirtualSiteTable listVirtualSites(int expected_vsite_count, const std::string &file_name,
                                  const std::vector<double> &tmp_masses,
                                  const BasicValenceTable &bvt,
                                  const std::vector<double> &tmp_bond_equilibria,
                                  const std::vector<double> &tmp_angl_equilibria,
                                  const std::vector<char4> &tmp_atom_types,
                                  const std::vector<char4> &tmp_atom_names,
                                  const std::vector<int> &vsite_custom_frames,
                                  const std::vector<double> &vsite_custom_details);

/// \brief Accumulate the bounds of an exclusion array, based on prior accumulations and lists of
///        previously excluded atoms.
///
/// \param current_bounds      The exclusion bounds array to (further) accumulate
/// \param atom_a              First atoms of each pair
/// \param atom_b              Second atoms of each pair
/// \param vsite_child_bounds  Bounds array of each parent atom's virtual site children
/// \param vsite_child_list    List of virtual site children for all atoms
void accumulateExclusionBounds(std::vector<int> *current_bounds, const int atom_a,
                               const int atom_b, const std::vector<int> &vsite_child_bounds,
                               const std::vector<int> &vsite_child_list);

/// \brief Accumulate two atoms that exclude one another in their respective exclusion lists.
///
/// \param excl_list  The exclusion list
/// \param counters   Increments for each atom beyond its bounds, marching through the list
/// \param bounds     Bounds array for the exclusion list
/// \param atom_i     The first atom
/// \param atom_j     The second atom
void markPairExclusion(std::vector<int> *excl_list, std::vector<int> *counters,
                       const std::vector<int> &bounds, const int atom_i, const int atom_j);

/// \brief Mark all exclusions for a pair of real atoms as well as any virtual site children that
///        either atom may have.  This routine does not consider duplicates or exclusions that
///        have already been marked in another list.  The exclusion list is assumed to have enough
///        space and will be culled later.
///
/// \param excl_list           Developing list of exclusions for all atoms
/// \param excl_counters       Counters for each atom's developing list of exclusions
/// \param excl_bounds         Boundaries of each atom's exclusions in the developing list
/// \param atom_a              The first atom of the interacting pair
/// \param atom_b              The second atom of the interacting pair
/// \param vsite_child_bounds  Boundaries of this atom's virtual site children in the list
/// \param vsite_child_list    List of virtual site children for all atoms
void markExclusions(std::vector<int> *excl_list, std::vector<int> *excl_counters,
                    const std::vector<int> &excl_bounds, const int atom_a, const int atom_b,
                    const std::vector<int> &vsite_child_bounds,
                    const std::vector<int> &vsite_child_list);

/// \brief Cull duplicate exclusions found within the same list.
///
/// \param excl_list    The list of exclusions (pointer, will change upon exit)
/// \param excl_bounds  The bounds of exclusions for each atom within the list (pointer, will
///                     change upon exit)
void cullDuplicateExclusions(std::vector<int> *excl_list, std::vector<int> *excl_bounds);

/// \brief Cull exclusions from one list based on coverage in some prior list
///
/// \param excl_list    The list of exclusions (pointer, will change upon exit)
/// \param excl_bounds  The bounds of exclusions for each atom within the list (pointer, will
///                      change upon exit)
/// \param prior_list    The list of exclusions (pointer, will change upon exit)
/// \param prior_bounds  The bounds of exclusions for each atom within the list (pointer, will
///                      change upon exit)
void cullPriorExclusions(std::vector<int> *excl_list, std::vector<int> *excl_bounds,
                         const std::vector<int> &prior_list,
                         const std::vector<int> &prior_bounds);

/// \brief Create a map of non-bonded exclusions of each order (1:1, 1:2, 1:3, and 1:4) for every
///        atom in the topology.  The map is reflexive: if atom I excludes atom J in some capacity,
///        then J likewise excludes I.  This can all be done based on information about bonds and
///        virtual sites: although some virtual sites may be implicit, they will be detected by
///        virtue of the fact that they have no mass, so that their parent atoms and other features
///        can be detected appropriately.
///
/// Overloaded:
///   - Provide the minimal arrays for the computation
///   - Provide some scratch work objects from which the necessary arrays can be extracted
///
/// \param atom_count       Total number of atoms in the system (for allocating exclusion bounds
///                         arrays, each atom has a given number of 1:1, 1:2, 1:3, or 1:4
///                         exclusions)
/// \param vs_particles     List of all virtual sites in the system
/// \param vs_parent_atoms  List of all virtual sites in the system
/// \param bond_i_atoms     List of I atoms in all bonds I -> J
/// \param bond_j_atoms     List of J atoms in all bonds I -> J
/// \param bvt              Table of basic valence terms and their indexing into the original
///                         topology
/// \param vst              Table of virtual sites and their frames
/// \{
Map1234 mapExclusions(int atom_count, const std::vector<int> &vs_particles,
                      const std::vector<int> &vs_parent_atoms,
                      const std::vector<int> &bond_i_atoms, const std::vector<int> &bond_j_atoms);
  
Map1234 mapExclusions(int atom_count, const BasicValenceTable &bvt, const VirtualSiteTable &vst);
/// \}

/// \brief Check over each atom's exclusions one by one, searching for its match within a range,
///        incrementing an array to indicate coverage if the matching exclusion is found.
///
/// \param coverage   Array of counters indicating how many times a particular exclusion in the
///                   struct developed from the original topology's exclusion list is addressed by
///                   exclusions in the de-novo lists created by this reader
/// \param bounds     Bounds array for one of the de-novo exclusion lists
/// \param excl_list  One of the de-novo exclusion lists
/// \param ce         Condensed array of exclusions created from lists in the original topology
void exclusionSearchLoop(std::vector<int> &coverage, const std::vector<int> &bounds,
                         const std::vector<int> &excl_list, const CondensedExclusions &ce);

/// \brief Check that the exclusions computed for an Amber topology match those in the original
///        exclusion list.
///
/// \param ce           The condensed exclusions list
/// \param all_nb_excl  The list of all 1:1, 1:2, 1:3, and 1:4 exclusions, double-counting all
///                     A -> B and B -> A pairs
/// \param file_name    Name of the topology source file, for error reporting purposes
void checkExclusions(const CondensedExclusions &ce, const Map1234 &all_nb_excl,
                     const std::string &file_name);

/// \brief Detecting hydrogen mass repartitioning is a tricky thing to do, but it must be done
///        when a topology does not have atomic numbers of its own.  What are the correct masses
///        for each element, to begin with?  If repartitioning has occurred, it must have moved
///        mass around within a molecule, so having a clear definition of where each molecule
///        begins and ends is critical.  First, assume that each atom takes the average mass of
///        isotopes in its natural abundance.  If that can explain the masses of every atom in a
///        given molecule, good.  If not, try explaining the masses of some of the lighter atoms
///        using the masses of less common isotopes.  If the masses still cannot be explained,
///        look to repartitioning, positing that there is one repartitioned mass for hydrogen
///        atoms throughout the simulation, starting from a natural 1.008 amu baseline and pulling
///        mass from atoms to which the hydrogens are bonded.
///
/// \param masses       The masses of all atoms in the system
/// \param atom_names   The names of all atoms, so that "EP" and "LP"-containing names may be
///                     flagged as virtual sites
/// \param nb12_list    Table of 1:2 non-bonded exclusions, to indicate bonding patterns
/// \param nb12_bounds  Limits for each atom in nb12_list
/// \param source       Name of the original topology file (for error reporting purposes)
/// \param policy       Direction for what to do in the event of bad input
std::vector<int> atomicNumbersFromMasses(const std::vector<double> &masses,
                                         const std::vector<char4> &atom_names,
                                         const std::vector<int> &nb12_list,
                                         const std::vector<int> &nb12_bounds,
                                         const std::string &source,
                                         const ExceptionResponse policy);

/// \brief Trace the bonding patterns to determine all of the atoms in each molecule.  Use a tree
///        search strategy until all branches are exhausted.  Return a vector with the index of
///        the molecule to which each atom belongs.  For most topologies, this list should place
///        series of contiguous atoms for each molecule in order and follow the molecule sizes, if
///        available, laid out elsewhere in the topology.  But check that.
///
/// \param all_nb_excl  A collection of all 1:1, 1:2, 1:3, and 1:4 interactions
std::vector<int> traceBondedPatterns(const Map1234 &all_nb_excl);

/// \brief Encapsulate the protocol to go from atoms connected by known bonds into a series of
///        defined molecules.  This function is an ugly mess of pointer arguments, but the code
///        would need to be replicated in other contexts, to reading other topology formats,
///        otherwise.
///
/// \param atom_count           The number of atoms in the system
/// \param molecule_count       Number of detected, individual molecules (returned)
/// \param all_nb_excl          Collection of all nonbonded exclusions (for molecule tracing)
/// \param molecule_membership  Membership of each atom in its home molecule (reallocated and
///                             returned)
/// \param molecule_limits      Ascertained boundaries of each molecule in molecule_contents
/// \param molecule_contents    Contents of each molecule, presented as a list of atom indices
///                             that could, in theory, hop around in the topology (in practice the
///                             contents of each molecule will be a continuous stream of numbers
///                             increasing monotonically)
void mapMolecules(const int atom_count, int *molecule_count, const Map1234 &all_nb_excl,
                  std::vector<int> *molecule_membership, std::vector<int> *molecule_limits,
                  std::vector<int> *molecule_contents);

/// \brief Compute a stencil for obtaining the first partial derivative of a cubic spline on a
///        regular grid.  The first derivatives at any given grid point are a function of all the
///        neighboring points, in both directions, along the axis of interest.  Mixed partial
///        derivatives are obtained by taking the outer product of the one-dimensional stencil.
///
/// \param npts  The number of grid points on the periodic grid (the non-periodic case can be
///              obtained by more a procedure incorporating the boundary conditions)
std::vector<double> cubicSplineDerivativeStencil(const int npts);

/// \brief Compute derivatives and tabulate other quantities for evaluating CMAP surfaces at
///        arbitrary points.
///
/// \param cmap_surf_count       A count of the number of CMAPs (fed in from the count in the
///                              function assembling the associated AtomGraph, should equal the
///                              length of tmp_cmap_surf_dims)
/// \param tmp_cmap_surf_dims    Dimensions of each CMAP surface
/// \param tmp_cmap_surf_bounds  Bounds of each surface in tmp_cmap_surfaces
/// \param tmp_cmap_surfaces     Values of CMAP surfaces at the actual grid points
CmapAccessories computeCmapDerivatives(int cmap_surf_count,
                                       const std::vector<int> tmp_cmap_surf_dims,
                                       const std::vector<int> tmp_cmap_surf_bounds,
                                       const std::vector<double> tmp_cmap_surfaces);

/// \brief Review the dihedrals and 1:4 exclusions, to verify that each exclusion is covered by
///        some dihedral which can assign it a parameter.  If there are 1:4 terms outside of the
///        dihedral coverage, they must be entered into a separate list and processed accordingly.
///        Return a list of any 1:4 interactions that still need coverage.
///
/// \param atom_count            The total number of atoms in the topology
/// \param atomic_numbers        List of atomic numbers for all atoms in the topology
/// \param bvt                   Table of basic valence terms (needed for dihedral information)
/// \param all_nb_excl           Lists of non-bonded exclusions (for 1:4 attenuated interactions)
/// \param vs_tbl                Virtual sites staging table
/// \param dihe_lj_screenings    Screening factors (van-der Waals / Lennard-Jones) for each
///                              dihedral parameter set (not indexed by dihedral terms, but by the
///                              unique dihedral parameter set)
/// \param dihe_elec_screenings  Electrostatic creening factors for each dihedral parameter set
/// \param policy                Indicator of what to do if a 1:4 interaction is double-covered
std::vector<int3> checkDihedral14Coverage(int atom_count, const std::vector<int> &atomic_numbers,
                                          const BasicValenceTable &bvt, const Map1234 &all_nb_excl,
                                          const VirtualSiteTable &vs_tbl,
                                          const AttenuationParameterSet &attn_parm,
                                          ExceptionResponse policy);

/// \brief Check that the largest residue is indeed the largest according to the residue limits
///        array.  Amber topologies explicitly record the number of atoms in this residue and the
///        number could be wrong if the residues have been altered for any reason.
///
/// \param residue_limits   The (authoritative) limits array
/// \param current_lr_size  The current guess as to the number of atoms in the largest residue
/// \param policy           Indicator of what to do if there is disagreement
int reviewLargestResidue(const std::vector<int> residue_limits, int current_lr_size,
                         ExceptionResponse policy);

/// \brief Match a long atom type name with its char4 codification.  The list of extended names
///        is all char4, three char4 devoted to names of up to 12 characters each, followed by
///        one char4 denoting the codified name found in the corresponding standard array.
///        As with other searches, this returns an index to the matching [ (char4 x 3) char4 ]
///        index, or the length of the list if no such entry can be found.
///
/// \param overflow_names  The list of extended names (this should be obtained from some AtomGraph
///                        abstract)
/// \param n_overflow      The number of overflow parameters to search (must not exceed the size
///                        of overflow_names
/// \param query           The extended name to match, with possible wildcards
/// \param wildcards       Series of wildcard indicators, equal in length to long_name or else no
///                        character in the vector is assumed to be a wildcard
std::vector<int> matchExtendedName(const char4* overflow_names, int n_overflow,
                                   const std::string &query,
                                   const std::vector<WildCardKind> &wildcards = {});

/// \brief Obtain SETTLE parameters for a particular water arrangement.  This will be repeated for
///        all SETTLE candidates, to collect an array of parameters for the SETTLE procedure if
///        necessary and treat different groups with the appropriate parameter set.  The SETTLE
///        parameters are returned as a double-precision real tuple of four numbers: the mass
///        of the oxygen (x member), mass of either hydrogen (y member), distance between the
///        oxygen atom and the center of mass of the constrained system (z member), and the
///        distance of either hydrogen to the center of mass (w member, the equilibrium geometry
///        must be an Isosceles triangle).
///
/// \param ox_idx           Topological index of the oxygen atom
/// \param h1_idx           Topological index of the first hydrogen atom
/// \param h2_idx           Topological index of the second hydrogen atom
/// \param atomic_mass      Masses of all atoms in the topology
/// \param bvt              Connectivity information from which to obtain parameter indices into
///                         the bond and angle equilibrium value arrays
/// \param bond_equilibria  Equilibrium lengths of all harmonic bond parameters
/// \param angl_equilibria  Equilibrium lengths of all harmonic angle parameters
SettleParm getSettleParameters(int ox_idx, int h1_idx, int h2_idx,
                               const std::vector<double> &atomic_masses,
                               const BasicValenceTable &bvt,
                               const std::vector<double> &bond_equilibria,
                               const std::vector<double> &angl_equilibria);

} // namespace topology
} // namespace stormm

#include "atomgraph_refinement.tpp"

#endif
