// -*-c++-*-
#ifndef STORMM_PUREMESH_H
#define STORMM_PUREMESH_H

#include <algorithm>
#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/tricubic_cell.h"
#include "Math/vector_ops.h"
#include "Math/rounding.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/polynumeric.h"
#include "Potential/energy_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "structure_enumerators.h"
#include "mesh_parameters.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using data_types::getStormmScalarTypeName;
using data_types::isScalarType;
using data_types::isFloatingPointScalarType;
using constants::CartesianDimension;
using constants::UnitCellAxis;
using energy::NonbondedPotential;
using energy::VdwCombiningRule;
using math::addScalarToVector;
using math::elementwiseMultiply;
using math::getTricubicMatrix;
using math::roundUp;
using math::TricubicCell;
using parse::NumberFormat;
using numerics::default_globalpos_scale_bits;
using numerics::fixedPrecisionGrid;
using topology::AtomGraph;
using topology::NonbondedKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::PhaseSpace;

/// \brief The templated, writeable abstract of a BackgroundMesh object.
template <typename Txfrm, typename Tdata> struct BackgroundMeshReader {

  /// \brief The constructor takes arguments for all member variables.
  BackgroundMeshReader(const MeshParamKit<Txfrm> &dims_in, GridDetail kind_in,
                       NonbondedPotential field_in, const llint* avec_x_in, const llint* avec_y_in,
                       const llint* avec_z_in, const llint* bvec_x_in, const llint* bvec_y_in,
                       const llint* bvec_z_in, const llint* cvec_x_in, const llint* cvec_y_in,
                       const llint* cvec_z_in, const llint* avec_abs_x_in,
                       const llint* avec_abs_y_in, const llint* avec_abs_z_in,
                       const int* avec_x_ovrf_in, const int* avec_y_ovrf_in,
                       const int* avec_z_ovrf_in, const int* bvec_x_ovrf_in,
                       const int* bvec_y_ovrf_in, const int* bvec_z_ovrf_in,
                       const int* cvec_x_ovrf_in, const int* cvec_y_ovrf_in,
                       const int* cvec_z_ovrf_in, const int* avec_abs_x_ovrf_in,
                       const int* avec_abs_y_ovrf_in, const int* avec_abs_z_ovrf_in,
                       const Tdata* coeffs_in, const int* ngbr_in, const size_t* ngbr_bounds_in);

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  BackgroundMeshReader(const BackgroundMeshReader<Txfrm, Tdata> &original) = default;
  BackgroundMeshReader(BackgroundMeshReader<Txfrm, Tdata> &&original) = default;
  /// \}

  const MeshParamKit<Txfrm> dims;  ///< Dimensions of the mesh.  These are pre-established.
  const GridDetail kind;           ///< The type of mesh, also pre-established.
  const NonbondedPotential field;  ///< The field described by the mesh, also pre-established.

  // The following are Cartesian X, Y, or Z fixed-precision coordinates of grid element origins.
  const llint* avec_x;         ///< Relative X coordinates of grid element origins stepping along
                               ///<   the mesh's "a" axis
  const llint* avec_y;         ///< Relative Y coordinates of "a" axis grid element origins
  const llint* avec_z;         ///< Relative Z coordinates of "a" axis grid element origins
  const llint* bvec_x;         ///< Relative X coordinates of "b" axis grid element origins
  const llint* bvec_y;         ///< Relative Y coordinates of "b" axis grid element origins
  const llint* bvec_z;         ///< Relative Z coordinates of "b" axis grid element origins
  const llint* cvec_x;         ///< Relative X coordinates of "c" axis grid element origins
  const llint* cvec_y;         ///< Relative Y coordinates of "c" axis grid element origins
  const llint* cvec_z;         ///< Relative Z coordinates of "c" axis grid element origins
  const llint* avec_abs_x;     ///< Absolute X coordinates of "a" axis grid element origins
  const llint* avec_abs_y;     ///< Absolute Y coordinates of "a" axis grid element origins
  const llint* avec_abs_z;     ///< Absolute Z coordinates of "a" axis grid element origins
  const int* avec_x_ovrf;      ///< Relative X overflow bits for "a" axis grid element origins
  const int* avec_y_ovrf;      ///< Relative Y overflow bits for "a" axis grid element origins
  const int* avec_z_ovrf;      ///< Relative Z overflow bits for "a" axis grid element origins
  const int* bvec_x_ovrf;      ///< Relative X overflow bits for "b" axis grid element origins
  const int* bvec_y_ovrf;      ///< Relative Y overflow bits for "b" axis grid element origins
  const int* bvec_z_ovrf;      ///< Relative Z overflow bits for "b" axis grid element origins
  const int* cvec_x_ovrf;      ///< Relative X overflow bits for "c" axis grid element origins
  const int* cvec_y_ovrf;      ///< Relative Y overflow bits for "c" axis grid element origins
  const int* cvec_z_ovrf;      ///< Relative Z overflow bits for "c" axis grid element origins
  const int* avec_abs_x_ovrf;  ///< Absolute X overflow bits for "a" axis grid element origins
  const int* avec_abs_y_ovrf;  ///< Absolute Y overflow bits for "a" axis grid element origins
  const int* avec_abs_z_ovrf;  ///< Absolute Z overflow bits for "a" axis grid element origins

  // The content of the mesh is what is actually writeable.
  const Tdata* coeffs;        ///< Coefficients for all mesh elements.  In an OCCLUSION mesh, these
                              ///<   are the bit-packed masks for each cubelet, 64 cubelets making
                              ///<   one element.  In a NONBONDED_FIELD or NONBONDED_ATOMIC mesh,
                              ///<   the coefficients are tricubic splines for each element.
  const int* ngbr;            ///< Concatenated lists of neighbor atoms
  const size_t* ngbr_bounds;  ///< Bounds array for the neighbor lists in ngbr
};
  
/// \brief The templated, writeable abstract of a BackgroundMesh object.
template <typename Txfrm, typename Tdata> struct BackgroundMeshWriter {

  /// \brief The constructor takes arguments for all member variables.
  BackgroundMeshWriter(const MeshParamKit<Txfrm> &measurements, GridDetail kind_in,
                       NonbondedPotential field_in, const llint* avec_x_in, const llint* avec_y_in,
                       const llint* avec_z_in, const llint* bvec_x_in, const llint* bvec_y_in,
                       const llint* bvec_z_in, const llint* cvec_x_in, const llint* cvec_y_in,
                       const llint* cvec_z_in, const llint* avec_abs_x_in,
                       const llint* avec_abs_y_in, const llint* avec_abs_z_in,
                       const int* avec_x_ovrf_in, const int* avec_y_ovrf_in,
                       const int* avec_z_ovrf_in, const int* bvec_x_ovrf_in,
                       const int* bvec_y_ovrf_in, const int* bvec_z_ovrf_in,
                       const int* cvec_x_ovrf_in, const int* cvec_y_ovrf_in,
                       const int* cvec_z_ovrf_in, const int* avec_abs_x_ovrf_in,
                       const int* avec_abs_y_ovrf_in, const int* avec_abs_z_ovrf_in,
                       Tdata* coeffs_in, int* ngbr_in, size_t* ngbr_bounds_in);

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  BackgroundMeshWriter(const BackgroundMeshWriter<Txfrm, Tdata> &original) = default;
  BackgroundMeshWriter(BackgroundMeshWriter<Txfrm, Tdata> &&original) = default;
  /// \}

  const MeshParamKit<Txfrm> dims;  ///< Dimensions of the mesh.  These are pre-established.
  const GridDetail kind;           ///< The type of mesh, also pre-established.
  const NonbondedPotential field;  ///< The field described by the mesh, also pre-established.

  // The following are Cartesian X, Y, or Z fixed-precision coordinates of grid element origins.
  const llint* avec_x;         ///< Relative X coordinates of grid element origins stepping along
                               ///<   the mesh's "a" axis
  const llint* avec_y;         ///< Relative Y coordinates of "a" axis grid element origins
  const llint* avec_z;         ///< Relative Z coordinates of "a" axis grid element origins
  const llint* bvec_x;         ///< Relative X coordinates of "b" axis grid element origins
  const llint* bvec_y;         ///< Relative Y coordinates of "b" axis grid element origins
  const llint* bvec_z;         ///< Relative Z coordinates of "b" axis grid element origins
  const llint* cvec_x;         ///< Relative X coordinates of "c" axis grid element origins
  const llint* cvec_y;         ///< Relative Y coordinates of "c" axis grid element origins
  const llint* cvec_z;         ///< Relative Z coordinates of "c" axis grid element origins
  const llint* avec_abs_x;     ///< Absolute X coordinates of "a" axis grid element origins
  const llint* avec_abs_y;     ///< Absolute Y coordinates of "a" axis grid element origins
  const llint* avec_abs_z;     ///< Absolute Z coordinates of "a" axis grid element origins
  const int* avec_x_ovrf;      ///< Relative X overflow bits for "a" axis grid element origins
  const int* avec_y_ovrf;      ///< Relative Y overflow bits for "a" axis grid element origins
  const int* avec_z_ovrf;      ///< Relative Z overflow bits for "a" axis grid element origins
  const int* bvec_x_ovrf;      ///< Relative X overflow bits for "b" axis grid element origins
  const int* bvec_y_ovrf;      ///< Relative Y overflow bits for "b" axis grid element origins
  const int* bvec_z_ovrf;      ///< Relative Z overflow bits for "b" axis grid element origins
  const int* cvec_x_ovrf;      ///< Relative X overflow bits for "c" axis grid element origins
  const int* cvec_y_ovrf;      ///< Relative Y overflow bits for "c" axis grid element origins
  const int* cvec_z_ovrf;      ///< Relative Z overflow bits for "c" axis grid element origins
  const int* avec_abs_x_ovrf;  ///< Absolute X overflow bits for "a" axis grid element origins
  const int* avec_abs_y_ovrf;  ///< Absolute Y overflow bits for "a" axis grid element origins
  const int* avec_abs_z_ovrf;  ///< Absolute Z overflow bits for "a" axis grid element origins

  // The content of the mesh is what is actually writeable.
  Tdata* coeffs;        ///< Coefficients for all mesh elements.  In an OCCLUSION mesh, these are
                        ///<   the bit-packed masks for each cubelet, 64 cubelets making one
                        ///<   element.  In a NONBONDED_FIELD or NONBONDED_ATOMIC mesh, the
                        ///<   coefficients are tricubic splines for each element.
  int* ngbr;            ///< Concatenated lists of neighbor atoms
  size_t* ngbr_bounds;  ///< Bounds array for the neighbor lists in ngbr
};
  
/// \brief A workspace for constructing a pure potential mesh based on the frozen atoms of a
///        large molecule.  If the large molecule has nonrigid components, they must be excluded
///        from contributing to the grid.  In addition, any atoms up to 1:4 (connected by three
///        bonds or less) must also be excluded from the grid-based potential.  Computations on
///        these atoms will not be accurate off the grid, but since they are frozen the
///        consequences are mitigated.
template <typename T> class BackgroundMesh {
public:

  /// \brief The constructor takes all dimension parameters plus an indication of what type of
  ///        potential, the molecular system, and what mask of atoms is to be mapped.  Variants
  ///        include different ways to define the limits of the mesh.  If a GPU is available, it
  ///        will be used to compute the mesh.
  ///
  /// \param ag               System topology (markings in its mobile_atoms array will be used to
  ///                         determine which atoms to map)
  /// \param cf               Cartesian coordinates of all particles
  /// \param scale_bits       Number of bits after the decimal when taking coordinates in units of
  ///                         Angstroms into the fixed-precision format of the mesh framework
  /// \param buffer           Breadth around the molecule drawing an orthhombic region in which to
  ///                         map the mesh
  /// \param mesh_bounds      Boundaries of the mesh, a six-element vector describing the lower
  ///                         and upper Cartesian X, Y, and Z limits of an orthorhombic region
  ///                         in which to define a rectilinear mesh.
  /// \param spacing          Grid spacings for the mesh cells.  Provide either a single number or
  ///                         three dimensions for the Cartesian X, Y, and Z widths of rectilinear
  ///                         elements.
  /// \param measurements_in  A full description of the mesh parameters.  This provides a means to
  ///                         define non-orthorhombic meshes.
  /// \param gpu              Details of the available GPU
  /// \{
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const MeshParameters &measurements_in = MeshParameters());

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                 double probe_radius_in = 0.0, double well_depth_in = 0.0,
                 VdwCombiningRule mixing_protocol_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                 const MeshParameters &measurements_in = MeshParameters(),
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Copy and move constructors must be defined explicitly to accommodate POINTER-kind
  ///        Hybrid objects.
  ///
  /// \param original  The original object to copy or move
  /// \{
  BackgroundMesh(const BackgroundMesh<T> &original);
  BackgroundMesh(BackgroundMesh<T> &&original);
  /// \}

  /// \brief Copy and move assignment operators must be defined explicitly to accommodate
  ///        POINTER-kind Hybrid objects.
  ///
  /// \param other  Assignments are made based on this pre-existing object.
  /// \{
  BackgroundMesh& operator=(const BackgroundMesh<T> &original);
  BackgroundMesh& operator=(BackgroundMesh<T> &&original);
  /// \}
  
  /// \brief Get a const pointer to the topology responsible for creating this mesh.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a const pointer to the coordinates responsible for creating this mesh.
  const CoordinateFrame* getCoordinatePointer() const;

  /// \brief Get an abstract of the mesh for performing calculations in double-precision.
  ///
  /// Overloaded:
  ///   - Get the reader for a const object
  ///   - Get the writer for a mutable object
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  /// \{
  BackgroundMeshReader<double, T>
  dpData(const HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  BackgroundMeshWriter<double, T> dpData(const HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get an abstract of the mesh for performing calculations in single-precision.
  ///
  /// Overloaded:
  ///   - Get the reader for a const object
  ///   - Get the writer for a const object
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  /// \{
  BackgroundMeshReader<float, T>
  spData(const HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  BackgroundMeshWriter<float, T> spData(const HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Set the topology that the mesh shall use.
  ///
  /// \param ag_in  Pointer to the topology of interest
  void setTopologyPointer(const AtomGraph *ag_in);

  /// \brief Set the coordinates that the mesh shall use.
  ///
  /// \param ag_in  Pointer to the coordinates of interest
  void setCoordinatePointer(const CoordinateFrame *ag_in);

  /// \brief Set the dimensions of the mesh.  Variants of this member function that do not include
  ///        a coordinate pointer or topology may assume that one has already been set.  Variants
  ///        that do accept a coordinate pointer will set that as the mesh's coordinate set.
  ///
  /// Overloaded:
  ///   - Take coordinate and topology pointers with a single parameter to define the mesh limits.
  ///   - Take a single parameter to define the mesh limits around a set of coordinates defined by
  ///     a particular topology, and assume that the mesh already has these things.
  ///   - Take the minimum and maximum Cartesian X, Y, and Z limits of a rectilinear mesh (this
  ///     will not require a system to implement)
  ///   - Take a single parameter to define isotropic mesh elements
  ///   - Take three or nine parameters to define anisotropic rectilinear or triclinic mesh
  ///     elements
  ///   - Indicate the scaling bits to be used in fixed-precision arithmetic when determining
  ///     positions on the mesh, or leave that parameter unchanged
  ///
  /// \param ag_in          New topology to build the mesh around (this will be incorporated into
  ///                       the object)
  /// \param cf_in          New coordinate set to build the mesh around (this will be incorporated
  ///                       into the object)
  /// \param padding        Region around the molecule of interest to spread the mesh
  /// \param mesh_bounds    Minimum and maximum Cartesian limits of the mesh
  /// \param spacing        Width (and length, and height, or bounding vectors) of mesh elements
  /// \param scale_bits_in  The number of bis after the decimal to use in fixed-precision
  ///                       computations of particle positions on the grid.  The default value of
  ///                       -100 is nonsensical and will be interpreted as "leave the current
  ///                       setting as it is." Limiting the number of function overloads in this
  ///                       way is a step towards controlling code bloat due to this large,
  ///                       templated class.
  /// \{
  void setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in, double padding,
                         double spacing, int scale_bits_in = -100);

  void setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in, double padding,
                         const std::vector<double> &spacing, int scale_bits_in = -100);

  void setMeshParameters(double padding, double spacing, int scale_bits_in = -100);

  void setMeshParameters(double padding, const std::vector<double> &spacing,
                         int scale_bits_in = -100);

  void setMeshParameters(const std::vector<double> &mesh_bounds, double spacing,
                         int scale_bits_in = -100);

  void setMeshParameters(const std::vector<double> &mesh_bounds,
                         const std::vector<double> &spacing, int scale_bits_in = -100);
  /// \}

  /// \brief Set the probe radius, meaning either the Lennard-Jones potential sigma radius or
  ///        the clash probe radius, depending on the mesh type.  This includes a validity check.
  ///
  /// \param probe_radius_in  The probe radius to set
  void setProbeRadius(double probe_radius_in);

  /// \brief Set the Lennard-Jones well depth for the probe that will generate the potential.  This
  ///        includes a validity check.
  ///
  /// \param well_depth_in  The well depth to set
  void setWellDepth(double well_depth_in);

  /// \brief Set the combining rule that will be used to make the probe interact with the receptor
  ///        on any mesh (with the exception of an electrostatic field).  With geometric combining
  ///        rules, it is possible to tailor a single mesh for all particles that might interact
  ///        with the receptor.  However, with Lorentz-Berthelot rules or any case of
  ///        non-conformant pair rules, new grids are required for each particle type that might
  ///        interact with the mesh potential.
  ///
  /// \param mixing_protocol_in  The method to use
  void setCombiningRule(VdwCombiningRule mixing_protocol_in);

  /// \brief Compute the appropriate field for the mesh.  This is called automatically by the
  ///        constructor if enough information is provided.
  ///
  /// \param gpu  Details of any available GPU
  /// 
  void computeField(const GpuDetails &gpu = null_gpu);

  /// \brief Compute neighbor lists for each mesh element, if appropriate for the mesh type.
  ///
  /// \param gpu  Details of any available GPU
  void computeNeighborLists(const GpuDetails &gpu = null_gpu);
  
private:

  /// Measurements for the mesh, including various transformation matrices and fixed-precision
  /// representations.
  MeshParameters measurements;

  /// The detail level of the mesh--this carries implications for the allowed data types.  An
  /// "OCCLUSION" mesh must be of type ullint (64-bit unsigned integer, to encode 4 x 4 x 4
  /// cubelets of bitwise information for "occluded" or "non-occluded." Such meshes can be very
  /// dense, as every mesh element constitutes 4 x 4 x 4 cubelets, or a total of 512 bytes of
  /// information, as much as 20k mesh points on each side.  In contrast, a "NONBONDED_FIELD" type
  /// grid must have data type float or double to encode the appropriate non-bonded potential or
  /// influence field (each Lennard-Jones type and Generalized Born parameter combination will get
  /// its own field).  A "NONBONDED_ATOMIC" type mesh stores field data in floating-point scalar
  /// data as well as local neighbor lists in an additional array of integer data, with a bounds
  /// array.
  GridDetail kind;

  /// The type of non-bonded potential to map to the mesh.  This is irrelevant for "OCCLUSION"-kind
  /// meshes (the effective potential is a clash based on the Lennard-Jones sigma radii and a probe
  /// width).
  NonbondedPotential field;

  /// The probe radius to use in computing clashes for an OCCLUSION-type mesh.
  double probe_radius;

  /// Lennard-Jones well depth of the probe particle's self interaction.  If the mesh is of type
  /// NONBONDED_FIELD or NONBONDED_ATOMIC, the probe_radius will be taken to mean the particle's
  /// self-interacting sigma parameter, and this is its epsilon.  For electrostatic fields, neither
  /// this number nor probe_radius has any relevance.
  double well_depth;
  
  /// The combining rule used (or to use) in computing clashes and Lennard-Jones interactions.
  VdwCombiningRule mixing_protocol;
  
  /// Fixed-precision Cartesian coordinates of the mesh grid points stepping along the lines
  /// [ i = 0...nx, 0, 0 ] (the "a" box vector), [ 0, j = 0...ny, 0 ] (the "b" box vector),
  /// [ 0, 0, k = 0...nz ] (the "c" box vector), and a fourth vector with absolute coordinates of
  /// the "a" bo vector points, [ i = 0...nx, 0, 0 ] + [ orig_x, orig_y, orig_z ].  Storing these
  /// values obviates the need to do expensive and complicated multiplications of fixed-precision
  /// numbers during grid-based particle calculations.
  /// \{
  Hybrid<llint> a_line_x;      //
  Hybrid<llint> a_line_y;      // Positions of "a" vector tick marks relative to the mesh origin
  Hybrid<llint> a_line_z;      //
  Hybrid<llint> b_line_x;      //
  Hybrid<llint> b_line_y;      // Positions of "b" vector tick marks relative to the mesh origin
  Hybrid<llint> b_line_z;      //
  Hybrid<llint> c_line_x;      //
  Hybrid<llint> c_line_y;      // Positions of "c" vector tick marks relative to the mesh origin
  Hybrid<llint> c_line_z;      //
  Hybrid<llint> a_abs_line_x;  //
  Hybrid<llint> a_abs_line_y;  // Absolute positions of the "a" vector tick marks
  Hybrid<llint> a_abs_line_z;  //
  /// \}

  /// Overflow for fixed-precision Cartesian coordinates of the mesh grid points stepping along
  /// the a, b, and c box vectors.
  /// \{
  Hybrid<int> a_line_x_overflow;      // Overflow bits for positions of the "a" vector tick marks
  Hybrid<int> a_line_y_overflow;      //   relative to the mesh origin.  These engage when the mesh
  Hybrid<int> a_line_z_overflow;      //   serves double-precision coefficients and calculatons.
  Hybrid<int> b_line_x_overflow;      //
  Hybrid<int> b_line_y_overflow;      // Overflow bits for "b" vector tick mark relative positions
  Hybrid<int> b_line_z_overflow;      //
  Hybrid<int> c_line_x_overflow;      //
  Hybrid<int> c_line_y_overflow;      // Overflow bits for "c" vector tick mark relative positions
  Hybrid<int> c_line_z_overflow;      //
  Hybrid<int> a_abs_line_x_overflow;  //
  Hybrid<int> a_abs_line_y_overflow;  // Overflow bits for "a" vector tick mark absolute positions
  Hybrid<int> a_abs_line_z_overflow;  //
  /// \}

  /// Coefficients for tricubic spline functions spanning each grid element.  Coefficients for
  /// element {Ea, Eb, Ec} traversing the "a", "b", and "c" axes, respectively, are held
  /// consecutively in a 64-element series with order {i, j, k} = (0, 0, 0), (0, 0, 1), ...,
  /// (0, 0, 3), (1, 0, 0), (1, 0, 1), ..., (1, 0, 3), ..., (3, 3, 3).  The 64-element series for
  /// element {Ea, Eb, Ec} will be found at 64 * ((((Ea * nb) + Eb) * nc) + Ec), where {na, nb, nc}
  /// are the mesh dimensions along each axis.
  Hybrid<T> coefficients;

  /// Pointer to the original topology
  AtomGraph *ag_pointer;
  
  /// Pointer to the coordinates responsible for creating the potential
  CoordinateFrame *cf_pointer;
  
  /// Snapshot of the frozen atoms in the topology (in case the topology is modified later)
  Hybrid<uint> frozen_atoms;

  /// Element-wise neighbor listts for nearby frozen atoms
  Hybrid<int> neighbor_list;

  /// Bounds array for the neighbor list atoms array above
  Hybrid<size_t> neighbor_list_bounds;
  
  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> int_data;

  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<llint> llint_data;

  /// \brief Obtain bounds for the mesh based on coordinates of frozen atoms.
  ///
  /// Overloaded:
  ///   - Provide a single buffer argument (in Angstroms), indicating the region to map around all
  ///     frozen atoms
  ///   - Provide explicit Cartesian minimum and maximum limits for the mapping (a different
  ///     overload of the constructor, one which does not call any form of getMeasurements(), must
  ///     be used to create a triclinic mesh for something like a crystallographic unit cell)
  ///   - Provide a single number for the length, width, and height of a rectilinear (orthorhombic)
  ///     mesh element.
  ///   - Provide three values for the length, width, and height of an anisotropic (but still
  ///     rectilinear) mesh element.
  ///
  /// \param ag           System topology, containing the list of frozen atoms
  /// \param cf           Coordinates of the system
  /// \param padding      Length to extend the mesh outside the extrema of the frozen atoms
  /// \param mesh_bounds  Six-element vector of minimum and maximum Cartesian X, Y, and Z limits
  ///                     for the mesh
  /// \param spacing      The length, width, and height of the mesh element (a real-valued scalar,
  ///                     or a three-element vector).  Units of Angstroms.
  /// \{
  MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf, double padding,
                                 double spacing, int scale_bits_in) const;

  MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                                 const std::vector<double> &mesh_bounds, double spacing,
                                 int scale_bits_in) const;

  MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf, double padding,
                                 const std::vector<double> &spacing, int scale_bits_in) const;

  MeshParameters getMeasurements(const std::vector<double> &mesh_bounds,
                                 const std::vector<double> &spacing, int scale_bits_in) const;
  /// \}
  
  /// \brief Allocate memory for the mesh coefficients, axis coordinates, and frozen atoms mask.
  ///        This does not allocate space for the concatenated mesh element neighbor lists or their
  ///        bounds array.
  void allocate();

  /// \brief Repair the POINTER-kind Hybrid objects in a newly copied or moved object.
  void rebase_pointers();

  /// \brief Certain types of meshes are restricted to certain data types.  This function will
  ///        ensure the correct relationships.
  void validateMeshKind() const;

  /// \brief Validate the fixed-precision scaling based on the mesh's data type.  Single-precision
  ///        nonbonded field meshes and clash maps must not overflow 64-bit positional
  ///        representations.
  void validateScalingBits() const;
  
  /// \brief Map out the axes of the mesh, drawing three lines radiating from the mesh origin along
  ///        its "a", "b", and "c" axes and marking the Cartesian coordinates of tick marks for the
  ///        origins of each grid point in lines [ i_a, 0, 0 ], [ 0, j_b, 0 ], or [ 0, 0, k_c ] on
  ///        these axes.  All locations are given relative to the origin, such that adding the
  ///        coordinates of the origin to the coordinates of point i_a from the "a" axis line, the
  ///        coordinates of point j_b from the "b" axis line, and point k_c from the "c" axis line
  ///        will give the location of the origin for grid element { i_a, j_b, k_c }.  A fourth
  ///        axis, the "a" axis points with the origin coordinates "baked in", is provided for
  ///        convenience.  Given that that the majority of transactions with the mesh will involve
  ///        a whole warp computing the 64-term tricubic polynomial for one particle after just one
  ///        thread computes the proper element and broadcasts the relative displacments within
  ///        that element, this savings may be significant.
  void computeMeshAxisCoordinates();

  /// \brief Color the exclusion mesh based on a particular set of coordinates.
  ///
  /// \param gpu  Details of any available GPU
  void colorExclusionMesh(const GpuDetails &gpu);

  /// \brief Map the electrostatics of the receptor's rigid atoms.
  ///
  /// \param gpu  Details of any available GPU
  void mapElectrostatics(const GpuDetails &gpu);

  /// \brief Common function for various NONBONDED_FIELD meshes making use of a series of grids
  ///        for the function values plus partial and mixed partial derivatives.  This works for
  ///        CPU-based code.  For an M x N x P mesh, annd grids measure (M + 1) x (N + 1) x (P + 1)
  ///        points.  The grids are populated with Cartesian X, Y, and Z derivatives and will feed
  ///        into functions for computing spline coefficients for any shape of mesh.
  ///
  /// \param u_grid       Grid of function values
  /// \param dudx_grid    Grid of Cartesian X partial derivatives
  /// \param dudy_grid    Grid of Cartesian Y partial derivatives
  /// \param dudz_grid    Grid of Cartesian Z partial derivatives
  /// \param dudxy_grid   Grid of Cartesian X and Y double derivatives
  /// \param dudxz_grid   Grid of Cartesian X and Z double derivatives
  /// \param dudyz_grid   Grid of Cartesian Y and Z double derivatives
  /// \param dudxyz_grid  Grid of triple derivatives
  void computeCoefficients(const std::vector<double> &u_grid,
                           const std::vector<double> &dudx_grid,
                           const std::vector<double> &dudy_grid,
                           const std::vector<double> &dudz_grid,
                           const std::vector<double> &dudxy_grid,
                           const std::vector<double> &dudxz_grid,
                           const std::vector<double> &dudyz_grid,
                           const std::vector<double> &dudxyz_grid);
};

} // namespace structure
} // namespace stormm

#include "background_mesh.tpp"

#endif
