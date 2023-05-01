// -*-c++-*-
#ifndef STORMM_NML_RECEPTOR_H
#define STORMM_NML_RECEPTOR_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Parsing/textfile.h"
#include "Structure/structure_enumerators.h"
#include "namelist_element.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::CartesianDimension;
using constants::ExceptionResponse;
using constants::UnitCellAxis;
using parse::TextFile;
using parse::WrapTextSearch;
using structure::BoundaryCondition;
using structure::GridDetail;
using structure::MeshPosition;
using structure::RMSDMethod;

/// \brief Default parameters for the &receptor namelist
constexpr int default_receptor_grid_dim = -1;
constexpr double default_receptor_grid_angle = 0.5 * symbols::pi;
constexpr double default_receptor_grid_origin = 0.0;
constexpr double default_receptor_grid_spacing = 1.0;
constexpr char default_receptor_grid_potential[] = "occlusion";
constexpr char default_receptor_grid_alignment[] = "molecule";
constexpr char default_receptor_grid_boundary[] = "isolated";
constexpr char default_receptor_alignment_method[] = "align_mass";

/// \brief Encapsulate the data extracted from a &receptor namelist to define a grid-mapped
///        representation of a rigid macromolecular structure.
///
class ReceptorControls {
public:
  
  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  ReceptorControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  ReceptorControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  ReceptorControls(const ReceptorControls &original) = default;
  ReceptorControls(ReceptorControls &&original) = default;
  ReceptorControls& operator=(const ReceptorControls &original) = default;
  ReceptorControls& operator=(ReceptorControls &&original) = default;
  /// \}

  /// \brief Get the label group in which to find all systems composing the receptor grid.  This
  ///        will correspond to a label from the &files namelist.
  const std::string& getLabelGroup() const;

  /// \brief Get the alignment mask showing which atoms and residues are most critical to the
  ///        alignment of multiple structures if more than one receptor structure is at hand.
  const std::string& getAlignmentMask() const;

  /// \brief Get the alignment method.
  RMSDMethod getAlignmentMethod() const;
  
  /// \brief Get the number of grid points in a particular dimension.  This can serve to query a
  ///        particular piece of information without invoking the more complex but comprehensive
  ///        MeshParameters object.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param dim  The axis along which to measure the number of mesh points
  /// \{
  int getAxisElementCount(UnitCellAxis dim) const;
  int getAxisElementCount(CartesianDimension dim) const;
  /// \}

  /// \brief Get the mesh spacing along a particular axis.  Overloads and input parameters to this
  ///        function follow from getAxisElementCount() above.
  /// \{
  double getMeshSpacing(UnitCellAxis dim) const;
  double getMeshSpacing(CartesianDimension dim) const;
  /// \}

  /// \brief Get the mesh alpha angle, between the b and c unit cell axes.  The value is returned
  ///        in units of radians.
  double getMeshAlpha() const;

  /// \brief Get the mesh beta angle, between the a and c unit cell axes, in units of radians.
  double getMeshBeta() const;

  /// \brief Get the mesh gamma angle, between the a and b unit cell axes, in units of radians.
  double getMeshGamma() const;
  
  /// \brief Get the origin of the mesh in one Cartesian dimension.
  double getMeshOrigin(CartesianDimension dim) const;

  /// \brief Get buffer widths between mesh boundaries and the receptor structure or structures.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param dim  The Cartesian or unit cell axis normal to the face of interest (along which the
  ///             width is relevant)
  /// \{
  double getBufferWidth(UnitCellAxis dim) const; 
  double getBufferWidth(CartesianDimension dim) const;
  /// \}

  /// \brief Get the type of content (potential energy field) that the mesh will represent.
  GridDetail getMeshContent() const;
  
  /// \brief Get the boundary conditions of the mesh
  BoundaryCondition getMeshBoundaries() const;
  
  /// \brief Get the manner in which the mesh is aligned to the rigid molecule it represents.
  MeshPosition getMeshPosition() const;

  /// \brief Get the number of bits after the decimal to be used in composing the mesh-based
  ///        field as well as the positions of its vertices.
  int getMeshScalingBits() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;

  /// \brief Set the label group for structures to be used in composing the mesh.
  void setLabelGroup(const std::string &label_group_in);

  /// \brief Set the number of mesh points along a particular axis.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param mesh_points_in  The desired number of mesh points
  /// \param dim             The unit cell axis in question
  /// \{
  void setMeshElementCount(int mesh_points_in, UnitCellAxis dim);
  void setMeshElementCount(int mesh_points_in, CartesianDimension dim);
  /// \}

  /// \brief Set the mesh spacing along a particular axis.  Overloading of this function follows
  ///        from setMeshElementCount() above.
  ///
  /// \param mesh_spacing_in  The desired distance between mesh points, in Angstroms
  /// \param dim              The unit cell axis in question
  /// \{
  void setMeshSpacing(double mesh_spacing_in, UnitCellAxis dim);
  void setMeshSpacing(double mesh_spacing_in, CartesianDimension dim);
  /// \}

  /// \brief Set the mesh alpha angle.
  ///
  /// \param alpha_in  The angle to set, in units of radians.
  void setMeshAlphaAngle(double alpha_in);

  /// \brief Set the mesh beta angle.
  ///
  /// \param beta_in  The angle to set, in units of radians.
  void setMeshBetaAngle(double beta_in);

  /// \brief Set the mesh gamma angle.
  ///
  /// \param gamma_in  The angle to set, in units of radians.
  void setMeshGammaAngle(double gamma_in);

  /// \brief Set the mesh origin along one Cartesian axis.
  ///
  /// \param mesh_origin_in  The desired mesh origin coordinate
  /// \param dim             The Cartesian axis of interest
  void setMeshOrigin(double mesh_origin_in, CartesianDimension dim);

  /// \brief Set the buffer distance between the mesh boundary and the nearest van-der Waals
  ///        sphere.
  ///
  /// Overloaded:
  ///   - Set a single parameter for all three unit cell or Cartesian axes
  ///   - Set separate parameters for secific unit cell or Cartesian axes
  ///
  /// \param buffer_width_in  The buffer width to set, in units of Angstroms
  /// \param dim              The unit cell or Cartesian axis normal to the face of interest
  /// \{
  void setBufferWidth(double buffer_width_in);
  void setBufferWidth(double buffer_width_in, UnitCellAxis dim);
  void setBufferWidth(double buffer_width_in, CartesianDimension dim);
  /// \}
  
  /// \brief Set the potential energy field that the mesh will represent.
  ///
  /// Overloaded:
  ///   - Provide the enumerated value explicitly
  ///   - Translate a string into the appropriate enumeration
  ///
  /// \param potential_in  The type of energetic field
  /// \{
  void setMeshPotential(const std::string &potential_in);
  void setMeshPotential(GridDetail potential_in);
  /// \}

  /// \brief Set the boundary conditions for the mesh.  Overloading in this function follows from
  ///        setMeshPotential() above.
  ///
  /// \param boundaries_in  The chosen boundary conditions
  /// \{
  void setMeshBoundaries(const std::string &boundaries_in);
  void setMeshBoundaries(BoundaryCondition boundaries_in);
  /// \}

  /// \brief Set the positioning of the mesh relative to the receptor molecule is describes.
  ///        Overloading in this function follows from setMeshPotential() above.
  ///
  /// \param alignment_in  The chosen mesh positioning
  /// \{
  void setMeshPosition(const std::string &alignment_in);
  void setMeshPosition(MeshPosition alignment_in);
  /// \}

  /// \brief Set the number of bits after the decimal to be used in mesh calculations.
  ///
  /// \param scaling_bits_in  The bit count
  void setMeshScalingBits(int mesh_scaling_bits_in);
  
private:
  ExceptionResponse policy;  ///< The course to take when encountering bad input
  std::string label_group;   ///< The label group from which to draw structures for the mesh.
                             ///<   Meshes expressing OCCLUSION potentials can draw from more than
                             ///<   one structure, but others may require single structures or more
                             ///<   strenuous approximations in order to draw meaningful maps of
                             ///<   receptor ensembles.
  std::string align_mask;    ///< Atom mask for aligning multiple structures
  std::string align_method;  ///< Chosen method for aligning multiple receptor structures according
                             ///<   to the atom mask given in align_mask
  int mesh_points_a;         ///< Number of mesh points along the mesh (or unit cell) c vector
  int mesh_points_b;         ///< Number of mesh points along the mesh (or unit cell) c vector
  int mesh_points_c;         ///< Number of mesh points along the mesh (or unit cell) c vector
  double mesh_spacing_a;     ///< Regular spacing between mesh points along the unit cell a vector
  double mesh_spacing_b;     ///< Regular spacing between mesh points along the unit cell b vector
  double mesh_spacing_c;     ///< Regular spacing between mesh points along the unit cell c vector
  double mesh_alpha;         ///< Unit cell angle between the b and c unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_beta;          ///< Unit cell angle between the a and c unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_gamma;         ///< Unit cell angle between the a and b unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_origin_x;      ///< Mesh origin along the Cartesian X axis
  double mesh_origin_y;      ///< Mesh origin along the Cartesian Y axis
  double mesh_origin_z;      ///< Mesh origin along the Cartesian Z axis
  double buffer_width_a;     ///< In isolated boundary conditions, this is the minimum distance
                             ///<   between any van-der Waals sphere surface and the nearest (b x c
                             ///<   plane) mesh boundary face.  In periodic boundary conditions,
                             ///<   the same concept applies although the distances technically
                             ///<   refer to the unit cell.
  double buffer_width_b;     ///< Minimum distance between any van-der Waals sphere of the receptor
                             ///<   and the plane of one of the mesh boundary's a x c faces
  double buffer_width_c;     ///< Minimum distance between any van-der Waals sphere of the receptor
                             ///<   and the plane of one of the mesh boundary's a x b faces
  std::string potential;     ///< The potential field that the mesh will display
  std::string boundaries;    ///< Boundary conditions on the mesh
  std::string mesh_position; ///< The way in which the mesh aligns to the molecule it represents
  int mesh_scaling_bits;     ///< Number of bits after the decimal to use when computing the
                             ///<   positions of mesh vertices and the potential expressed on it

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};

/// \brief Free function to read the &receptor namelist.
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator receptorInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
