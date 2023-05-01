#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "Structure/mesh_parameters.h"
#include "input.h"
#include "nml_receptor.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using parse::minimalRealFormat;
using structure::default_mesh_scaling_bits;
using structure::getEnumerationName;
using structure::translateBoundaryCondition;
using structure::translateGridDetail;
using structure::translateMeshPosition;
using structure::translateRMSDMethod;

//-------------------------------------------------------------------------------------------------
ReceptorControls::ReceptorControls(const ExceptionResponse policy_in) :
  policy{policy_in},
  label_group{std::string("")}, align_mask{std::string("")},
  align_method{std::string(default_receptor_alignment_method)},
  mesh_points_a{default_receptor_grid_dim},
  mesh_points_b{default_receptor_grid_dim},
  mesh_points_c{default_receptor_grid_dim},
  mesh_spacing_a{default_receptor_grid_spacing},
  mesh_spacing_b{default_receptor_grid_spacing},
  mesh_spacing_c{default_receptor_grid_spacing},
  mesh_alpha{default_receptor_grid_angle},
  mesh_beta{default_receptor_grid_angle},
  mesh_gamma{default_receptor_grid_angle},
  mesh_origin_x{0.0}, mesh_origin_y{0.0}, mesh_origin_z{0.0},
  buffer_width_a{0.0}, buffer_width_b{0.0}, buffer_width_c{0.0},
  potential{std::string(default_receptor_grid_potential)},
  boundaries{std::string(default_receptor_grid_boundary)},
  mesh_position{std::string(default_receptor_grid_alignment)},
  mesh_scaling_bits{default_mesh_scaling_bits},
  nml_transcript{"receptor"}
{}

//-------------------------------------------------------------------------------------------------
ReceptorControls::ReceptorControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  ReceptorControls(policy_in)
{
  const NamelistEmulator t_nml = receptorInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  if (t_nml.getKeywordStatus("label_group") == InputStatus::MISSING) {

    // Always fatal: the receptor grid must have a structure to operate on, and a distinct label
    // is the only way to specify that structure.
    rtErr("A label group must be provided in order to construct a mesh based on a rigid or "
          "semi-rigid receptor molecule.", "ReceptorControls");
  }
  else {
    label_group = t_nml.getStringValue("label_group");
  }
  t_nml.assignVariable(&align_mask, "alignment_mask");
  align_method = t_nml.getStringValue("alignment_method");

  // Mesh dimensions are set to default values during the object's construction, but can be
  // replaced with other values if the associated keywords are specified in the user input.
  t_nml.assignVariable(&mesh_points_a, &mesh_points_b, &mesh_points_c, "mesh_dim");
  t_nml.assignVariable(&mesh_points_a, "mesh_dim_a");
  t_nml.assignVariable(&mesh_points_b, "mesh_dim_b");
  t_nml.assignVariable(&mesh_points_c, "mesh_dim_c");
  t_nml.assignVariable(&mesh_spacing_a, &mesh_spacing_b, &mesh_spacing_c, "mesh_spacing");
  t_nml.assignVariable(&mesh_spacing_a, "mesh_spacing_a");
  t_nml.assignVariable(&mesh_spacing_b, "mesh_spacing_b");
  t_nml.assignVariable(&mesh_spacing_c, "mesh_spacing_c");
  t_nml.assignVariable(&mesh_alpha, symbols::pi / 180.0, "mesh_alpha");
  t_nml.assignVariable(&mesh_beta,  symbols::pi / 180.0, "mesh_beta");
  t_nml.assignVariable(&mesh_gamma, symbols::pi / 180.0, "mesh_gamma");
  t_nml.assignVariable(&mesh_origin_x, "mesh_origin_x");
  t_nml.assignVariable(&mesh_origin_y, "mesh_origin_y");
  t_nml.assignVariable(&mesh_origin_z, "mesh_origin_z");
  t_nml.assignVariable(&buffer_width_a, &buffer_width_b, &buffer_width_c, "buffer_width");
  t_nml.assignVariable(&buffer_width_a, "buffer_width_a");
  t_nml.assignVariable(&buffer_width_b, "buffer_width_b");
  t_nml.assignVariable(&buffer_width_c, "buffer_width_c");
  setMeshPotential(t_nml.getStringValue("potential"));
  setMeshBoundaries(t_nml.getStringValue("boundary"));
  setMeshPosition(t_nml.getStringValue("mesh_position"));
  mesh_scaling_bits = t_nml.getIntValue("scaling_bits");
}

//-------------------------------------------------------------------------------------------------
const std::string& ReceptorControls::getLabelGroup() const {
  return label_group;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReceptorControls::getAlignmentMask() const {
  return align_mask;
}

//-------------------------------------------------------------------------------------------------
RMSDMethod ReceptorControls::getAlignmentMethod() const {
  return translateRMSDMethod(align_method);
}
  
//-------------------------------------------------------------------------------------------------
int ReceptorControls::getAxisElementCount(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mesh_points_a;
  case UnitCellAxis::B:
    return mesh_points_b;
  case UnitCellAxis::C:
    return mesh_points_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ReceptorControls::getAxisElementCount(const CartesianDimension dim) const {
  return getAxisElementCount(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshSpacing(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mesh_spacing_a;
  case UnitCellAxis::B:
    return mesh_spacing_b;
  case UnitCellAxis::C:
    return mesh_spacing_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshSpacing(const CartesianDimension dim) const {
  return getMeshSpacing(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshAlpha() const {
  return mesh_alpha;
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshBeta() const {
  return mesh_beta;
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshGamma() const {
  return mesh_gamma;
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getMeshOrigin(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return mesh_origin_x;
  case CartesianDimension::Y:
    return mesh_origin_y;
  case CartesianDimension::Z:
    return mesh_origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getBufferWidth(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return buffer_width_a;
  case UnitCellAxis::B:
    return buffer_width_b;
  case UnitCellAxis::C:
    return buffer_width_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ReceptorControls::getBufferWidth(const CartesianDimension dim) const {
  return getBufferWidth(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
GridDetail ReceptorControls::getMeshContent() const {
  return translateGridDetail(potential);
}
  
//-------------------------------------------------------------------------------------------------
BoundaryCondition ReceptorControls::getMeshBoundaries() const {
  return translateBoundaryCondition(boundaries);
}

//-------------------------------------------------------------------------------------------------
MeshPosition ReceptorControls::getMeshPosition() const {
  return translateMeshPosition(mesh_position);
}

//-------------------------------------------------------------------------------------------------
int ReceptorControls::getMeshScalingBits() const {
  return mesh_scaling_bits;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& ReceptorControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setLabelGroup(const std::string &label_group_in) {
  label_group = label_group_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshElementCount(const int mesh_points_in, const UnitCellAxis dim) {
  if (mesh_points_in <= 0) {
    rtErr("Invalid mesh dimension " + std::to_string(mesh_points_in) + " along the " +
          getEnumerationName(dim) + " axis.", "ReceptorControls", "setMeshElementCount");
  }
  switch (dim) {
  case UnitCellAxis::A:
    mesh_points_a = mesh_points_in;
    break;
  case UnitCellAxis::B:
    mesh_points_b = mesh_points_in;
    break;
  case UnitCellAxis::C:
    mesh_points_c = mesh_points_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshElementCount(const int mesh_points_in,
                                           const CartesianDimension dim) {
  setMeshElementCount(mesh_points_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshSpacing(const double mesh_spacing_in, const UnitCellAxis dim) {
  if (mesh_spacing_in <= 0.0) {
    rtErr("Invalid mesh spacing " + minimalRealFormat(mesh_spacing_in, 1.0e-4) + ".",
          "ReceptorControls", "setMeshSpacing");
  }
  switch (dim) {
  case UnitCellAxis::A:
    mesh_spacing_a = mesh_spacing_in;
    break;
  case UnitCellAxis::B:
    mesh_spacing_b = mesh_spacing_in;
    break;
  case UnitCellAxis::C:
    mesh_spacing_c = mesh_spacing_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshSpacing(const double mesh_spacing_in, const CartesianDimension dim) {
  setMeshSpacing(mesh_spacing_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshAlphaAngle(const double alpha_in) {
  mesh_alpha = alpha_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshBetaAngle(const double beta_in) {
  mesh_beta = beta_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshGammaAngle(const double gamma_in) {
  mesh_gamma = gamma_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshOrigin(const double mesh_origin_in, const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    mesh_origin_x = mesh_origin_in;
    break;
  case CartesianDimension::Y:
    mesh_origin_y = mesh_origin_in;
    break;
  case CartesianDimension::Z:
    mesh_origin_z = mesh_origin_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setBufferWidth(const double buffer_width_in) {
  buffer_width_a = buffer_width_in;
  buffer_width_b = buffer_width_in;
  buffer_width_c = buffer_width_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setBufferWidth(const double buffer_width_in, const UnitCellAxis dim) {
  switch (dim) {
  case UnitCellAxis::A:
    buffer_width_a = buffer_width_in;
  case UnitCellAxis::B:
    buffer_width_b = buffer_width_in;
  case UnitCellAxis::C:
    buffer_width_c = buffer_width_in;
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setBufferWidth(const double buffer_width_in, const CartesianDimension dim) {
  setBufferWidth(buffer_width_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPotential(const std::string &potential_in) {
  potential = potential_in;
  try {
    const GridDetail interp = translateGridDetail(potential);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid potential form " + potential + " was specified.", "ReceptorControls",
            "setMeshPotential");
    case ExceptionResponse::WARN:
      rtWarn("An invalid potential form " + potential + " was specified.  The default of " +
             std::string(default_receptor_grid_potential) + " will be restored.",
             "ReceptorControls", "setMeshPotential");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    potential = std::string(default_receptor_grid_potential);
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPotential(const GridDetail potential_in) {
  potential = getEnumerationName(potential_in);
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshBoundaries(const std::string &boundaries_in) {
  boundaries = boundaries_in;
  try {
    const BoundaryCondition interp = translateBoundaryCondition(boundaries);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid boundary condition " + boundaries + " was specified.", "ReceptorControls",
            "setMeshBoundaries");
    case ExceptionResponse::WARN:
      rtWarn("An invalid boundary condition " + boundaries + " was specified.  The default of " +
             std::string(default_receptor_grid_boundary) + " will be restored.",
             "ReceptorControls", "setMeshBoundaries");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    boundaries = std::string(default_receptor_grid_boundary);
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshBoundaries(const BoundaryCondition boundaries_in) {
  boundaries = getEnumerationName(boundaries_in);
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPosition(const std::string &mesh_position_in) {
  mesh_position = mesh_position_in;
  try {
    const MeshPosition interp = translateMeshPosition(mesh_position);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid mesh positioning " + mesh_position + " was specified.", "ReceptorControls",
            "setMeshPosition");
    case ExceptionResponse::WARN:
      rtWarn("An invalid mesh positioning " + mesh_position + " was specified.  The default of " +
             std::string(default_receptor_grid_alignment) + " will be restored.",
             "ReceptorControls", "setMeshPosition");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    mesh_position = std::string(default_receptor_grid_alignment);
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPosition(const MeshPosition mesh_position_in) {
  mesh_position = getEnumerationName(mesh_position_in);
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshScalingBits(const int mesh_scaling_bits_in) {
  mesh_scaling_bits = mesh_scaling_bits_in;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator receptorInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("receptor", CaseSensitivity::AUTOMATIC, policy, "Collects details of a "
                         "receptor density field for docking and conformer selection in STORMM.");

  // Keyword: label_group
  t_nml.addKeyword("label_group", NamelistType::STRING, std::string(""));
  t_nml.addHelp("label_group", "Files management is confined to the &files namelist by specifying "
                "a label group to define the coordinates and topology (or topologies) composing "
                "the receptor.");

  // Keyword: alignment_mask
  t_nml.addKeyword("alignment_mask", NamelistType::STRING, std::string(""));
  t_nml.addHelp("alignment_mask", "When multiple structures are provided to compose the mesh, a "
                "method for aligning them is critical.  This atom mask indicates which atoms and "
                "residues of each molecule shall be compared in order to bring all structures of "
                "the receptor into the most relevant alignment.");

  // Keyword: alignment_method
  t_nml.addKeyword("alignment_method", NamelistType::STRING,
                   std::string(default_receptor_alignment_method));
  t_nml.addHelp("alignment_method", "Multiple receptor structures should be aligned to one "
                "another in order to provide the best focus on a particular binding site.  Use "
                "this keyword to indicate how to position each structure relative to an average "
                "set of positions for atoms described by alignment_mask.");
  
  // Keyword: mesh_dim{_a, _b, _c}
  t_nml.addKeyword("mesh_dim_a", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_a", "The number of mesh points to space at regular intervals along the "
                "Cartesian X axis (mesh a vector).  Negative values in this keyword imply that "
                "the dimension shall be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim_b", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_b", "The number of mesh points to space at regular intervals along the "
                "mesh b vector.  Negative values in this keyword imply that the dimension shall "
                "be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim_c", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_c", "The number of mesh points to space at regular intervals along the "
                "mesh c vector.  Negative values in this keyword imply that the dimension shall "
                "be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim", "General mesh dimension to apply along all three mesh axes.  "
                "Specifying a dimension for a particular axis will override this general "
                "directive.");

  // Keyword: mesh_spacing{_a, _b, _c}
  t_nml.addKeyword("mesh_spacing_a", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_a", "The regular spacing between mesh points along the mesh a "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing_b", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_b", "The regular spacing between mesh points along the mesh b "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing_c", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_c", "The regular spacing between mesh points along the mesh c "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing", "General mesh spacing to apply along all three mesh axes.  "
                "Specifying the spacing for a particular axis will override this general "
                "directive.");

  // Keyword: mesh_alpha
  t_nml.addKeyword("mesh_alpha", NamelistType::REAL);
  t_nml.addHelp("mesh_alpha", "Angle between the mesh b and c vectors, in units of radians.");

  // Keyword: mesh_beta
  t_nml.addKeyword("mesh_beta", NamelistType::REAL);
  t_nml.addHelp("mesh_beta", "Angle between the mesh a and c vectors, in units of radians.");

  // Keyword: mesh_gamma
  t_nml.addKeyword("mesh_gamma", NamelistType::REAL);
  t_nml.addHelp("mesh_gamma", "Angle between the mesh a and b vectors, in units of radians.");

  // Keyword: mesh_origin{_x, _y, _z}
  t_nml.addKeyword("mesh_origin_x", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_x", "Origin of the mesh along the Cartesian X axis, in units of "
                "Angstroms.");
  t_nml.addKeyword("mesh_origin_y", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_y", "Origin of the mesh along the Cartesian Y axis, in units of "
                "Angstroms.");
  t_nml.addKeyword("mesh_origin_z", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_z", "Origin of the mesh along the Cartesian Z axis, in units of "
                "Angstroms.");

  // Keyword: buffer_width{_a, _b, _c}
  t_nml.addKeyword("buffer_width_a", NamelistType::REAL);
  t_nml.addHelp("buffer_width_a", "The minimum distance between van-der Waals spheres of the "
                "receptor and the mesh boundary faces parallel to the b x c plane.");
  t_nml.addKeyword("buffer_width_b", NamelistType::REAL);
  t_nml.addHelp("buffer_width_b", "The minimum distance between van-der Waals spheres of the "
                "receptor and the mesh boundary faces parallel to the a x c plane.");
  t_nml.addKeyword("buffer_width_c", NamelistType::REAL);
  t_nml.addHelp("buffer_width_c", "The minimum distance between van-der Waals spheres of the "
                "receptor and the mesh boundary faces parallel to the a x b plane.");
  t_nml.addKeyword("buffer_width", NamelistType::REAL);
  t_nml.addHelp("buffer_width", "The minimum distance between van-der Waals spheres of the "
                "receptor and any mesh boundary faces.");
  
  // Keyword: potential
  t_nml.addKeyword("potential", NamelistType::STRING,
                   std::string(default_receptor_grid_potential));
  t_nml.addHelp("potential", "The type of non-bonded receptor potential to map onto the mesh (the "
                "mesh will express this potential in units of kcal/mol).");

  // Keyword: boundary
  t_nml.addKeyword("boundary", NamelistType::STRING,
                   std::string(default_receptor_grid_boundary));
  t_nml.addHelp("boundary", "The type of boundary conditions in which the mesh exists.  Applying "
                "periodic boundary conditions here implies that the receptor molecule will be "
                "treated as periodic in the given mesh dimensions.  Options include ISOLATED, "
                "NONPERIODIC, or NON_PERIODIC (isolated boundary conditions), and PERIODIC "
                "(periodic boundary conditions).");
  
  // Keyword: mesh_position
  t_nml.addKeyword("mesh_position", NamelistType::STRING,
                   std::string(default_receptor_grid_alignment));
  t_nml.addHelp("mesh_position", "The manner in which the receptor's potential mesh shall be "
                "aligned to the receptor molecule itself.  Aligning to the molecule centers the "
                "mesh on the receptor with equal space (or overhang) between the molecule's atoms "
                "and the planar faces of the mesh.  Options include MOLECULE or STRUCTURE (center "
                "on the molecule), ORIGIN (set the origin of the mesh to the origin of the "
                "receptor's Cartesian coordinate system), and ARBITRARY (set the origin of the "
                "mesh to some user-specified point in space).");

  // Keyword: scaling_bits
  t_nml.addKeyword("scaling_bits", NamelistType::INTEGER,
                   std::to_string(default_mesh_scaling_bits));
  t_nml.addHelp("scaling_bits", "The number of bits after the decimal to use in setting mesh "
                "vertex positions as well as computing potential quantities.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}
  
} // namespace namelist
} // namespace stormm
