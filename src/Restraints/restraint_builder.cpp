#include "bounded_restraint.h"
#include "restraint_apparatus.h"
#include "restraint_builder.h"

namespace omni {
namespace restraints {

//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph &ag, const CoordinateFrameReader &cframe,
                          const CoordinateFrameReader &reference_cframe, const AtomMask &mask,
                          double displacement_penalty, double displacement_onset,
                          double displacement_plateau, double proximity_penalty,
                          double proximity_onset, double proximity_plateau) {

}
  
//-------------------------------------------------------------------------------------------------
std::vector<BoundedRestraint>
applyChiralityInversionRestraints(const AtomGraph &ag, const CoordinateFrameReader &cframe,
                                  const AtomMask &chiral_mask) {

}

} // namespace restraints
} // namespace omni
