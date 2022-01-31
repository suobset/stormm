// -*-c++-*-
#ifndef OMNI_CHEMISTRY_ENUMERATORS_H
#define OMNI_CHEMISTRY_ENUMERATORS_H

namespace omni {
namespace chemistry {

// Enumerate the chiral orientations of a center with four unique bonded groups
enum class ChiralOrientation {
  RECTUS,   // R- or D-chiral
  SINISTER, // S- or L-chiral
  NONE      // No chirality, or no preference
};

} // namespace chemistry
} // namespace omni

#endif

