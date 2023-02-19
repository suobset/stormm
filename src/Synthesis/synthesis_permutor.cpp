#include "copyright.h"
#include "Chemistry/chemistry_enumerators.h"
#include "FileManagement/file_listing.h"
#include "Namelists/nml_conformer.h"
#include "synthesis_permutor.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using chemistry::ChiralInversionProtocol;
using chemistry::MapRotatableGroups;
using namelist::default_conf_rotation_samples;
using namelist::default_conf_cis_trans_samples;
using namelist::default_conf_rotation_set_zero;
using namelist::default_conf_rotation_set_one;
using namelist::default_conf_rotation_set_two;
using namelist::default_conf_cis_trans_set_zero;
using namelist::default_conf_cis_trans_set_one;
using diskutil::getBaseName;

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor() :
    system_count{0}, permutor_map_count{0},
    rotatable_bond_samples{default_conf_rotation_samples},
    cis_trans_bond_samples{default_conf_cis_trans_samples},
    features_objects_accessible{false},
    permutor_map_indices{HybridKind::POINTER, "synp_perm_map_idx"},
    permutor_elements{HybridKind::ARRAY, "synp_perm_elements"},
    permutor_element_bounds{HybridKind::POINTER, "synp_perm_ele_bnds"},
    system_settings{HybridKind::POINTER, "synp_sys_sett"},
    synthesis_data{HybridKind::ARRAY, "synp_synth_data"},
    rotatable_group_atoms{HybridKind::POINTER, "synp_rot_grp_atoms"},
    rotatable_group_bounds{HybridKind::POINTER, "synp_rot_grp_bnds"},
    permutor_rotatable_group_bounds{HybridKind::POINTER, "synp_prot_grp_bnds"},
    cis_trans_group_atoms{HybridKind::POINTER, "synp_ctx_grp_atoms"},
    cis_trans_group_bounds{HybridKind::POINTER, "synp_ctx_grp_bnds"},
    permutor_cis_trans_group_bounds{HybridKind::POINTER, "synp_pctx_grp_bnds"},
    invertible_group_atoms{HybridKind::POINTER, "synp_inv_grp_atoms"},
    invertible_group_bounds{HybridKind::POINTER, "synp_inv_grp_bnds"},
    permutor_invertible_group_bounds{HybridKind::POINTER, "synp_pinv_grp_bnds"},
    rotatable_bond_markers{HybridKind::POINTER, "synp_rot_bond_mark"},
    cis_trans_bond_markers{HybridKind::POINTER, "synp_ctx_bond_mark"},
    chiral_markers{HybridKind::POINTER, "synp_chiral_mark"},
    rotatable_bond_settings{HybridKind::ARRAY, "synp_rot_bond_sett"},
    cis_trans_bond_settings{HybridKind::ARRAY, "synp_ctx_bond_sett"},
    sp_rotatable_bond_settings{HybridKind::ARRAY, "synp_rot_bond_sett_sp"},
    sp_cis_trans_bond_settings{HybridKind::ARRAY, "synp_ctx_bond_sett_sp"},
    rotatable_bond_settings_bounds{HybridKind::POINTER, "synp_rot_bond_sbnd"},
    cis_trans_bond_settings_bounds{HybridKind::POINTER, "synp_ctx_bond_sbnd"},
    group_data{HybridKind::ARRAY, "synp_group_data"},
    marker_data{HybridKind::ARRAY, "synp_marker_data"},
    chiral_settings{HybridKind::ARRAY, "synp_chiral_sett"},
    chiral_settings_bounds{HybridKind::ARRAY, "synp_chiral_sbnd"},
    general_rotation_settings{}, general_cis_trans_settings{}, rotatable_groups{},
    cis_trans_groups{}, invertible_groups{}, state_trackers{}, topologies{}, features{}
{
  // Load the default settings for rotatable bonds
  general_rotation_settings.resize(rotatable_bond_samples);
  general_rotation_settings[0] = stod(std::string(default_conf_rotation_set_zero));
  general_rotation_settings[1] = stod(std::string(default_conf_rotation_set_one));
  general_rotation_settings[2] = stod(std::string(default_conf_rotation_set_two));
  general_cis_trans_settings.resize(cis_trans_bond_samples);
  general_cis_trans_settings[0] = stod(std::string(default_conf_cis_trans_set_zero));
  general_cis_trans_settings[1] = stod(std::string(default_conf_cis_trans_set_one));
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis *poly_ag, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ag->getSystemCount();
  permutor_map_count = poly_ag->getUniqueTopologyCount();
  topologies = poly_ag->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);  
  fillPermutorMaps(all_chemfe_ptr, timer);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis &poly_ag, StopWatch *timer) :
    SynthesisPermutor(poly_ag.getSelfPointer(), timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ps->getSystemCount();
  permutor_map_count = poly_ps->getUniqueTopologyCount();
  topologies = poly_ps->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);
  fillPermutorMaps(all_chemfe_ptr, timer);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps, StopWatch *timer) :
    SynthesisPermutor(poly_ps.getSelfPointer(), timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                                     const bool retain_pointers, StopWatch *timer) :
    SynthesisPermutor()
{
  permutor_map_count = chemfe_in.size();
  features_objects_accessible = retain_pointers;
  fillPermutorMaps(chemfe_in, timer);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in,
                                     const bool retain_pointers, StopWatch *timer) :
    SynthesisPermutor(temporaryFeaturesPointers(chemfe_in), retain_pointers)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const SynthesisPermutor &original) :
  system_count{original.system_count}, permutor_map_count{original.permutor_map_count},
  rotatable_bond_samples{original.rotatable_bond_samples},
  cis_trans_bond_samples{original.cis_trans_bond_samples},
  features_objects_accessible{original.features_objects_accessible},
  permutor_map_indices{original.permutor_map_indices},
  permutor_elements{original.permutor_elements},
  permutor_element_bounds{original.permutor_element_bounds},
  system_settings{original.system_settings},
  synthesis_data{original.synthesis_data},
  rotatable_group_atoms{original.rotatable_group_atoms},
  rotatable_group_bounds{original.rotatable_group_bounds},
  permutor_rotatable_group_bounds{original.permutor_rotatable_group_bounds},
  cis_trans_group_atoms{original.cis_trans_group_atoms},
  cis_trans_group_bounds{original.cis_trans_group_bounds},
  permutor_cis_trans_group_bounds{original.permutor_cis_trans_group_bounds},
  invertible_group_atoms{original.invertible_group_atoms},
  invertible_group_bounds{original.invertible_group_bounds},
  permutor_invertible_group_bounds{original.permutor_invertible_group_bounds},
  rotatable_bond_markers{original.rotatable_bond_markers},
  cis_trans_bond_markers{original.cis_trans_bond_markers},
  chiral_markers{original.chiral_markers},
  rotatable_bond_settings{original.rotatable_bond_settings},
  cis_trans_bond_settings{original.cis_trans_bond_settings},
  sp_rotatable_bond_settings{original.sp_rotatable_bond_settings},
  sp_cis_trans_bond_settings{original.sp_cis_trans_bond_settings},
  rotatable_bond_settings_bounds{original.rotatable_bond_settings_bounds},
  cis_trans_bond_settings_bounds{original.cis_trans_bond_settings_bounds},
  group_data{original.group_data},
  marker_data{original.marker_data},
  chiral_settings{original.chiral_settings},
  chiral_settings_bounds{original.chiral_settings_bounds},
  general_rotation_settings{original.general_rotation_settings},
  general_cis_trans_settings{original.general_cis_trans_settings},
  rotatable_groups{original.rotatable_groups},
  cis_trans_groups{original.cis_trans_groups},
  invertible_groups{original.invertible_groups},
  state_trackers{original.state_trackers},
  topologies{original.topologies},
  features{original.features}
{
  // Rebase the POINTER-kind Hybrid objects 
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(SynthesisPermutor &&original) :
  system_count{original.system_count}, permutor_map_count{original.permutor_map_count},
  rotatable_bond_samples{original.rotatable_bond_samples},
  cis_trans_bond_samples{original.cis_trans_bond_samples},
  features_objects_accessible{original.features_objects_accessible},
  permutor_map_indices{std::move(original.permutor_map_indices)},
  permutor_elements{std::move(original.permutor_elements)},
  permutor_element_bounds{std::move(original.permutor_element_bounds)},
  system_settings{std::move(original.system_settings)},
  synthesis_data{std::move(original.synthesis_data)},
  rotatable_group_atoms{std::move(original.rotatable_group_atoms)},
  rotatable_group_bounds{std::move(original.rotatable_group_bounds)},
  permutor_rotatable_group_bounds{std::move(original.permutor_rotatable_group_bounds)},
  cis_trans_group_atoms{std::move(original.cis_trans_group_atoms)},
  cis_trans_group_bounds{std::move(original.cis_trans_group_bounds)},
  permutor_cis_trans_group_bounds{std::move(original.permutor_cis_trans_group_bounds)},
  invertible_group_atoms{std::move(original.invertible_group_atoms)},
  invertible_group_bounds{std::move(original.invertible_group_bounds)},
  permutor_invertible_group_bounds{std::move(original.permutor_invertible_group_bounds)},
  rotatable_bond_markers{std::move(original.rotatable_bond_markers)},
  cis_trans_bond_markers{std::move(original.cis_trans_bond_markers)},
  chiral_markers{std::move(original.chiral_markers)},
  rotatable_bond_settings{std::move(original.rotatable_bond_settings)},
  cis_trans_bond_settings{std::move(original.cis_trans_bond_settings)},
  sp_rotatable_bond_settings{std::move(original.sp_rotatable_bond_settings)},
  sp_cis_trans_bond_settings{std::move(original.sp_cis_trans_bond_settings)},
  rotatable_bond_settings_bounds{std::move(original.rotatable_bond_settings_bounds)},
  cis_trans_bond_settings_bounds{std::move(original.cis_trans_bond_settings_bounds)},
  group_data{std::move(original.group_data)},
  marker_data{std::move(original.marker_data)},
  chiral_settings{std::move(original.chiral_settings)},
  chiral_settings_bounds{std::move(original.chiral_settings_bounds)},
  general_rotation_settings{std::move(original.general_rotation_settings)},
  general_cis_trans_settings{std::move(original.general_cis_trans_settings)},
  rotatable_groups{std::move(original.rotatable_groups)},
  cis_trans_groups{std::move(original.cis_trans_groups)},
  invertible_groups{std::move(original.invertible_groups)},
  state_trackers{std::move(original.state_trackers)},
  topologies{std::move(original.topologies)},
  features{std::move(original.features)}
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor& SynthesisPermutor::operator=(const SynthesisPermutor &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all member variables
  system_count = other.system_count;
  permutor_map_count = other.permutor_map_count;
  rotatable_bond_samples = other.rotatable_bond_samples;
  cis_trans_bond_samples = other.cis_trans_bond_samples;
  features_objects_accessible = other.features_objects_accessible;
  permutor_map_indices = other.permutor_map_indices;
  permutor_elements = other.permutor_elements;
  permutor_element_bounds = other.permutor_element_bounds;
  system_settings = other.system_settings;
  synthesis_data = other.synthesis_data;
  rotatable_group_atoms = other.rotatable_group_atoms;
  rotatable_group_bounds = other.rotatable_group_bounds;
  permutor_rotatable_group_bounds = other.permutor_rotatable_group_bounds;
  cis_trans_group_atoms = other.cis_trans_group_atoms;
  cis_trans_group_bounds = other.cis_trans_group_bounds;
  permutor_cis_trans_group_bounds = other.permutor_cis_trans_group_bounds;
  invertible_group_atoms = other.invertible_group_atoms;
  invertible_group_bounds = other.invertible_group_bounds;
  permutor_invertible_group_bounds = other.permutor_invertible_group_bounds;
  rotatable_bond_markers = other.rotatable_bond_markers;
  cis_trans_bond_markers = other.cis_trans_bond_markers;
  chiral_markers = other.chiral_markers;
  rotatable_bond_settings = other.rotatable_bond_settings;
  cis_trans_bond_settings = other.cis_trans_bond_settings;
  sp_rotatable_bond_settings = other.sp_rotatable_bond_settings;
  sp_cis_trans_bond_settings = other.sp_cis_trans_bond_settings;
  rotatable_bond_settings_bounds = other.rotatable_bond_settings_bounds;
  cis_trans_bond_settings_bounds = other.cis_trans_bond_settings_bounds;
  group_data = other.group_data;
  marker_data = other.marker_data;
  chiral_settings = other.chiral_settings;
  chiral_settings_bounds = other.chiral_settings_bounds;
  general_rotation_settings = other.general_rotation_settings;
  general_cis_trans_settings = other.general_cis_trans_settings;
  rotatable_groups = other.rotatable_groups;
  cis_trans_groups = other.cis_trans_groups;
  invertible_groups = other.invertible_groups;
  state_trackers = other.state_trackers;
  topologies = other.topologies;
  features = other.features;

  // Rebase pointers and return
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor& SynthesisPermutor::operator=(SynthesisPermutor &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all member variables
  system_count = other.system_count;
  permutor_map_count = other.permutor_map_count;
  rotatable_bond_samples = other.rotatable_bond_samples;
  cis_trans_bond_samples = other.cis_trans_bond_samples;
  features_objects_accessible = other.features_objects_accessible;
  permutor_map_indices = std::move(other.permutor_map_indices);
  permutor_elements = std::move(other.permutor_elements);
  permutor_element_bounds = std::move(other.permutor_element_bounds);
  system_settings = std::move(other.system_settings);
  synthesis_data = std::move(other.synthesis_data);
  rotatable_group_atoms = std::move(other.rotatable_group_atoms);
  rotatable_group_bounds = std::move(other.rotatable_group_bounds);
  permutor_rotatable_group_bounds = std::move(other.permutor_rotatable_group_bounds);
  cis_trans_group_atoms = std::move(other.cis_trans_group_atoms);
  cis_trans_group_bounds = std::move(other.cis_trans_group_bounds);
  permutor_cis_trans_group_bounds = std::move(other.permutor_cis_trans_group_bounds);
  invertible_group_atoms = std::move(other.invertible_group_atoms);
  invertible_group_bounds = std::move(other.invertible_group_bounds);
  permutor_invertible_group_bounds = std::move(other.permutor_invertible_group_bounds);
  rotatable_bond_markers = std::move(other.rotatable_bond_markers);
  cis_trans_bond_markers = std::move(other.cis_trans_bond_markers);
  chiral_markers = std::move(other.chiral_markers);
  rotatable_bond_settings = std::move(other.rotatable_bond_settings);
  cis_trans_bond_settings = std::move(other.cis_trans_bond_settings);
  sp_rotatable_bond_settings = std::move(other.sp_rotatable_bond_settings);
  sp_cis_trans_bond_settings = std::move(other.sp_cis_trans_bond_settings);
  rotatable_bond_settings_bounds = std::move(other.rotatable_bond_settings_bounds);
  cis_trans_bond_settings_bounds = std::move(other.cis_trans_bond_settings_bounds);
  group_data = std::move(other.group_data);
  marker_data = std::move(other.marker_data);
  chiral_settings = std::move(other.chiral_settings);
  chiral_settings_bounds = std::move(other.chiral_settings_bounds);
  general_rotation_settings = std::move(other.general_rotation_settings);
  general_cis_trans_settings = std::move(other.general_cis_trans_settings);
  rotatable_groups = std::move(other.rotatable_groups);
  cis_trans_groups = std::move(other.cis_trans_groups);
  invertible_groups = std::move(other.invertible_groups);
  state_trackers = std::move(other.state_trackers);
  topologies = std::move(other.topologies);
  features = std::move(other.features);

  // As with other move assignment operations, rebasing is not needed.  Simply return the result.
  return *this;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorSetCount() const {
  return permutor_map_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getRotatableBondSampleCount() const {
  return rotatable_bond_samples;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getCisTransBondSampleCount() const {
  return rotatable_bond_samples;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisPermutor::getElementSampleCounts(int system_index) const {
  validateSystemIndex(system_index);
  const int map_idx = permutor_map_indices.readHost(system_index);
  return state_trackers[map_idx].getStateLimits();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisPermutor::getElementSampleCounts(const AtomGraph *query_ag) const {
  for (int i = 0; i < permutor_map_count; i++) {
    if (topologies[i] == query_ag) {
      return state_trackers[i].getStateLimits();
    }
  }
  rtErr("No topology match was found for an object originating in file " +
        getBaseName(query_ag->getFileName()) + ".", "SynthesisPermutor", "getElementSampleCounts");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisPermutor::getElementSampleCounts(const AtomGraph &query_ag) const {
  return getElementSampleCounts(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
const SyPermutorKit<double> SynthesisPermutor::dpData(const HybridTargetLevel tier) const {
  return SyPermutorKit<double>(system_count, permutor_map_count, permutor_map_indices.data(tier),
                               permutor_elements.data(tier), permutor_element_bounds.data(tier),
                               system_settings.data(tier), rotatable_group_atoms.data(tier),
                               rotatable_group_bounds.data(tier),
                               permutor_rotatable_group_bounds.data(tier),
                               cis_trans_group_atoms.data(tier), cis_trans_group_bounds.data(tier),
                               permutor_cis_trans_group_bounds.data(tier),
                               invertible_group_atoms.data(tier),
                               invertible_group_bounds.data(tier),
                               permutor_invertible_group_bounds.data(tier),
                               rotatable_bond_markers.data(tier),
                               cis_trans_bond_markers.data(tier), chiral_markers.data(tier),
                               rotatable_bond_settings.data(tier),
                               cis_trans_bond_settings.data(tier), 
                               rotatable_bond_settings_bounds.data(tier),
                               cis_trans_bond_settings_bounds.data(tier),
                               chiral_settings.data(tier), chiral_settings_bounds.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::upload() {
  permutor_elements.upload();
  synthesis_data.upload();
  rotatable_bond_settings.upload();
  cis_trans_bond_settings.upload();
  sp_rotatable_bond_settings.upload();
  sp_cis_trans_bond_settings.upload();
  group_data.upload();
  marker_data.upload();
  chiral_settings.upload();
  chiral_settings_bounds.upload();
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::download() {
  permutor_elements.download();
  synthesis_data.download();
  rotatable_bond_settings.download();
  cis_trans_bond_settings.download();
  sp_rotatable_bond_settings.download();
  sp_cis_trans_bond_settings.download();
  group_data.download();
  marker_data.download();
  chiral_settings.download();
  chiral_settings_bounds.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::rebasePointers() {
  permutor_map_indices.swapTarget(&synthesis_data);
  permutor_element_bounds.swapTarget(&synthesis_data);
  system_settings.swapTarget(&synthesis_data);
  rotatable_group_atoms.swapTarget(&group_data);
  rotatable_group_bounds.swapTarget(&group_data);
  permutor_rotatable_group_bounds.swapTarget(&group_data);  
  cis_trans_group_atoms.swapTarget(&group_data);
  cis_trans_group_bounds.swapTarget(&group_data);
  permutor_cis_trans_group_bounds.swapTarget(&group_data);  
  invertible_group_atoms.swapTarget(&group_data);
  invertible_group_bounds.swapTarget(&group_data);
  permutor_invertible_group_bounds.swapTarget(&group_data);  
  rotatable_bond_markers.swapTarget(&marker_data);
  cis_trans_bond_markers.swapTarget(&marker_data);
  chiral_markers.swapTarget(&marker_data);
  rotatable_bond_settings_bounds.swapTarget(&group_data);
  cis_trans_bond_settings_bounds.swapTarget(&group_data);
}

//------------------------------------------------------------------------------------------------
std::vector<ChemicalFeatures> SynthesisPermutor::temporaryFeatures(StopWatch *timer) const {
  std::vector<ChemicalFeatures> result(permutor_map_count);
  result.reserve(topologies.size());
  for (int i = 0; i < permutor_map_count; i++) {
    result.emplace_back(ChemicalFeatures(topologies[i], MapRotatableGroups::YES, 300.0, timer));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<ChemicalFeatures*>
SynthesisPermutor::temporaryFeaturesPointers(const std::vector<ChemicalFeatures> &feat_in) const {
  std::vector<ChemicalFeatures*> result(permutor_map_count);
  for (int i = 0; i < permutor_map_count; i++) {

    // Even passing the vector by reference, pointers to the individual elements will be valid.
    // A pointer to the vector itself would not.
    result[i] = const_cast<ChemicalFeatures*>(&feat_in[i]);
  }  
  return result;
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::fillPermutorMaps(const std::vector<ChemicalFeatures*> &chemfe_in,
                                         StopWatch *timer) {

  // Check that, for each set of features, rotatable bond groups have been mapped.  If not, create
  // a new ChemicalFeatures object which does map the rotating groups (the input features are taken
  // as const).  
  features.resize(permutor_map_count);
  std::vector<ChemicalFeatures> fallback_features;
  int n_unmapped = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    n_unmapped += (chemfe_in[i]->rotatableGroupsMapped() == false);
  }
  fallback_features.reserve(n_unmapped);
  n_unmapped = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    if (chemfe_in[i]->rotatableGroupsMapped() == false) {
      fallback_features.emplace_back(chemfe_in[i]->getTopologyPointer(), MapRotatableGroups::YES,
                                     300.0, timer);
      features[i] = &fallback_features[n_unmapped];
      n_unmapped++;
    }
    else {
      features[i] = chemfe_in[i];
    }
  }

  // If any features had to be recalculated in order to get the maps, the temporary objects will
  // not persist and the feature pointers as a whole should be considered unreliable.
  if (n_unmapped > 0) {
    features_objects_accessible = false;
  }

  // Obtain the isomerization groups
  std::vector<int> tmp_permutor_rotatable_group_bounds(permutor_map_count + 1);
  std::vector<int> tmp_permutor_cis_trans_group_bounds(permutor_map_count + 1);
  std::vector<int> tmp_permutor_invertible_group_bounds(permutor_map_count + 1);
  int n_rotg = 0;
  int n_ctxg = 0;
  int n_invg = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    tmp_permutor_rotatable_group_bounds[i] = n_rotg;
    tmp_permutor_cis_trans_group_bounds[i] = n_ctxg;
    tmp_permutor_invertible_group_bounds[i] = n_invg;
    n_rotg += features[i]->getRotatableBondCount();
    n_ctxg += features[i]->getCisTransBondCount();
    n_invg += features[i]->getChiralCenterCount();
  }
  tmp_permutor_rotatable_group_bounds[permutor_map_count] = n_rotg;
  tmp_permutor_cis_trans_group_bounds[permutor_map_count] = n_ctxg;
  tmp_permutor_invertible_group_bounds[permutor_map_count] = n_invg;
  rotatable_groups.reserve(n_rotg);
  cis_trans_groups.reserve(n_ctxg);
  invertible_groups.reserve(n_invg);
  std::vector<int4> tmp_rotatable_bond_markers(n_rotg);
  std::vector<int4> tmp_cis_trans_bond_markers(n_ctxg);
  std::vector<int4> tmp_chiral_markers(n_invg);
  n_rotg = 0;
  n_ctxg = 0;
  n_invg = 0;  
  for (int i = 0; i < permutor_map_count; i++) {
    const int ni_rotg = features[i]->getRotatableBondCount();
    const int ni_ctxg = features[i]->getCisTransBondCount();
    const int ni_invg = features[i]->getChiralCenterCount();
    const std::vector<IsomerPlan> rotg_vec = features[i]->getRotatableBondGroups();
    for (int j = 0; j < ni_rotg; j++) {
      rotatable_groups.push_back(rotg_vec[j]);
      tmp_rotatable_bond_markers[n_rotg] = { rotg_vec[j].getRootHandle(),
                                             rotg_vec[j].getRootAtom(), rotg_vec[j].getPivotAtom(),
                                             rotg_vec[j].getPivotHandle() };
      n_rotg++;
    }
    const std::vector<IsomerPlan> ctxg_vec = features[i]->getCisTransIsomerizationGroups();
    for (int j = 0; j < ni_ctxg; j++) {
      cis_trans_groups.push_back(ctxg_vec[j]);
      tmp_cis_trans_bond_markers[n_rotg] = { rotg_vec[j].getRootHandle(),
                                             rotg_vec[j].getRootAtom(), rotg_vec[j].getPivotAtom(),
                                             rotg_vec[j].getPivotHandle() };
      n_ctxg++;
    }
    const std::vector<IsomerPlan> invg_vec = features[i]->getChiralInversionGroups();
    const std::vector<int4> chiral_base_atoms = features[i]->getChiralArmBaseAtoms();
    for (int j = 0; j < ni_invg; j++) {
      invertible_groups.push_back(invg_vec[j]);
      tmp_chiral_markers[n_rotg] = chiral_base_atoms[j];
      n_invg++;
    }
  }
  
  // Count the numbers of atoms in each distinct system's rotatable bond groups.  Track to ensure
  // that the group mappings do not exceed a maximum of two billion (2^31) atoms, and raise an
  // exception if they do.  The int-format mapping vectors are limited, but any problem set that
  // hits these limits is almost certain to run up against other barriers like the memory
  // available to the card.
  bool problem = false;
  std::vector<int> tmp_rotatable_group_bounds(n_rotg + 1);
  std::vector<int> tmp_cis_trans_group_bounds(n_ctxg + 1);
  std::vector<int> tmp_invertible_group_bounds(n_invg + 1);
  int rotg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_rotatable_group_bounds[permutor_map_count]; i++) {
    tmp_rotatable_group_bounds[i] = rotg_atom_count;
    rotg_atom_count += rotatable_groups[i].getMovingAtomCount() + 2;
    problem = (problem || rotg_atom_count < 0);
  }
  int ctxg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_cis_trans_group_bounds[permutor_map_count]; i++) {
    tmp_cis_trans_group_bounds[i] = ctxg_atom_count;
    ctxg_atom_count += cis_trans_groups[i].getMovingAtomCount() + 2;
    problem = (problem || ctxg_atom_count < 0);
  }
  int invg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_invertible_group_bounds[permutor_map_count]; i++) {
    tmp_invertible_group_bounds[i] = invg_atom_count;

    // For chiral inversions, rather than specifying the root and atom and pivot atom, specify
    // the central atom and cast the Chiral inversion protocol into an int right behind the central
    // atom index.
    switch (invertible_groups[i].getChiralPlan()) {
    case ChiralInversionProtocol::ROTATE:
      invg_atom_count += invertible_groups[i].getMovingAtomCount() + 2;
      break;
    case ChiralInversionProtocol::REFLECT:
    case ChiralInversionProtocol::DO_NOT_INVERT:
      invg_atom_count += 2;
      break;
    }
    problem = (problem || invg_atom_count < 0);
  }
  tmp_rotatable_group_bounds[n_rotg] = rotg_atom_count;
  tmp_cis_trans_group_bounds[n_ctxg] = ctxg_atom_count;
  tmp_invertible_group_bounds[n_invg] = invg_atom_count;
  if (problem) {
    rtErr("Rotatable, cis-trans, or invertible atom groups specify too many atoms for available "
          "indexing.  Reduce the number of systems if possible.", "SynthesisPermutor",
          "FillPermutorMaps");
  }
  
  // Allocate temporary storage arrays for the moving atom groups, then fill them.
  std::vector<int> tmp_rotatable_group_atoms(rotg_atom_count);
  std::vector<int> tmp_cis_trans_group_atoms(ctxg_atom_count);
  std::vector<int> tmp_invertible_group_atoms(invg_atom_count);
  rotg_atom_count = 0;
  ctxg_atom_count = 0;
  invg_atom_count = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    const std::vector<IsomerPlan> rotg_vec = features[i]->getRotatableBondGroups();
    const int ni_rotg = rotg_vec.size();
    for (int j = 0; j < ni_rotg; j++) {
      tmp_rotatable_group_atoms[rotg_atom_count] = rotg_vec[j].getRootAtom();
      rotg_atom_count++;
      tmp_rotatable_group_atoms[rotg_atom_count] = rotg_vec[j].getPivotAtom();
      rotg_atom_count++;
      const std::vector<int>& mv_atoms = rotg_vec[j].getMovingAtoms();
      const int nij_atoms = rotg_vec[j].getMovingAtomCount();
      for (int k = 0; k < nij_atoms; k++) {
        tmp_rotatable_group_atoms[rotg_atom_count] = mv_atoms[k];
        rotg_atom_count++;
      }
      rotg_atom_count += nij_atoms;
    }
  }
  for (int i = 0; i < permutor_map_count; i++) {
    const std::vector<IsomerPlan> ctxg_vec = features[i]->getCisTransIsomerizationGroups();
    const int ni_ctxg = ctxg_vec.size();
    for (int j = 0; j < ni_ctxg; j++) {
      tmp_cis_trans_group_atoms[ctxg_atom_count] = ctxg_vec[j].getRootAtom();
      ctxg_atom_count++;
      tmp_cis_trans_group_atoms[ctxg_atom_count] = ctxg_vec[j].getPivotAtom();
      ctxg_atom_count++;
      const std::vector<int>& mv_atoms = ctxg_vec[j].getMovingAtoms();
      const int nij_atoms = ctxg_vec[j].getMovingAtomCount();
      for (int k = 0; k < nij_atoms; k++) {
        tmp_cis_trans_group_atoms[ctxg_atom_count] = mv_atoms[k];
        ctxg_atom_count++;
      }
      ctxg_atom_count += nij_atoms;
    }
  }
  for (int i = 0; i < permutor_map_count; i++) {
    const std::vector<IsomerPlan> invg_vec = features[i]->getChiralInversionGroups();    
    const int ni_invg = invg_vec.size();
    for (int j = 0; j < ni_invg; j++) {
      tmp_invertible_group_atoms[invg_atom_count] = invg_vec[j].getRootAtom();
      invg_atom_count++;
      tmp_invertible_group_atoms[invg_atom_count] = static_cast<int>(invg_vec[j].getChiralPlan());
      invg_atom_count++;
      switch (invg_vec[j].getChiralPlan()) {
      case ChiralInversionProtocol::ROTATE:
        {
          const std::vector<int>& mv_atoms = invg_vec[j].getMovingAtoms();
          const int nij_atoms = invg_vec[j].getMovingAtomCount();
          for (int k = 0; k < nij_atoms; k++) {
            tmp_invertible_group_atoms[invg_atom_count] = mv_atoms[k];
            invg_atom_count++;
          }
          invg_atom_count += nij_atoms;
        }
        break;
      case ChiralInversionProtocol::REFLECT:
      case ChiralInversionProtocol::DO_NOT_INVERT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::setVariableRanges(const ConformerControls &confcon) {

  // Obtain rotatable and cis-trans bond settings.  Set pointers to host data for convenience.
  const std::vector<double>& rotb_nml_settings = confcon.getRotationSampleValues();
  const std::vector<double>& ctxb_nml_settings = confcon.getCisTransSampleValues();
  double* rotb_settings_ptr = rotatable_bond_settings.data();
  double* ctxb_settings_ptr = cis_trans_bond_settings.data();
  float* sp_rotb_settings_ptr = sp_rotatable_bond_settings.data();
  float* sp_ctxb_settings_ptr = sp_cis_trans_bond_settings.data();
  
  // Lay out the state trackers
  for (int i = 0; i < permutor_map_count; i++) {
    int llim = permutor_rotatable_group_bounds.readHost(i);
    int hlim = permutor_rotatable_group_bounds.readHost(i + 1);
    for (int j = llim; j < hlim; j++) {
      rotb_settings_ptr[j]    = rotb_nml_settings[j - llim];
      sp_rotb_settings_ptr[j] = rotb_nml_settings[j - llim];
    }
    llim = permutor_cis_trans_group_bounds.readHost(i);
    hlim = permutor_cis_trans_group_bounds.readHost(i + 1);
    for (int j = llim; j < hlim; j++) {
      ctxb_settings_ptr[j]    = ctxb_nml_settings[j - llim];
      sp_ctxb_settings_ptr[j] = ctxb_nml_settings[j - llim];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::validateSystemIndex(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("System index " + std::to_string(system_index) + " is out of bounds for a collection "
          "of " + std::to_string(system_count) + " systems.", "SynthesisPermutor",
          "getElementSampleCounts");
  }
}

} // namespace synthesis
} // namespace stormm
