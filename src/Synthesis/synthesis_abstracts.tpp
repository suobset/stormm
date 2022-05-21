// -*-c++-*-
namespace omni {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
SyValenceKit<T>::SyValenceKit(const int nvwu_in, const T coulomb_in, const T* bond_keq_in,
                              const T* bond_leq_in, const T* angl_keq_in, const T* angl_theta_in,
                              const T* dihe_amp_in, const T* dihe_freq_in, const T* dihe_phi_in,
                              const T* attn14_elec_in, const T* attn14_vdw_in, const T* charges_in,
                              const T* lja_14_coeff_in, const T* ljb_14_coeff_in,
                              const T* ljc_14_coeff_in, const int* lj_idx_in,
                              const int* n_lj_types_in, const int* ljabc_offsets_in,
                              const T* ubrd_keq_in, const T* ubrd_leq_in, const T* cimp_keq_in,
                              const T* cimp_phi_in, const int* cmap_dim_in,
                              const T* cmap_patches_in, const int* cmap_patch_bounds_in,
                              const int2* vwu_abstracts_in, const int* vwu_imports_in,
                              const uint2* vwu_manip_in, const uint2* cbnd_insr_in,
                              const uint2* angl_insr_in, const uint2* cdhe_insr_in,
                              const uint* cdhe_ovrt_insr_in, const uint2* cmap_insr_in,
                              const uint* infr14_insr_in, const uint* cbnd_acc_in,
                              const uint* angl_acc_in, const uint* cdhe_acc_in,
                              const uint* cmap_acc_in, const uint* infr14_acc_in) :
    nvwu{nvwu_in}, coulomb{coulomb_in}, bond_keq{bond_keq_in}, bond_leq{bond_leq_in},
    angl_keq{angl_keq_in}, angl_theta{angl_theta_in}, dihe_amp{dihe_amp_in},
    dihe_freq{dihe_freq_in}, dihe_phi{dihe_phi_in}, attn14_elec{attn14_elec_in},
    attn14_vdw{attn14_vdw_in}, charges{charges_in}, lja_14_coeff{lja_14_coeff_in},
    ljb_14_coeff{ljb_14_coeff_in}, ljc_14_coeff{ljc_14_coeff_in}, lj_idx{lj_idx_in},
    n_lj_types{n_lj_types_in}, ljabc_offsets{ljabc_offsets_in}, ubrd_keq{ubrd_keq_in},
    ubrd_leq{ubrd_leq_in}, cimp_keq{cimp_keq_in}, cimp_phi{cimp_phi_in}, cmap_dim{cmap_dim_in},
    cmap_patches{cmap_patches_in}, cmap_patch_bounds{cmap_patch_bounds_in},
    vwu_abstracts{vwu_abstracts_in}, vwu_imports{vwu_imports_in}, vwu_manip{vwu_manip_in},
    cbnd_insr{cbnd_insr_in}, angl_insr{angl_insr_in}, cdhe_insr{cdhe_insr_in},
    cdhe_ovrt_insr{cdhe_ovrt_insr_in}, cmap_insr{cmap_insr_in}, infr14_insr{infr14_insr_in},
    cbnd_acc{cbnd_acc_in}, angl_acc{angl_acc_in}, cdhe_acc{cdhe_acc_in}, cmap_acc{cmap_acc_in},
    infr14_acc{infr14_acc_in}
{}

} // namespace synthesis
} // namespace omni
