// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T>
ImplicitSolventRecipe<T>::ImplicitSolventRecipe(const ImplicitSolventKit<T> &isk,
                                                const NeckGeneralizedBornKit<T> &ngbk) :
    natom{isk.natom},
    igb{isk.igb},
    table_size{ngbk.table_size},
    dielectric{isk.dielectric},
    kappa{(isk.saltcon > constants::tiny) ? sqrt(0.10806 * isk.saltcon) : 0.0},
    gb_offset{(igb == ImplicitSolventModel::NECK_GB_II) ? 0.195141 : 0.09},
    gb_neckscale{(igb == ImplicitSolventModel::NECK_GB_II) ? 0.826836 : 0.361825},
    gb_neckcut{ngbk.neck_cut},
    neck_gb_idx{isk.neck_gb_idx},
    pb_radii{isk.pb_radii},
    gb_screen{isk.gb_screen},
    gb_alpha{isk.gb_alpha},
    gb_beta{isk.gb_beta},
    gb_gamma{isk.gb_gamma},
    neck_max_sep{ngbk.max_separation},
    neck_max_val{ngbk.max_value}
{}

} // namespace energy
} // namespace omni
