// -*-c++-*-
namespace omni {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T>
CacheResourceKit<T>::CacheResourceKit(const int max_blocks_in, const int max_atoms_in,
                                      llint* xcrd_in, llint* ycrd_in, llint* zcrd_in,
                                      llint* xvel_in, llint* yvel_in, llint* zvel_in,
                                      int* xcrd_overflow_in, int* ycrd_overflow_in,
                                      int* zcrd_overflow_in, int* xvel_overflow_in,
                                      int* yval_overflow_in, int* zvel_overflow_in,
                                      int* xfrc_overflow_in, int* yfrc_overflow_in,
                                      int* zfrc_overflow_in, T* charges_in, int* lj_idx_in) :
    max_blocks{max_blocks_in}, max_atoms{max_atoms_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xcrd_overflow{xfrc_overflow_in},
    ycrd_overflow{yfrc_overflow_in}, zcrd_overflow{zfrc_overflow_in},
    xvel_overflow{xfrc_overflow_in}, yvel_overflow{yfrc_overflow_in},
    zvel_overflow{zfrc_overflow_in}, xfrc_overflow{xfrc_overflow_in},
    yfrc_overflow{yfrc_overflow_in}, zfrc_overflow{zfrc_overflow_in}, charges{charges_in},
    lj_idx{lj_idx_in}
{}

} // namespace energy
} // namespace omni
