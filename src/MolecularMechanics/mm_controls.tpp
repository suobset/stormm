// -*-c++-*-
namespace omni {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename T>
MMControlKit<T>::MMControlKit(const int step_in, const T dt_in, const T rattle_tol_in,
                              int* vwu_progress_in, int* nbwu_progress_in,
                              int* pmewu_progress_in) :
    step{step_in}, dt{dt_in}, rattle_tol{rattle_tol_in}, vwu_progress{vwu_progress_in},
    nbwu_progress{nbwu_progress_in}, pmewu_progress{pmewu_progress_in}
{}

} // namespace mm
} // namespace omni
