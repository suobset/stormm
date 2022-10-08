include_directories(../)

add_library(line_minimization line_minimization.cpp)
add_library(minimization minimization.cpp)
add_library(mm_controls mm_controls.cpp)
add_library(mm_evaluation mm_evaluation.cpp)

install(TARGETS line_minimization minimization mm_controls mm_evaluation DESTINATION lib)
install(FILES line_minimization.h minimization.h minimization.tpp mm_controls.h mm_controls.tpp
        mm_enumerators.h mm_evaluation.h mm_evaluation.tpp DESTINATION include)

