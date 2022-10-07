include_directories(../)

add_library(behavior behavior.cpp)
add_library(generalized_born generalized_born.cpp)

install(TARGETS behavior generalized_born DESTINATION lib)
install(FILES behavior.h fixed_precision.h generalized_born.h hpc_bounds.h scaling.h
        symbol_values.h DESTINATION include)
