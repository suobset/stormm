add_library(matrix_ops matrix_ops.cpp)
add_library(reduction reduction.cpp)
add_library(rounding rounding.cpp)
add_library(series_ops series_ops.cpp)
add_library(statistics statistics.cpp)
add_library(summation summation.cpp)
add_library(tickcounter tickcounter.cpp)
add_library(vector_ops vector_ops.cpp)

install(TARGETS matrix_ops reduction rounding series_ops statistics summation tickcounter
        vector_ops DESTINATION lib)
install(FILES matrix.h matrix_ops.h matrix_ops.tpp matrix.tpp multiplication.h multiplication.tpp
        reduction_abstracts.h reduction_abstracts.tpp reduction_bridge.h reduction_enumerators.h
        reduction.h reduction.tpp reduction_workunit.h rounding.h rounding.tpp series_ops.h
        series_ops.tpp sorting.h statistics.h summation.h summation.tpp tickcounter.h
        tickcounter.tpp vector_ops.h vector_ops.tpp DESTINATION include)
