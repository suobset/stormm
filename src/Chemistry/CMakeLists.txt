include_directories(../)

add_library(atommask atommask.cpp)
add_library(chemical_features chemical_features.cpp)
add_library(indigo indigo.cpp)
add_library(znumber znumber.cpp)

install(TARGETS atommask chemical_features indigo znumber DESTINATION lib)
install(FILES atommask.h chemistry_enumerators.h chemical_features.h indigo.h periodic_table.h
        znumber.h DESTINATION include)
