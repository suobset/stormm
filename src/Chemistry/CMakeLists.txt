add_library(atommask atommask.cpp)
add_library(chirality chirality.cpp)
add_library(znumber znumber.cpp)

install(TARGETS atommask chirality znumber DESTINATION lib)
install(FILES atommask.h chirality.h znumber.h DESTINATION include)
