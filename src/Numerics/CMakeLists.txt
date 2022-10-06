add_library(split_fixed_precision split_fixed_precision.cpp)

install(TARGETS split_fixed_precision DESTINATION lib)
install(FILES split_fixed_precision.h DESTINATION include)
