add_library(random random.cpp)

install(TARGETS random DESTINATION lib)
install(FILES random.h random.tpp DESTINATION include)
