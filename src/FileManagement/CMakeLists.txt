add_library(directory_util directory_util.cpp)
add_library(file_listing file_listing.cpp)
add_library(file_util file_util.cpp)

install(TARGETS directory_util file_listing file_util DESTINATION lib)
install(FILES directory_util.h file_listing.h file_util.h DESTINATION include)
