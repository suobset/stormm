add_library(code_dox code_dox.cpp)
add_library(display display.cpp)
add_library(error_format error_format.cpp)
add_library(summary_file summary_file.cpp)

install(TARGETS code_dox display error_format summary_file DESTINATION lib)
install(FILES code_dox.h display.h error_format.h summary_file.h DESTINATION include)
