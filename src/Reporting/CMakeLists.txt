include_directories(../)

add_library(code_dox code_dox.cpp)
add_library(display display.cpp)
add_library(error_format error_format.cpp)
add_library(reporting_enumerators reporting_enumerators.cpp)
add_library(summary_file summary_file.cpp)

install(TARGETS code_dox display error_format reporting_enumerators summary_file DESTINATION lib)
install(FILES code_dox.h display.h error_format.h reporting_enumerators.h summary_file.h
        DESTINATION include)
