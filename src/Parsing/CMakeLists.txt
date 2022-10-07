include_directories(../)

add_library(ascii_numbers ascii_numbers.cpp)
add_library(citation citation.cpp)
add_library(parse parse.cpp)
add_library(polynumeric polynumeric.cpp)
add_library(tabulation tabulation.cpp)
add_library(textfile textfile.cpp)
add_library(textguard textguard.cpp)

install(TARGETS ascii_numbers citation parse polynumeric tabulation textfile textguard
        DESTINATION lib)
install(FILES ascii_numbers.h citation.h parse.h polynumeric.h polynumeric.tpp tabulation.h
        textfile.h textguard.h DESTINATION include)
