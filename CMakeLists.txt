# Minimum version required
cmake_minimum_required(VERSION 3.1)

# Set the project name
project(STORMM VERSION 0.1)

# Definitions
add_definitions(-std=c++17)

# Specify the C++ 11 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the src subdirectory, and all of its own subdirectories
message(STATUS "C++ standard is: ${CMAKE_CXX_STANDARD}")

add_subdirectory(src)
