cmake_minimum_required(VERSION 2.6)
project(Nanopore1)

SET(CMAKE_VERBOSE_MAKEFILE ON)

# Compiler options.
add_definitions(-std=c++11)
add_definitions(-ggdb3)
add_definitions(-O3 -msse4.2)
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Define our shared library.
file(GLOB SOURCES *.cpp)
add_library(Nanopore2 SHARED ${SOURCES})

# To make sure the library is named Nanopore2.so,
# get rid of the "lib" prefix.
set_target_properties(Nanopore2 PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(Nanopore2 PROPERTIES  DEFINE_SYMBOL "")

# Include directories.
include_directories(/usr/include/python3.5m)
include_directories(/usr/local/include/python3.5)

# Libraries to link with.
# Currently none needed. Uncomment the next line and 
# add them here, before the closed parenthesis.
# target_link_libraries(Nanopore2)



