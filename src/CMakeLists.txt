cmake_minimum_required(VERSION 2.6)
project(Nanopore1)

SET(CMAKE_VERBOSE_MAKEFILE ON)

# Compiler options.
add_definitions(-std=c++11)
add_definitions(-ggdb3)
add_definitions(-O3 -msse4.2)
add_definitions(-Wall -Wconversion -Wno-unused-result)

# This is required by dset64.hpp.
# It allows std::atomic<__uint128_t>::is_lock_free()
# to return true.
add_definitions(-mcx16)

# Definitions needed to eliminate runtime dependency 
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

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
target_link_libraries(Nanopore2 png)



