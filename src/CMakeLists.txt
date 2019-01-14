cmake_minimum_required(VERSION 3.0)
project(shasta-library)

# SET(CMAKE_VERBOSE_MAKEFILE ON)

# Compiler options.
add_definitions(-std=c++14)
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
add_library(shasta SHARED ${SOURCES})

# To make sure the library is named shasta.so,
# get rid of the "lib" prefix.
set_target_properties(shasta PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(shasta PROPERTIES  DEFINE_SYMBOL "")

# Python 3 and pybind11.
# I was not able to get find_package to work with pybind11 as installed by pip3.
# The code below creates include paths for both Python 3 and pybind11.
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA_PYTHON_INCLUDES)
add_definitions(${SHASTA_PYTHON_INCLUDES})

# Libraries to link with.
target_link_libraries(shasta atomic png pthread z spoa MarginCore)
# target_link_libraries(shasta MarginCore)

# The shared library goes to the bin directory.
install(TARGETS shasta DESTINATION shasta-install/bin)




