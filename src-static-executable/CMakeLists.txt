# This creates a Shasta static executable that provides
# basic functionality and reduced performance. 
# For full functionality use the shared library built
# under directory src.

cmake_minimum_required(VERSION 3.0)
project(shasta-static-executable)

file(GLOB SOURCES main.cpp ../src/*.cpp)
include_directories(../src)

add_executable(shasta-static-executable ${SOURCES})
set_target_properties(shasta-static-executable PROPERTIES LINK_FLAGS "-static" )
set_target_properties(shasta-static-executable PROPERTIES OUTPUT_NAME "shasta")


# Compiler options.
add_definitions(-std=c++14)
add_definitions(-ggdb3)
add_definitions(-O3 -msse4.2)
add_definitions(-Wall -Wconversion -Wno-unused-result)

# This is required by dset64.hpp.
# It allows std::atomic<__uint128_t>::is_lock_free()
# to return true. However this does not happen in all gcc versions.
add_definitions(-mcx16)

# This is used in the C++ code to ifdef out
# code that should not be included in the static executable.
add_definitions(-DSHASTA_STATIC_EXECUTABLE)

# Specify the libraries.
target_link_libraries(shasta-static-executable pthread boost_system spoa)

# The static executable goes to the bin directory.
install(TARGETS shasta-static-executable DESTINATION shasta-install/bin)


