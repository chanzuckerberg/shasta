cmake_minimum_required(VERSION 3.0)

project(shastaDynamicLibrary)

# C++ dialect.
add_definitions(-std=c++14)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
    add_definitions(-fsanitize=address)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is required by dset64.hpp/dset64-gccAtomic.hpp".
add_definitions(-mcx16)

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Definitions that control what components get built.
add_definitions(-DSHASTA_HTTP_SERVER)
add_definitions(-DSHASTA_PYTHON_API)

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Include directory.
include_directories(../src)

# Define our library.
add_library(shastaDynamicLibrary SHARED ${SOURCES})

# Make sure the library is named shasta.so.
set_target_properties(shastaDynamicLibrary PROPERTIES OUTPUT_NAME "shasta")
set_target_properties(shastaDynamicLibrary PROPERTIES PREFIX "")
set_target_properties(shastaDynamicLibrary PROPERTIES DEFINE_SYMBOL "")

# Python 3 and pybind11.
# I was not able to get find_package to work with pybind11 as installed by pip3.
# The code below creates include paths for both Python 3 and pybind11.

# Python 3.8 requires the `--embed` flag to be passed in.
# find_python cmake module was introduced in v3.12, which isn't available for installation
# via `apt` on Ubuntu 18.04 LTS. So I am not using it. 
execute_process(COMMAND /usr/bin/python3-config --embed --libs OUTPUT_VARIABLE
SHASTA_PYTHON_LIBRARIES RESULT_VARIABLE COMMAND_RESULT)
if(${COMMAND_RESULT} EQUAL 0)
    # Python3 version is >= 3.8.0. So the --embed flag worked. Otherwise the command
    # would have failed with a non-zero exit code.
else()
    execute_process(COMMAND /usr/bin/python3-config --libs OUTPUT_VARIABLE SHASTA_PYTHON_LIBRARIES)
endif()
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA_PYTHON_INCLUDES)
add_definitions(${SHASTA_PYTHON_INCLUDES})
string(STRIP ${SHASTA_PYTHON_LIBRARIES} SHASTA_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${SHASTA_PYTHON_LIBRARIES}")

# Libraries to link with.
target_link_libraries(
    shastaDynamicLibrary 
     atomic png boost_program_options pthread z spoa ${SHASTA_PYTHON_LIBRARIES})

# Install the shared library into the bin directory.
install(TARGETS shastaDynamicLibrary DESTINATION shasta-install/bin)




