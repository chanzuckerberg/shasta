cmake_minimum_required(VERSION 3.16)

project(shastaDynamicLibrary)

# C++ dialect.
add_definitions(-std=c++20)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result -Wno-trigraphs)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
    # Address sanitizer does not work well with the Python API,
    # so it is commented out here. But it could be useful
    # to debug the dynamic executable. 
    # add_definitions(-fsanitize=address)
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
    # NDEBUG is required to turn off SeqAn debug code.
    add_definitions(-DNDEBUG)
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is recommended for dset64.hpp/dset64-gccAtomic.hpp".
# It's available only on x86 architectures.
if(X86_64)
  add_definitions(-mcx16)
endif(X86_64)

# Native build.
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# This is needed to avoid some Boost warnings about deprecated messages.
# We canot fix this in Shasta as it is a Boost problem.
# It can be removed if Boost is fixed.
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

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
execute_process(COMMAND /usr/bin/python3-config --embed --libs OUTPUT_VARIABLE SHASTA_PYTHON_LIBRARIES)
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA_PYTHON_INCLUDES)
add_definitions(${SHASTA_PYTHON_INCLUDES})
string(STRIP ${SHASTA_PYTHON_LIBRARIES} SHASTA_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${SHASTA_PYTHON_LIBRARIES}")



# Libraries to link with.
target_link_libraries(
    shastaDynamicLibrary 
     atomic png boost_program_options pthread z spoa lapack blas ${SHASTA_PYTHON_LIBRARIES})

# Install the shared library into the bin directory.
install(TARGETS shastaDynamicLibrary DESTINATION shasta-install/bin)




