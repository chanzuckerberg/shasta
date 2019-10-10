cmake_minimum_required(VERSION 3.15)
project(shastaDynamicLibrary)

# Use nvcc for GPU compilation
if(BUILD_FOR_GPU)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
    find_package(CUDA 10.0 REQUIRED)
    set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
    set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")
    project(shastaDynamicLibrary LANGUAGES C CXX CUDA)
endif(BUILD_FOR_GPU)

# C++ dialect.
#add_definitions(-std=c++14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Compilation warnings.
#add_definitions(-Wall -Wconversion -Wno-unused-result)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wconversion -Wno-unused-result")

# Optimization and debug options.
if(BUILD_DEBUG)
    #add_definitions(-ggdb3)
    #add_definitions(-O0)
    #add_definitions(-fsanitize=address)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3 -O0 -fsanitize=address")
else(BUILD_DEBUG)
    #add_definitions(-g0)
    #add_definitions(-O3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O3")
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is required by dset64.hpp/dset64-gccAtomic.hpp".
#add_definitions(-mcx16)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcx16")

# Native build.
if(BUILD_NATIVE)
    #add_definitions(-march=native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif(BUILD_NATIVE)

# NVCC/CXX flags and include dirs for GPU
if(BUILD_FOR_GPU)
    set(CMAKE_NVCC_FLAGS "${CMAKE_NVCC_FLAGS} -std=c++14 -c -Xcompiler -fPIC -shared")
    include_directories("${CUDA_INCLUDE_DIRS}")
endif(BUILD_FOR_GPU)

# Build id.
add_definitions(-DBUILD_ID=${BUILD_ID})

# Definitions needed to eliminate dependency on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Definitions that control what components get built.
add_definitions(-DSHASTA_HTTP_SERVER)
add_definitions(-DSHASTA_PYTHON_API)
add_definitions(-DSHASTA_MARGINPHASE)
if(BUILD_FOR_GPU)
    add_definitions(-DSHASTA_BUILD_FOR_GPU)
endif(BUILD_FOR_GPU)

# Sources files.
if(BUILD_FOR_GPU)
    file(GLOB SOURCES ../src/*.cpp ../src/gpu/*.cpp ../src/gpu/*.cu)
else(BUILD_FOR_GPU)
    file(GLOB SOURCES ../src/*.cpp)
endif(BUILD_FOR_GPU)

# Include directory.
include_directories(../src)

# Define our library.
add_library(shastaDynamicLibrary SHARED ${SOURCES})

# Make sure the library is named shasta.so.
set_target_properties(shastaDynamicLibrary PROPERTIES OUTPUT_NAME "shasta")
set_target_properties(shastaDynamicLibrary PROPERTIES PREFIX "")
set_target_properties(shastaDynamicLibrary PROPERTIES DEFINE_SYMBOL "")

if(BUILD_FOR_GPU)
    set_property(TARGET shastaDynamicLibrary PROPERTY CUDA_STANDARD 11)
    set_target_properties(shastaDynamicLibrary PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif(BUILD_FOR_GPU)

# Python 3 and pybind11.
# I was not able to get find_package to work with pybind11 as installed by pip3.
# The code below creates include paths for both Python 3 and pybind11.
execute_process(COMMAND python3 -m pybind11 --includes OUTPUT_VARIABLE SHASTA_PYTHON_INCLUDES)
add_definitions(${SHASTA_PYTHON_INCLUDES})
execute_process(COMMAND /usr/bin/python3-config --libs OUTPUT_VARIABLE SHASTA_PYTHON_LIBRARIES)
string(STRIP ${SHASTA_PYTHON_LIBRARIES} SHASTA_PYTHON_LIBRARIES)
SET(CMAKE_LINKER_FLAGS  "${CMAKE_LINKER_FLAGS} ${SHASTA_PYTHON_LIBRARIES}")

# Libraries to link with.
target_link_libraries(
    shastaDynamicLibrary 
    atomic png pthread z spoa MarginCore ${SHASTA_PYTHON_LIBRARIES})

# Install the shared library into the bin directory.
install(TARGETS shastaDynamicLibrary DESTINATION shasta-install/bin)




