cmake_minimum_required(VERSION 3.0)
project(shastaDynamicExecutable)

if(BUILD_FOR_GPU)
    find_package(CUDA 10.0 REQUIRED)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
    project(shastaDynamicExecutable LANGUAGES C CXX CUDA)
endif(BUILD_FOR_GPU)

# C++ dialect.
add_definitions(-std=c++14)

# Compilation warnings.
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
    add_definitions(-fsanitize=address)
    SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
else(BUILD_DEBUG)
    add_definitions(-g0)
    add_definitions(-O3)
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
add_definitions(-DSHASTA_MARGINPHASE)
if(BUILD_FOR_GPU)
    add_definitions(-DSHASTA_BUILD_FOR_GPU)
endif(BUILD_FOR_GPU)

# Source files
file(GLOB SOURCES ../srcMain/*.cpp)

# Include directory.
include_directories(../src)

# Define our executable.
add_executable(shastaDynamicExecutable ${SOURCES})
set_target_properties(shastaDynamicExecutable PROPERTIES OUTPUT_NAME "shastaDynamic")

# This is needed to make sure the executable looks for the shared
# library in the same directory.
set_target_properties(shastaDynamicExecutable PROPERTIES INSTALL_RPATH "\$ORIGIN")

# Libraries to link with.
target_link_libraries(
    shastaDynamicExecutable
    shastaDynamicLibrary
    atomic boost_system boost_program_options boost_chrono spoa png z pthread)

# Install the dynamic executable into the bin directory.
install(TARGETS shastaDynamicExecutable DESTINATION shasta-install/bin)

# When install is complete, create the AppImage.
if(BUILD_APPIMAGE)
install(CODE "execute_process(COMMAND ${CMAKE_SOURCE_DIR}/AppImage/CreateAppImage.py shasta-install)")
endif(BUILD_APPIMAGE)

