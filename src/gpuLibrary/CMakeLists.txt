cmake_minimum_required(VERSION 3.10)

find_package(CUDA 10.0 REQUIRED)

# Use nvcc for GPU compilation
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

add_definitions(-Xcompiler -static-libgcc)
add_definitions(-Xcompiler -static-libstdc++)
add_definitions(-dlink)

project(shastaGpuLibrary LANGUAGES C CXX CUDA)

# Source files
add_library(shastaGpuLibrary STATIC ../gpu/GPU.cu)

set_property(TARGET shastaGpuLibrary PROPERTY CUDA_STANDARD 11)
set_property(TARGET shastaGpuLibrary PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
set_property(TARGET shastaGpuLibrary PROPERTY CUDA_SEPARABLE_COMPILATION ON)
