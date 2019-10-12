cmake_minimum_required(VERSION 3.10)

find_package(CUDA 10.0 REQUIRED)

# Use nvcc for GPU compilation
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(shastaGPULibrary LANGUAGES C CXX CUDA)

# Source files
add_library(shastaGPULibrary STATIC GPU.cu)

set_property(TARGET shastaGPULibrary PROPERTY CUDA_STANDARD 11)
