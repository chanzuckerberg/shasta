cmake_minimum_required(VERSION 3.10)
find_package(CUDA 10.0)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")

project(shastaGpu)

set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost REQUIRED COMPONENTS system program_options chrono)

# C++ dialect.
add_definitions(-std=c++14)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
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
add_definitions(-DSHASTA_BUILD_FOR_GPU)

# Sources files.
file(GLOB GPUSRCS ../src/gpu/*.cu)
cuda_add_library(gpuLib ${GPUSRCS} STATIC)
target_include_directories(gpuLib PUBLIC ../src/gpu)

# Sources files.
file(GLOB LIBSRCS ../src/*.cpp ../src/gpu/*.cpp ../src/gpu/*.cu)
add_library(shastaStaticGpuLibrary STATIC ${LIBSRCS})
target_include_directories(shastaStaticGpuLibrary PUBLIC ../src ../src/gpu)
target_link_libraries(shastaStaticGpuLibrary
    gpuLib)

# Sources files.
file(GLOB SOURCES ../srcMain/*.cpp)

# Define our executable.
add_executable(shastaGpu ${SOURCES})

set_target_properties(shastaGpu PROPERTIES OUTPUT_NAME "shastaGpu")

# Libraries to link with.
target_link_libraries(shastaGpu
    shastaStaticGpuLibrary
    ${Boost_LIBRARIES}
    atomic /usr/local/lib/libspoa.a png z
    -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

# The static executable goes to the bin directory.
install(TARGETS shastaGpu DESTINATION shasta-install/bin)


