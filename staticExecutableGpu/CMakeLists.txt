cmake_minimum_required(VERSION 3.10)

set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

project(shastaGpu LANGUAGES CXX CUDA)

set(Boost_USE_STATIC_LIBS   ON)
find_package(Boost REQUIRED COMPONENTS system program_options chrono)

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11")
# C++ dialect.
SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# Optimization and debug options.
if(BUILD_DEBUG)
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb3")
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
else(BUILD_DEBUG)
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0")
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(BUILD_DEBUG)

# 16-byte compare and swap.
# This is required by dset64.hpp/dset64-gccAtomic.hpp".
SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcx16")

# Native build.
if(BUILD_NATIVE)
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
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
add_library(gpuLib STATIC ${GPUSRCS})
target_include_directories(gpuLib PUBLIC ../src/gpu)

# Sources files.
file(GLOB LIBSRCS ../src/*.cpp ../src/gpu/*.cpp ../src/gpu/*.cu)
add_library(shastaStaticGpuLibrary STATIC ${LIBSRCS})
set_property(TARGET shastaStaticGpuLibrary PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES "CUDA")
target_include_directories(shastaStaticGpuLibrary PUBLIC ../src ../src/gpu)
target_link_libraries(shastaStaticGpuLibrary
    gpuLib)

# Sources files.
file(GLOB SOURCES ../srcMain/*.cpp)

# Define our executable.
add_executable(shastaGpu ${SOURCES})

set_target_properties(shastaGpu PROPERTIES OUTPUT_NAME "shastaGpu")

#if(NOT CMAKE_CUDA_DEVICE_LINK_EXECUTABLE)
#  set(CMAKE_CUDA_DEVICE_LINK_EXECUTABLE
#    "<CMAKE_CUDA_COMPILER> ${CMAKE_CUDA_HOST_FLAGS} <FLAGS> ${CMAKE_CUDA_COMPILE_OPTIONS_PIC} ${_CMAKE_CUDA_EXTRA_DEVICE_LINK_FLAGS} -shared -dlink <OBJECTS> -o <TARGET> <LINK_LIBRARIES>${__IMPLICT_DLINK_FLAGS}")
#endif()

# Libraries to link with.
target_link_libraries(shastaGpu
    shastaStaticGpuLibrary
    ${Boost_LIBRARIES}
    atomic /usr/local/lib/libspoa.a png z
    -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)

# The static executable goes to the bin directory.
if(NOT SYSTEM_INSTALL)
  install(TARGETS shastaGpu DESTINATION shasta-install/bin)
else()
  install(TARGETS shastaGpu DESTINATION bin)
endif()


