cmake_minimum_required(VERSION 3.0)
project(shastaStaticLibrary)

# C++ dialect.
add_definitions(-std=c++14)

# Compilation warnings.
if(NOT MACOS)
    add_definitions(-Wall -Wconversion -Wno-unused-result)
endif(NOT MACOS)

# Optimization and debug options.
if(BUILD_DEBUG)
    add_definitions(-ggdb3)
    add_definitions(-O0)
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

# Definitions that control what components get built.
add_definitions(-DSHASTA_HTTP_SERVER)

# Sources files.
file(GLOB SOURCES ../src/*.cpp)

# Include directory.
include_directories(../src)

# Include directories required on macOS.
if(MACOS)
    set(Boost_USE_STATIC_LIBS ON)
    find_package(Boost REQUIRED COMPONENTS program_options chrono)
    include_directories(/usr/local/include/)
    include_directories(${Boost_INCLUDE_DIRS})
endif(MACOS)

# Define our library.
add_library(shastaStaticLibrary STATIC ${SOURCES})

# Make sure the library is named shasta.a.
set_target_properties(shastaStaticLibrary PROPERTIES OUTPUT_NAME "shasta")
set_target_properties(shastaStaticLibrary PROPERTIES PREFIX "")
set_target_properties(shastaStaticLibrary PROPERTIES DEFINE_SYMBOL "")

# Install the static library into the bin directory.
install(TARGETS shastaStaticLibrary DESTINATION shasta-install/bin)




