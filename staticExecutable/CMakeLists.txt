cmake_minimum_required(VERSION 3.0)
project(shastaStaticExecutable)

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

# Source files
file(GLOB SOURCES ../srcMain/*.cpp)

# Include directory.
include_directories(../src)

# Include directories required on macOS.
if(MACOS)
    set(Boost_USE_STATIC_LIBS ON)
    find_package(Boost REQUIRED COMPONENTS program_options chrono)
    include_directories(/usr/local/include/)
    include_directories(${Boost_INCLUDE_DIRS})
endif(MACOS)

# Define our executable.
add_executable(shastaStaticExecutable ${SOURCES})
set_target_properties(shastaStaticExecutable PROPERTIES OUTPUT_NAME "shasta")

# Request a static executable.
# For macOS, this is not supported, and we only link some of the libraries statically.
if(NOT MACOS)
    set_target_properties(shastaStaticExecutable PROPERTIES LINK_FLAGS "-static" )
endif(NOT MACOS)

# Libraries to link with.
if(MACOS)
    find_library(PNG_LIBRARY libpng.a REQUIRED)
    #find_library(Z_LIBRARY libz.a REQUIRED)
    target_link_libraries(
    shastaStaticExecutable 
    shastaStaticLibrary
    ${Boost_LIBRARIES} /usr/local/lib/libspoa.a ${PNG_LIBRARY} z pthread)
else(MACOS)
    # For arcane reasons, statically linking with the pthread
    # library on Linux requires "--whole-archive".
    if(X86_64 AND BUILD_USE_SPOA_WITH_CPU_DISPATCH)
        target_link_libraries(
            shastaStaticExecutable
            shastaStaticLibrary
            atomic boost_system boost_program_options boost_chrono spoa cpu_features png z
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
    else(X86_64 AND BUILD_USE_SPOA_WITH_CPU_DISPATCH)
        target_link_libraries(
            shastaStaticExecutable
            shastaStaticLibrary
            atomic boost_system boost_program_options boost_chrono spoa png z
            -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
    endif(X86_64 AND BUILD_USE_SPOA_WITH_CPU_DISPATCH) 
    
endif(MACOS)
  
# Also request the glibc library to be linked statically. 
# This is not supported on macOS.   
if(NOT MACOS)
    SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
endif(NOT MACOS)

# The static executable goes to the bin directory.
install(TARGETS shastaStaticExecutable DESTINATION shasta-install/bin)


