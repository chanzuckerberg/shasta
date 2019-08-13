# This creates a Shasta static executable that provides
# basic functionality and reduced performance. 
# For full functionality use the shared library built
# under directory src.

cmake_minimum_required(VERSION 3.0)
project(shasta-static-executable)

# Define the source files that go into our executable.
file(GLOB SOURCES *.cpp ../src/*.cpp)
include_directories(../src)

# Define our executable.
add_executable(shasta-static-executable ${SOURCES})
set_target_properties(shasta-static-executable PROPERTIES OUTPUT_NAME "shasta")

# For Linux, link statically.
# For macOS, this is not supported, and we only link some of the libraries statically.
if(NOT MACOS)
    set_target_properties(shasta-static-executable PROPERTIES LINK_FLAGS "-static" )
endif(NOT MACOS)

# Compiler options.
add_definitions(-std=c++14)
add_definitions(-ggdb3)
add_definitions(-O3)
if(BUILD_NATIVE)
    add_definitions(-march=native)
endif(BUILD_NATIVE)

# Turn compilation warnings on, but not on macOS.
if(NOT MACOS)
    add_definitions(-Wall -Wconversion -Wno-unused-result)
endif(NOT MACOS)

# Include directories required on macOS.
if(MACOS)
    set(Boost_USE_STATIC_LIBS ON)
    find_package(Boost REQUIRED COMPONENTS program_options chrono)
    include_directories(/usr/local/include/)
    include_directories(${Boost_INCLUDE_DIRS})
endif(MACOS)

# Pass the build id to the compile command line.
add_definitions(-DBUILD_ID=${BUILD_ID})

# This is required by dset64.hpp.
# It allows std::atomic<__uint128_t>::is_lock_free()
# to return true. However this does not happen in all gcc versions.
add_definitions(-mcx16)



# Definitions that control what components get built.
add_definitions(-DSHASTA_HTTP_SERVER)



# Specify the libraries to link with.
if(MACOS)
    target_link_libraries(shasta-static-executable ${Boost_LIBRARIES} /usr/local/lib/libspoa.a pthread)
else(MACOS)
    # For arcane reasons, statically linking with the pthread
    # library on Linux requires "--whole-archive". 
    target_link_libraries(
        shasta-static-executable
        atomic boost_system boost_program_options boost_chrono spoa png z
        -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
endif(MACOS)
  
  
    
# Also request the glibc library to be linked statically. 
# This is not supported on macOS.   
if(NOT MACOS)
    SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
endif(NOT MACOS)

# The static executable goes to the bin directory.
install(TARGETS shasta-static-executable DESTINATION shasta-install/bin)


