# This creates a Shasta executable that is
# highly instrumented for debugging.
# It is not built during a Shasta build.

cmake_minimum_required(VERSION 3.0)
project(shasta-static-executable)

# Define the source files that go into our executable.
file(GLOB SOURCES ../src-static-executable/*.cpp ../src/*.cpp)
include_directories(../src)

# Define our executable.
add_executable(shasta-dynamic-executable-debug ${SOURCES})
set_target_properties(shasta-dynamic-executable-debug PROPERTIES OUTPUT_NAME "shasta")


# Compiler options.
add_definitions(-std=c++14)
add_definitions(-ggdb3)
add_definitions(-O3)

# This traces lots of memory errors.
add_definitions(-fsanitize=address)
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")

# Turn compilation warnings on, but not on macOS.
add_definitions(-Wall -Wconversion -Wno-unused-result)

# Definitions that control what components get built.
add_definitions(-DSHASTA_HTTP_SERVER)

# Specify the libraries to link with.
target_link_libraries(
    shasta-dynamic-executable-debug
    atomic boost_system boost_program_options boost_chrono spoa png z
    pthread)



