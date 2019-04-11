cmake_minimum_required(VERSION 2.6)
project(shasta)

# Add the directory that build the shared library.
option(BUILD_LIBRARY "Build the shared library." ON)
message(STATUS "BUILD_LIBRARY is " ${BUILD_LIBRARY})
if(BUILD_LIBRARY)
    add_subdirectory(src)
endif(BUILD_LIBRARY)

# Add the directory that builds the static executable.
option(BUILD_EXECUTABLE "Build the static executable." OFF)
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})
if(BUILD_EXECUTABLE)
    add_subdirectory(src-static-executable)
endif(BUILD_EXECUTABLE)

# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta-install USE_SOURCE_PERMISSIONS)

# Install the docs directory.
install(DIRECTORY docs DESTINATION shasta-install)

# The shared library is installed by the cmake file in the src directory.
# The static executable is installed by the cmake file in the src-static-executable directory.




