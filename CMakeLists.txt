cmake_minimum_required(VERSION 2.6)
project(shasta)

# Figure out if we are on macOS.
set(MACOS OFF)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOS ON)
endif()
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})
message(STATUS "MACOS is " ${MACOS})



# Figure out what we want to build.
# By default, we build only the executable on macOS,
# and both the executable and the library on Linux.
if(MACOS)
    option(BUILD_LIBRARY "Build the shared library." OFF)
    option(BUILD_EXECUTABLE "Build the static executable." ON)
else(MACOS)
    option(BUILD_LIBRARY "Build the shared library." ON)
    option(BUILD_EXECUTABLE "Build the static executable." ON)
endif(MACOS)



# Option to build the shared library with -march=native.
# The static executable is never built with -march=native,
# regardless of how this is set.
option(BUILD_NATIVE "Build with -march=native (shared library only)." ON)
message(STATUS "BUILD_NATIVE is " ${BUILD_NATIVE})


# The BUILD_ID can be specified to identify the build
# This is normally used only when building a new GitHub release,
# in which case we use the following option when running Cmake:
# -DBUILD_ID="Shasta Release X.Y.Z"
if(NOT DEFINED BUILD_ID)
    set(BUILD_ID "Shasta development build. This is not a released version.")
endif(NOT DEFINED BUILD_ID)
message(STATUS "BUILD_ID is: " ${BUILD_ID})



# Add the subdirectories.
if(BUILD_LIBRARY)
    add_subdirectory(src)
endif(BUILD_LIBRARY)
if(BUILD_EXECUTABLE)
    add_subdirectory(src-static-executable)
endif(BUILD_EXECUTABLE)
message(STATUS "BUILD_LIBRARY is " ${BUILD_LIBRARY})
message(STATUS "BUILD_EXECUTABLE is " ${BUILD_EXECUTABLE})



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




