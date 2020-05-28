cmake_minimum_required(VERSION 3.0)
project(shasta)

# This is the top level cmake file for a Shasta build. 
# Use "cmake -D" to modify these option defined below 
# in each "option" line. For example, to turn off
# building the dynamic executable use
# cmake -DBUILD_DYNAMIC_EXECUTABLE=OFF



# Figure out if we are on macOS.
set(MACOS OFF)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOS ON)
endif()
message(STATUS "CMAKE_SYSTEM_NAME is " ${CMAKE_SYSTEM_NAME})
message(STATUS "MACOS is " ${MACOS})



# Decide what we want to build.
if(MACOS)
    # For MacOS, only do the static builds.
    option(BUILD_STATIC_LIBRARY "Build the static library." ON)
    option(BUILD_STATIC_EXECUTABLE "Build the static executable." ON)
    option(BUILD_DYNAMIC_LIBRARY "Build the shared library." OFF)
    option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." OFF)
    option(BUILD_APPIMAGE "Build the AppImage." OFF)
else(MACOS)
    # For Linux, build everything.
    option(BUILD_USE_SPOA_WITH_CPU_DISPATCH "Link to SPOA library which uses cpu dispatch" ON)
    option(BUILD_STATIC_LIBRARY "Build the static library." ON)
    option(BUILD_STATIC_EXECUTABLE "Build the static executable." ON)
    option(BUILD_DYNAMIC_LIBRARY "Build the shared library." ON)
    option(BUILD_DYNAMIC_EXECUTABLE "Build the dynamic executable." ON)
    option(BUILD_APPIMAGE "Build the AppImage." OFF)
endif(MACOS)



# To build the static executable, we also need the static library.
if(BUILD_STATIC_EXECUTABLE AND NOT BUILD_STATIC_LIBRARY)
    set(BUILD_STATIC_LIBRARY ON)
    message(STATUS "Turned on BUILD_STATIC_LIBRARY because it is needed for BUILD_STATIC_EXECUTABLE.")
endif(BUILD_STATIC_EXECUTABLE AND NOT BUILD_STATIC_LIBRARY)



# To build the AppImage, we also need the dynamic executable.
if(BUILD_APPIMAGE AND NOT BUILD_DYNAMIC_EXECUTABLE)
    set(BUILD_DYNAMIC_EXECUTABLE ON)
    message(STATUS "Turned on BUILD_DYNAMIC_EXECUTABLE because it is needed for BUILD_APPIMAGE.")
endif(BUILD_APPIMAGE AND NOT BUILD_DYNAMIC_EXECUTABLE)



# To build the dynamic executable, we also need the dynamic library.
if(BUILD_DYNAMIC_EXECUTABLE AND NOT BUILD_DYNAMIC_LIBRARY)
    set(BUILD_DYNAMIC_LIBRARY ON)
    message(STATUS "Turned on BUILD_DYNAMIC_LIBRARY because it is needed for BUILD_DYNAMIC_EXECUTABLE.")
endif(BUILD_DYNAMIC_EXECUTABLE AND NOT BUILD_DYNAMIC_LIBRARY)



# Write out what we will actually build.
message(STATUS "BUILD_STATIC_LIBRARY is " ${BUILD_STATIC_LIBRARY})
message(STATUS "BUILD_STATIC_EXECUTABLE is " ${BUILD_STATIC_EXECUTABLE})
message(STATUS "BUILD_DYNAMIC_LIBRARY is " ${BUILD_DYNAMIC_LIBRARY})
message(STATUS "BUILD_DYNAMIC_EXECUTABLE is " ${BUILD_DYNAMIC_EXECUTABLE})
message(STATUS "BUILD_APPIMAGE is " ${BUILD_APPIMAGE})



# Option to build with -march=native.
option(BUILD_NATIVE "Build with -march=native." OFF)
message(STATUS "BUILD_NATIVE is " ${BUILD_NATIVE})



# Option to request a debug build.
option(BUILD_DEBUG "Make a debuggable build." OFF)
message(STATUS "BUILD_DEBUG is " ${BUILD_DEBUG})



# The BUILD_ID can be specified to identify the build
# This is normally used only when building a new GitHub release,
# in which case we use the following option when running Cmake:
# -DBUILD_ID="Shasta Release X.Y.Z"
# BUILD_ID should not contain a comma, otherwise the C preprocessor does not handle it correctly.
if(NOT DEFINED BUILD_ID)
    set(BUILD_ID "Shasta development build. This is not a released version.")
endif(NOT DEFINED BUILD_ID)
message(STATUS "BUILD_ID is: " ${BUILD_ID})



# Add the subdirectories we need.
if(BUILD_STATIC_LIBRARY)
    add_subdirectory(staticLibrary)
endif(BUILD_STATIC_LIBRARY)
if(BUILD_STATIC_EXECUTABLE)
    add_subdirectory(staticExecutable)
endif(BUILD_STATIC_EXECUTABLE)
if(BUILD_DYNAMIC_LIBRARY)
    add_subdirectory(dynamicLibrary)
endif(BUILD_DYNAMIC_LIBRARY)
if(BUILD_DYNAMIC_EXECUTABLE)
    add_subdirectory(dynamicExecutable)
endif(BUILD_DYNAMIC_EXECUTABLE)



# Install to the shasta-install directory.
set(CMAKE_INSTALL_PREFIX .)

# Install the scripts.
file(GLOB SCRIPTS scripts/*.py scripts/*.sh)
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta-install USE_SOURCE_PERMISSIONS)

# Install the docs directory.
install(DIRECTORY docs DESTINATION shasta-install)

# Install the docker directory.
install(DIRECTORY docker DESTINATION shasta-install)

# The targets built in each subdirectory are
# installed by the cmake file of each subdirectory.




