cmake_minimum_required(VERSION 2.6)
project(shasta)

set(CMAKE_INSTALL_PREFIX .)

add_subdirectory(src)

# Install the scripts.
file(GLOB SCRIPTS "scripts/*.py")
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)

# Install the configuration files.
install(DIRECTORY conf DESTINATION shasta-install USE_SOURCE_PERMISSIONS)

# Install the docs directory.
install(DIRECTORY docs DESTINATION shasta-install)

# The shared library is installed by the cmake file in the src directory.




