cmake_minimum_required(VERSION 2.6)
project(shasta)

set(CMAKE_INSTALL_PREFIX .)

add_subdirectory(src)

file(GLOB SCRIPTS "scripts/*.py")
install(PROGRAMS ${SCRIPTS} DESTINATION shasta-install/bin)




