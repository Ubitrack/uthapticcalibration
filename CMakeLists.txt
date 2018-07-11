cmake_minimum_required(VERSION 3.4.0)

project(ubitrack_hapticcalibration VERSION 1.3.0 LANGUAGES CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

include(GNUInstallDirs)
include(cmake/UbitrackConfig)

# External Libraries (could this be solved with conan dependencies ??)

# library
add_subdirectory(src)
