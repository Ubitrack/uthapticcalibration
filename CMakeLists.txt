project(ubitrack_hapticcalibration)

cmake_minimum_required(VERSION 3.4.0)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(KEEP_RPATHS)

include(GNUInstallDirs)
include(cmake/UbitrackConfig)

# External Libraries (could this be solved with conan dependencies ??)

# library
add_subdirectory(src)
