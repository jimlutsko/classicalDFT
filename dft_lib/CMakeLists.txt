cmake_minimum_required(VERSION 3.10)

# Project definition and description (DFTLib)
project(classicalDFTv2
    VERSION 2.0
    DESCRIPTION "classical DFT library"
    LANGUAGES CXX)

#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_COLOR_MAKEFILE   ON)

# Setting up the cxx compiler to
set(CMAKE_CXX_STANDARD 14)
set(dft_lib_include_path ${PROJECT_SOURCE_DIR}/core)
set(dft_lib_src_path ${PROJECT_SOURCE_DIR}/core/src)

# Main include folder
include_directories(${dft_lib_include_path})

# DFTLib sub-projects.
add_subdirectory(core)
add_subdirectory(tests)