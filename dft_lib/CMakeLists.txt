cmake_minimum_required(VERSION 3.10)

# Project definition and description
# Name: classicalDFTv2 to avoid overlapping with the already existing
project(classicalDFTv2
        VERSION 2.0
        DESCRIPTION "classical DFT library"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Finding dependencies:
include(FindDependencies)
find_library_or_set_default(GRACE_NP grace_np /usr/local)
find_library_or_set_default(GTEST gtest /usr/local)

# Setting up the cxx compiler
set(CMAKE_CXX_STANDARD 14)
set(dft_lib_include_path ${PROJECT_SOURCE_DIR}/core)
set(dft_lib_src_path ${PROJECT_SOURCE_DIR}/core/src)

if (LINUX AND NOT APPLE)
    # Main include folder
    message(STATUS "Linux system detected...")
    include_directories(
        ${dft_lib_include_path})
elseif(WIN32)
    message(STATUS "Windows system detected...")
    include_directories(
        ${dft_lib_include_path})
elseif(APPLE)
    # Check whether we are linux or macOS
    set(brew_include_path /usr/local/include)
    include_directories(
            ${dft_lib_include_path}
            ${brew_include_path})
endif()

# dft_lib sub-projects:
add_subdirectory(core)

# dft_lib tests:
enable_testing()
add_subdirectory(tests)

# dft_lib examples:
add_subdirectory(examples/console)