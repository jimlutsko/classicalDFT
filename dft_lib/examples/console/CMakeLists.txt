cmake_minimum_required(VERSION 3.10)

project(examples_console LANGUAGES CXX)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)
set(ALL_SOURCES ${TEST_SOURCES})
set(ENV{DISPLAY} ":0")

include_directories(${classical_dft_include_dirs})

add_executable(${PROJECT_NAME} ${ALL_SOURCES})
target_link_libraries(
        "${PROJECT_NAME}"
        "classical_dft"
        "${GRACE_NP}"
        "${GTEST}"
        "pthread")