cmake_minimum_required(VERSION 3.10)

project(examples_geometry LANGUAGES CXX)

file(GLOB_RECURSE ALL_SOURCES LIST_DIRECTORIES false *.h *.cpp)
set(ENV{DISPLAY} ":0")

include_directories(${classical_dft_include_dirs})

add_executable(${PROJECT_NAME} ${ALL_SOURCES})
target_link_libraries(
        "${PROJECT_NAME}"
        "classical_dft"
        "${GRACE_NP}"
        "pthread")