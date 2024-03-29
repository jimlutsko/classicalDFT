cmake_minimum_required(VERSION 3.10)

project(classical_dft_tests LANGUAGES CXX)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp ./**.cpp)
set(ALL_SOURCES ${TEST_SOURCES})
set(ENV{DISPLAY} ":0")

include_directories(${classical_dft_include_dirs})

add_executable(${PROJECT_NAME} ${ALL_SOURCES})
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
target_link_libraries(
        "${PROJECT_NAME}"
        "classical_dft"
        "${GSLCBLAS}"
        "${GSLLIB}"
        "${GRACE_NP}"
        "${GTEST}"
        "pthread")

add_custom_target(gtests ${CMAKE_CTEST_COMMAND} -V)

configure_file(utils/config_files/config.ini config.ini COPYONLY)
configure_file(utils/config_files/config.json config.json COPYONLY)
configure_file(utils/config_files/config.xml config.xml COPYONLY)
configure_file(utils/config_files/config.info config.info COPYONLY)