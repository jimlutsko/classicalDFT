cmake_minimum_required(VERSION 3.10)

#This allows us to include the git SHA to identify exactly which code is being used
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_stuff/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/GitSHA1.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/src/GitSHA1.cpp" @ONLY)

# Project name and a few useful settings. Other commands can pick up the results
project(classicalDFT VERSION 1.0 DESCRIPTION "classical DFT library and apps" LANGUAGES CXX)

#Version information
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_SOURCE_DIR}/include/config.h"
  )

include_directories(${PROJECT_SOURCE_DIR}/include) 

file(GLOB all_SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library(${PROJECT_NAME} STATIC ${all_SRCS})

# All users of this library will need at least C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME} PUBLIC -O3 -MMD)


# enable omp if it exists
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_compile_options(${PROJECT_NAME} PUBLIC -fopenmp)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD  COMMAND mv libclassicalDFT* ../)



