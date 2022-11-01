cmake_minimum_required(VERSION 3.10)

set(USE_OMP false)

#This allows us to include the git SHA to identify exactly which code is being used
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_stuff/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/GitSHA1.cpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/src/GitSHA1.cpp" @ONLY)

# Project name and a few useful settings. Other commands can pick up the results
project(classicalDFT VERSION 1.0 DESCRIPTION "classical DFT library and apps" LANGUAGES CXX)

#Required
find_package(Boost REQUIRED REQUIRED COMPONENTS serialization)

#Version information
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_SOURCE_DIR}/include/config.h"
  )

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        "/usr/local/include"
)

message("Project source directory is " ${PROJECT_SOURCE_DIR})

#file(GLOB all_SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
set(all_SRCS
  ${PROJECT_SOURCE_DIR}/src/Dynamical_Matrix.cpp
  ${PROJECT_SOURCE_DIR}/src/Density.cpp
  ${PROJECT_SOURCE_DIR}/src/DFT.cpp
  ${PROJECT_SOURCE_DIR}/src/DFT_Coex.cpp
  ${PROJECT_SOURCE_DIR}/src/Species.cpp
  ${PROJECT_SOURCE_DIR}/src/FMT_Species.cpp
  ${PROJECT_SOURCE_DIR}/src/FMT.cpp
  ${PROJECT_SOURCE_DIR}/src/Minimizer.cpp
  ${PROJECT_SOURCE_DIR}/src/DDFT.cpp
  ${PROJECT_SOURCE_DIR}/src/Arnoldi.cpp
  ${PROJECT_SOURCE_DIR}/src/Eigenvalues.cpp
  ${PROJECT_SOURCE_DIR}/src/Log_Det.cpp
  ${PROJECT_SOURCE_DIR}/src/Interaction.cpp
  ${PROJECT_SOURCE_DIR}/src/GitSHA1.cpp
  ${PROJECT_SOURCE_DIR}/src/myColor.cpp
  ${PROJECT_SOURCE_DIR}/src/options.cpp
  ${PROJECT_SOURCE_DIR}/src/Table.cpp
  ${PROJECT_SOURCE_DIR}/src/poly34.cpp
  ${PROJECT_SOURCE_DIR}/src/visit_writer.cpp
)

add_library(${PROJECT_NAME} STATIC ${all_SRCS})

# All users of this library will need at least C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME}  PUBLIC -fconcepts)

# find gsl
link_directories("/usr/local/lib")

target_link_libraries(${PROJECT_NAME} PUBLIC Boost::serialization)

# enable omp if it exists

message("USE_OMP = " ${USE_OMP})

if(USE_OMP)
  find_package(OpenMP REQUIRED)
  
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD  COMMAND mv libclassicalDFT* ../)



