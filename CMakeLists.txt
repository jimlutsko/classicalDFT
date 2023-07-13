cmake_minimum_required(VERSION 3.10)

set(USE_OMP true)
set(USE_MPI true)
set(USE_SLEPC true)

set(USE_FMT_WEIGHTS_BEFORE_JUN_2021 false)
set(USE_FMT_WEIGHTS_BEFORE_JAN_2023 false)

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
  ${PROJECT_SOURCE_DIR}/src/FMT_Species_EOS.cpp
  ${PROJECT_SOURCE_DIR}/src/FMT.cpp
  ${PROJECT_SOURCE_DIR}/src/Potential.cpp
  ${PROJECT_SOURCE_DIR}/src/Minimizer.cpp
  ${PROJECT_SOURCE_DIR}/src/DDFT.cpp
  ${PROJECT_SOURCE_DIR}/src/Arnoldi.cpp
  ${PROJECT_SOURCE_DIR}/src/Eigenvalues.cpp
  ${PROJECT_SOURCE_DIR}/src/DFT_Petsc.cpp
  ${PROJECT_SOURCE_DIR}/src/DFT_Slepc.cpp
  ${PROJECT_SOURCE_DIR}/src/Log_Det.cpp
  ${PROJECT_SOURCE_DIR}/src/Interaction.cpp
  ${PROJECT_SOURCE_DIR}/src/GitSHA1.cpp
  ${PROJECT_SOURCE_DIR}/src/myColor.cpp
  ${PROJECT_SOURCE_DIR}/src/options.cpp
  ${PROJECT_SOURCE_DIR}/src/Table.cpp
#  ${PROJECT_SOURCE_DIR}/src/poly34.cpp
  ${PROJECT_SOURCE_DIR}/src/visit_writer.cpp
)

add_library(${PROJECT_NAME} STATIC ${all_SRCS})

# All users of this library will need at least C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME}  PUBLIC -fconcepts)

# find gsl
link_directories("/usr/local/lib")

target_link_libraries(${PROJECT_NAME} PUBLIC Boost::serialization)

message("USE_SLEPC = " ${USE_SLEPC})

if(USE_SLEPC)
  find_package(PkgConfig REQUIRED) 
  set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  pkg_check_modules(PETSC_PKG REQUIRED IMPORTED_TARGET $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig/PETSc.pc)
  pkg_check_modules(SLEPC_PKG REQUIRED IMPORTED_TARGET $ENV{SLEPC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig/slepc.pc)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_SLEPC)
  target_link_libraries(${PROJECT_NAME} PUBLIC PkgConfig::PETSC_PKG PkgConfig::SLEPC_PKG)
endif()

message("USE_MPI = " ${USE_MPI})

if(USE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_INCLUDE_PATH})
  
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
    set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
  
  #target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MPI)
  target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX)
endif()

# enable omp if it exists

message("USE_OMP = " ${USE_OMP})

if(USE_OMP)
  find_package(OpenMP REQUIRED)
  
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD  COMMAND mv libclassicalDFT* ../)

## Uncomment the next line for debug information on the FMT weights calculation
#target_compile_definitions(${PROJECT_NAME} PUBLIC DEBUG_FMT_WEIGHTS)

message("USE_FMT_WEIGHTS_BEFORE_JUN_2021 = " ${USE_FMT_WEIGHTS_BEFORE_JUN_2021})
message("USE_FMT_WEIGHTS_BEFORE_JAN_2023 = " ${USE_FMT_WEIGHTS_BEFORE_JAN_2023})

if(USE_FMT_WEIGHTS_BEFORE_JUN_2021)
  target_compile_definitions(${PROJECT_NAME} PUBLIC FMT_WEIGHTS_BEFORE_JUN_2021)
endif()

if(USE_FMT_WEIGHTS_BEFORE_JAN_2023)
  target_compile_definitions(${PROJECT_NAME} PUBLIC FMT_WEIGHTS_BEFORE_JAN_2023)
endif()

## Enable/Disable Arnoldi (can pose problems as the current implementation only works with Armadillo)
#target_compile_definitions(${PROJECT_NAME} PUBLIC ARNOLDI)


