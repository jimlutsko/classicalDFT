cmake_minimum_required(VERSION 3.10)

project(Wall_Interaction)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

## Set up the compilation
add_executable(${PROJECT_NAME}  ${PROJECT_SOURCE_DIR}/wall.cpp)
set(CMAKE_CXX_EXTENSIONS ON)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

## Find packages needed
find_package(GSL REQUIRED)
find_package(Armadillo  REQUIRED)
find_package(OpenMP)
#find_package(GRACE_NP REQUIRED)

## Get the classicalDFT library and headers
find_library(CLASSICALDFT classicalDFT  PATHS ${PROJECT_SOURCE_DIR}/../Lib)
include_directories(${PROJECT_SOURCE_DIR}/../Lib/include)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()
      
target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CLASSICALDFT}
  "grace_np"
  "gsl"
  "armadillo"
  "mgl"
  "fftw3_omp"
  "fftw3"
  )

