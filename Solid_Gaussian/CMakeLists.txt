cmake_minimum_required(VERSION 3.10)

# MODIFY: Project name: used to name the executable 
project(Solid)

# MODIFY: This needs to point to the Lib directory
set(PATH_TO_DFT_LIB ${PROJECT_SOURCE_DIR}/../Lib)

# MODIFY: this is the app file 
add_executable(${PROJECT_NAME}  ${PROJECT_SOURCE_DIR}/solid_gaussian.cpp ../Lib/src/VDW1.cpp)
add_executable(LIQVAP  ${PROJECT_SOURCE_DIR}/LiquidVapor.cpp ../Lib/src/VDW1.cpp)

#MODIFY: enable these lines if you want to use Grace (Grace.h, Grace.cpp)
#set(USE_GRACE true)

#MODIFY: enable these lines if you want to use Display object (Display.h: creates png images)
#set(USE_MGL true)

################### End of customizations

## Set up the compilation
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_EXTENSIONS ON)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_features(LIQVAP PUBLIC cxx_std_11)

## Find other packages needed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../scripts)

find_package(GSL REQUIRED)
find_package(OpenMP)
find_package(Armadillo  REQUIRED)

if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

## Get the classicalDFT library and headers
find_library(CLASSICALDFT classicalDFT  PATHS ${PATH_TO_DFT_LIB})
include_directories(${PATH_TO_DFT_LIB}/include)

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CLASSICALDFT}
  "gsl"
  "gslcblas"
  "armadillo"
  "fftw3_omp"
  "fftw3"
  )

target_link_libraries(LIQVAP PUBLIC
  ${CLASSICALDFT}
  "gsl"
  "gslcblas"
  "armadillo"
  "fftw3_omp"
  "fftw3"
  )


if(USE_MGL)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_MGL)
  target_include_directories(${PROJECT_NAME} PUBLIC ${PATH_TO_DFT_LIB}/additionalCode)
  find_package(MathGL2 REQUIRED)
  target_link_libraries(${PROJECT_NAME} PUBLIC "mgl")
endif()

if(USE_GRACE)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_GRACE)
  target_sources(${PROJECT_NAME} PUBLIC ${PATH_TO_DFT_LIB}/additionalCode/Grace.cpp)
  target_include_directories(${PROJECT_NAME} PUBLIC ${PATH_TO_DFT_LIB}/additionalCode)
  target_link_libraries(${PROJECT_NAME} PUBLIC "grace_np")
endif()
