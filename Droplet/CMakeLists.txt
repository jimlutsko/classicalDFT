cmake_minimum_required(VERSION 3.10)

# MODIFY: Project name is only for the log files generated when running the program
project(Droplet)

# MODIFY: name of the application file
add_executable(${PROJECT_NAME}  ${PROJECT_SOURCE_DIR}/droplet.cpp)

#MODIFY: enable this line if you want to use Grace
target_compile_definitions(${PROJECT_NAME} PUBLIC USE_GRACE)

## Set up the compilation
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_EXTENSIONS ON)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

## Find other packages needed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../scripts)

find_package(GSL REQUIRED)
find_package(OpenMP)
find_package(MathGL2 REQUIRED)
find_package(Armadillo  REQUIRED)


## Get the classicalDFT library and headers
find_library(CLASSICALDFT classicalDFT  PATHS ${PROJECT_SOURCE_DIR}/../Lib)
include_directories(${PROJECT_SOURCE_DIR}/../Lib/include)

if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

#get_target_property(XX  ${PROJECT_NAME} COMPILE_DEFINITIONS)
#message("target property is : ${XX}")


target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CLASSICALDFT}
  "gsl"
  "armadillo"
  "mgl"
  "fftw3_omp"
  "fftw3"
  )

get_target_property(XX  ${PROJECT_NAME} COMPILE_DEFINITIONS)
message("variable xx is ${XX}")


if(XX MATCHES "USE_GRACE")
  target_link_libraries(${PROJECT_NAME} PUBLIC "grace_np")
endif()

