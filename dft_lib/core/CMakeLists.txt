cmake_minimum_required(VERSION 3.10)

project(classical_dft LANGUAGES CXX)
include_directories(${dft_lib_include_path})

# Set up the source files required for the compilation to succeed
set(ALL_SOURCES
        # Graphics utils:
        ${dft_lib_src_path}/graph/grace.cpp
        # Library exceptions:
        ${dft_lib_src_path}/exceptions/grace_exceptions.cpp
        ${dft_lib_src_path}/exceptions/parameter_exceptions.cpp
        # Utils:
        ${dft_lib_src_path}/utils/config_parser.cpp
        # Numeric:
        ${dft_lib_src_path}/numerics/integration.cpp
        ${dft_lib_src_path}/numerics/arithmetic.cpp
        # Physics:
        ${dft_lib_src_path}/physics/potentials/intermolecular/potential.cpp
        # Geometry:
        ${dft_lib_src_path}/geometry/vertex.cpp
        )

# Set up the dependencies and naming of the library
set(BINARY ${PROJECT_NAME})
add_library(${BINARY} STATIC ${ALL_SOURCES})
target_link_libraries(${BINARY})

# Expose project's public includes to other subprojects through cache variable.
set(${PROJECT_NAME}_include_dirs ${dft_lib_include_path}
        CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

# Move the static library file lib{project_name}.a to the root folder
# This does not work OK with tests, as they need to link "-L dft_core"
# add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD  COMMAND mv lib${PROJECT_NAME}* ../)