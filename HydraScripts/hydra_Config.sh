#!/bin/bash

module purge
module load CMake/3.10.2-GCCcore-6.4.0

# GCCcore-6.4.0 is compatible with foss-2018a
module load  GSL/2.4-GCCcore-6.4.0
module load  Armadillo/8.400.0-foss-2018a
module load MathGL/2.4.1-foss-2018a

../Config.sh "-DARMADILLO_LIBRARY=$EBROOTARMADILLO/lib64 -DARMADILLO_INCLUDE_DIR=$EBROOTARMADILLO/include -DMATHGL2_LIBRARY=$EBROOTMATHGL/lib64 -DMATHGL2_INCLUDE_DIR=$EBROOTMATHGL/include"

