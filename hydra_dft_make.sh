#!/bin/bash

module purge
# GCCcore-6.4.0 is compatible with foss-2018a
module load  GSL/2.4-GCCcore-6.4.0
module load  Armadillo/8.400.0-foss-2018a

../dft_make.sh $@


	
