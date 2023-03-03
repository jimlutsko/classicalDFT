#! /bin/bash

host=$(cat /etc/hostname)
echo "Hostname: $host"

if [ $host == "hercules" ]
then
    	echo "Loading modules for hercules"

	module load releases/2019b
	module load Boost/1.71.0-gompi-2019b
	module load CMake/3.15.3-GCCcore-8.3.0
	module load Eigen/3.3.7
	module load GSL/2.6-GCC-8.3.0
	module load FFTW/3.3.8-gompi-2019b

elif [ $host == "dragon1-h0.umons.ac.be" ]
then
    	echo "Loading modules for dragon1"

	module load releases/2019b
	module load Boost/1.71.0-gompi-2019b
	module load CMake/3.15.3-GCCcore-8.3.0
	module load Eigen/3.3.7
	module load GSL/2.6-GCC-8.3.0
	module load FFTW/3.3.8-gompi-2019b

elif [ $host == "dragon2-ctrl0.umons.ac.be" ]
then
	echo "Loading modules for dragon2"

	module load releases/2020b
	module load Armadillo/10.5.3-foss-2020b
	module load Boost/1.74.0-GCC-10.2.0
	module load CMake/3.18.4-GCCcore-10.2.0
	module load Eigen/3.3.9-GCCcore-10.2.0
	module load GSL/2.6-GCC-10.2.0

#	module load releases/2018b
#	module load Armadillo/9.600.5-foss-2018b
#	module load Boost/1.67.0-foss-2018b
#	module load CMake/3.12.1-GCCcore-7.3.0
#	module load GSL/2.5-GCC-7.3.0-2.30

else
    	echo "Unknown host. Cannot load modules"
fi

