## Note: I cannot make a script directly loading the correct modules based on machine hostname
## because the hostname on worker machines where the jobs are executed is not always the same.
## This is why I save the module list in ~/ceci_modules, because the home directory is not shared
## between clusters but is accessible from the worker machines.

#! /bin/bash

host=$(cat /etc/hostname)
echo "Hostname: $host"
host=${host%-*}
echo "Truncated: $host"

destination=~/ceci_modules
echo "## File generated on $(date)" > $destination
echo "## Hostname: $(cat /etc/hostname)" >> $destination
echo " " >> $destination

if [ $host == "hercules" ]
then
    	echo "Creating file $destination on hercules"

	echo "module load releases/2019b" >> $destination
	echo "module load Boost/1.71.0-gompi-2019b" >> $destination
	echo "module load CMake/3.15.3-GCCcore-8.3.0" >> $destination
	echo "module load Eigen/3.3.7" >> $destination
	echo "module load GSL/2.6-GCC-8.3.0" >> $destination
	echo "module load FFTW/3.3.8-gompi-2019b" >> $destination

elif [ $host == "dragon1" ]
then
    	echo "Creating file $destination on dragon1"
	
	echo "module load releases/2019b" >> $destination
	echo "module load Boost/1.71.0-gompi-2019b" >> $destination
	echo "module load CMake/3.15.3-GCCcore-8.3.0" >> $destination
	echo "module load Eigen/3.3.7" >> $destination
	echo "module load GSL/2.6-GCC-8.3.0" >> $destination
	echo "module load FFTW/3.3.8-gompi-2019b" >> $destination

elif [ $host == "dragon2" ]
then
	echo "Creating file $destination on dragon2"

	echo "module load releases/2020b" >> $destination
	echo "module load Armadillo/10.5.3-foss-2020b" >> $destination
	echo "module load Boost/1.74.0-GCC-10.2.0" >> $destination
	echo "module load CMake/3.18.4-GCCcore-10.2.0" >> $destination
	echo "module load Eigen/3.3.9-GCCcore-10.2.0" >> $destination
	echo "module load GSL/2.6-GCC-10.2.0" >> $destination

#	echo "module load releases/2018b" >> $destination
#	echo "module load Armadillo/9.600.5-foss-2018b" >> $destination
#	echo "module load Boost/1.67.0-foss-2018b" >> $destination
#	echo "module load CMake/3.12.1-GCCcore-7.3.0" >> $destination
#	echo "module load GSL/2.5-GCC-7.3.0-2.30" >> $destination

else
    	echo "Unknown host. Cannot load modules"
fi

