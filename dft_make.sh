#!/bin/bash

DIR1="build"
FLAG1=" "
Lib1=0

while test $# -gt 0
do
    case "$1" in
        clean) FLAG1="--clean-first"
            ;;
        debug) DIR1="debug"
               ;;
	lib) Lib1=1
            ;;
        ceci)
            module load releases/2018b
            module load Boost/1.67.0-foss-2018b
            module load CMake/3.12.1-GCCcore-7.3.0
            module load GSL/2.5-GCC-7.3.0-2.30
            ;;
        *) echo "bad option $1"
            ;;
    esac
    shift
done

if [ $Lib1 == '1' ]
then
    here=$(pwd)
    echo $here
    cd ../legacy_lib
    echo "Making Library"
    echo cmake --build $DIR1 $FLAG1
    cmake --build $DIR1 $FLAG1
    cd $here
fi

echo cmake --build $DIR1 $FLAG1 
cmake --build $DIR1 $FLAG1
