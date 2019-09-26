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
        *) echo "bad option $1"
            ;;
    esac
    shift
done

if [ $Lib1 == '1' ]
then
    here=$(pwd)
    echo $here
    cd ../Lib
    echo "Making Library"
    echo cmake --build $DIR1 $FLAG1
    cmake --build $DIR1 $FLAG1
    cd $here
fi

echo cmake --build $DIR1 $FLAG1 
cmake --build $DIR1 $FLAG1


	
