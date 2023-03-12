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
        ceci) source load_ceci_modules.sh
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
