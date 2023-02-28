while test $# -gt 0
do
    case "$1" in
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

#! /bin/sh
root_dir=$(pwd)
parent_dir="$(dirname "$root_dir")"

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release $1 ..
cd $root_dir

mkdir -p debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug $1 ..
cd $root_dir
