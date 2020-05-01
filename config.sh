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