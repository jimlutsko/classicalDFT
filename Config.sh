#! /bin/sh
mkdir build
cd build 
cmake -DCMAKE_BUILD_TYPE=Release $1 ..
cd ..
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug $1 ..
cd ..
