#! /bin/sh
mkdir build
cd build 
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..
