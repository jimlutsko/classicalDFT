#! /bin/sh
mkdir build
cd build
cmake ..
cd ..
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..
