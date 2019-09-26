#! /bin/sh
mkdir build
cd build -DCMAKE_BUILD_TYPE=Release ..
cmake ..
cd ..
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..
