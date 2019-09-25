## Welcome to classicalDFT

This is a suite of code for doing classical DFT calculations.

# Installation

1. Switch to the branch MultiSpecies (command "git checkout MultiSpecies").

2. erform the following to setup and compile the library:

cd Lib
../Config.sh
../dft_make.sh

3. To build the model application do the following

cd ../Droplet
../Config.sh
../dft_make.sh


You can cd into TEST and run using ../Droplet input.dat.

# Notes:
1. The Config.sh command is only run the first time the application is created. It configures things for cmake.

2. Rebuild using ../dft_make.sh. Note that this also takes three possible arguments in any order: "clean", "debug" and "lib". Their effects are as follows:

..1 "clean" causes a clean build (all objects are first deleted)

..2 "debug" performs a debug build

..3 "lib" causes the libraries to be rebuilt also.

3. Note that "dft_make.sh debug lib" causes both the app and the libraries to be built in debug mode. Similarly, "dft_make.sh lib" causes both to be built in release mode.

4. The first few lines of CMakeLists.txt contain information that you will want to modify if this is used as a model to create another application.




