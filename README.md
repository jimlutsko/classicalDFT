# Welcome to classicalDFT

This is a suite of code for doing classical DFT calculations. Documentation will eventually appear here.  

The main branch of this repository is no longer supported: all functionality will eventually be restored in a new branch called "MultiSpecies".

## Getting Started

TBD

### Prerequisites

What things you need to install the software and how to install them

TBD


## Installation

1. Perform the following to setup and compile the library:

>cd Lib

>../Config.sh

>../dft_make.sh

2. To build the model application do the following

>cd ../Droplet

>../Config.sh

>../dft_make.sh'


You can cd into TEST and run using ../Droplet input.dat.

## Notes:
1. The Config.sh command is only run the first time the application is created. It configures things for cmake.
2. Rebuild using ../dft_make.sh. Note that this also takes three possible arguments in any order: "clean", "debug" and "lib". Their effects are as follows:
   * "clean" causes a clean build (all objects are first deleted)
   * "debug" performs a debug build
   * "lib" causes the libraries to be rebuilt also.
   * Note that "dft_make.sh debug lib" causes both the app and the libraries to be built in debug mode. Similarly, "dft_make.sh lib" causes both to be built in release mode.

4. The first few lines of CMakeLists.txt contain information that you will want to modify if this is used as a model to create another application.



## Authors

* **James F. Lutsko** - *Initial work* - [lutsko.com](http://lutsko.com)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Development was funded by ESA




