[![Unit tests](https://github.com/jimlutsko/classicalDFT/workflows/Unit%20tests/badge.svg)](https://github.com/jimlutsko/classicalDFT/actions?query=workflow%3A%22Unit%20tests%22)

# classicalDFT

## Welcome to classicalDFT

*classicalDFT* aims at making possible to do [classical density functional theory](https://en.wikipedia.org/wiki/Density_functional_theory#Classical_Density_Functional_Theory) (DFT) calculations in a smooth manner.

**classicalDFT** is a repository which started as a suite of code for doing advanced calculations in various research projects on the broad field of statistical physics. This suite of code has matured over the years and it successfully served its purpose. The evolution of this repository into a standard library (*classicalDFT*) was motivated by the apparent lack of a standard open-source repository for classical DFT calculations in a compiled and robust language, as is the case of `C++`.

In the future, a mailing list will be advertised here for questions, discussions, and other topics which might be worth channeling in a unified way.

### Getting started

The standard rules for installing a `C++` library from an external repository apply in this case too. However, such standard steps which might be obvious for experienced developers could become overwhelming for developers or researchers who just have some basic understanding of software development with `CMake` and `C++`. Considering the broad audience this project might be subject to, we recommend the following documentation links for the different users:

* Beginners: [classicalDFT Primer](README.md)
* Experienced programmers: [How to install classicalDFT](documentation/installation/README.md)

## Requirements

**classicalDFT** is thought to keep requirements for building at minimum. However, there are some requirements which 
need to be satisfied for the correct funcitoning of the library. Such dependencies are:

* C++14-standard-compliant compiler, we recommend to check [GNU Compiler Collection](https://gcc.gnu.org/) 
* [CMake](https://cmake.org/download/) (version >= 3.8)
* [Boost](https://www.boost.org/)
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)
* [FFTW3](http://www.fftw.org/) a library designed to compute discrete Fourier transforms
* [Grace](http://plasma-gate.weizmann.ac.il/Grace/) plotting tool for the X Window System
* [Google Test](https://github.com/google/googletest) a testing framework developed by Google's Testing Technology team

If you notice any problems on your platform, please notify [classicaldft@classicaldftgroup.com](). Patches for fixing them are more than welcome!

## Features

In this section we will list the main features implemented in **classicalDFT**, e.g. the different DFT energy models, differentiation methods, etc. The library is currently being refactored with the aim of achieving a better harmonisation between the different modules implemented. For this reason, we'll be listing here only the modules which are ready in the new refactored version (though you can always check the original version in [here](legacy_lib/)):

1. [Console](dft_lib/examples/console/README.md): A namespace with some convenient functionality when interacting with the console/terminal output
2. [Grace](dft_lib/examples/graphs/README.md): A wrapper class which allows us to produce live graphs by interacting with xmgrace

## Platforms

**classicalDFT** has been [used](dft_lib/docs/setup-env.md) on a variety of platforms:

- Linux
- Mac OSX
- Windows (by using [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux))

## Authors

* **James F. Lutsko** - *Initial work* - [lutsko.com](http://lutsko.com)

## Contributing

If you feel like contributing to this project, please read the [CONTRIBUTING.md](CONTRIBUTING.md) where you will find the basic details on how to contribute to **classicalDFT**.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details

## Acknowledgments

Development was funded by ESA