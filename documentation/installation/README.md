## Generic build instructions

### Building classicalDFT v2.0

The new refactored version of the library [dft_lib](../../dft_lib/) can be easily set up and tested by utilising the bash script [buildlib.sh](../../buildlib.sh). This script has the same modes as the [dft_make.sh](../../dft_make.sh) which was used to build the version 1.0. By placing ourselves at the root folder of the repository, we just need to do:

```
$ chmod a+x buildlib.sh
$ ./buildlib.sh [MODE]
```

Where `[MODE]` can be substituted by any of the following words:

* **clean** causes a clean build (all objects are first deleted)
* **debug** performs a debug build
* **lib** causes the libraries to be rebuilt also
* **test** causes the rebuild and testing of the library

As an example, running tests will produces something like:

```bash
$ ./buildlib.sh test
======================================
  Starting folder configuration
======================================
release
...
debug
...
======================================
  classicalDFT library
======================================
[Step 1/2] Making Library...
cmake --build build/release
Scanning dependencies of target classical_dft
...
[Step 2/2] Running tests...
...

Start 1: classical_dft_tests

1: Test command: /Users/mduranol/Dev/github/classicalDFT/dft_lib/build/release/tests/classical_dft_tests
1: Test timeout computed to be: 10000000
1: [==========] Running tests
....
1: [  PASSED  ] N tests.
1/1 Test #1: classical_dft_tests ..............   Passed    0.01 sec
100% tests passed, 0 tests failed out of N
Total Test time (real) =   0.01 sec
Built target gtests
```

### Building classicalDFT v1.0

To build and run an example that use it, you need to go through the following steps.

1. Perform the following to setup and compile the library:

   ```bash
   $ cd legacy_lib
   $ ../Config.sh
   $ ../dft_make.sh
   ```

2. To build the model application do the following:

   ```bash
   $ cd ../Examples/Droplet
   $ ./../Config.sh
   $ ./../dft_make.sh
   ```

You can cd into TEST and run using:

```bash
$ cd TEST
$ ../Droplet input.dat
```

#### Remarks

1. The [Config.sh](../../config.sh) command is only run the first time the application is created. It configures things for cmake.
2. Rebuild using [dft_make.sh](../../dft_make.sh). Note that this also takes three possible arguments in any order: **clean**, **debug** and **lib**. Their effects are as follows:
   * **clean** causes a clean build (all objects are first deleted)
   * **debug** performs a debug build
   * **lib** causes the libraries to be rebuilt also.
   * Note that `dft_make.sh debug lib` causes both the app and the libraries to be built in debug mode. Similarly, "dft_make.sh lib" causes both to be built in release mode.

3. The first few lines of [CMakeLists.txt](../../legacy_lib/CMakeLists.txt) contain information that you will want to modify if this is used as a model to create another application.
