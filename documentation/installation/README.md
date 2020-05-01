## Generic build instructions

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

### Building classicalDFT v2.0

This section will be populated shortly as we continue with the development of a refactored version of **classicalDFT** (see [dft_lib](../../dft_lib/)) which we hope will prepare classicalDFT for better collaboration.