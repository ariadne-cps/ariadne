

# Ariadne <img align="right" src="http://www.ariadne-cps.org/img/ariadne-transparent.png" alt="Ariadne" width="80"/> 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://github.com/ariadne-cps/ariadne/workflows/Continuous%20Integration/badge.svg)](https://github.com/ariadne-cps/ariadne/actions) [![codecov](https://codecov.io/gh/ariadne-cps/ariadne/branch/master/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/ariadne)

Ariadne is a tool for reachability analysis and model checking of hybrid systems. Additionally, it is a framework for rigorous computation featuring arithmetic, linear algebra, calculus, geometry, algebraic and differential equations, and optimization solvers.

## Installation ##

The installation instructions are presented for Ubuntu Linux systems and derivatives (using Aptitude) and macOS systems (using Homebrew). However, openSUSE and Fedora are known to be working when using their own package managers. Windows installations are not supported yet.

Official packages are available for Ubuntu derivatives and macOS, but if your architecture or particular setup necessarily requires compilation from sources, instructions are provided. The build system used is CMake. The library is tested for compilation using gcc (minimum required: 10.2) and clang (minimum required: 11.0). AppleClang currently does not support C++20 Concepts yet and therefore is not usable on macOS at the moment.

### Official packages

Supplied packages are published using the official Launchpad platform for Ubuntu, and a custom Homebrew tap repository for macOS. Packages simply require all the dependencies, namely: MPFR, Cairo and Python 3.

#### Ubuntu

In order to install the Aptitude package, first you need to import the ppa repository to the ppa list:

```
sudo add-apt-repository ppa:ariadne-cps/ariadne
```

Then you can install Ariadne (along with any missing dependencies):

```
sudo apt-get install ariadne
```

which will be updated with the latest release as with other Aptitude packages.

#### macOS

In order to install the Homebrew package, you just need to

```
brew install ariadne-cps/tap/ariadne
```

which in one line both sets up the "tap" for Ariadne and installs the package along with any missing dependencies. The package will be upgraded with any new versions after a `brew upgrade` is issued. The package currently supports x86-64 architectures on macOS Big Sur. Other configurations trigger an automatic build, therefore you should not need to deal with sources any time.

### Dependencies

If installed from sources, the only required library dependency is MPFR. To enable the graphical output you will require Cairo in order to save into png files. Finally, the Python bindings require the Python headers (version 3 is only supported, since version 2 is discontinued). In particular for Python, there is an internal Git submodule dependency on the header-only [pybind11](https://github.com/pybind/pybind11) library. Therefore in order to build the Python interface, Git must be installed even if Ariadne has been downloaded as an archive. Download of the dependency is automatic though.

Finally, if you want to build the documentation, you need Doxygen and a working Latex distribution (including the Math packages).

Please note that adding new dependencies after preparing the build environment requires to re-run the CMake command.

Specific instructions for Ubuntu and macOS follow, starting from installation from pre-compiled packages.

#### Ubuntu

Aptitude packages: `cmake pkg-config git libmpfr-dev libcairo2-dev` and either `clang-11` or `g++-10` for the compiler toolchain.

Additional package required for the Python interface: `python3-dev`.

Additional packages required for documentation: `doxygen doxygen-latex` 

#### macOS

Homebrew packages: `cmake git mpfr cairo` and `gcc@10` if using GCC.

For Cairo support, you may need to set up a permanent variable for the path of pkgconfig by adding the following line in your `~\.bash_profile`:

```
export PKG_CONFIG_PATH=/usr/local/opt/libffi/lib/pkgconfig
```

To allow building the documentation: `brew cask install mactex-no-gui` and `brew install doxygen`.

### Building

To build the library from sources in a clean way, it is preferable that you set up a build subdirectory, say:

```
$ mkdir build
 $ cd build
```

Then you can prepare the build environment, choosing a Release build for maximum performance:

```
$ cmake .. -DCMAKE_BUILD_TYPE=Release
```

At this point, if no error arises, you can build the library itself:

```
$ cmake --build . --target ariadne --parallel
```

If you prefer to use the Python interface over the C++ library, you should build the Python module with:

```
$ cmake --build . --target pyariadne --parallel
```

Optionally, you can also build and run the test suite for the library:

```
$ cmake --build . --target tests --parallel
 $ ctest
```

where no error should appear.

To build libraries, tests, examples and tutorials, simply type:

```
$ cmake --build . --target everything --parallel
```

To build the documentation, instead use:

```
$ cmake --build . --target doc --parallel
```

You can access the built documentation from the `docs/html/index.html` file in the build directory.


### Installing globally

To install the library globally from built sources, you must do

```
$ cmake --build . --target install --parallel
```

using `sudo` if you require administrator privileges for a Linux installation. Please note that the installation will build the whole distribution beforehand, hence it is preferable that you first build the other targets without administrator privileges, build the install target.

To find the installed library under Ubuntu, you may need to set the LD\_LIBRARY\_PATH in the .bashrc file of your home directory:

```
export LD_LIBRARY_PATH=/usr/local/lib
```

### Building executables using Ariadne

The tutorials directory contains two CMake projects that rely on a correct installation of Ariadne, either by using a package or by building the sources. You can copy a project directory in any place on your file system and follow the instructions on the README file inside to check that your installation was successful.

Due to limitations of the C++ standard library on macOS since C++11, you won't be able to build an executable with GCC if the Ariadne library has been built using Clang, and viceversa. Hence on macOS you shall use the same compiler for both Ariadne and any projects that depend on it. If Ariadne comes from the Homebrew package, then it has been built using g++ 10.

## Contribution guidelines ##

If you would like to contribute to Ariadne, please contact the developers. We are especially interested to hear how the documentation and user interface could be improved.

* Luca Geretti <luca.geretti@univr.it>
* Pieter Collins <pieter.collins@maastrichtuniversity.nl>
