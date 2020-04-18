

# Ariadne <img align="right" src="http://www.ariadne-cps.org/img/ariadne-transparent.png" alt="Ariadne" width="80"/> 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/ariadne-cps/ariadne.svg?branch=master)](https://travis-ci.org/ariadne-cps/ariadne) [![codecov](https://codecov.io/gh/ariadne-cps/ariadne/branch/master/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/ariadne)

Ariadne is a tool for reachability analysis and model checking of hybrid systems. Additionally, it is a framework for rigorous computation featuring arithmetic, linear algebra, calculus, geometry, algebraic and differential equations, and optimization solvers.

### Installation ###

The installation instructions are presented for Ubuntu 20.04 and macOS 10.15 only. However, openSUSE and Fedora are known to be working when using their own package managers. Windows installations are not supported yet.

For the Ubuntu installation, we will refer to packages available on Aptitude. The macOS installation instead will assume you are using the Brew package manager.

The build system is CMake. The library is tested for compilation using gcc and clang.

#### Dependencies

The only required library dependencies of Ariadne are GMP and MPFR. If you want to enable the graphical output you will require Cairo in order to save into png files. Finally, the Python bindings require the Python headers (either Python 2 or 3 are supported). In particular for Python, there is an internal Git submodule dependency on the header-only [pybind11](https://github.com/pybind/pybind11) library. Therefore in order to build the Python interface, Git must be installed even if Ariadne has been downloaded as an archive.

Finally, if you want to build the documentation, you need Doxygen and a working Latex distribution (including the Math packages).

Please note that adding new dependencies after preparing the build environment requires to re-run the CMake command.

Specific instructions for Ubuntu and macOS follow.

##### Ubuntu
Aptitude packages: `cmake pkg-config git libgmp-dev libmpfr-dev libcairo2-dev`

Additional package required for the Python interface: `python3-dev` or `python-dev`.

Additional packages required for documentation: `doxygen doxygen-latex` 

##### macOS
1. Install the Command Line Developer Tools (will also be asked when installing Homebrew) from the Apple Store

2. Install Homebrew from http://brew.sh/ . Homebrew packages required: `cmake git mpfr gmp cairo`

For Cairo support, you may need to set up a permanent variable for the path of pkgconfig by adding the following line in your `~\.bash_profile`:

```
export PKG_CONFIG_PATH=/usr/local/opt/libffi/lib/pkgconfig
```

To allow building the documentation: `brew cask install mactex-no-gui` and `brew install doxygen`.

#### Building

To build the library in a clean way, it is preferable that you set up a build subdirectory:

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
$ make ariadne
```

If you prefer to use the Python interface over the C++ library, you should build the Python module with:


```
$ make pyariadne
```


Optionally, you can also build and run the test suite for the library:

```
$ make tests
$ make test
```

where no error should appear.

To build libraries, tests, examples and tutorials, simply type:

```
$ make
```

To build the documentation, instead use:

```
$ make doc
```

You can access the built documentation from the `docs/html/index.html` file in the build directory.


### Installing globally

To install the library globally, you must do

```
$ make install
```

or

```
$ sudo make install
```

if you require administrator privileges, in particular for a Linux installation. Please note that the installation will build the whole distribution beforehand, hence it is preferable that you first build the binaries without administrator privileges, then install.

To find the installed library under Ubuntu, you may need to set the LD\_LIBRARY\_PATH in the .bashrc file of your home directory:

```
export LD_LIBRARY_PATH=/usr/local/lib
```

### Building executables using Ariadne

The tutorials directory contains two CMake projects that rely on a correct installation of Ariadne. You can copy a project directory in any place on your file system and follow the instructions on the README file inside to check that your installation was successful.

Due to limitations of the C++ standard library on macOS since C++11, you won't be able to build an executable with GCC if the Ariadne library has been built using Clang, and viceversa. Hence on macOS you shall use the same compiler for both Ariadne and any projects that depend on it.

### Contribution guidelines ###

If you would like to contribute to Ariadne, please contact the developers. We are especially interested to hear how the documentation and user interface could be improved.

* Pieter Collins <pieter.collins@maastrichtuniversity.nl>
* Luca Geretti <luca.geretti@univr.it>
