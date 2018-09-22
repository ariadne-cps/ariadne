

# Ariadne <img align="right" src="http://www.ariadne-cps.org/img/ariadne-transparent.png" alt="Ariadne" width="80"/> 

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Build Status](https://travis-ci.org/ariadne-cps/development.svg?branch=master)](https://travis-ci.org/ariadne-cps/development) [![codecov](https://codecov.io/gh/ariadne-cps/development/branch/master/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/development)

Ariadne is a tool for reachability analysis and model checking of hybrid systems. Additionally, it is a framework for rigorous computation featuring arithmetic, linear algebra, calculus, geometry, algebraic and differential equations, and optimization solvers.

* This repository contains the main development fork of the tool. For a more stable version with a less sophisticated user interface, see the *release-1.0* repository.
* The latest stable tagged version is internal_release-1.9.0. However, the code in the master branch should always be usable.

### Installation ###

These installation instructions have been tested on Ubuntu 18.04 and macOS 10.13.

For the Ubuntu installation, we will refer to packages available on Aptitude. The macOS installation instead will assume you are using the Brew package manager.

The build system is CMake. The library is tested for compilation using gcc and clang.

#### Dependencies

The only required library dependencies of Ariadne are GMP and MPFR. If you want to enable the graphical output you will require Cairo (to save into png files) and GTK2 (for window display). Finally, the Python bindings require the Python headers (Python 2 or 3 are supported). In particular for Python, there is an internal Git submodule dependency on the header-only [pybind11](https://github.com/pybind/pybind11) library. Therefore in order to fetch the dependency, Git must be installed.

Finally, if you want to build the documentation, you need Doxygen and a working Latex distribution (including the Math packages).

Specific instructions for Ubuntu and macOS follow (documentation packages are excluded).

##### Ubuntu
Aptitude packages: `cmake git libgmp-dev libmpfr-dev libgtk2.0-dev libcairo2-dev`

Additional Aptitude package required for the Python interface: `python3-dev` or `python-dev`.

##### OSX
1. Install the Command Line Developer Tools (will also be asked when installing Homebrew) from the Apple Store

2. Install Homebrew from http://brew.sh/ . Homebrew packages required: `cmake git mpfr gmp gtk cairo`

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


### Installing globally

To install the library globally, you must do

```
$ make install
```

or

```
$ sudo make install
```

if you require administrator privileges, in particular for a Linux installation. Please note that the installation will build the whole distribution beforehand.

To find the installed library under Ubuntu, you may need to set the LD_LIBRARY_PATH in the .bashrc file of your home directory:

```
export LD_LIBRARY_PATH=/usr/local/lib
```

### Building executables using Ariadne

The tutorials directory contains three CMake projects that rely on a correct installation of Ariadne. You can copy a project directory in any place on your file system and follow the instructions on the README file inside to check that your installation was successful.

Due to limitations of the C++ standard library on macOS since C++11, you won't be able to build an executable with GCC if the Ariadne library has been built using Clang, and viceversa. Hence on macOS you shall use the same compiler for both Ariadne and any projects that depend on it.

### Contribution guidelines ###

* If you would like to contribute to Ariadne, please contact the developers. We are especially interested to hear how the documentation and user interface could be improved.

* Pieter Collins <pieter.collins@maastrichtuniversity.nl>
* Luca Geretti <luca.geretti@univr.it>
