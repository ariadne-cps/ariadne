# README #

Ariadne is a tool for reachability analysis and model checking of hybrid systems. Additionally, it is a framework for rigorous computation featuring arithmetic, linear algebra, calculus, geometry, algebraic and differential equations, and optimization solvers.

* This repository contains the main development fork of the tool. For a more stable version with a less sophisticated user interface, see the *stable* repository
* The latest semi-stable tagged version is internal_release-1.1.0. However, the code in the master branch should always be usable.

### Installation ###

These installation instructions have been tested on Ubuntu 18.04 and macOS 10.13.

For the Ubuntu installation, we will refer to packages available on Aptitude. The macOS installation instead will assume you are using the Brew package manager.

The build system is CMake. The library is tested for compilation using gcc and clang.

#### Dependencies

The library dependencies of ARIADNE are the following:

##### Ubuntu
Aptitude packages required: `cmake libgmp-dev libmpfr-dev libgtk2.0-dev libcairo2-dev`

Additional Aptitude packages required for the Python interface: `python2.7-dev libboost-python-dev`

##### OSX
1. Install the Command Line Developer Tools (will also be asked when installing Homebrew) from the Apple Store

2. Install Homebrew from http://brew.sh/ . Homebrew packages required: `cmake boost gtk cairo`

Optionally, if you want to build the documentation, you need Doxygen and a working Latex distribution (including the Math packages).

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

At this point, if no error arises, you can build the library `libariadne.so` itself:

```
$ make ariadne
```

If you prefer to use the Python interface over the C++ library, you should build the file `ariadne.so` with:


```
$ make pyariadne
```


Optionally, you can also build and run the test suite for the library:

```
$ make tests; make test
```

where no error should appear.

To build libraries, tests and examples, simply type:

```
$ make
```
or

```
$ make all
```

To build the documentation, use:

```
$ make doc
```


### Installing globally

To install the library globally, you must do

```
$ make install
```

To find the installed library under Ubuntu, you may need to set the LD_LIBRARY_PATH in the .bashrc file:

```
export LD_LIBRARY_PATH=/usr/local/lib
```

The tutorials directory contains two CMake projects that rely on the installation of Ariadne. You can copy a project directory in any place on your file system and follow the instructions on the README file inside.

### Contribution guidelines ###

* If you would like to contribute to Ariadne, please contact the developers. We are especially interested to hear how the documentation and user interface could be improved.

* Pieter Collins <pieter.collins@maastrichtuniversity.nl>
* Luca Geretti <luca.geretti@univr.it>
