# Ariadne

[![Ubuntu](https://github.com/ariadne-cps/ariadne/actions/workflows/ubuntu.yml/badge.svg?branch=main)](https://github.com/ariadne-cps/ariadne/actions/workflows/ubuntu.yml?query=branch%3Amain) [![macOS](https://github.com/ariadne-cps/ariadne/actions/workflows/macos.yml/badge.svg?branch=main)](https://github.com/ariadne-cps/ariadne/actions/workflows/macos.yml?query=branch%3Amain) [![Windows](https://github.com/ariadne-cps/ariadne/actions/workflows/win.yml/badge.svg?branch=main)](https://github.com/ariadne-cps/ariadne/actions/workflows/win.yml?query=branch%3Amain)

Ariadne provides rigorous computation, numerical evolution, reachability
analysis and verification of dynamical systems.

This repository builds the top-level `ariadne` library from
[Reasoning](https://github.com/ariadne-cps/reasoning) and its transitive
dependencies. When a compatible Python development environment is available,
it also builds the combined `pyariadne` extension from the Reasoning Python
bindings.

## Dependencies

The direct Git submodules are [Build](https://github.com/ariadne-cps/build),
which bootstraps CMake, and [Reasoning](https://github.com/ariadne-cps/reasoning).
Other dependencies, including [Ariadne Testing](https://github.com/ariadne-cps/testing)
for the upstream C++ suites, are supplied transitively by Reasoning.
This top-level aggregator defines no separate C++ test executables.
CMake verifies that the direct and transitive Build commits agree.

## Build

```bash
git clone --branch main --recurse-submodules https://github.com/ariadne-cps/ariadne.git
cd ariadne
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target ariadne
```

If Python 3.11 or newer with development headers is available, build the
optional extension with `cmake --build build --target pyariadne`.
To disable Python bindings explicitly, configure with
`-DBUILD_PYTHON_BINDINGS=OFF`.

## Install and uninstall

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel
cmake --install build
cmake --build build --target uninstall
```

The uninstall target removes files listed in the CMake installation manifest.

## License

Ariadne is released under the GNU General Public License v3.0.
