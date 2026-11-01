# Ariadne

Ariadne provides rigorous computation, numerical evolution, reachability
analysis and verification of dynamical systems.

This repository builds the top-level `ariadne` library from
[Reasoning](https://github.com/ariadne-cps/reasoning) and its transitive
dependencies. When a compatible Python development environment is available,
it also builds the combined `pyariadne` extension from the Reasoning Python
bindings.

## Dependencies

The direct Git submodules are [Configuration](https://github.com/ariadne-cps/configuration),
required to bootstrap CMake, and Reasoning. Other dependencies are
provided transitively by Reasoning.

## Build

```bash
git clone --branch decoupling --recurse-submodules https://github.com/ariadne-cps/ariadne.git
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
