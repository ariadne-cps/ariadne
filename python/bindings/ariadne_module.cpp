#include "pybind11.hpp"

void reasoning_submodule(pybind11::module&);

void ariadne_module(pybind11::module& m) {
    reasoning_submodule(m);
}
