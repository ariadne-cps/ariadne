#include "pybind11.hpp"

void evolution_submodule(pybind11::module&);
void reasoning_submodule(pybind11::module&);

void ariadne_module(pybind11::module& m) {
    evolution_submodule(m);
    reasoning_submodule(m);
}
