#include "pybind11.hpp"

void evolution_submodule(pybind11::module&);
void analysis_submodule(pybind11::module&);

void ariadne_module(pybind11::module& m) {
    evolution_submodule(m);
    analysis_submodule(m);
}
