#include "pybind11.hpp"

void evolution_submodule(pybind11::module& module);
void analysis_submodule(pybind11::module& module);

PYBIND11_MODULE(pyariadne, module) {
    evolution_submodule(module);
    analysis_submodule(module);
}
