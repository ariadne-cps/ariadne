
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(pyariadne-numeric-src OBJECT
    numeric_module.cpp
    numeric_submodule.cpp
)

add_library(pyariadne-algebra-src OBJECT
    utilities.cpp
    algebra_module.cpp
    numeric_submodule.cpp
    linear_algebra_submodule.cpp
    differentiation_submodule.cpp
    calculus_submodule.cpp
)


add_library(pyariadne-src OBJECT
    utilities.cpp
    ariadne_module.cpp
    numeric_submodule.cpp
    symbolic_submodule.cpp
    linear_algebra_submodule.cpp
    optimization_submodule.cpp
    differentiation_submodule.cpp
    function_submodule.cpp
    calculus_submodule.cpp
    geometry_submodule.cpp
    solver_submodule.cpp
    storage_submodule.cpp
    system_submodule.cpp
    evolution_submodule.cpp
    graphics_submodule.cpp
)

cotire(pyariadne-numeric-src)
cotire(pyariadne-src)
