find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS python)

find_package(PythonLibs 2.4)

if(PYTHONLIBS_FOUND)
  message("Found Python headers")
  message("  PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
  message("  PYTHON_LIBRARY_DIRS: ${PYTHON_LIBRARY_DIRS}")
  message("  PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")
  include_directories(${PYTHON_INCLUDE_DIRS})
else()
  message(WARNING "Unable to find Python headers. You will not be able to use the Python interface.")
endif()


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
    expression_submodule.cpp
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
