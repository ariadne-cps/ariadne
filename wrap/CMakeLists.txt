find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS python)

find_package(PythonLibs 2.4)
message("  PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
message("  PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

IF(PYTHONLIBS_FOUND)
  include_directories(${PYTHON_INCLUDE_DIRS})
ELSE()
  message(FATAL_ERROR "Unable to find Python headers")
ENDIF()


add_library(pynumeric-src OBJECT
    numeric_module.cc
    numeric_submodule.cc
)

add_library(pyalgebra-src OBJECT
    utilities.cc
    algebra_module.cc
    numeric_submodule.cc
    linear_algebra_submodule.cc
    differentiation_submodule.cc
    calculus_submodule.cc
)


add_library(pyariadne-src OBJECT
    utilities.cc
    ariadne_module.cc
    numeric_submodule.cc
    linear_algebra_submodule.cc
    optimization_submodule.cc
    differentiation_submodule.cc
    function_submodule.cc
    calculus_submodule.cc
    geometry_submodule.cc
    solver_submodule.cc
    storage_submodule.cc
    system_submodule.cc
    evolution_submodule.cc
    graphics_submodule.cc
)

cotire(pynumeric-src)
cotire(pyalgebra-src)
cotire(pyariadne-src)
