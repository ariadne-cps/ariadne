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
    numeric_module.cpp
    numeric_submodule.cpp
)

add_library(pyalgebra-src OBJECT
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

cotire(pynumeric-src)
cotire(pyalgebra-src)
cotire(pyariadne-src)
