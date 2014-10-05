find_package(Boost ${Boost_REQUIRED_VERSION} COMPONENTS python)

find_package(PythonLibs 2.4)
message("  PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")
message("  PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

IF(PYTHONLIBS_FOUND)
  include_directories(${PYTHON_INCLUDE_DIRS})
ELSE()
  message(FATAL_ERROR "Unable to find Python headers")
ENDIF()


add_library(pyinterval
    interval_module.cc numeric_submodule.cc linear_algebra_submodule.cc
)
set_target_properties(pyinterval PROPERTIES PREFIX "" OUTPUT_NAME "interval")
target_link_libraries (pyinterval ariadne boost_python)

add_library(pydifferentiation
    differentiation_module.cc differentiation_submodule.cc
)
set_target_properties(pydifferentiation PROPERTIES PREFIX "" OUTPUT_NAME "differentiation")
target_link_libraries (pydifferentiation ariadne boost_python)

add_library(pyariadne
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
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne")
target_link_libraries (pyariadne ariadne-extra ariadne boost_python)

