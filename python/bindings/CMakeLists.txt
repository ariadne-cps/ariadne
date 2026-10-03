include_directories(${Python3_INCLUDE_DIRS})

get_target_property(PYTHON_COMMON_SOURCE_DIR pyariadne-python-common SOURCE_DIR)
get_target_property(NUMERIC_SOURCE_DIR numeric SOURCE_DIR)
get_target_property(INTERVAL_SOURCE_DIR interval SOURCE_DIR)
include_directories(
    "${PYTHON_COMMON_SOURCE_DIR}/include"
    "${NUMERIC_SOURCE_DIR}/python/bindings"
    "${INTERVAL_SOURCE_DIR}/python/bindings"
)

message(STATUS "Please note that you need to use the Python 3 interpreter to import pyariadne correctly.")

add_library(pyariadne-algebra-submodule-obj OBJECT
    utilities.cpp
    linear_algebra_submodule.cpp
    differentiation_submodule.cpp
    calculus_submodule.cpp
)

add_library(pyariadne-extra-submodule-obj OBJECT
    symbolic_submodule.cpp
    optimization_submodule.cpp
    function_submodule.cpp
    geometry_submodule.cpp
    io_submodule.cpp
    solver_submodule.cpp
    storage_submodule.cpp
    system_submodule.cpp
    evolution_submodule.cpp
    graphics_submodule.cpp
    hybrid_submodule.cpp
)

add_library(pyariadne-module-obj OBJECT
    ariadne_module.cpp
)

foreach(target
    pyariadne-algebra-submodule-obj
    pyariadne-extra-submodule-obj
    pyariadne-module-obj
)
    target_link_libraries(${target} PRIVATE
        ariadne
        pyariadne-python-common
        pybind11::module
    )
    set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endforeach()
