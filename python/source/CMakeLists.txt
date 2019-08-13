include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(SYSTEM ${PYBIND11_INCLUDE_DIR})

message(STATUS "You will have to use the Python ${PYTHON_VERSION_MAJOR} interpreter to import Ariadne correctly.")

add_library(pyariadne-numeric-submodule-obj OBJECT
    numeric_submodule.cpp
)

add_library(pyariadne-numeric-module-obj OBJECT
    numeric_module.cpp
)

add_library(pyariadne-algebra-submodule-obj OBJECT
    utilities.cpp
    linear_algebra_submodule.cpp
    differentiation_submodule.cpp
    calculus_submodule.cpp
)

add_library(pyariadne-algebra-module-obj OBJECT
    algebra_module.cpp
)

add_library(pyariadne-extra-submodule-obj OBJECT  
    symbolic_submodule.cpp
    optimization_submodule.cpp
    function_submodule.cpp
    geometry_submodule.cpp
    solver_submodule.cpp
    storage_submodule.cpp
    system_submodule.cpp
    evolution_submodule.cpp
    graphics_submodule.cpp
)

add_library(pyariadne-module-obj OBJECT
    ariadne_module.cpp
)

cotire(pyariadne-numeric-submodule-obj CONFIGURATIONS None)
cotire(pyariadne-algebra-submodule-obj CONFIGURATIONS None)
cotire(pyariadne-extra-submodule-obj CONFIGURATIONS None)

if(COVERAGE)
    target_link_libraries(pyariadne-numeric-submodule-obj PUBLIC coverage_config)
    target_link_libraries(pyariadne-numeric-module-obj PUBLIC coverage_config)
    target_link_libraries(pyariadne-algebra-submodule-obj PUBLIC coverage_config)
    target_link_libraries(pyariadne-algebra-module-obj PUBLIC coverage_config)
    target_link_libraries(pyariadne-extra-submodule-obj PUBLIC coverage_config)
    target_link_libraries(pyariadne-module-obj PUBLIC coverage_config)
endif()