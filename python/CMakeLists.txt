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
cotire(pyariadne-algebra-src)
cotire(pyariadne-src)

if(APPLE)
    set(PYTHON_LIBRARIES_TO_USE "-undefined dynamic_lookup")
else()
    set(PYTHON_LIBRARIES_TO_USE ${PYTHON_LIBRARIES})
endif()

add_library(pyariadne $<TARGET_OBJECTS:pyariadne-src>)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne" SUFFIX ".so")
target_link_libraries (pyariadne ariadne ${PYTHON_LIBRARIES_TO_USE} ${BOOST_PYTHON_LIBRARIES})

add_library(pyariadne-numeric SHARED $<TARGET_OBJECTS:pyariadne-numeric-src> $<TARGET_OBJECTS:ariadne-utility> $<TARGET_OBJECTS:ariadne-numeric>)
set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric" SUFFIX ".so")
target_link_libraries (pyariadne-numeric ${PYTHON_LIBRARIES_TO_USE} ${BOOST_PYTHON_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

add_library(pyariadne-algebra $<TARGET_OBJECTS:pyariadne-algebra-src>)
set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra" SUFFIX ".so")
target_link_libraries (pyariadne-algebra ariadne ${PYTHON_LIBRARIES_TO_USE} ${BOOST_PYTHON_LIBRARIES})

add_subdirectory(tests)

