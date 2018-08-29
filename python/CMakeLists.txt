include(FindPythonInterp)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pytest --version RESULT_VARIABLE PYTEST_result)
if(${PYTEST_result} EQUAL 0)
  set(PYTHON_TESTS ON)
  message(STATUS "pytest package found, Python tests will be performed.")
else()
  message(WARNING "pytest package not available, will not perform Python tests.")
endif()

if(PYTHON_TESTS)
    file(GLOB PYTHON_TEST_FILES ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.py)
    file(COPY ${PYTHON_TEST_FILES} DESTINATION ${PROJECT_BINARY_DIR}/python)
    add_test(NAME test_python COMMAND ${PYTHON_EXECUTABLE} -m pytest WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/python)
endif()

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

add_library(pyariadne $<TARGET_OBJECTS:pyariadne-src>)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne" SUFFIX ".so")
target_link_libraries (pyariadne ariadne ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIBRARIES})

add_library(pyariadne-numeric SHARED $<TARGET_OBJECTS:pyariadne-numeric-src> $<TARGET_OBJECTS:ariadne-utility> $<TARGET_OBJECTS:ariadne-numeric>)
set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric" SUFFIX ".so")
target_link_libraries (pyariadne-numeric ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

add_library(pyariadne-algebra $<TARGET_OBJECTS:pyariadne-algebra-src>)
set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra" SUFFIX ".so")
target_link_libraries (pyariadne-algebra ariadne ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIBRARIES})


