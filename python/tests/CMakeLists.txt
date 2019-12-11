include(FindPythonInterp)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pytest --version RESULT_VARIABLE PYTEST_result OUTPUT_QUIET ERROR_QUIET)
if(${PYTEST_result} EQUAL 0)
  set(PYTHON_TESTS ON)
  message(STATUS "pytest package found, Python tests will be performed.")
else()
  message(WARNING "pytest package not available, will not perform Python tests.")
endif()

if(PYTHON_TESTS)

    set(TESTS_PYTHON
        import
        numeric
        linear_algebra
        function
        calculus
        solvers
        symbolic
    )

    foreach(TESTNAME ${TESTS_PYTHON})
        add_custom_target(test_python_${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_${TESTNAME}.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        add_test(NAME test_python_${TESTNAME} COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_${TESTNAME}.py WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        set_property(TEST test_python_${TESTNAME} PROPERTY ENVIRONMENT "PYTHONPATH=${PROJECT_BINARY_DIR}/python")
        set_property(TEST test_python_${TESTNAME} PROPERTY LABELS "python;primary")
    endforeach()

    add_dependencies(tests pyariadne)

endif()
