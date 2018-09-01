include(FindPythonInterp)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pytest --version RESULT_VARIABLE PYTEST_result OUTPUT_QUIET ERROR_QUIET)
if(${PYTEST_result} EQUAL 0)
  set(PYTHON_TESTS ON)
  message(STATUS "pytest package found, Python tests will be performed.")
else()
  message(WARNING "pytest package not available, will not perform Python tests.")
endif()

if(PYTHON_TESTS)
    add_test(NAME test_python COMMAND ${PYTHON_EXECUTABLE} -m pytest WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    SET_PROPERTY(TEST test_python PROPERTY ENV "PYTHONPATH=${PYTHONPATH}:${PROJECT_BINARY_DIR}/python" )
endif()
