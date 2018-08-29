set(PYTHON_TESTS
    test_import
    test_numeric
    test_linear_algebra
    test_function
    test_calculus
    test_geometry
)

#TODO: Python tests don't always run; need to set correct PYTHONPATH
foreach(test ${PYTHON_TESTS})
    add_test(${test} ./${test}.py)
endforeach()

