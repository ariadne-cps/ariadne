set(PYTHON_TESTS 
    test_function
    test_linear_algebra
    test_geometry  
    test_numeric
    test_calculus   
    test_import
)

#TODO: Python tests don't always run; need to set correct PYTHONPATH
foreach(test ${PYTHON_TESTS})
    add_test(${test} ./${test}.py)
endforeach()

