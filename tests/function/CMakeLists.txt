set(UNIT_TESTS
    test_function
    test_polynomial
    test_procedure
    test_taylor_function
    test_taylor_model
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne-kernel)
endforeach()

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})