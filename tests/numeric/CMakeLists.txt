set(UNIT_TESTS
    test_dyadic
    test_float
    test_integer
    test_logical
    test_number
    test_rational
    test_real    
    test_rounding
    test_validated_float
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne-core)    
endforeach()

add_dependencies(tests-core ${UNIT_TESTS})
add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})