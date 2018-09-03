set(UNIT_TESTS
    test_differential
    test_expansion
    test_matrix
    test_multi_index
    test_series
    test_vector
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne-core)
endforeach()

add_dependencies(tests-core ${UNIT_TESTS})
add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})