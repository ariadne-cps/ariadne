set(UNIT_TESTS
    test_continuous_evolution
    test_differential_inclusion
    test_discrete_evolution
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    add_test(${TEST} ${TEST})
    target_link_libraries(${TEST} ariadne)
endforeach()

add_dependencies(tests ${UNIT_TESTS})