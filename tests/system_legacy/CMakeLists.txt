# Legacy model-only regressions retained until equivalent coverage is verified
# in the standalone system repository. No changes to ss/sys are required.
set(UNIT_TESTS
    test_differential_inclusion
    test_iterated_map
    test_vector_field
)

foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "system")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
