set(UNIT_TESTS
    test_logical
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne-foundations utility)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "foundations;core;kernel")
endforeach()

add_dependencies(tests-core ${UNIT_TESTS})
add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
