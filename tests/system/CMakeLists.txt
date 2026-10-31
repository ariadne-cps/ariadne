# system layer behavioral tests (moved without changing test logic).
set(UNIT_TESTS
    test_hybrid_automaton
    test_first_order_pde_model_header
)

foreach(TEST IN LISTS UNIT_TESTS)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "hybrid")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
