set(UNIT_TESTS
    test_enclosure
    test_iterated_map
    test_iterated_map_evolver
    test_vector_field
    test_vector_field_evolver
    test_vector_field_simulator
    test_differential_inclusion
    test_differential_inclusion_evolver
    test_finite_time_reachability
    test_infinite_time_reachability
    test_verify_safety
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "dynamics;primary")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
